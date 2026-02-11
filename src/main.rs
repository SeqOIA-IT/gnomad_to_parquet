use tikv_jemallocator::Jemalloc;

#[global_allocator]
static GLOBAL: Jemalloc = Jemalloc;


use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::ops::Index;
use std::time::{Instant};

use std::sync::Arc;

use anyhow::anyhow;
use parquet::file::properties::WriterProperties;
use parquet::file::properties::WriterVersion;
use arrow::record_batch::RecordBatch;
use parquet::arrow::arrow_writer::ArrowWriter;
use arrow::datatypes::*;
use arrow::array::*;
use parquet::basic::Compression;

use flate2::bufread::MultiGzDecoder;

use std::env;
use std::path::PathBuf;
use log::{error,info,warn};

fn main() -> anyhow::Result<()> {
    env_logger::init();

    let args: Vec<String> = env::args().collect();
    let input_vcf_file = args.get(1).unwrap();
    // ex : /annotations/Human/GRCh38/DB/public/gnomad/2.1.1/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz
    process(        
        
        input_vcf_file
        )    
}



fn process(
    input_vcf:&str    
) -> anyhow::Result<()> {
        
        
        
        
        let reader : Box<dyn BufRead> = if input_vcf.ends_with("gz") {
            let f = BufReader::new(File::open(&input_vcf).map_err(|e| anyhow::anyhow!("erreur fichier {} : {} ",input_vcf,e))?);
            let gz = MultiGzDecoder::new(f);
            Box::new(BufReader::new(gz))                        
        } else {
                Box::new(BufReader::with_capacity(
                100 * 128 * 1024,
                File::open(&input_vcf).map_err(|e| anyhow::anyhow!("erreur fichier {} : {} ",input_vcf,e))?,
            ))
        };

        

        let mut count:u32 = 0;

        let mut it = reader.lines();
        
        
        let mut infos = HashMap::new();
        

        while let Some(Ok(l2)) = it.next() {
            if l2.starts_with("#CHROM") {
                break;
            }

            if l2.starts_with("##INFO") {
                info!("{}",l2);
                let info_f = parse_info(&l2).unwrap();
                infos.insert(info_f.name.to_string(),info_f);    
            }
            
        }
        let mut infos_keys_x = infos.keys().collect::<Vec<_>>();
        infos_keys_x.sort();
        /* 
        infos_keys = infos_keys.into_iter().filter(|x| 
            if x.starts_with("VRS") || *x=="vep" { // on ne parse pas
            false
        } else {
            true
        }
        ).collect();
        */

        // on ne prends que les clés {annot}_{pop}_{sexe}
        let mut special = Vec::new();
        let annots = vec!["AC","AF","AN","nhomalt"];
        let pops = vec!["","_afr","_ami","_amr","_asj","_eas","_fin","_nfe","_remaining","_sas","_mid"];        
        let sexes = vec!["","_female","_male"];

        for annot in annots.iter() {
            for pop in pops.iter() {
                for sexe in sexes.iter() {
                    special.push(format!("{}{}{}",annot,pop,sexe))
                }
            }
        }

        let mut special_parquet_name = Vec::new();
        let annots_parquet_name = vec!["AC","AF","AN","nhomalt"];
        let pops_parquet_name = vec!["","_afr","_ami","_amr","_asj","_eas","_fin","_nfe","_remaining","_sas","_mid"];        
        let sexes_parquet_name = vec!["","_XX","_XY"];

        for annot_parquet_name in annots_parquet_name.iter() {
            for pop_parquet_name in pops_parquet_name.iter() {
                for sexe_parquet_name in sexes_parquet_name.iter() {
                    special_parquet_name.push(format!("{}{}{}",annot_parquet_name,pop_parquet_name,sexe_parquet_name))
                }
            }
        }



        let infos_keys:Vec<(&String,&String)> = infos_keys_x.into_iter().filter(|x| {
            
            if special.contains(x) {
                true
            } else {                
                println!("info : key {} in vcf but not selected for parquet",x);
                false
            }
            
        }).map(|x| (x,special_parquet_name.get(special.iter().position(|n| n == x).unwrap()).unwrap())).collect();

        println!("--------");
        
        special.iter().for_each(|target| {
            
            if infos_keys.iter().find(|(x,y)| target==*y).is_none() {
                println!("warning : desired key {} not in vcf",target);
            }

        });



        
        let mut current_chrom=1;
        let mut file = std::fs::File::create(format!("gnomad_v4_chrom{}.parquet",current_chrom)).unwrap();
        let batch_size = 100_000;
        // Default writer properties
        let props = WriterProperties::builder()
                    //.set_write_batch_size(batch_size)
                    .set_writer_version(WriterVersion::PARQUET_1_0)
                    .set_statistics_enabled(true)
                    .set_max_row_group_size(batch_size)
                    .set_compression(Compression::SNAPPY)
                    .build();

        let mut fields = Vec::new();        
        fields.push(Field::new("chrom", DataType::Int32, false));        
        fields.push(Field::new("pos", DataType::Int32, false));
        fields.push(Field::new("ref_alt_hash", DataType::Int64, false));        
        fields.push(Field::new("a_ref", DataType::Utf8, false));
        fields.push(Field::new("a_alt", DataType::Utf8, false));

        for (vcf_key,parquet_key) in &infos_keys {


            let info_f = infos.get(*vcf_key).unwrap();
            match info_f.value_type {
                ValueType::Flag => fields.push(Field::new(parquet_key, DataType::Boolean, true)),
                ValueType::Float => fields.push(Field::new(parquet_key, DataType::Float64, true)),
                ValueType::String => fields.push(Field::new(parquet_key, DataType::Utf8, true)),
                ValueType::Integer => fields.push(Field::new(parquet_key, DataType::Int32, true)),
            }
        
            
        }

        let mut contig_builder = Int32Builder::new(batch_size);        
        let mut pos_builder = Int32Builder::new(batch_size);
        let mut ref_alt_hash_builder = Int64Builder::new(batch_size);
        
        let mut a_ref_builder = StringBuilder::new(batch_size);
        let mut a_alt_builder = StringBuilder::new(batch_size);

        let mut fields_buffer = Vec::new();

        for (vcf_key,parquet_key) in infos_keys.iter() {
            let info = infos.get(*vcf_key).unwrap();

            let builder = match info.value_type {
                
                ValueType::String =>  Box::new(StringBuilder::new(batch_size)) as Box<dyn ArrayBuilder>,
                ValueType::Float =>  Box::new(Float64Builder::new(batch_size)) as Box<dyn ArrayBuilder>,
                ValueType::Integer =>  Box::new(Int32Builder::new(batch_size)) as Box<dyn ArrayBuilder>,
                ValueType::Flag =>  Box::new(BooleanBuilder::new(batch_size)) as Box<dyn ArrayBuilder>,
                
            };
           
            
            fields_buffer.push(builder);

        }

        let schema = Arc::new(Schema::new(fields));
        
        let mut writer = ArrowWriter::try_new(file, schema.clone(), Some(props.clone())).unwrap();


        
        
        
        let mut last_pos = 1;


        while let Some(Ok(l2)) = it.next() {
            let one_line = l2.split("\t").collect::<Vec<_>>();
            //let info_fields=;
            let chrom_r= compute_contig_num(one_line[0]);

            if chrom_r.is_err() {
                warn!("can't parse chrom {} ",one_line[0]);
                continue;
            }

            let chrom = chrom_r?;

            if chrom!=current_chrom {
                // on change de chrom
                if contig_builder.len()>0 {
                    // write to parquet
                    let mut rb = Vec::new();                        
                    rb.push(Arc::new(contig_builder.finish()) as ArrayRef);        
                    rb.push(Arc::new(pos_builder.finish()) as ArrayRef);
                    rb.push(Arc::new(ref_alt_hash_builder.finish()) as ArrayRef);
                    rb.push(Arc::new(a_ref_builder.finish()) as ArrayRef);
                    rb.push( Arc::new(a_alt_builder.finish()) as ArrayRef);
                    
                    for mut f in fields_buffer.into_iter() {                   
                        rb.push( f.finish() as ArrayRef);                   
                    }
                    let batch = RecordBatch::try_new(schema.clone(),rb).unwrap();
        
                    
                    writer.write(&batch).expect("Writing batch");
        
                    
                }
        
                writer.close()?;

                current_chrom = chrom;
                file = std::fs::File::create(format!("gnomad_v4_chrom{}.parquet",current_chrom)).unwrap();
                writer = ArrowWriter::try_new(file, schema.clone(), Some(props.clone())).unwrap();
                last_pos = 1;
                contig_builder = Int32Builder::new(batch_size);        
                pos_builder = Int32Builder::new(batch_size);
                ref_alt_hash_builder = Int64Builder::new(batch_size);
                a_ref_builder = StringBuilder::new(batch_size);
                a_alt_builder = StringBuilder::new(batch_size);
                fields_buffer = Vec::new();

                for (vcf_key,parquet_key) in infos_keys.iter() {
                    let info = infos.get(*vcf_key).unwrap();
        
                    let builder = match info.value_type {
                        
                        ValueType::String =>  Box::new(StringBuilder::new(batch_size)) as Box<dyn ArrayBuilder>,
                        ValueType::Float =>  Box::new(Float64Builder::new(batch_size)) as Box<dyn ArrayBuilder>,
                        ValueType::Integer =>  Box::new(Int32Builder::new(batch_size)) as Box<dyn ArrayBuilder>,
                        ValueType::Flag =>  Box::new(BooleanBuilder::new(batch_size)) as Box<dyn ArrayBuilder>,
                        
                    };
                   
                    
                    fields_buffer.push(builder);
        
                }

            }
            let pos = one_line[1].parse::<u32>()?;
            let a_ref = one_line[3].to_string();
            let a_alt = one_line[4].to_string();

            
            if pos<last_pos {
                return Err(anyhow!("current pos<last_pos, donc non trié !"))
            }

            
            last_pos = pos;

            let h = parse_all_info(&one_line.get(7).unwrap(),&infos,&infos_keys)?;
            
            
            contig_builder.append_value(chrom as i32).unwrap();                
            pos_builder.append_value(pos as i32).unwrap();
            use xxhash_rust::xxh64::xxh64;
            let hash = xxh64(format!("{}-{}",a_ref,a_alt).as_bytes(),0);
            ref_alt_hash_builder.append_value(hash as i64).unwrap();
            
            a_ref_builder.append_value(a_ref).unwrap();
            a_alt_builder.append_value(a_alt).unwrap();
            
            for (index,(vcf_key,parquet_key)) in infos_keys.iter().enumerate() {
                let  b = &mut fields_buffer.get_mut(index).unwrap();

                let info = infos.get(*vcf_key).unwrap();
                let value = h.get(*vcf_key);
    
                match info.value_type {
                    
                    
                    ValueType::Float =>  {
                        let builder : &mut arrow::array::Float64Builder=  b.as_any_mut().downcast_mut::<Float64Builder>().unwrap();
                        match value {
                            Some(ValueParsed::Float(n)) => builder.append_value(*n).unwrap(),
                            None | Some(ValueParsed::Null) => builder.append_null().unwrap(),
                            _=> return Err(anyhow!("value shoud be float for field {}",vcf_key))
                        }
                    }
                    ValueType::Integer =>  {
                        let builder : &mut arrow::array::Int32Builder=  b.as_any_mut().downcast_mut::<Int32Builder>().unwrap();
                        match value {
                            Some(ValueParsed::Integer(n)) => builder.append_value(*n).unwrap(),
                            None | Some(ValueParsed::Null)  => builder.append_null().unwrap(),
                            _=> return Err(anyhow!("value shoud be integer for field {}",vcf_key))
                        }
                    },
                    ValueType::String =>  {
                        let builder : &mut arrow::array::StringBuilder=  b.as_any_mut().downcast_mut::<StringBuilder>().unwrap();
                        match value {
                            Some(ValueParsed::String(n)) => builder.append_value(n).unwrap(),
                            None | Some(ValueParsed::Null)  => builder.append_null().unwrap(),
                            _=> return Err(anyhow!("value shoud be string for field {}",vcf_key))
                        }
                    },
                    ValueType::Flag =>  {
                        let builder : &mut arrow::array::BooleanBuilder=  b.as_any_mut().downcast_mut::<BooleanBuilder>().unwrap();
                        match value {
                            Some(ValueParsed::Flag(n)) => builder.append_value(*n).unwrap(),
                            None | Some(ValueParsed::Null)  => builder.append_null().unwrap(),
                            _=> return Err(anyhow!("value shoud be flag for field {}",vcf_key))
                        }
                    },

                    
                    
                };
               
                
              
            }
            
            if contig_builder.len()==batch_size {
                // write to parquet
                let mut rb = Vec::new();                        
                rb.push(Arc::new(contig_builder.finish()) as ArrayRef);        
                rb.push(Arc::new(pos_builder.finish()) as ArrayRef);
                rb.push(Arc::new(ref_alt_hash_builder.finish()) as ArrayRef);
                rb.push(Arc::new(a_ref_builder.finish()) as ArrayRef);
                rb.push( Arc::new(a_alt_builder.finish()) as ArrayRef);
                
                for mut f in fields_buffer.into_iter() {                   
                    rb.push( f.finish() as ArrayRef);                   
                }
                let batch = RecordBatch::try_new(schema.clone(),rb).unwrap();

                
                writer.write(&batch).expect("Writing batch");


                contig_builder = Int32Builder::new(batch_size);        
                pos_builder = Int32Builder::new(batch_size);
                ref_alt_hash_builder = Int64Builder::new(batch_size);
                a_ref_builder = StringBuilder::new(batch_size);
                a_alt_builder = StringBuilder::new(batch_size);
                fields_buffer = Vec::new();

                for (vcf_key,parquet_key) in infos_keys.iter() {
                    let info = infos.get(*vcf_key).unwrap();
        
                    let builder = match info.value_type {
                        
                        ValueType::String =>  Box::new(StringBuilder::new(batch_size)) as Box<dyn ArrayBuilder>,
                        ValueType::Float =>  Box::new(Float64Builder::new(batch_size)) as Box<dyn ArrayBuilder>,
                        ValueType::Integer =>  Box::new(Int32Builder::new(batch_size)) as Box<dyn ArrayBuilder>,
                        ValueType::Flag =>  Box::new(BooleanBuilder::new(batch_size)) as Box<dyn ArrayBuilder>,
                        
                    };
                   
                    
                    fields_buffer.push(builder);
        
                }
            }
            count+=1;
            
        }

        if contig_builder.len()>0 {
            // write to parquet
            let mut rb = Vec::new();                        
            rb.push(Arc::new(contig_builder.finish()) as ArrayRef);        
            rb.push(Arc::new(pos_builder.finish()) as ArrayRef);
            rb.push(Arc::new(ref_alt_hash_builder.finish()) as ArrayRef);
            rb.push(Arc::new(a_ref_builder.finish()) as ArrayRef);
            rb.push( Arc::new(a_alt_builder.finish()) as ArrayRef);
            
            for mut f in fields_buffer.into_iter() {                   
                rb.push( f.finish() as ArrayRef);                   
            }
            let batch = RecordBatch::try_new(schema.clone(),rb).unwrap();

            
            writer.write(&batch).expect("Writing batch");

            
        }

        writer.close()?;

        println!("count={}",count);
        Ok(())
        
    }

    

    pub fn parse_all_info(info_fields:&str,infos:&HashMap<String,InfoFormatConfig>,infos_keys:&Vec<(&String,&String)>) -> anyhow::Result<HashMap<String,ValueParsed>>{
        let mut h = HashMap::new();
        // maintenant nous parsons les lignes
        let k_vs = info_fields.split(";").collect::<Vec<_>>().iter().map(|kv| {
            let ar = kv.split("=").collect::<Vec<_>>();
            match ar.as_slice() {
                [a,b] => Ok((a.to_string(),b.to_string())),
                [a] => Ok((a.to_string(),"true".to_string())), // une seul valeur donc Flag
                _ => Ok(("dummy".to_string(),kv.to_string())) // return Err(anyhow!("")) // cas vep=xxx|yyy=u| .. debilité !
            }
            
        }
        ).collect::<anyhow::Result<HashMap<_,_>>>()?;
        // liste des [k,v]
        for (k,_x) in infos_keys {
        //for k_v in k_vs {
            if let Some(v) = k_vs.get(*k) {
                h.insert(k.to_string(),convert(k,v,infos)?);
                //let x = k_v.get(0).unwrap();
                /* 
                if x.starts_with("VRS") { // on ne parse pas
                    continue;
                }
                if *x == "vep" { // on ne parse pas
                    continue;
                }
                
            
                match k_v.as_slice() {
                    [k,v] => { h.insert(k.to_string(),convert(k,v,infos)?); () },
                    [k] => { h.insert(k.to_string(),ValueParsed::Flag(true)); () },
                    _ => return Err(anyhow!("not parsable : {:?}",k_v))
                }
                */
            } else {
                // on mets à null la valeur
                h.insert(k.to_string(),ValueParsed::Null);
            }
        }

        Ok(h)
    }
    pub fn convert(k:&str,v:&str,infos:&HashMap<String,InfoFormatConfig>)-> anyhow::Result<ValueParsed> {
        if v=="." {
            return Ok(ValueParsed::Null)
        }
        if let Some(info) = infos.get(k) {
            match info.value_type {
                ValueType::String => return Ok(ValueParsed::String(v.to_string())),
                ValueType::Integer => return Ok(ValueParsed::Integer(v.parse::<i32>()?)),
                ValueType::Float => return Ok(ValueParsed::Float(v.parse::<f64>()?)),
                ValueType::Flag => return Ok(ValueParsed::Flag(true)),                
            }
        } else {
            return Err(anyhow!("info header not found for key {}",k))
        }
    }
    pub fn parse_info(info:&str) -> anyhow::Result<InfoFormatConfig> {
        let parsed : Vec<&str>= info[8..].split(",").collect(); // ID=xx  Number=XX Type=xx
        let mut h = HashMap::new();

        for keyval in parsed {
            match keyval.split("=").collect::<Vec<_>>().as_slice() {
                [k,v] => { h.insert(k.to_string(),v.to_string()); () },
                _ => ()
            }
        }

        let value_type = match h.get("Type").unwrap().to_string().as_str() {
            "Flag" => ValueType::Flag,
            "String" => ValueType::String,
            "Integer" => ValueType::Integer,
            "Float" => ValueType::Float,
            x=> return Err(anyhow!("info type unknown : {}",x))

        };


        Ok(InfoFormatConfig {
            name:h.get("ID").unwrap().to_string(),
            value_number:h.get("Number").unwrap().to_string(),
            value_type:value_type,
        })

        
    }

    #[derive(Debug)]
    pub struct InfoFormatConfig {
        pub name: String,
        pub value_type: ValueType,
        pub value_number: String       
    }

    #[derive(Debug)]
        pub enum ValueParsed {
        String(String),
        Integer(i32),
        Float(f64),        
        Flag(bool),
        Null
    }

    #[derive(Debug)]
    pub enum ValueType {        
        String,
        Integer,
        Float,        
        Flag
    }
    


    pub fn compute_contig_num(chrom_str_p: &str) -> anyhow::Result<u8> {
        
            let chrom_str = if chrom_str_p.starts_with("chr") {
                &chrom_str_p[3..]
            } else {
                chrom_str_p
            };
            let c = chrom_str.parse::<u8>();
            match c {
                Ok(num) if num >= 1 && num <= 22 => {
                    Ok(num)
                },
                Ok(_) => Err(anyhow!("contig not reconize {}",chrom_str)) ,
                Err(_) => match chrom_str.to_uppercase().as_str() {
                    "X" => Ok(23),
                    "Y" => Ok(24),
                    "MT" => Ok(25),
                    "M" => Ok(25),
                    _ => Err(anyhow!("contig not reconize {}",chrom_str)) 
                }
    
            }
           
        
    }


    