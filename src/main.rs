use anyhow::anyhow;
use arrow::array::*;
use arrow::datatypes::*;
use arrow::record_batch::RecordBatch;
use parquet::arrow::arrow_writer::ArrowWriter;
use parquet::basic::Compression;
use parquet::file::properties::WriterVersion;
use parquet::file::properties::{EnabledStatistics, WriterProperties};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::sync::Arc;

use flate2::bufread::MultiGzDecoder;

use tikv_jemallocator::Jemalloc;
#[global_allocator]
static GLOBAL: Jemalloc = Jemalloc;

use clap::Parser;

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Cli {
    #[arg(long, default_value = "gnomad_")]
    prefix_output_filename: String,

    #[arg(long, default_value = ".")]
    output_dir: String,

    #[arg(long)]
    input_fullname: String,

    #[arg(long)]
    only_chrom: Option<String>, // ex: 1
}

fn main() -> anyhow::Result<()> {
    env_logger::init();

    let cli = Cli::parse();

    let ouput_fullname_prefix = format!("{}/{}", cli.output_dir, cli.prefix_output_filename);

    process(&cli.input_fullname, &ouput_fullname_prefix, cli.only_chrom)
}

fn process(
    input_fullname: &str,
    ouput_fullname_prefix: &str,
    only_chrom: Option<String>,
) -> anyhow::Result<()> {
    let reader: Box<dyn BufRead> = if input_fullname.ends_with("gz") {
        let f = BufReader::new(
            File::open(input_fullname)
                .map_err(|e| anyhow::anyhow!("error file {} : {} ", input_fullname, e))?,
        );
        let gz = MultiGzDecoder::new(f);
        Box::new(BufReader::new(gz))
    } else {
        Box::new(BufReader::with_capacity(
            100 * 128 * 1024,
            File::open(input_fullname)
                .map_err(|e| anyhow::anyhow!("error file {} : {} ", input_fullname, e))?,
        ))
    };

    let mut count: u32 = 0;

    let mut it = reader.lines();

    let mut infos = HashMap::new();

    while let Some(Ok(l2)) = it.next() {
        if l2.starts_with("#CHROM") {
            // stop reading the header
            break;
        }

        if l2.starts_with("##INFO") {
            // reader an INFO line
            // reading the annotation format
            let info_f = parse_info(&l2).unwrap();
            infos.insert(info_f.name.to_string(), info_f);
        }
    }

    let infos_keys = get_infos_keys(&infos)?;

    let batch_size = 100_000;

    // construct arrow builders
    let mut fields_buffer = get_arrow_builders(&infos, &infos_keys, batch_size)?;

    let schema = get_schema(&infos, &infos_keys)?;
    let mut current_chrom = if let Some(c) = &only_chrom {
        compute_contig_num(c)?
    } else {
        1
    };
    let mut last_pos = 1;

    let mut writer = get_writer(
        schema.clone(),
        format!("{}{}.parquet", ouput_fullname_prefix, current_chrom),
        batch_size,
    )?;

    while let Some(Ok(l2)) = it.next() {
        let one_line = l2.split("\t").collect::<Vec<_>>();

        let chrom = compute_contig_num(one_line[0])?;
        if only_chrom.is_some() && chrom != current_chrom {
            last_pos = 1; // reset last_pos
            continue; // only_chrom is given, but chrom in line doesn't match, bypass
        }

        if only_chrom.is_none() && chrom != current_chrom {
            // new chrom, write last data
            if !fields_buffer[0].is_empty() {
                // write to parquet
                let mut rb = Vec::new();
                for mut f in fields_buffer.into_iter() {
                    rb.push(f.finish() as ArrayRef);
                }
                let batch = RecordBatch::try_new(schema.clone(), rb).unwrap();

                writer.write(&batch).expect("Writing batch");
            }

            writer.close()?;

            // create new chrom file
            current_chrom = chrom;
            last_pos = 1;
            writer = get_writer(
                schema.clone(),
                format!("{}{}.parquet", ouput_fullname_prefix, current_chrom),
                batch_size,
            )?;
            fields_buffer = get_arrow_builders(&infos, &infos_keys, batch_size)?;
        }

        let pos = one_line[1].parse::<u32>()?;
        let a_ref = one_line[3].to_string();
        let a_alt = one_line[4].to_string();

        // check if gnomad file is well ordered
        if pos < last_pos {
            return Err(anyhow!("current pos<last_pos, donc non trié !"));
        }

        last_pos = pos;

        let h = parse_all_info(one_line.get(7).unwrap(), &infos, &infos_keys)?;

        fields_buffer[0]
            .as_any_mut()
            .downcast_mut::<Int32Builder>()
            .unwrap()
            .append_value(chrom as i32);
        fields_buffer[1]
            .as_any_mut()
            .downcast_mut::<Int32Builder>()
            .unwrap()
            .append_value(pos as i32);
        use xxhash_rust::xxh64::xxh64;
        let hash = xxh64(format!("{}-{}", a_ref, a_alt).as_bytes(), 0);
        fields_buffer[2]
            .as_any_mut()
            .downcast_mut::<Int64Builder>()
            .unwrap()
            .append_value(hash as i64);
        fields_buffer[3]
            .as_any_mut()
            .downcast_mut::<StringBuilder>()
            .unwrap()
            .append_value(a_ref);
        fields_buffer[4]
            .as_any_mut()
            .downcast_mut::<StringBuilder>()
            .unwrap()
            .append_value(a_alt);

        for (index, info_key) in infos_keys.iter().enumerate() {
            let b = &mut fields_buffer.get_mut(index + 4).unwrap();

            let info = infos.get(*info_key).unwrap();
            let value = h.get(*info_key);

            match info.value_type {
                ValueType::Float => {
                    let builder: &mut arrow::array::Float64Builder =
                        b.as_any_mut().downcast_mut::<Float64Builder>().unwrap();
                    match value {
                        Some(ValueParsed::Float(n)) => builder.append_value(*n),
                        None | Some(ValueParsed::Null) => builder.append_null(),
                        _ => return Err(anyhow!("value shoud be float for field {}", info_key)),
                    }
                }
                ValueType::Integer => {
                    let builder: &mut arrow::array::Int32Builder =
                        b.as_any_mut().downcast_mut::<Int32Builder>().unwrap();
                    match value {
                        Some(ValueParsed::Integer(n)) => builder.append_value(*n),
                        None | Some(ValueParsed::Null) => builder.append_null(),
                        _ => return Err(anyhow!("value shoud be integer for field {}", info_key)),
                    }
                }
                ValueType::String => {
                    let builder: &mut arrow::array::StringBuilder =
                        b.as_any_mut().downcast_mut::<StringBuilder>().unwrap();
                    match value {
                        Some(ValueParsed::String(n)) => builder.append_value(n),
                        None | Some(ValueParsed::Null) => builder.append_null(),
                        _ => return Err(anyhow!("value shoud be string for field {}", info_key)),
                    }
                }
                ValueType::Flag => {
                    let builder: &mut arrow::array::BooleanBuilder =
                        b.as_any_mut().downcast_mut::<BooleanBuilder>().unwrap();
                    match value {
                        Some(ValueParsed::Flag(n)) => builder.append_value(*n),
                        None | Some(ValueParsed::Null) => builder.append_null(),
                        _ => return Err(anyhow!("value shoud be flag for field {}", info_key)),
                    }
                }
            };
        }

        if fields_buffer[0].len() == batch_size {
            // write to parquet
            let mut rb = Vec::new();
            for mut f in fields_buffer.into_iter() {
                rb.push(f.finish() as ArrayRef);
            }
            let batch = RecordBatch::try_new(schema.clone(), rb).unwrap();

            writer.write(&batch).expect("Writing batch");

            fields_buffer = get_arrow_builders(&infos, &infos_keys, batch_size)?;
        }
        count += 1;
    }

    if !fields_buffer[0].is_empty() {
        // write to parquet
        let mut rb = Vec::new();
        for mut f in fields_buffer.into_iter() {
            rb.push(f.finish() as ArrayRef);
        }
        let batch = RecordBatch::try_new(schema.clone(), rb).unwrap();

        writer.write(&batch).expect("Writing batch");
    }

    writer.close()?;

    log::info!("nb variants converted {}", count);
    Ok(())
}

pub fn get_arrow_builders(
    infos: &HashMap<String, InfoFormatConfig>,
    infos_keys: &Vec<&String>,
    batch_size: usize,
) -> anyhow::Result<Vec<Box<dyn ArrayBuilder>>> {
    let mut fields_buffer = Vec::new();

    let contig_builder = Box::new(Int32Builder::with_capacity(batch_size)) as Box<dyn ArrayBuilder>;
    let pos_builder = Box::new(Int32Builder::with_capacity(batch_size)) as Box<dyn ArrayBuilder>;
    let ref_alt_hash_builder =
        Box::new(Int64Builder::with_capacity(batch_size)) as Box<dyn ArrayBuilder>;

    let a_ref_builder = Box::new(StringBuilder::new()) as Box<dyn ArrayBuilder>;
    let a_alt_builder = Box::new(StringBuilder::new()) as Box<dyn ArrayBuilder>;

    fields_buffer.push(contig_builder);
    fields_buffer.push(pos_builder);
    fields_buffer.push(ref_alt_hash_builder);
    fields_buffer.push(a_ref_builder);
    fields_buffer.push(a_alt_builder);

    for info_key in infos_keys.iter() {
        let info = infos.get(*info_key).unwrap();

        let builder = match info.value_type {
            ValueType::String => Box::new(StringBuilder::new()) as Box<dyn ArrayBuilder>,
            ValueType::Float => {
                Box::new(Float64Builder::with_capacity(batch_size)) as Box<dyn ArrayBuilder>
            }
            ValueType::Integer => {
                Box::new(Int32Builder::with_capacity(batch_size)) as Box<dyn ArrayBuilder>
            }
            ValueType::Flag => {
                Box::new(BooleanBuilder::with_capacity(batch_size)) as Box<dyn ArrayBuilder>
            }
        };

        fields_buffer.push(builder);
    }
    Ok(fields_buffer)
}
pub fn get_infos_keys(infos: &HashMap<String, InfoFormatConfig>) -> anyhow::Result<Vec<&String>> {
    let mut infos_keys = infos.keys().collect::<Vec<_>>();
    infos_keys.sort();

    // take only gnomad annotation with format : {annot}_{pop}_{sex}
    let mut special = Vec::new();
    let annots = ["AC", "AF", "AN", "nhomalt"];
    let pops = [
        "",
        "_afr",
        "_ami",
        "_amr",
        "_asj",
        "_eas",
        "_fin",
        "_nfe",
        "_remaining",
        "_sas",
        "_mid",
    ];
    let sexes = ["", "_XX", "_XY"];

    for annot in annots.iter() {
        for pop in pops.iter() {
            for sexe in sexes.iter() {
                special.push(format!("{}{}{}", annot, pop, sexe))
            }
        }
    }

    infos_keys.retain(|x| special.contains(x));

    Ok(infos_keys)
}
pub fn get_schema(
    infos: &HashMap<String, InfoFormatConfig>,
    infos_keys: &Vec<&String>,
) -> anyhow::Result<Arc<Schema>> {
    let mut fields = vec![
        Field::new("chrom", DataType::Int32, false),
        Field::new("pos", DataType::Int32, false),
        Field::new("ref_alt_hash", DataType::Int64, false),
        Field::new("a_ref", DataType::Utf8, false),
        Field::new("a_alt", DataType::Utf8, false),
    ];

    for key in infos_keys.iter() {
        let info_f = infos.get(*key).unwrap();
        match info_f.value_type {
            ValueType::Flag => fields.push(Field::new(*key, DataType::Boolean, true)),
            ValueType::Float => fields.push(Field::new(*key, DataType::Float64, true)),
            ValueType::String => fields.push(Field::new(*key, DataType::Utf8, true)),
            ValueType::Integer => fields.push(Field::new(*key, DataType::Int32, true)),
        }
    }

    let schema = Arc::new(Schema::new(fields));
    Ok(schema)
}
pub fn get_writer(
    schema: Arc<Schema>,
    output_file: String,
    batch_size: usize,
) -> anyhow::Result<ArrowWriter<File>> {
    // Default writer properties
    let props = WriterProperties::builder()
        //.set_write_batch_size(batch_size)
        .set_writer_version(WriterVersion::PARQUET_1_0)
        .set_statistics_enabled(EnabledStatistics::Page)
        .set_max_row_group_size(batch_size)
        .set_compression(Compression::SNAPPY)
        .build();

    let file = std::fs::File::create(output_file).unwrap();
    let writer = ArrowWriter::try_new(file, schema, Some(props)).unwrap();
    Ok(writer)
}
pub fn parse_all_info(
    info_fields: &str,
    infos: &HashMap<String, InfoFormatConfig>,
    infos_keys: &Vec<&String>,
) -> anyhow::Result<HashMap<String, ValueParsed>> {
    let mut h = HashMap::new();
    // maintenant nous parsons les lignes
    let k_vs = info_fields
        .split(";")
        .collect::<Vec<_>>()
        .iter()
        .map(|kv| {
            let ar = kv.split("=").collect::<Vec<_>>();
            match ar.as_slice() {
                [a, b] => Ok((a.to_string(), b.to_string())),
                [a] => Ok((a.to_string(), "true".to_string())), // une seul valeur donc Flag
                _ => Ok(("dummy".to_string(), kv.to_string())), // return Err(anyhow!("")) // cas vep=xxx|yyy=u| .. debilité !
            }
        })
        .collect::<anyhow::Result<HashMap<_, _>>>()?;
    // liste des [k,v]
    for k in infos_keys {
        //for k_v in k_vs {
        if let Some(v) = k_vs.get(*k) {
            h.insert(k.to_string(), convert(k, v, infos)?);
        } else {
            h.insert(k.to_string(), ValueParsed::Null);
        }
    }

    Ok(h)
}

pub fn convert(
    k: &str,
    v: &str,
    infos: &HashMap<String, InfoFormatConfig>,
) -> anyhow::Result<ValueParsed> {
    if v == "." {
        return Ok(ValueParsed::Null);
    }
    if let Some(info) = infos.get(k) {
        match info.value_type {
            ValueType::String => Ok(ValueParsed::String(v.to_string())),
            ValueType::Integer => Ok(ValueParsed::Integer(v.parse::<i32>()?)),
            ValueType::Float => Ok(ValueParsed::Float(v.parse::<f64>()?)),
            ValueType::Flag => Ok(ValueParsed::Flag(true)),
        }
    } else {
        Err(anyhow!("info header not found for key {}", k))
    }
}
pub fn parse_info(info: &str) -> anyhow::Result<InfoFormatConfig> {
    let parsed: Vec<&str> = info[8..].split(",").collect(); // ID=xx  Number=XX Type=xx
    let mut h = HashMap::new();

    for keyval in parsed {
        if let [k, v] = keyval.split("=").collect::<Vec<_>>().as_slice() {
            h.insert(k.to_string(), v.to_string());
        }
    }

    let value_type = match h.get("Type").unwrap().to_string().as_str() {
        "Flag" => ValueType::Flag,
        "String" => ValueType::String,
        "Integer" => ValueType::Integer,
        "Float" => ValueType::Float,
        x => return Err(anyhow!("info type unknown : {}", x)),
    };

    Ok(InfoFormatConfig {
        name: h.get("ID").unwrap().to_string(),
        value_number: h.get("Number").unwrap().to_string(),
        value_type,
    })
}

#[derive(Debug)]
pub struct InfoFormatConfig {
    pub name: String,
    pub value_type: ValueType,
    pub value_number: String,
}

#[derive(Debug)]
pub enum ValueParsed {
    String(String),
    Integer(i32),
    Float(f64),
    Flag(bool),
    Null,
}

#[derive(Debug)]
pub enum ValueType {
    String,
    Integer,
    Float,
    Flag,
}

pub fn compute_contig_num(chrom_str_p: &str) -> anyhow::Result<u8> {
    let chrom_str = if let Some(p) = chrom_str_p.strip_prefix("chr") {
        p
    } else {
        chrom_str_p
    };

    let c = chrom_str.parse::<u8>();
    match c {
        Ok(num) if (1..=22).contains(&num) => Ok(num),
        Ok(_) => Err(anyhow!("contig not reconize {}", chrom_str)),
        Err(_) => match chrom_str {
            "X" => Ok(23),
            "Y" => Ok(24),
            "MT" => Ok(25),
            "M" => Ok(25),
            _ => Err(anyhow!("contig not reconize {}", chrom_str)),
        },
    }
}
