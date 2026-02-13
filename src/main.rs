use anyhow::anyhow;
use anyhow::Context;
use arrow::array::*;
use arrow::datatypes::*;
use arrow::record_batch::RecordBatch;
use clap::Parser;
use flate2::bufread::MultiGzDecoder;
use parquet::arrow::arrow_writer::ArrowWriter;
use parquet::basic::Compression;
use parquet::file::properties::WriterVersion;
use parquet::file::properties::{EnabledStatistics, WriterProperties};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::sync::Arc;
use tikv_jemallocator::Jemalloc;
use xxhash_rust::xxh64::xxh64;

/// always better performance with an alternative memory allocator than default
#[global_allocator]
static GLOBAL: Jemalloc = Jemalloc;

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

/// CLI program to convert a gnomad vcf file to a gLeaves compatible parquet annotation file
/// ```
/// gnomad_to_parquet --prefix-output-filename my_gnomad_name_ --output-dir /tmp --input-fullname /data/gnomad.bgz --only-chrom 1
/// ```
fn main() -> anyhow::Result<()> {
    env_logger::init();

    let cli = Cli::parse();

    let ouput_fullname_prefix = format!("{}/{}", cli.output_dir, cli.prefix_output_filename);

    process(&cli.input_fullname, &ouput_fullname_prefix, cli.only_chrom)
}

/// function to convert a gnomad vcf file to a gLeaves compatible parquet annotation file
/// * `input_fullname` - full path name to the gnomad vcf file to process ex : /xxx/gnomad.vcf.gz
/// * `ouput_fullname_prefix` - prefix of the final output file name generated : ex /yyy/gnomad_gleaves_ (will create file /yyy/gnomad_gleaves_CHROM.parquet)
/// * `only_chrom` - if given, generate parquet file only for the specified chrom (in case of one gnomad file containing all chrom)
fn process(
    input_fullname: &str,
    ouput_fullname_prefix: &str,
    only_chrom: Option<String>,
) -> anyhow::Result<()> {
    // check if file is compressed
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
        // read an INFO line
        if l2.starts_with("##INFO") {
            // reading the annotation format
            let info_f = parse_info(&l2)?;
            infos.insert(info_f.name.to_string(), info_f);
        }

        // stop reading the header
        if l2.starts_with("#CHROM") {
            break;
        }
    }

    if infos.is_empty() {
        return Err(anyhow::anyhow!("no info lines found"));
    }

    // keep in infos only desire annotations
    filter_infos(&mut infos);

    // sort keys (order is important to parquet writing)
    let mut infos_keys = infos.keys().collect::<Vec<_>>();
    infos_keys.sort();

    let batch_size = 100_000;

    // construct parquet schema from infos
    let schema = get_schema(&infos, &infos_keys)?;

    let mut current_chrom = if let Some(c) = &only_chrom {
        compute_contig_num(c)?
    } else {
        0 // mean : "read the first chrom from the first line of the given gnomad file"
    };

    // last_pos in the chrom (to check if ordering in the gnomad file is ok)
    let mut last_pos = 1;

    // current chrom parquet writer
    let mut writer: Option<ArrowWriter<File>> = None;

    // construct arrow builders
    let mut fields_buffer = get_arrow_builders(&infos, &infos_keys, batch_size)?;

    // process each "data" line of the gnomad input
    while let Some(Ok(l2)) = it.next() {
        let one_line = l2.split("\t").collect::<Vec<_>>();

        let chrom = compute_contig_num(one_line[0])?;

        // current_chrom read from the first gnomad data line
        if current_chrom == 0 {
            current_chrom = chrom;
        }

        // "only_chrom" cli param is given, but chrom in current line doesn't match, bypass the line
        if only_chrom.is_some() && chrom != current_chrom {
            last_pos = 1; // reset last_pos
            continue;
        }

        if only_chrom.is_none() && chrom != current_chrom {
            // new chrom, write last data
            if !fields_buffer[0].is_empty() {
                // convert arrow buffer (fields_buffer) to parquet recordbatch
                let rb = fields_buffer
                    .into_iter()
                    .map(|mut buf| buf.finish() as ArrayRef)
                    .collect::<Vec<_>>();
                let batch = RecordBatch::try_new(schema.clone(), rb)?;

                // no writer currently define, create it
                if writer.is_none() {
                    writer = Some(get_writer(
                        schema.clone(),
                        format!("{}{}.parquet", ouput_fullname_prefix, current_chrom),
                        batch_size,
                    )?);
                }
                writer
                    .as_mut()
                    .unwrap()
                    .write(&batch)
                    .expect("Writing batch");
            }

            writer.take().unwrap().close()?;

            // chrom have changed
            current_chrom = chrom;
            last_pos = 1;
            fields_buffer = get_arrow_builders(&infos, &infos_keys, batch_size)?;
        }

        let pos = one_line[1].parse::<u32>()?;
        let a_ref = one_line[3].to_string();
        let a_alt = one_line[4].to_string();

        // check if gnomad file is well ordered
        if pos < last_pos {
            return Err(anyhow!(
                "current pos<last_pos, gnomad not sorted by pos abort !"
            ));
        }

        last_pos = pos;

        // annotation position in the VCF line is 7, parse it, and get only annotations presents in the infos variable
        let h = parse_all_info(
            one_line.get(7).context("one_line problem")?,
            &infos,
            &infos_keys,
        )?;

        // now fill the arrow buffer with data
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

        // fill arrow buffer with annotations (ordering by the infos_keys)
        for (index, info_key) in infos_keys.iter().enumerate() {
            let b = &mut fields_buffer.get_mut(index + 5).unwrap();

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

        // write batch of data if buffer reach batch_size
        if fields_buffer[0].len() == batch_size {
            // write to parquet

            // convert arrow buffer (fields_buffer) to parquet recordbatch
            let rb = fields_buffer
                .into_iter()
                .map(|mut buf| buf.finish() as ArrayRef)
                .collect::<Vec<_>>();
            let batch = RecordBatch::try_new(schema.clone(), rb)?;

            if writer.is_none() {
                writer = Some(get_writer(
                    schema.clone(),
                    format!("{}{}.parquet", ouput_fullname_prefix, current_chrom),
                    batch_size,
                )?);
            }
            writer
                .as_mut()
                .unwrap()
                .write(&batch)
                .expect("Writing batch");

            // create new empty arrow buffer
            fields_buffer = get_arrow_builders(&infos, &infos_keys, batch_size)?;
        }
        count += 1;
    }
    // end of file reached

    // check if we still have data in buffer
    if !fields_buffer[0].is_empty() {
        // write to parquet

        // convert arrow buffer (fields_buffer) to parquet recordbatch
        let rb = fields_buffer
            .into_iter()
            .map(|mut buf| buf.finish() as ArrayRef)
            .collect::<Vec<_>>();
        let batch = RecordBatch::try_new(schema.clone(), rb)?;

        if writer.is_none() {
            writer = Some(get_writer(
                schema.clone(),
                format!("{}{}.parquet", ouput_fullname_prefix, current_chrom),
                batch_size,
            )?);
        }
        writer
            .as_mut()
            .unwrap()
            .write(&batch)
            .expect("Writing batch");
    }

    writer.take().unwrap().close()?;

    log::info!("nb gnomad variants converted {}", count);
    Ok(())
}

/// return arrow array builders from the infos annotations format
/// * `infos` - info annotations format
/// * `infos_keys` - desired order of the keys
/// * `batch_size` - used to preallocate array
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

/// remove from the list of info lines (=annots) all non population reference
/// * `infos` - info annotations format
pub fn filter_infos(infos: &mut HashMap<String, InfoFormatConfig>) {
    // take only gnomad annotation with format : {annot}_{pop}_{sex}
    let mut elected_annots = Vec::new();
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
                elected_annots.push(format!("{}{}{}", annot, pop, sexe))
            }
        }
    }

    infos.retain(|k, _v| elected_annots.contains(k));
}

/// return parquet schema from the infos annotations format
/// * `infos` - info annotations format
/// * `infos_keys` - desired order of the keys
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

/// return a parquet writer corresponding to the infos annotations format
/// * `schema` - schema of the parquet file return by the get_schema function
/// * `output_file` - full path name of the target parquet file
/// * `batch_size` - the parquet rowgroup size
pub fn get_writer(
    schema: Arc<Schema>,
    output_file: String,
    batch_size: usize,
) -> anyhow::Result<ArrowWriter<File>> {
    // Default writer properties
    let props = WriterProperties::builder()
        .set_writer_version(WriterVersion::PARQUET_1_0)
        .set_statistics_enabled(EnabledStatistics::Page)
        .set_max_row_group_size(batch_size)
        .set_compression(Compression::SNAPPY)
        .build();

    let file = std::fs::File::create(output_file)?;
    let writer = ArrowWriter::try_new(file, schema, Some(props))?;
    Ok(writer)
}

/// parse all field in the INFO part of a VCF line
/// * `info_fields` - the INFO part
/// * `infos` - info annotations format
/// * `infos_keys` - desired order of the keys
pub fn parse_all_info(
    info_fields: &str,
    infos: &HashMap<String, InfoFormatConfig>,
    infos_keys: &Vec<&String>,
) -> anyhow::Result<HashMap<String, ValueParsed>> {
    let mut h = HashMap::new();

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

    // keep only annot from the infos_keys
    for k in infos_keys {
        if let Some(v) = k_vs.get(*k) {
            h.insert(k.to_string(), convert(k, v, infos)?);
        } else {
            h.insert(k.to_string(), ValueParsed::Null);
        }
    }

    Ok(h)
}

/// convert a string annotation into a ValueParsed using the corresponding InfoFormatConfig
/// * `k` - annotation name (= index in the infos datastructure)
/// * `v` - the value to parse
/// * `infos` - annot -> InfoFormatConfig mapping
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

/// parse a vcf header INFO line
/// * `info` - the header info
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

/// store a vcf INFO line specification
#[derive(Debug)]
pub struct InfoFormatConfig {
    /// name of the annotation
    pub name: String,

    /// VCF type format
    pub value_type: ValueType,

    /// VCF number format
    pub value_number: String,
}

/// VCF field annotation enum
#[derive(Debug)]
pub enum ValueParsed {
    String(String),
    Integer(i32),
    Float(f64),
    Flag(bool),
    Null,
}

/// VCF annotation datatype
#[derive(Debug)]
pub enum ValueType {
    String,
    Integer,
    Float,
    Flag,
}

/// convert a chrom in string format to numeric (ex : "1" -> 1, ... "X"=> 23, "Y"=> 24, "MT"=>25)
/// * `chrom_str_p` - the chrom in string format
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

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_input_file_not_found() {
        let result = process("dummy_file_name", "/tmp", None);
        assert_eq!(result.is_err(), true)
    }

    #[test]
    fn test_output_dir_notfound() {
        let result = process(
            "./test_data/gnomad-4.1-chr1-1000.vcf.bgz",
            "/tmp/xx/toto",
            None,
        );
        assert_eq!(result.is_err(), true)
    }

    #[test]
    fn test_no_info() {
        let result = process(
            "./test_data/gnomad-4.1-chr1-without-info.vcf.bgz",
            "/tmp/toto_",
            None,
        );
        assert_eq!(result.is_err(), true)
    }

    #[test]
    fn test_not_sorted() {
        let result = process(
            "./test_data/gnomad-4.1-chr1-1000-not-sorted.vcf.gz",
            "/tmp/toto_",
            None,
        );
        assert_eq!(result.is_err(), true)
    }

    #[test]
    fn test_convert_monochrom() {
        let prefix = "/tmp/test_convert_monochrom_gnomad_";
        let p1_f = format!("{prefix}1.parquet");
        let p2_f = format!("{prefix}2.parquet");
        let p1 = std::path::Path::new(&p1_f);
        let p2 = std::path::Path::new(&p2_f);

        if p1.exists() {
            std::fs::remove_file(p1)
                .context("can't delete tmp file")
                .unwrap();
        }
        if p2.exists() {
            std::fs::remove_file(p2)
                .context("can't delete tmp file")
                .unwrap();
        }

        let result1 = process("./test_data/gnomad-4.1-chr1-1000.vcf.bgz", prefix, None);
        let result2 = process("./test_data/gnomad-4.1-chr2-1000.vcf.bgz", prefix, None);

        let p1_exists = p1.exists();
        let p2_exists = p2.exists();

        if p1.exists() {
            std::fs::remove_file(p1)
                .context("can't delete tmp file")
                .unwrap();
        }
        if p2.exists() {
            std::fs::remove_file(p2)
                .context("can't delete tmp file")
                .unwrap();
        }

        assert_eq!(result1.is_ok(), true);
        assert_eq!(p1_exists, true);
        assert_eq!(result2.is_ok(), true);
        assert_eq!(p2_exists, true);
    }

    #[test]
    fn test_convert_multichrom_only_one_chrom() {
        let prefix = "/tmp/test_convert_multichrom_one_chrom_gnomad_";
        let p1_f = format!("{prefix}1.parquet");
        let p2_f = format!("{prefix}2.parquet");
        let p1 = std::path::Path::new(&p1_f);
        let p2 = std::path::Path::new(&p2_f);

        if p1.exists() {
            std::fs::remove_file(p1)
                .context("can't delete tmp file")
                .unwrap();
        }
        if p2.exists() {
            std::fs::remove_file(p2)
                .context("can't delete tmp file")
                .unwrap();
        }

        let result2 = process(
            "./test_data/gnomad-4.1-multi_chrom.vcf.gz",
            prefix,
            Some("2".to_string()),
        );

        let p1_exists = p1.exists();
        let p2_exists = p2.exists();

        if p1.exists() {
            std::fs::remove_file(p1)
                .context("can't delete tmp file")
                .unwrap();
        }
        if p2.exists() {
            std::fs::remove_file(p2)
                .context("can't delete tmp file")
                .unwrap();
        }

        assert_eq!(result2.is_ok(), true);
        assert_eq!(p1_exists, false);
        assert_eq!(p2_exists, true);
    }

    #[test]
    fn test_convert_multichrom_all_chrom() {
        let prefix = "/tmp/test_convert_multichrom_all_chrom_gnomad_";
        let p1_f = format!("{prefix}1.parquet");
        let p2_f = format!("{prefix}2.parquet");
        let p1 = std::path::Path::new(&p1_f);
        let p2 = std::path::Path::new(&p2_f);

        if p1.exists() {
            std::fs::remove_file(p1)
                .context("can't delete tmp file")
                .unwrap();
        }
        if p2.exists() {
            std::fs::remove_file(p2)
                .context("can't delete tmp file")
                .unwrap();
        }

        let result = process("./test_data/gnomad-4.1-multi_chrom.vcf.gz", prefix, None);

        let p1_exists = p1.exists();
        let p2_exists = p2.exists();

        if p1.exists() {
            std::fs::remove_file(p1)
                .context("can't delete tmp file")
                .unwrap();
        }
        if p2.exists() {
            std::fs::remove_file(p2)
                .context("can't delete tmp file")
                .unwrap();
        }

        assert_eq!(result.is_ok(), true);
        assert_eq!(p1_exists, true);
        assert_eq!(p2_exists, true);
    }
}
