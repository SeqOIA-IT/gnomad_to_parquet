### CLI tool : gnomad_to_parquet
convert a gnomad vcf file into parquet file format to be used with gLeaves dynamic annotation

## compilation
to be cross linux version compatible :
cargo build --release --target=x86_64-unknown-linux-musl

## usage
gnomad_to_parquet --prefix-output-filename my_gnomad_name_ --output-dir /tmp --input-fullname /data/gnomad.bgz --only-chrom 1
=> creation du fichier gnomad_v4_chrom1.parquet 

## documentation
cargo doc --no-deps