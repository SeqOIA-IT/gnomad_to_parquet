### CLI tool : gnomad_to_parquet
Convert a GnomAD VCF file into a Parquet file format to be used with gLeaves dynamic annotation

## compilation
To be cross linux version compatible :
``cargo build --release --target=x86_64-unknown-linux-musl``

## usage
``gnomad_to_parquet --prefix-output-filename my_gnomad_name_ --output-dir /tmp --input-fullname /data/gnomad.bgz --only-chrom 1``

=> create the file my_gnomad_name_1.parquet into /tmp

## documentation
``cargo doc --no-deps``
