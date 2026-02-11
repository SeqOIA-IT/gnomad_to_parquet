### CLI tool : gnomad_to_parquet
convert a gnomad vcf file into parquet file format to be used with gLeaves dynamic annotation

## compilation
to be cross linux version compatible :
cargo build --release --target=x86_64-unknown-linux-musl

## usage
gnomad_to_parquet gnomad.genomes.v4.0.sites.chr1.vcf.bgz 1
=> creation du fichier gnomad_v4_chrom1.parquet 