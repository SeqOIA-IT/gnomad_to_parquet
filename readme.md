### CLI tool : gnomad_to_parquet

Convert a gnomad vcf file into parquet file format to be used with gLeaves dynamic annotation

## compilation
to be cross linux version compatible :
```
cargo build --release --target=x86_64-unknown-linux-musl
```

## usage

Create a file gnomad_v4_chrom1.parquet
```
gnomad_to_parquet --prefix-output-filename my_gnomad_name_ --output-dir /tmp --input-fullname /data/gnomad.bgz --only-chrom 1
```

## documentation

Generate documentation:
```
cargo doc --no-deps
```
