# CLI Usage

## Create an Index

Create a queryname-sorted BAM:

```sh
samtools sort -N -o reads.name.bam reads.bam
```

Build the index:

```sh
bni index reads.name.bam
```

This writes:

```text
reads.name.bam.bni
```

BNIv2 stores one entry per BGZF block that contains BAM record starts. Each
entry remembers the first and last QNAME in that block-like range.

## Extract Records

Extract one read name as BAM:

```sh
bni get reads.name.bam 'READ_ID' > one_read.bam
```

Extract as SAM:

```sh
bni get -O sam reads.name.bam 'READ_ID' > one_read.sam
```

Extract many names:

```sh
bni get -f names.txt -o subset.bam reads.name.bam
```

For name-file queries, `bni get` batches requested names by indexed BGZF block
to reduce random seeks. Output records are emitted in BAM/index order, and
duplicate requested names are fetched once.

Show index metadata:

```sh
bni stats reads.name.bam
```

Check BAM/index consistency:

```sh
bni check reads.name.bam
bni check --full reads.name.bam
```

## `bni index`

```text
Usage:
  bni index [options] <in.name.bam>

Options:
  -o, --output FILE          output index file [default: <in.bam>.bni]
  -f, --force                overwrite an existing index
  -@, --threads INT          decompression threads
      --no-header-check      do not require @HD SO/SS tags
  -v, --verbose              print progress summary
  -h, --help                 show help
```

By default, `bni index` expects:

```text
@HD SO:queryname SS:queryname:lexicographical
```

If `SS` is absent, `bni` still validates actual lexicographic `QNAME` order
while indexing. If the BAM is naturally queryname-sorted with `samtools sort -n`,
indexing will normally fail when a natural-order pair disagrees with `strcmp`
order. This is deliberate: `bni` is a `samtools sort -N`-oriented index.

## `bni get`

```text
Usage:
  bni get [options] <in.name.bam> <read-name>
  bni get [options] -f names.txt <in.name.bam>

Options:
  -i, --index FILE           index file [default: <in.bam>.bni]
  -o, --output FILE          output file [default: stdout]
  -O, --output-format FMT    bam|sam|cram [default: bam, or inferred from -o]
  -f, --name-file FILE       read names, one per line
  -@, --threads INT          input/output threads
      --no-header            do not write SAM header; valid only with -O sam
      --with-header          write header [default]
      --missing-ok           exit 0 even if some names are absent
      --list-missing         print missing names to stderr
      --ignore-metadata      do not require BAM size/mtime/header hash to match index
  -h, --help                 show help
```

## `bni stats`

```text
Usage:
  bni stats [options] <in.name.bam>

Options:
  -i, --index FILE           index file [default: <in.bam>.bni]
  -h, --help                 show help
```

## `bni check`

```text
Usage:
  bni check [options] <in.name.bam>

Options:
  -i, --index FILE           index file [default: <in.bam>.bni]
      --quick                check metadata only [default]
      --full                 seek every indexed BGZF-block range and verify boundaries/counts
  -@, --threads INT          decompression threads
  -h, --help                 show help
```
