# bni

`bni` is a BAM Name Index for read-centric random access in BAM files sorted with:

```sh
samtools sort -N -o reads.name.bam reads.bam
```

It builds a sidecar index:

```text
reads.name.bam
reads.name.bam.bni
```

The index maps each read name / SAM `QNAME` to one contiguous BGZF virtual-offset range:

```text
QNAME -> beg_voff, end_voff, n_records
```

This is intentionally different from a record-offset index for arbitrary BAM order. `bni` assumes that all records with the same `QNAME` are contiguous in the BAM. That makes extracting a read's primary, supplementary, secondary, and unmapped records a single BGZF seek followed by a short sequential read.

## Status

Prototype C implementation, format version `BNIv1`.

The code is written from scratch. Existing BAM read-index projects are useful context, but this repository does not copy their code.

## Requirements

- C compiler
- htslib development headers and library
- zlib / bzip2 / lzma as required by your htslib build

## Build

With a system htslib available via `pkg-config`:

```sh
make
```

With a local htslib source/build tree:

```sh
make HTSLIB_DIR=/path/to/htslib
```

Install:

```sh
make install PREFIX=$HOME/.local
```

## Usage

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

Show index metadata:

```sh
bni stats reads.name.bam
```

Check BAM/index consistency:

```sh
bni check reads.name.bam
bni check --full reads.name.bam
```

## Commands

### `bni index`

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

If `SS` is absent, `bni` still validates actual lexicographic `QNAME` order while indexing. If the BAM is naturally queryname-sorted with `samtools sort -n`, indexing will normally fail when a natural-order pair disagrees with `strcmp` order. This is deliberate: `bni` is a `samtools sort -N`-oriented index.

### `bni get`

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

### `bni stats`

```text
Usage:
  bni stats [options] <in.name.bam>

Options:
  -i, --index FILE           index file [default: <in.bam>.bni]
  -h, --help                 show help
```

### `bni check`

```text
Usage:
  bni check [options] <in.name.bam>

Options:
  -i, --index FILE           index file [default: <in.bam>.bni]
      --quick                check metadata only [default]
      --full                 also seek every indexed range and verify QNAME/count
  -@, --threads INT          decompression threads
  -h, --help                 show help
```

## BNIv1 file format

`.bni` is an uncompressed binary sidecar file. BAM remains BGZF-compressed; BNI itself is not compressed.

```text
BNI file
├── 128-byte header
├── fixed-width entry table
└── NUL-terminated string table
```

Header fields are little-endian:

```text
magic                  4 bytes  "BNI\1"
version                u32      1
header_size            u32      128
flags                  u32      BNI_FLAG_QUERYNAME_GROUPED
n_names                u64
n_records              u64
entries_offset         u64
strings_offset         u64
strings_size           u64
bam_size               u64
bam_mtime              i64
header_hash            u64      FNV-1a 64-bit over SAM header text
sort_order             u32      1 = queryname:lexicographical
entry_size             u32      32
reserved               zero-filled to 128 bytes
```

Each entry is 32 bytes, little-endian:

```text
name_offset            u64      offset into string table
beg_voff               u64      BGZF virtual offset of first record
end_voff               u64      BGZF virtual offset just after last record
n_records              u32      records with this QNAME
reserved               u32
```

Entries are sorted lexicographically by QNAME, and lookup is a binary search over the entry table.

## Limitations

- Input is BAM only.
- Input must be BGZF-compressed and seekable.
- The intended sort order is `samtools sort -N` / `queryname:lexicographical`.
- CRAM input is not supported in BNIv1.
- Output CRAM mode is passed through htslib and may require htslib-specific reference configuration.
- The implementation currently loads the full `.bni` into memory. The binary layout is mmap-friendly, but mmap is not yet implemented.

## Why this exists

Coordinate indexes such as `.bai` and `.csi` answer:

```text
reference interval -> BAM offsets
```

`bni` answers:

```text
read name -> BAM BGZF offset range
```

For long-read work, this is useful when you want to inspect all alignments belonging to one molecule/read without scanning the whole BAM.
