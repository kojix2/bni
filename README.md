# bni

[![CI](https://github.com/kojix2/bni/actions/workflows/ci.yml/badge.svg)](https://github.com/kojix2/bni/actions/workflows/ci.yml)
[![Lines of Code](https://img.shields.io/endpoint?url=https%3A%2F%2Ftokei.kojix2.net%2Fbadge%2Fgithub%2Fkojix2%2Fbni%2Flines)](https://tokei.kojix2.net/github/kojix2/bni)

`bni` builds and queries a read-name index for BAM files sorted with
`samtools sort -N`.

## Status

- Prototype C implementation
- File format: [BNIv2](docs/format.md)
- Input: BGZF-compressed BAM
- Required sort order: `queryname:lexicographical`

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

Build outputs:

```text
bni
libbni.a
libbni.dylib    # macOS
libbni.so       # Linux
```

Install:

```sh
make install PREFIX=$HOME/.local
```

## Usage

```sh
samtools sort -N -o reads.name.bam reads.bam
bni index reads.name.bam
bni get reads.name.bam 'READ_ID' > one_read.bam
```

Additional checks and metadata:

```sh
bni stats reads.name.bam
bni check --full reads.name.bam
```

Default index path:

```text
reads.name.bam.bni
```

## Format Summary

BNIv2 stores one entry per BGZF block that contains BAM record starts.

Each entry stores:

```text
first_qname
last_qname
beg_voff
end_voff
n_records
```

Lookup uses the first entry whose `last_qname >= target_qname`, seeks to
`beg_voff`, and scans forward until the target QNAME is found and passed.

## Documentation

- [CLI usage](docs/cli.md)
- [Library API](docs/library-api.md)
- [BNIv2 file format](docs/format.md)

## Limitations

- Input is BAM only.
- Input must be BGZF-compressed and seekable.
- Input must be sorted with `samtools sort -N`.
- CRAM input is not supported in BNIv2.
- Index entries are loaded into memory; the string table is mmap-backed when available.

## Acknowledgements

`bni` is inspired by Jared Simpson's [`bri`](https://github.com/jts/bri).
