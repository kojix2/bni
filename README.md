# bni

[![CI](https://github.com/kojix2/bni/actions/workflows/ci.yml/badge.svg)](https://github.com/kojix2/bni/actions/workflows/ci.yml)
[![Lines of Code](https://img.shields.io/endpoint?url=https%3A%2F%2Ftokei.kojix2.net%2Fbadge%2Fgithub%2Fkojix2%2Fbni%2Flines)](https://tokei.kojix2.net/github/kojix2/bni)

`bni` is a BAM Name Index for read-centric random access in BAM files sorted with:

```sh
samtools sort -N -o reads.name.bam reads.bam
```

It builds a sidecar index:

```text
reads.name.bam
reads.name.bam.bni
```

The index maps each SAM `QNAME` to one contiguous BGZF virtual-offset range:

```text
QNAME -> beg_voff, end_voff, n_records
```

Unlike coordinate indexes such as `.bai` and `.csi`, `bni` answers "read name -> BAM offsets". This is useful for long-read work where all alignments for one molecule/read need to be fetched without scanning the whole BAM.

## Status

Prototype C implementation, format version `BNIv1`.

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

This builds the CLI plus libraries:

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

## Quick Start

```sh
samtools sort -N -o reads.name.bam reads.bam
bni index reads.name.bam
bni get reads.name.bam 'READ_ID' > one_read.bam
bni stats reads.name.bam
bni check --full reads.name.bam
```

## Documentation

- [CLI usage](docs/cli.md)
- [Library API](docs/library-api.md)
- [BNIv1 file format](docs/format.md)

## Limitations

- Input is BAM only.
- Input must be BGZF-compressed and seekable.
- The intended sort order is `samtools sort -N` / `queryname:lexicographical`.
- CRAM input is not supported in BNIv1.
- The implementation currently loads the full `.bni` into memory. The binary layout is mmap-friendly, but mmap is not yet implemented.

## Acknowledgements

`bni` is inspired by Jared Simpson's [`bri`](https://github.com/jts/bri), but
targets queryname-sorted BAMs rather than coordinate-sorted BAMs.
