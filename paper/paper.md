---
title: 'bni: A BAM name index for read-centric random access'
tags:
  - C
  - bioinformatics
  - genomics
  - BAM
  - long-read sequencing
authors:
  - name: Kojix2
    affiliation: 1
affiliations:
  - name: Independent researcher
    index: 1
date: 11 June 2026
bibliography: references.bib
---

# Summary

`bni` is a C command-line tool and library for read-centric random access in
queryname-sorted BAM files. It builds a binary sidecar index for BAM files
sorted with `samtools sort -N`, allowing records for a read to be found by a
local scan rather than by reading the whole file. The implementation uses the
BGZF virtual-offset model defined by the SAM/BAM specification [@samv1] and
htslib [@htslib].

The current `BNIv2` format is a BGZF-block name-range index. It stores one entry
for each BGZF block that contains the start of at least one BAM record. Each
entry records the block's first and last `QNAME`, virtual-offset range, and
record count. Lookup finds the earliest candidate block by its last `QNAME`,
seeks to the first record boundary, and scans forward until the requested name
is found and passed.

The package provides a CLI for index construction, read extraction, metadata
inspection, and index validation, plus a C API for embedding the same index and
fetch operations in other software.

# Statement of need

BAM files are commonly indexed by genomic coordinate using formats such as BAI
and CSI in the SAM/BAM ecosystem [@samv1; @samtools2009]. These indexes are well
suited to queries of the form "which alignments overlap this reference
interval?" but they do not answer the complementary question "where are all
records for this read name?" In long-read sequencing workflows, read-centric
inspection is often needed to investigate a molecule across primary,
supplementary, secondary, and unmapped alignments. Existing tools can filter by
query name, but without a read-name index this usually requires a full file
scan.

Earlier software such as `bri` demonstrated the utility of a BAM read-name
index for extracting alignments by read name in long-read sequencing data
[@bri], and `bripy` exposed related functionality through a Python interface
[@bripy]. `bni` focuses on the queryname-sorted case. Instead of storing one
index entry per distinct `QNAME`, it stores the name range covered by each BGZF
block. This keeps the index tied to compressed-block structure while using
lexicographic `QNAME` order to locate candidate records.

# Design and implementation

The `BNIv2` file format is an uncompressed binary sidecar consisting of a
128-byte header, a fixed-width entry table, and a NUL-terminated string table.
The header records the format version, table offsets, entry and record counts,
BAM metadata, a 64-bit FNV-1a hash of the SAM header, sort order, and entry
size. Each 40-byte entry stores offsets for the first and last boundary names,
the beginning and ending BGZF virtual offsets, and the number of assigned
records.

Index construction streams through a BGZF-compressed BAM file with htslib,
records the BGZF virtual offset before and after each BAM record, validates
lexicographic query-name order, and starts a new index entry when the compressed
BGZF block address changes. Because a BGZF block can begin in the middle of a
BAM record, `bni` stores the first BAM record boundary inside the block rather
than assuming an uncompressed offset of zero.

Lookup loads the `.bni` file into memory and performs binary search over entry
last names. Fetching records seeks to the candidate entry's beginning virtual
offset and reads sequentially with htslib until the requested `QNAME` is found
and then passed. This scan may continue beyond the first candidate entry,
because the same read name can span adjacent BGZF block entries. Validation
checks BAM metadata and can also verify every indexed block range.

# Functionality

The command-line interface includes:

- `bni index` to build a `.bam.bni` sidecar index.
- `bni get` to extract one or many read names as BAM, SAM, or CRAM output.
- `bni stats` to inspect index metadata.
- `bni check` to verify metadata and, optionally, every indexed range.

The library interface exposes index construction, index loading, direct entry
lookup, and a reusable reader object for fetching records by `QNAME`.

# Availability

The source code is released under the MIT license and is available at
<https://github.com/kojix2/bni>. The package builds with a C compiler and
htslib, and includes smoke tests covering CLI extraction and the library API.

# Acknowledgements

`bni` builds on the BAM, SAM, and BGZF ecosystem implemented by htslib and
SAMtools.

# AI usage disclosure

The initial draft of this paper was prepared with assistance from OpenAI Codex.
The author reviewed and edited the generated text and remains responsible for
the manuscript content.

# References
