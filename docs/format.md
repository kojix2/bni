# BNIv2 File Format

`.bni` is an uncompressed binary sidecar file for BAM files sorted with
`samtools sort -N`. BAM remains BGZF-compressed; BNI itself is not compressed.

BNIv2 is a **BGZF-block name-range index**. It does not store one entry per
read name. Instead, it stores one entry for each BGZF block that contains the
start of at least one BAM record.

```text
BNI file
├── 128-byte header
├── fixed-width BGZF-block entry table
└── NUL-terminated string table
```

Header fields are little-endian:

```text
magic                  4 bytes  "BNI\1"
version                u32      2
header_size            u32      128
flags                  u32      BNI_FLAG_BGZF_BLOCKS
n_blocks               u64      number of BGZF-block entries
n_records              u64      total BAM records indexed
entries_offset         u64      byte offset of entry table
strings_offset         u64      byte offset of string table
strings_size           u64      byte size of string table
bam_size               u64      source BAM file size
bam_mtime              i64      source BAM mtime
header_hash            u64      FNV-1a 64-bit over SAM header text
sort_order             u32      1 = queryname:lexicographical
entry_size             u32      40
reserved               zero-filled to 128 bytes
```

Each entry is 40 bytes, little-endian:

```text
first_name_offset      u64      offset of first QNAME in string table
last_name_offset       u64      offset of last QNAME in string table
beg_voff               u64      virtual offset of first record start in this BGZF block
end_voff               u64      virtual offset just after this entry's last assigned record
n_records              u32      records whose starts are assigned to this entry
reserved               u32
```

The string table contains only boundary names:

```text
first0\0last0\0first1\0last1\0...
```

## Semantics

For each BGZF block that contains BAM record starts, BNI records:

```text
BGZF-record-start block -> first_qname, last_qname, beg_voff, end_voff
```

`beg_voff` is a BAM record boundary. It is usually inside the corresponding
BGZF block, not necessarily at uncompressed offset zero. This matters because a
BGZF block can begin in the middle of a previous BAM record. Seeking to
`beg_voff` is therefore safer than blindly seeking to the raw BGZF block start.

Entries are ordered as they occur in the queryname-sorted BAM. Because the BAM
is lexicographically queryname-sorted, `first_qname` and `last_qname` are
non-decreasing across entries. A QNAME may span adjacent entries, so equality is
allowed.

Lookup finds the earliest entry whose `last_qname >= target_qname`, seeks to
that entry's `beg_voff`, and scans forward with `sam_read1()` until the target
QNAME is found and then passed.

```text
find first entry where last_qname >= target
seek to entry.beg_voff
while sam_read1():
    if qname < target: continue
    if qname == target: emit
    if qname > target: stop
```

This is intentionally different from record-offset indexes for arbitrary BAM
order. BNIv2 uses the queryname sort order and BGZF block structure to keep the
index small while still avoiding a full BAM scan.
