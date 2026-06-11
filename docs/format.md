# BNIv1 File Format

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

## Semantics

The index maps each read name / SAM `QNAME` to one contiguous BGZF virtual-offset range:

```text
QNAME -> beg_voff, end_voff, n_records
```

This requires all records for the same `QNAME` to be contiguous in the BAM. The intended sort order is `samtools sort -N` / `queryname:lexicographical`.
