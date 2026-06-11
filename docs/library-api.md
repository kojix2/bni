# Library API

`bni` can be used as a C library as well as a CLI.

Build with:

```sh
make
```

This produces:

```text
libbni.a
libbni.dylib    # macOS
libbni.so       # Linux
```

Install with:

```sh
make install PREFIX=$HOME/.local
```

This installs `bni`, `include/bni.h`, and the static/shared libraries.

## Link

Include `bni.h` and link with `-lbni` plus the htslib linker flags.

Example:

```sh
cc example.c -I$HOME/.local/include -L$HOME/.local/lib -lbni $(pkg-config --libs htslib)
```

## Build an Index

```c
#include <bni.h>

int rc = bni_build_index("reads.name.bam", NULL, NULL, NULL);
```

Pass a `bni_build_options_t` to control threads, overwrite behavior, and header checking:

```c
bni_build_options_t opts = {0};
opts.threads = 4;
opts.force = 1;

bni_build_stats_t stats = {0};
int rc = bni_build_index("reads.name.bam", "reads.name.bam.bni", &opts, &stats);
```

## Search an Index

Load a `.bni` directly and search the BGZF-block entry table:

```c
bni_index_t *idx = bni_index_open("reads.name.bam.bni");
const bni_entry_t *entry = bni_find_entry(idx, "READ_ID");
if (entry != NULL) {
    const char *first = bni_entry_first_name(idx, entry);
    const char *last  = bni_entry_last_name(idx, entry);
    /* candidate block: first <= possible READ_ID <= last */
    /* entry->beg_voff, entry->end_voff, entry->n_records */
}
bni_index_close(idx);
```

`bni_find_entry()` returns the earliest entry whose `last_qname >= target`.
The target may still be absent; the reader must scan forward and stop when
`QNAME > target`. Usually you should use `bni_reader_fetch()` rather than
interpreting entries directly.

Use `bni_load_index_file()` and `bni_index_destroy()` if you want caller-owned storage:

```c
bni_index_t idx;
if (bni_load_index_file("reads.name.bam.bni", &idx) == 0) {
    const bni_entry_t *entry = bni_find_entry(&idx, "READ_ID");
    bni_index_destroy(&idx);
}
```

## Fetch BAM Records

```c
#include <stdio.h>
#include <bni.h>

static int on_record(const bam1_t *record, const sam_hdr_t *header, void *user) {
    (void)header;
    (void)user;
    puts(bam_get_qname(record));
    return 0;
}

int main(void) {
    bni_reader_t *reader = NULL;
    if (bni_reader_open("reads.name.bam", NULL, NULL, &reader) != 0) return 1;

    uint32_t n_records = 0;
    int rc = bni_reader_fetch(reader, "READ_ID", on_record, NULL, &n_records);
    if (rc == 1) {
        /* not found */
    }

    bni_reader_close(reader);
    return rc < 0 ? 1 : 0;
}
```

`bni_reader_fetch()` returns:

- `0` when the name is found
- `1` when the name is absent
- `-1` on error

The callback receives a reused htslib `bam1_t`; copy it inside the callback if you need to keep it after the callback returns.

## Reader Options

```c
bni_reader_options_t opts = {0};
opts.threads = 4;
opts.ignore_metadata = 1;

bni_reader_t *reader = NULL;
int rc = bni_reader_open("reads.name.bam", "custom.bni", &opts, &reader);
```

By default, `bni_reader_open()` checks BAM size, mtime, and SAM header hash against the index metadata.
