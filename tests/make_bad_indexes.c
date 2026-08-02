#include "bni.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int write_index(const char *path, const bni_file_header_t *header,
                       const bni_entry_t *entries, const char *strings) {
  if (bni_write_index_file(path, header, entries, strings) != 0) {
    fprintf(stderr, "failed writing %s\n", path);
    return -1;
  }
  return 0;
}

int main(int argc, char **argv) {
  if (argc != 7) {
    fprintf(stderr, "usage: %s IN EMPTY EMPTY_COUNT COUNT TRUNCATED GAP\n", argv[0]);
    return 2;
  }

  bni_index_t idx;
  if (bni_load_index_file(argv[1], &idx) != 0) {
    return 1;
  }
  if (idx.header.n_blocks < 2 || idx.header.n_records == UINT64_MAX) {
    fprintf(stderr, "test index needs at least two entries and a non-maximal record count\n");
    bni_index_destroy(&idx);
    return 1;
  }

  int status = 1;
  bni_file_header_t header = idx.header;
  header.n_blocks = 0;
  header.n_records = 0;
  header.strings_offset = header.entries_offset;
  header.strings_size = 0;
  if (write_index(argv[2], &header, NULL, "") != 0) {
    goto done;
  }

  header.n_records = 123;
  if (write_index(argv[3], &header, NULL, "") != 0) {
    goto done;
  }

  header = idx.header;
  header.n_records++;
  if (write_index(argv[4], &header, idx.entries, idx.strings) != 0) {
    goto done;
  }

  header = idx.header;
  header.n_blocks--;
  header.n_records -= idx.entries[idx.header.n_blocks - 1].n_records;
  header.strings_offset = header.entries_offset + (header.n_blocks * header.entry_size);
  if (write_index(argv[5], &header, idx.entries, idx.strings) != 0) {
    goto done;
  }

  bni_entry_t *entries = malloc((size_t)idx.header.n_blocks * sizeof(*entries));
  if (entries == NULL) {
    goto done;
  }
  memcpy(entries, idx.entries, (size_t)idx.header.n_blocks * sizeof(*entries));
  if (entries[1].beg_voff == UINT64_MAX || entries[1].beg_voff + 1 >= entries[1].end_voff) {
    fprintf(stderr, "second test entry is too short to introduce a gap\n");
    free(entries);
    goto done;
  }
  entries[1].beg_voff++;
  if (write_index(argv[6], &idx.header, entries, idx.strings) == 0) {
    status = 0;
  }
  free(entries);

done:
  bni_index_destroy(&idx);
  return status;
}
