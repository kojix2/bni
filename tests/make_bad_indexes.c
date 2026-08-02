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
  if (argc != 13) {
    fprintf(stderr,
            "usage: %s IN EMPTY EMPTY_COUNT COUNT TRUNCATED GAP ENTRY_ORDER RANGE_ORDER "
            "INSIDE_OFFSET UNTERMINATED UNKNOWN_SORT UNKNOWN_FLAGS\n",
            argv[0]);
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
    goto entries_done;
  }
  entries[1].beg_voff++;
  if (write_index(argv[6], &idx.header, entries, idx.strings) != 0) {
    goto entries_done;
  }

  memcpy(entries, idx.entries, (size_t)idx.header.n_blocks * sizeof(*entries));
  entries[0].first_name_offset = idx.entries[1].last_name_offset;
  entries[0].last_name_offset = idx.entries[0].first_name_offset;
  if (write_index(argv[7], &idx.header, entries, idx.strings) != 0) {
    goto entries_done;
  }

  memcpy(entries, idx.entries, (size_t)idx.header.n_blocks * sizeof(*entries));
  entries[1].first_name_offset = idx.entries[0].first_name_offset;
  entries[1].last_name_offset = idx.entries[0].first_name_offset;
  if (write_index(argv[8], &idx.header, entries, idx.strings) != 0) {
    goto entries_done;
  }

  memcpy(entries, idx.entries, (size_t)idx.header.n_blocks * sizeof(*entries));
  entries[0].first_name_offset++;
  if (write_index(argv[9], &idx.header, entries, idx.strings) != 0) {
    goto entries_done;
  }

  char *strings = malloc((size_t)idx.header.strings_size);
  if (strings == NULL) {
    goto entries_done;
  }
  memcpy(strings, idx.strings, (size_t)idx.header.strings_size);
  strings[idx.header.strings_size - 1] = 'X';
  if (write_index(argv[10], &idx.header, idx.entries, strings) != 0) {
    free(strings);
    goto entries_done;
  }
  free(strings);

  header = idx.header;
  header.sort_order = UINT32_MAX;
  if (write_index(argv[11], &header, idx.entries, idx.strings) != 0) {
    goto entries_done;
  }

  header = idx.header;
  header.flags |= UINT32_C(0x80000000);
  if (write_index(argv[12], &header, idx.entries, idx.strings) == 0) {
    status = 0;
  }

entries_done:
  free(entries);

done:
  bni_index_destroy(&idx);
  return status;
}
