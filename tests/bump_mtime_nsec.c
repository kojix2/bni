#define _POSIX_C_SOURCE 200809L

#include <fcntl.h>
#include <stdint.h>
#include <stdio.h>
#include <sys/stat.h>
#include <time.h>

int main(int argc, char **argv) {
  if (argc != 2) {
    fprintf(stderr, "usage: %s FILE\n", argv[0]);
    return 2;
  }

  struct stat before;
  if (stat(argv[1], &before) != 0) {
    perror("stat");
    return 1;
  }
#if defined(__APPLE__)
  struct timespec times[2] = {before.st_atimespec, before.st_mtimespec};
#else
  struct timespec times[2] = {before.st_atim, before.st_mtim};
#endif
  times[1].tv_nsec = times[1].tv_nsec == 999999999L ? 999999998L : times[1].tv_nsec + 1;
  if (utimensat(AT_FDCWD, argv[1], times, 0) != 0) {
    perror("utimensat");
    return 1;
  }

  struct stat after;
  if (stat(argv[1], &after) != 0) {
    perror("stat");
    return 1;
  }
#if defined(__APPLE__)
  long before_nsec = before.st_mtimespec.tv_nsec;
  long after_nsec = after.st_mtimespec.tv_nsec;
#else
  long before_nsec = before.st_mtim.tv_nsec;
  long after_nsec = after.st_mtim.tv_nsec;
#endif
  if (before.st_size != after.st_size || before.st_mtime != after.st_mtime ||
      before_nsec == after_nsec) {
    fprintf(stderr, "filesystem did not preserve a nanosecond-only mtime change\n");
    return 1;
  }
  return 0;
}
