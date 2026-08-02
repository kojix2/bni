#!/usr/bin/env bash
set -euo pipefail

BNI=${1:-./bni}
if ! command -v samtools >/dev/null 2>&1; then
  echo "smoke test skipped: samtools not found" >&2
  exit 0
fi

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

cat > "$tmp/in.sam" <<'SAM'
@HD	VN:1.6	SO:unsorted
@SQ	SN:chr1	LN:1000
read2	0	chr1	20	60	5M	*	0	0	ACGTA	IIIII
read1	0	chr1	10	60	5M	*	0	0	ACGTA	IIIII
read1	2048	chr1	30	60	5M	*	0	0	ACGTA	IIIII
read10	0	chr1	40	60	5M	*	0	0	ACGTA	IIIII
SAM

samtools view -b -o "$tmp/in.bam" "$tmp/in.sam"
samtools sort -N -o "$tmp/in.name.bam" "$tmp/in.bam"
"$BNI" index -f "$tmp/in.name.bam"
"$BNI" check --full "$tmp/in.name.bam"
printf '@HD\tVN:1.6\tSO:queryname\tSS:queryname:lexicographical\n@SQ\tSN:chr1\tLN:1000\n' \
  > "$tmp/empty.sam"
samtools view -b -o "$tmp/empty.bam" "$tmp/empty.sam"
"$BNI" index -f "$tmp/empty.bam"
"$BNI" check --full "$tmp/empty.bam" >/dev/null
"$BNI" get -O sam --no-header "$tmp/in.name.bam" read1 > "$tmp/read1.sam"
count=$(cut -f1 "$tmp/read1.sam" | grep -c '^read1$')
if [ "$count" -ne 2 ]; then
  echo "expected 2 read1 records, got $count" >&2
  exit 1
fi
printf 'read2\nread1\nread1\n' > "$tmp/names.txt"
"$BNI" get -O sam --no-header -f "$tmp/names.txt" "$tmp/in.name.bam" > "$tmp/many.sam"
actual=$(cut -f1 "$tmp/many.sam" | tr '\n' ' ')
if [ "$actual" != "read1 read1 read2 " ]; then
  echo "unexpected batched name-file output order/counts: $actual" >&2
  exit 1
fi
if "$BNI" get -O sam --no-header -f "$tmp" "$tmp/in.name.bam" \
  > /dev/null 2> "$tmp/name-read-error.err"; then
  echo "expected name-file read error to fail" >&2
  exit 1
fi
if ! grep -q 'failed reading name file' "$tmp/name-read-error.err"; then
  echo "expected name-file read error message" >&2
  exit 1
fi
printf 'read11\n' > "$tmp/missing_names.txt"
"$BNI" get -O sam --no-header --missing-ok --list-missing -f "$tmp/missing_names.txt" "$tmp/in.name.bam" > "$tmp/missing.sam" 2> "$tmp/missing.err"
if ! grep -qx 'read11' "$tmp/missing.err"; then
  echo "expected read11 to be reported missing" >&2
  exit 1
fi
"$BNI" stats "$tmp/in.name.bam" >/dev/null

expect_collision() {
  local description=$1
  local protected_file=$2
  local backup_file=$3
  shift 3
  if "$@" 2> "$tmp/collision.err"; then
    echo "expected $description collision to fail" >&2
    exit 1
  fi
  if ! grep -q 'would overwrite' "$tmp/collision.err"; then
    echo "expected $description collision error" >&2
    exit 1
  fi
  if ! cmp -s "$protected_file" "$backup_file"; then
    echo "$description collision modified the protected file" >&2
    exit 1
  fi
}

cp "$tmp/in.name.bam" "$tmp/protected.bam"
cp "$tmp/protected.bam" "$tmp/protected.bam.backup"
cp "$tmp/in.name.bam.bni" "$tmp/protected.bam.bni"
expect_collision "BAM output" "$tmp/protected.bam" "$tmp/protected.bam.backup" \
  "$BNI" get -O sam --no-header -o "$tmp/protected.bam" "$tmp/protected.bam" read1

cp "$tmp/protected.bam.bni" "$tmp/protected.bam.bni.backup"
expect_collision "index output" "$tmp/protected.bam.bni" "$tmp/protected.bam.bni.backup" \
  "$BNI" get -O sam --no-header -o "$tmp/protected.bam.bni" "$tmp/protected.bam" read1

printf 'read1\n' > "$tmp/protected.names"
cp "$tmp/protected.names" "$tmp/protected.names.backup"
expect_collision "name-file output" "$tmp/protected.names" "$tmp/protected.names.backup" \
  "$BNI" get -O sam --no-header -f "$tmp/protected.names" -o "$tmp/protected.names" \
  "$tmp/protected.bam"

ln "$tmp/protected.bam" "$tmp/protected.hardlink.bam"
expect_collision "hard-link output" "$tmp/protected.bam" "$tmp/protected.bam.backup" \
  "$BNI" get -O sam --no-header -o "$tmp/protected.hardlink.bam" "$tmp/protected.bam" read1

ln -s "$tmp/protected.bam" "$tmp/protected.symlink.bam"
expect_collision "symbolic-link output" "$tmp/protected.bam" "$tmp/protected.bam.backup" \
  "$BNI" get -O sam --no-header -o "$tmp/protected.symlink.bam" "$tmp/protected.bam" read1

cp "$tmp/in.name.bam" "$tmp/index-protected.bam"
cp "$tmp/index-protected.bam" "$tmp/index-protected.bam.backup"
expect_collision "index command output" "$tmp/index-protected.bam" \
  "$tmp/index-protected.bam.backup" "$BNI" index -f -o "$tmp/index-protected.bam" \
  "$tmp/index-protected.bam"

cp "$tmp/in.name.bam.bni" "$tmp/atomic-failure.bni"
cp "$tmp/atomic-failure.bni" "$tmp/atomic-failure.bni.backup"
if (trap '' XFSZ; ulimit -f 0; \
  "$BNI" index -f -o "$tmp/atomic-failure.bni" "$tmp/in.name.bam") 2>/dev/null; then
  echo "expected size-limited index write to fail" >&2
  exit 1
fi
if ! cmp -s "$tmp/atomic-failure.bni" "$tmp/atomic-failure.bni.backup"; then
  echo "failed atomic write modified the existing index" >&2
  exit 1
fi
if find "$tmp" -maxdepth 1 -name 'atomic-failure.bni.tmp.*' -print -quit | grep -q .; then
  echo "failed atomic write left a temporary index" >&2
  exit 1
fi

cc -Wall -Wextra -std=c11 tests/bump_mtime_nsec.c -o "$tmp/bump_mtime_nsec"
cp "$tmp/in.name.bam" "$tmp/nsec-stale.bam"
"$BNI" index -f "$tmp/nsec-stale.bam"
"$tmp/bump_mtime_nsec" "$tmp/nsec-stale.bam"
if "$BNI" check "$tmp/nsec-stale.bam" > /dev/null 2> "$tmp/nsec-check.err"; then
  echo "expected nanosecond-only mtime change to fail metadata check" >&2
  exit 1
fi
if ! grep -q 'mtime nanoseconds mismatch' "$tmp/nsec-check.err"; then
  echo "expected nanosecond mtime mismatch error from check" >&2
  exit 1
fi
if "$BNI" get -O sam --no-header "$tmp/nsec-stale.bam" read1 \
  > /dev/null 2> "$tmp/nsec-get.err"; then
  echo "expected nanosecond-only mtime change to fail get" >&2
  exit 1
fi
if ! grep -q 'mtime nanoseconds differ' "$tmp/nsec-get.err"; then
  echo "expected nanosecond mtime mismatch error from get" >&2
  exit 1
fi

if command -v pkg-config >/dev/null 2>&1 && pkg-config --exists htslib; then
  cc $(pkg-config --cflags htslib) -Iinclude tests/library_api.c libbni.a $(pkg-config --libs htslib) -o "$tmp/library_api"
  "$tmp/library_api" "$tmp/in.name.bam" read1 2

  awk 'BEGIN {
    print "@HD\tVN:1.6\tSO:queryname\tSS:queryname:lexicographical"
    print "@SQ\tSN:chr1\tLN:1000"
    for (i = 0; i < 3000; ++i) {
      printf "bulk%05d\t4\t*\t0\t0\t*\t*\t0\t0\tACGTA\tIIIII\n", i
    }
  }' > "$tmp/bulk.sam"
  samtools view -b -o "$tmp/bulk.bam" "$tmp/bulk.sam"
  "$BNI" index -f "$tmp/bulk.bam"
  "$BNI" check --full "$tmp/bulk.bam" >/dev/null

  cc $(pkg-config --cflags htslib) -Iinclude tests/make_bad_indexes.c libbni.a \
    $(pkg-config --libs htslib) -o "$tmp/make_bad_indexes"
  "$tmp/make_bad_indexes" "$tmp/bulk.bam.bni" "$tmp/empty.bni" \
    "$tmp/empty-count.bni" "$tmp/count.bni" "$tmp/truncated.bni" "$tmp/gap.bni" \
    "$tmp/entry-order.bni" "$tmp/range-order.bni" "$tmp/inside-offset.bni" \
    "$tmp/unterminated.bni"

  expect_bad_full_check() {
    local description=$1
    local index_file=$2
    local expected_error=$3
    if "$BNI" check --full -i "$index_file" "$tmp/bulk.bam" \
      > "$tmp/bad-check.out" 2> "$tmp/bad-check.err"; then
      echo "expected $description full check to fail" >&2
      exit 1
    fi
    if ! grep -q "$expected_error" "$tmp/bad-check.err"; then
      echo "unexpected $description full-check error" >&2
      exit 1
    fi
  }

  expect_bad_full_check "empty index" "$tmp/empty.bni" 'not covered by the index'
  expect_bad_full_check "empty index record count" "$tmp/empty-count.bni" \
    'index record count mismatch'
  expect_bad_full_check "header record count" "$tmp/count.bni" 'index record count mismatch'
  expect_bad_full_check "truncated index" "$tmp/truncated.bni" 'not covered by the index'
  expect_bad_full_check "entry gap" "$tmp/gap.bni" 'non-contiguous virtual-offset ranges'
  expect_bad_full_check "entry QNAME order" "$tmp/entry-order.bni" 'decreasing QNAME range'
  expect_bad_full_check "QNAME range order" "$tmp/range-order.bni" \
    'QNAME ranges decrease between entries'
  expect_bad_full_check "inside-string offset" "$tmp/inside-offset.bni" \
    'is not at a string boundary'
  expect_bad_full_check "unterminated string" "$tmp/unterminated.bni" \
    'string table is not NUL-terminated'
else
  echo "library API smoke test skipped: pkg-config htslib not found" >&2
fi

echo "smoke test OK" >&2
