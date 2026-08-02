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

if command -v pkg-config >/dev/null 2>&1 && pkg-config --exists htslib; then
  cc $(pkg-config --cflags htslib) -Iinclude tests/library_api.c libbni.a $(pkg-config --libs htslib) -o "$tmp/library_api"
  "$tmp/library_api" "$tmp/in.name.bam" read1 2
else
  echo "library API smoke test skipped: pkg-config htslib not found" >&2
fi

echo "smoke test OK" >&2
