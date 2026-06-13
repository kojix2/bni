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

if command -v pkg-config >/dev/null 2>&1 && pkg-config --exists htslib; then
  cc $(pkg-config --cflags htslib) -Iinclude tests/library_api.c libbni.a $(pkg-config --libs htslib) -o "$tmp/library_api"
  "$tmp/library_api" "$tmp/in.name.bam" read1 2
else
  echo "library API smoke test skipped: pkg-config htslib not found" >&2
fi

echo "smoke test OK" >&2
