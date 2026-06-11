# Examples

```sh
samtools sort -N -o reads.name.bam reads.bam
bni index -@ 4 reads.name.bam
bni stats reads.name.bam
bni get -O sam reads.name.bam 'read-id-here' | less
```

For a list of read names:

```sh
bni get -f names.txt -o subset.bam reads.name.bam
```
