# Preprocess sequence files
## Demultiplexing
If necessary (Bittleston), samples were demultiplexed with [https://github.com/yhwu/idemp](idemp):
```bash
idemp -b prok.bc -I1 *I1*.fastq -R1 *R1*.fastq -R2 *R2*.fastq -m n -o demux
```

...where prok.bc holds the barcodes with the affiliated samples in following, tab-delimited format:

```bash
Barcode	Sampleid
ACTGTTTACTGT	Sample01
CAGGCCACTCTC	Sample02
ACCCAAGCGTTA    Sample03
...	...
```

## SFF to FASTQ
If necessary (Boynton), SFF files were converted to FASTQ files  with [https://github.com/indraniel/sff2fastq](sff2fastq):
```bash
sff2fastq file.sff > file.fastq
```
