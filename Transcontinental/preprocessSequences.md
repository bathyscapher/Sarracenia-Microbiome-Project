# Preprocess sequence files
## Demultiplex Illumina FASTQ with R1, R2 and I1
If necessary (Bittleston), samples were demultiplexed with [mothur](https://mothur.org/):

```bash
fastq.info(file=lb.file, oligos=lb_16S_ext.oligo, bdiffs=1, fasta=t, qfile=t)
fastq.info(file=lb.file, oligos=lb_18S_ext_rc.oligo, bdiffs=1, fasta=t, qfile=t)
```

...where lb.file holds the tab separated names of R1, R2 and I1 filenames (`none` specifies the missing reverse index I2 file):
```bash
*R1*.fastq *R2_001*.fastq *I1*.fastq none
```

...and lb_16S_ext.oligo/lb_18S_ext_rc.oligo the tab separated barcodes with affiliated samples in the format (whereby for the 18S the reverse complements of the barcodes is needed and `NONE` specifies the missing reverse index file I2):
```bash
barcode	ACTGTTTACTGT	NONE	Sample01
barcode	CAGGCCACTCTC	NONE	Sample02
barcode	ACCCAAGCGTTA	NONE    Sample03
...	...	...	...
```


## Demultiplex 454 Roche SFF and convert to FASTQ
Demultiplex into SFFs and after, convert to FASTQ with [mothur](https://mothur.org/):
```mothur
sffinfo(sff=temporal.sff, flow=F)
trim.seqs(fasta=temporal.fasta, qfile=temporal.qual, oligos=temporal.oligo, bdiffs=1)
sort.seqs(fasta=temporal.trim.fasta, qfile=temporal.trim.qual)
make.fastq(fasta=temporal.trim.sorted.fasta, qfile=temporal.trim.sorted.qual)
split.groups(fastq=temporal.trim.sorted.fastq, group=temporal.groups)

sffinfo(sff=geographic.sff, flow=F)
trim.seqs(fasta=geographic.fasta, qfile=geographic.qual, oligos=geographic.oligo, bdiffs=1)
sort.seqs(fasta=geographic.trim.fasta, qfile=geographic.trim.qual)
make.fastq(fasta=geographic.trim.sorted.fasta, qfile=geographic.trim.sorted.qual)
split.groups(fastq=geographic.trim.sorted.fastq, group=geographic.groups)
```

...with the tab separated oligo file in the format:
```
barcode	ACTGTTTACTGT	Sample01
barcode	CAGGCCACTCTC	Sample02
barcode	ACCCAAGCGTTA    Sample03
...	...	...
```


## Simplify filenames
### Bittleston
```bash
gzip *fastq
rename -n 's/Undetermined_S0_L001_//' *q.gz
rename -n 's/001.//' *q.gz

mv R1_Sro16..fastq.gz R1_Sro162.fastq.gz
mv R2_Sro16..fastq.gz R2_Sro162.fastq.gz
```

## Boynton
```bash
gzip *fastq
rename -n 's/temporal.trim.sorted.//' *q.gz
rename -n 's/geographic.trim.sorted.//' *q.gz
```

## Young
```bash
rename -n 's/Young16S\d{1,2}_/16S-/' *fastq
rename -n 's/Young18S\d{1,2}_/18S-/' *fastq

rename -n 's/L001_//' *fastq
rename -n 's/_001//' *fastq
```


## Young Greenhouse
```bash
for run in {1..2}; do rename 's/_/-/' *fastq.gz; done
```
