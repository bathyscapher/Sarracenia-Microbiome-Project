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
sffinfo(sff=geographic.sff, oligos=its_geo.oligo, flow=F)
make.fastq(fasta=geographic.fasta)
fastq.info(fastq=geographic.fastq, oligos=its_geo.oligo)

sffinfo(sff=temporal.sff, oligos=its_temp.oligo, flow=F)
make.fastq(fasta=temporal.fasta)
fastq.info(fastq=temporal.fastq, oligos=its_temp.oligo)
```

...with the tab separated oligo file in the format:
```
linker	GTGTGYCAGCMGCCGCGGTAA
ACTGTTTACTGT	Sample01
CAGGCCACTCTC	Sample02
ACCCAAGCGTTA    Sample03
...	...
```


## Simplify filenames
### Bittleston
```bash
rename -n 's/Undetermined_S0_L001_//' *q.gz
rename -n 's/001.//' *q.gz
```

## Boynton
```bash
gzip *fastq
rename -n 's/temporal.//' *q.gz
rename -n 's/geographic.//' *q.gz
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
