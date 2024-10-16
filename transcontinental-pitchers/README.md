# Workflow
## Preprocess FASTQ files
### Demultiplex Illumina FASTQ with R1, R2 and I1
If necessary (Bittleston), samples were demultiplexed with [mothur](https://mothur.org/):

```bash
fastq.info(file=lb.file, oligos=lb_16S_ext.oligo, bdiffs=1, fasta=f, qfile=f)
fastq.info(file=lb.file, oligos=lb_18S_ext_rc.oligo, bdiffs=1, fasta=f, qfile=f)
```

...where lb.file holds the tab separated names of R1, R2 and I1 (`none` specifies the missing reverse index I2 file):
```bash
*R1*.fastq *R2_001*.fastq *I1*.fastq none
```

...and lb_16S_ext.oligo/lb_18S_ext_rc.oligo the tab separated barcodes with affiliated samples in the format (whereby for the 18S the reverse complements of the barcodes is needed and `NONE` specifies the missing reverse index file I2):
```bash
barcode	ACTGTTTACTGT	NONE	Sample01
barcode	CAGGCCACTCTC	NONE	Sample02
barcode	ACCCAAGCGTTA	NONE	Sample03
...	...	...	...
```

### Demultiplex Illumina FASTQ with R1 and R2 (without index files I1 or I2)
If necessary (Freedman), samples were demultiplexed with [mothur](https://mothur.org/):

```bash
set.current(processors=6)
fastq.info(file=zf.Samp1.file, oligos=zf.oligo, pdiffs=2, bdiffs=2, fasta=f, qfile=f)
fastq.info(file=zf.Sam55.file, oligos=zf.oligo, pdiffs=2, bdiffs=2, fasta=f, qfile=f)
```

...where `zf.\*.file` holds the tab separated names of R1 and R2 (`none none` specifies the missing index files):
```bash
*R1*.fastq *R2_001*.fastq none none
```

...with the tab separated `zf.oligo` in the format:
```bash
barcode	ACTGTTTACTGT	NONE	Sample01
barcode	CAGGCCACTCTC	NONE	Sample02
barcode	ACCCAAGCGTTA	NONE    Sample03
...	...	...	...
```


### Demultiplex 454 Roche SFF and convert to FASTQ
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

### Boynton
```bash
gzip *fastq
rename -n 's/temporal.trim.sorted.//' *q.gz
rename -n 's/geographic.trim.sorted.//' *q.gz
```

### Freedman
```bash
gzip *fastq

## For Samp1*
rename -n 's/Samp1-54_S3_L001_//' *q.gz
## For Sam55*
rename -n 's/Sam55-108_S4_L001_//' *q.gz

## For both
#FIX THIS: rename -n 's/{Samp1-54_S3|Sam55-108_S4}_L001_//' *q.gz

rename -n 's/001.//' *q.gz

rename -n 's/\.(?=[^.]*\.)/-/g' '{}' *q.gz # replace dots with dashes, skip file extension
rename -n 's/-fastq/.fastq/' *q.gz

```

### Korn
```bash
rename 's/MI.M03\d{3}[.|_]?0\d{3}.001.FLD0\d{3}.//' *fastq.gz
for run in {1..3}; do rename 's/_//' *fastq.gz; done
```

### Young
```bash
rename -n 's/Young16S\d{1,2}_/16S-/' *fastq
rename -n 's/Young18S\d{1,2}_/18S-/' *fastq

rename -n 's/L001_//' *fastq
rename -n 's/_001//' *fastq
```


### Young Greenhouse
```bash
for run in {1..2}; do rename 's/_/-/' *fastq.gz; done
```

## ASV construction
[ASVs](https://en.wikipedia.org/wiki/Amplicon_sequence_variant) were constructed from the raw reads mainly following the [dada2 tutorial](https://benjjneb.github.io/dada2/tutorial.html) for reads sequenced with 16S- and 18S primers and
the [ITS data2 tutorial](https://benjjneb.github.io/dada2/ITS_workflow.html) for reads sequenced with ITS primers.
ASVs were classified with [SILVA 138](https://www.arb-silva.de/download/arb-files/) and processed with the R package [phyloseq](https://joey711.github.io/phyloseq/index.html).
