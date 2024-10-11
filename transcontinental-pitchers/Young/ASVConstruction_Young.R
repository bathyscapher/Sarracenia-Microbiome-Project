# SMP ASV Young data ###########################################################

library("dada2")
library("ShortRead")
library("DECIPHER")
library("gridExtra")


rm(list = ls())


primer <- "16S"
# primer <- "18S"

if (primer == "16S") {
  wd <- "Young/16S"
} else if (primer == "18S") {
  wd <- "Young/18S"
}


ncore <- 6 # specify number of available cores


## Read fastq ##################################################################
list.files(path = wd, pattern = "fastq")

rF <- sort(list.files(path = wd, pattern = "R1", full.names = TRUE))
rR <- sort(list.files(path = wd, pattern = "R2", full.names = TRUE))


## Extract sample names
sample.names <- sapply(strsplit(basename(rF), "_"), `[`, 2)


## Check for primers ###########################################################
## Remove Ns from reads
rF.fN <- file.path(paste0(wd, "/filtN"), basename(rF))
rR.fN <- file.path(paste0(wd, "/filtN"), basename(rR))

filterAndTrim(rF, rF.fN, rR, rR.fN, maxN = 0, multithread = TRUE)


## Forward and reverse primer
if (primer == "16S") {
  FWD <- "ACTCCTACGGRAGGCAGCAG" # F338
  REV <- "TACNVGGGTATCTAATCC" # R802
} else if (primer == "18S") {
  FWD <- "TCCAAGGAAGGCAGCAGG" # F426
  REV <- "AGTCCTATTCCATTATTCCATG" # R853
}


## Compile all orientations of the primers
allOrients <- function(primer) {
  dna <- DNAString(primer)
  orients <- c(Forward = dna, Complement = complement(dna),
               Reverse = reverse(dna), RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)


## Count occurence of all primer orientations in the reads
primerHits <- function(primer, fn) {
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}


## Cut the primers
cutadapt <- "/usr/bin/cutadapt"
# system2(cutadapt, args = "--version")


path.cut <- file.path(wd, "cutPrimers")
if(!dir.exists(path.cut)) {
  dir.create(path.cut)
}

rF.cut <- file.path(path.cut, basename(rF.fN))
rR.cut <- file.path(path.cut, basename(rR.fN))


## Trim FWD and the reverse-complement of REV off of R1 &  REV and the
## reverse-complement of FWD off of R2
FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

R1.flags <- paste("-g", FWD, "-a", REV.RC)
R2.flags <- paste("-G", REV, "-A", FWD.RC)


## Run Cutadapt
for(i in seq_along(rF)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2,
                             "-o", rF.cut[i], "-p", rR.cut[i], # output
                             rF.fN[i], rR.fN[i])) # input
}


## Compare primers before and after trimming
(reads <- ceiling(runif(1, 1, length(rF.cut))))
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = rF.fN[[reads]]),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = rR.fN[[reads]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = rF.fN[[reads]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = rR.fN[[reads]]))
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = rF.cut[[reads]]),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = rR.cut[[reads]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = rF.cut[[reads]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = rR.cut[[reads]]))


### Filter and trim. Place filtered files in filtered/ subdirectory
if (primer == "16S") {
  rF.cut.f <- file.path(paste0(wd, "/filtered"),
                        paste0(sample.names, "_16S_R1_filt.fastq.gz"))
  rR.cut.f <- file.path(paste0(wd, "/filtered"),
                        paste0(sample.names, "_16S_R2_filt.fastq.gz"))
} else if (primer == "18S") {
  rF.cut.f <- file.path(paste0(wd, "/filtered"),
                        paste0(sample.names, "_18S_R1_filt.fastq.gz"))
  rR.cut.f <- file.path(paste0(wd, "/filtered"),
                        paste0(sample.names, "_18S_R2_filt.fastq.gz"))
}


names(rF.cut.f) <- sample.names
names(rR.cut.f) <- sample.names


## Quality filtering ###########################################################
out <- filterAndTrim(rF.fN, rF.cut.f, rR.fN, rR.cut.f,
                     maxN = 0, maxEE = c(2, 2), minLen = 100,
                     truncQ = 2, rm.phix = TRUE,
                     compress = TRUE, multithread = ncore, verbose = TRUE)


### Read quality profiles exemplarily
set.seed(74398)
startPlot <- ceiling(runif(1, 0, length(rF) - 3))

plot.rF <- plotQualityProfile(rF.fN[startPlot:(startPlot + 2)])
plot.rF.f <- plotQualityProfile(rF.cut.f[startPlot:(startPlot + 2)])

plot.rR <- plotQualityProfile(rR.fN[startPlot:(startPlot + 2)])
plot.rR.f <- plotQualityProfile(rR.cut.f[startPlot:(startPlot + 2)])

grid.arrange(plot.rF, plot.rR, plot.rF.f, plot.rR.f, ncol = 2)


## Estimate and plot the error rates ###########################################
errF <- learnErrors(rF.cut.f, multithread = ncore, verbose = TRUE)
errR <- learnErrors(rR.cut.f, multithread = ncore, verbose = TRUE)


plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)


## Core sample inference algorithm #############################################
dadaF <- dada(rF.cut.f, err = errF, multithread = ncore)
dadaR <- dada(rR.cut.f, err = errR, multithread = ncore)
# dadaFs[[1]]


if (primer == "16S") {
  saveRDS(dadaF, paste0(wd, "/dadaF_EY_16S.rds"))
  saveRDS(dadaR, paste0(wd, "/dadaR_EY_16S.rds"))
} else if (primer == "18S") {
  saveRDS(dadaF, paste0(wd, "/dadaF_EY_18S.rds"))
  saveRDS(dadaR, paste0(wd, "/dadaR_EY_18S.rds"))
}


### Merge paired reads
contigs <- mergePairs(dadaF, rF.cut.f, dadaR, rR.cut.f, verbose = TRUE)
head(contigs[[1]])


if (primer == "16S") {
  saveRDS(contigs, paste0(wd, "/contigs_EY_16S.rds"))
} else if (primer == "18S") {
  saveRDS(contigs, paste0(wd, "/contigs_EY_18S.rds"))
}


## Construct ASV table #########################################################
seqtab <- makeSequenceTable(contigs)
dim(seqtab)


### Distribution of sequence lengths
table(nchar(getSequences(seqtab)))


## Chimera detection ###########################################################
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                    multithread = ncore, verbose = TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim) / sum(seqtab)


if (primer == "16S") {
  saveRDS(seqtab.nochim, paste0(wd, "/seqtab.nochim_EY_16S.rds"))
  # seqtab.nochim <- readRDS(paste0(wd, "/seqtab.nochim_EY_16s.rds"))
} else if (primer == "18S") {
  saveRDS(seqtab.nochim, paste0(wd, "/seqtab.nochim_EY_18S.rds"))
  # seqtab.nochim <- readRDS(paste0(wd, "/seqtab.nochim_EY_18s.rds"))
}


### Track reads through the pipeline
getN <- function(x){
  sum(getUniques(x))
  }

track <- cbind(out, sapply(dadaF, getN), sapply(dadaR, getN),
               sapply(contigs, getN), rowSums(seqtab.nochim))

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged",
                     "nonchim")
rownames(track) <- sample.names
head(track)


## Assign taxonomy #############################################################
### RDP classifier #############################################################
taxa <- assignTaxonomy(seqtab.nochim,
                       "refs/silva_nr_v132_train_set.fa.gz",
                       multithread = TRUE, verbose = TRUE)

if (primer == "16S") {
  saveRDS(taxa, paste0(wd, "/taxa_EY_16S.rds"))
} else if (primer == "18S") {
  saveRDS(taxa, paste0(wd, "/taxa_EY_18S.rds"))
}


### IdTaxa and SILVA ###########################################################
dna <- DNAStringSet(getSequences(seqtab.nochim))

load("refs/SILVA_SSU_r138_2019.RData")


ids <- IdTaxa(dna, trainingSet, strand = "top", processors = ncore,
              verbose = TRUE)
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")


### Convert output object of class "Taxa"
taxa.id <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa
  }))

colnames(taxa.id) <- ranks
rownames(taxa.id) <- getSequences(seqtab.nochim)


if (primer == "16S") {
  saveRDS(taxa.id, paste0(wd, "/taxa.id_EY_16S.rds"))
} else if (primer == "18S") {
  saveRDS(taxa.id, paste0(wd, "/taxa.id_EY_18S.rds"))
}
