# SMP ASV Goddard data #########################################################

library("dada2")
library("ShortRead")
library("DECIPHER")
library("gridExtra")


rm(list = ls())
gc()

primer <- "16S"
# primer <- "18S"

if (primer == "16S") {
  wd <- "Goddard/16S"
} else if (primer == "18S") {
  wd <- "Goddard/18S"
}


ncore <- 6 # specify number of available cores


## Read fastq ##################################################################
list.files(path = wd, pattern = "fastq")

rF <- sort(list.files(path = wd, pattern = "R1", full.names = TRUE))
rR <- sort(list.files(path = wd, pattern = "R2", full.names = TRUE))


## Extract sample names
sample.names <- sapply(strsplit(basename(rF), "_"), `[`, 1)

sample.names <- gsub("_[A-Z]{6}-[A-Z]{6}_L001_R[1|2]_001.fastq.gz",
                          "",
                          basename(rF))


## Check for primers ###########################################################
## Remove Ns from reads
rF.fN <- file.path(paste0(wd, "/filtN"), basename(rF))
rR.fN <- file.path(paste0(wd, "/filtN"), basename(rR))

filterAndTrim(rF, rF.fN, rR, rR.fN, maxN = 0, multithread = TRUE)


## Forward and reverse primer ## ADD CORRECT PRIMERS ####
if (primer == "16S") {
  FWD <- "GTGYCAGCMGCCGCGGTAA" #
  REV <- "GGACTACNVGGGTWTCTAAT" #
} else if (primer == "18S") {
  FWD <- "CCCTGCCHTTTGTACACAC" #
  REV <- "CCTTCYGCAGGTTCACCTAC" #
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


path.cut <- file.path(".", "cutPrimers")
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


## Compare reads before and after trimming
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
  rF.cut.f <- file.path("filtered",
                        paste0(sample.names, "_16S_R1_filt.fastq.gz"))
  rR.cut.f <- file.path("filtered",
                        paste0(sample.names, "_16S_R2_filt.fastq.gz"))
} else if (primer == "18S") {
  rF.cut.f <- file.path("filtered",
                        paste0(sample.names, "_18S_R1_filt.fastq.gz"))
  rR.cut.f <- file.path("filtered",
                        paste0(sample.names, "_18S_R2_filt.fastq.gz"))
}

names(rF.cut.f) <- sample.names
names(rR.cut.f) <- sample.names


## Quality filtering ###########################################################
out <- filterAndTrim(rF.fN, rF.cut.f, rR.fN, rR.cut.f,
                     maxN = 0, maxEE = c(2, 2), minLen = 100,
                     truncQ = 2, rm.phix = TRUE,
                     compress = TRUE, multithread = ncore, verbose = TRUE)


### Plot quality profiles exemplarily
set.seed(74398)
startPlot <- ceiling(runif(1, 0, length(rF) - 3))

plot.rF <- plotQualityProfile(rF.fN[startPlot:(startPlot + 2)])
plot.rF.f <- plotQualityProfile(rF.cut.f[startPlot:(startPlot + 2)])

plot.rR <- plotQualityProfile(rR.fN[startPlot:(startPlot + 2)])
plot.rR.f <- plotQualityProfile(rR.cut.f[startPlot:(startPlot + 2)])

grid.arrange(plot.rF, plot.rR,
             plot.rF.f, plot.rR.f,
             ncol = 2)


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
  saveRDS(dadaF, "dadaF_RK_16S.rds")
  saveRDS(dadaR, "dadaR_RK_16S.rds")
} else if (primer == "18S") {
  saveRDS(dadaF, "dadaF_RK_18S.rds")
  saveRDS(dadaR, "dadaR_RK_18S.rds")
}


## Merge paired reads for 16S, concatenate them for 18S ########################
if (primer == "16S") {
  contigs <- mergePairs(dadaF, rF.cut.f, dadaR, rR.cut.f, verbose = TRUE)
} else if (primer == "18S") {
  contigs <- mergePairs(dadaF, rF.cut.f, dadaR, rR.cut.f, verbose = TRUE,
                        justConcatenate = TRUE)
}

head(contigs[[1]])


if (primer == "16S") {
  saveRDS(contigs, "contigs_RK_16S.rds")
} else if (primer == "18S") {
  saveRDS(contigs, "contigs_RK_18S.rds")
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
  saveRDS(seqtab.nochim, "seqtab.nochim_RK_16S.rds")
  # seqtab.nochim <- readRDS("seqtab.nochim_RK_16S.rds")
} else if (primer == "18S") {
  saveRDS(seqtab.nochim, "seqtab.nochim_RK_18S.rds")
  # seqtab.nochim <- readRDS("seqtab.nochim_RK_18S.rds")
}


## Track reads through the pipeline ############################################
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
                       "../../refs/",
                       multithread = TRUE, verbose = TRUE)


if (primer == "16S") {
  saveRDS(taxa, "taxa_RK_16S.rds")
} else if (primer == "18S") {
  saveRDS(taxa, "taxa_RK_18S.rds")
}


### IdTaxa and SILVA ###########################################################
dna <- DNAStringSet(getSequences(seqtab.nochim))

load("../../SILVA_SSU_r138_2019.RData")

ids <- IdTaxa(dna,
              trainingSet,
              strand = "top",
              processors = ncore,
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
  saveRDS(taxa.id, "taxa.id_RK_16S.rds")
} else if (primer == "18S") {
  saveRDS(taxa.id, "taxa.id_RK_18S.rds")
}
