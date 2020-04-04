################################################################################
################################################################################
################################################################################
################################################################################
### SMP ASV Young Greenhouse data
### Author: korn@cumulonimbus.at University of Fribourg 2020
################################################################################

rm(list = ls())
setwd("/scratch/pSMP/Young_Greenhouse/")


library("dada2")
library("DECIPHER")
library("gridExtra")


ncore <- 8 # specify number of available cores

################################################################################
### Read fastq
list.files(pattern = "fastq")

rF <- sort(list.files(pattern = "R1", full.names = TRUE))
rR <- sort(list.files(pattern = "R2", full.names = TRUE))


## Extract sample names
sample.names <- sapply(strsplit(basename(rF), "_"), `[`, 1)


### Filter and trim. Place filtered files in filtered/ subdirectory
rF.f <- file.path("filtered", paste0(sample.names, "_R1_filt.fastq.gz"))
rR.f <- file.path("filtered", paste0(sample.names, "_R2_filt.fastq.gz"))

names(rF.f) <- sample.names
names(rR.f) <- sample.names


### Filtering
out <- filterAndTrim(rF, rF.f, rR, rR.f,
                     maxN = 0, maxEE = c(2, 2),
                     truncQ = 2, rm.phix = TRUE,
                     compress = TRUE, multithread = ncore, verbose = TRUE)


### Read quality profiles exemplarily
set.seed(74398)
startPlot <- ceiling(runif(1, 0, length(rF) - 3))

plot.rF <- plotQualityProfile(rF[startPlot:(startPlot + 2)])
plot.rF.f <- plotQualityProfile(rF.f[startPlot:(startPlot + 2)])

plot.rR <- plotQualityProfile(rR[startPlot:(startPlot + 2)])
plot.rR.f <- plotQualityProfile(rR.f[startPlot:(startPlot + 2)])

grid.arrange(plot.rF, plot.rR, plot.rF.f, plot.rR.f, ncol = 2)


### Estimate and plot the error rates
errF <- learnErrors(rF.f, multithread = ncore, verbose = TRUE)
errR <- learnErrors(rR.f, multithread = ncore, verbose = TRUE)

plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)


### Core sample inference algorithm
dadaF <- dada(rF.f, err = errF, multithread = ncore)
dadaR <- dada(rR.f, err = errR, multithread = ncore)
# dadaFs[[1]]


# saveRDS(dadaF, "dadaF_EY_GH.rds")
# saveRDS(dadaR, "dadaR_EY_GH.rds")


### Merge paired reads
contigs <- mergePairs(dadaF, rF.f, dadaR, rR.f, verbose = TRUE)
head(contigs[[1]])


# saveRDS(contigs, "contigs_EY_GH.rds")


### Construct ASV table
seqtab <- makeSequenceTable(contigs)
dim(seqtab)


### Distribution of sequence lengths
table(nchar(getSequences(seqtab)))


### Chimera detection
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                    multithread = ncore, verbose = TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim) / sum(seqtab)


saveRDS(seqtab.nochim, "seqtab.nochim_EY_GH.rds")


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


### Assign taxonomy with IdTaxa and SILVA
dna <- DNAStringSet(getSequences(seqtab.nochim))

load("~/Desktop/SMP_unsynced/silva/SILVA_SSU_r138_2019.RData")


ids <- IdTaxa(dna, trainingSet, strand = "top", processors = ncore,
              verbose = TRUE)
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")


### Convert output object of class "Taxa"
taxa.id <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  # taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
  }))

colnames(taxa.id) <- ranks
rownames(taxa.id) <- getSequences(seqtab.nochim)

# saveRDS(taxa.id, "taxa.id_EY_GH.rds")


################################################################################
################################################################################
################################################################################
################################################################################
