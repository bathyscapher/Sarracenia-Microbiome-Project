################################################################################
################################################################################
################################################################################
################################################################################
### SMP ASV Bittleson data
### Author: korn@cumulonimbus.at University of Fribourg 2020
################################################################################

rm(list = ls())


# library("rSFFreader")
library("dada2")
library("DECIPHER")
library("gridExtra")


# setwd("/scratch/pSMP/Boynton/geo/")
setwd("/scratch/pSMP/Boynton/temp/")


################################################################################
## SFF read lenght
# sff.geo <- readSffGeometry("/scratch/pSMP/Boynton/geo/geographic.sff")
# sff.geo$nReads
#
# sff.temp <- readSffGeometry("/scratch/pSMP/Boynton/temp/temporal.sff")
# sff.temp$nReads
#
# par(mfrow = c(2,1))
# plot(sff.geo$Read_Widths)
# abline(h = mean(sff.geo$Read_Widths), col = "red")
#
# plot(sff.temp$Read_Widths)
# abline(h = mean(sff.temp$Read_Widths), col = "red")
# dev.off()


################################################################################
## Specify number of available cores
ncore <- 8


################################################################################
### Read fastq
list.files(pattern = "fastq.gz")

rF <- sort(list.files(pattern = "gz", full.names = TRUE))


## Extract sample names
sample.names <- sapply(strsplit(basename(rF), "\\."), `[`, 1)


### Filter and trim. Place filtered files in filtered/ subdirectory
rF.f <- file.path("filtered", paste0(sample.names, "_filt.fastq.gz"))

names(rF.f) <- sample.names


### Filtering, min and max values as in Boynton et al. 2019
filt <- filterAndTrim(rF, rF.f, minLen = 200, maxLen = 1000,
                      maxN = 0, rm.phix = TRUE, truncQ = 2, #maxEE = c(2, 2),
                      compress = TRUE, multithread = ncore, verbose = TRUE)


### Read quality profiles exemplarily
set.seed(74398)
startPlot <- ceiling(runif(1, 0, length(rF) - 3))

plot.rF <- plotQualityProfile(rF[startPlot:(startPlot + 2)])
plot.rF.f <- plotQualityProfile(rF.f[startPlot:(startPlot + 2)])

grid.arrange(plot.rF, plot.rF.f, nrow = 2)


### Estimate and plot the error rates
errF <- learnErrors(rF.f, multithread = ncore, verbose = TRUE)

plotErrors(errF, nominalQ = TRUE)


### Core sample inference algorithm
dadaF <- dada(rF.f, errF, HOMOPOLYMER_GAP_PENALTY = -1, BAND_SIZE = 32,
              multithread = ncore)
# dadaFs[[1]]


# saveRDS(dadaF, "dadaF_PB.rds")


### Merge paired reads
# contigs <- mergePairs(dadaF, rF.f, dadaR, rR.f, verbose = TRUE)
# head(contigs[[1]])


# saveRDS(contigs, "contigs_PB.rds")


### Construct ASV table
seqtab <- makeSequenceTable(dadaF)
dim(seqtab)


### Distribution of sequence lengths
table(nchar(getSequences(seqtab)))
plot(table(nchar(getSequences(seqtab))))


### Chimera detection
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                    multithread = ncore, verbose = TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim) / sum(seqtab)


saveRDS(seqtab.nochim, "seqtab.nochim_PB.rds")


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

saveRDS(taxa.id, "taxa.id_PB.rds")


################################################################################
################################################################################
################################################################################
################################################################################
