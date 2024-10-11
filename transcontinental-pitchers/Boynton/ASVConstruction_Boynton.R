# SMP ASV Boynton data #########################################################

rm(list = ls())
gc()


# library("rSFFreader")
library("dada2")
library("ShortRead")
library("Biostrings")
library("DECIPHER")
library("gridExtra")


# primer <- "geo"
primer <- "temp"

if (primer == "geo") {
  wd <- "Boynton/geographic"
} else if (primer == "temp") {
  wd <- "Boynton/temporal"
}


ncore <- 6 # specify number of available cores


## SFF read lenght #############################################################
# sff.geo <- readSffGeometry("geographic.sff")
# sff.geo$nReads
#
# sff.temp <- readSffGeometry("temporal.sff")
# sff.temp$nReads
#
# par(mfrow = c(2,1))
# plot(sff.geo$Read_Widths)
# abline(h = mean(sff.geo$Read_Widths), col = "red")
#
# plot(sff.temp$Read_Widths)
# abline(h = mean(sff.temp$Read_Widths), col = "red")
# dev.off()


## Read fastq ##################################################################
list.files(path = wd,
           pattern = "fastq.gz")

rF <- sort(list.files(path = wd, pattern = "gz", full.names = TRUE))


## Extract sample names
sample.names <- sapply(strsplit(basename(rF), "\\."), `[`, 1)


## Check for primers ###########################################################
## Remove Ns from reads
rF.fN <- file.path(paste0(wd, "/filtN"), basename(rF))
filterAndTrim(rF, rF.fN, maxN = 0, multithread = TRUE)


## Forward and reverse primer
FWD <- "CTTGGTCATTTAGAGGAAGTAA" # ITS1F
REV <- "TCCTCCGCTTATTGATATGC" # ITS4


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

(reads <- ceiling(runif(1, 1, length(rF.fN))))
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = rF.fN[[reads]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = rF.fN[[reads]]))


## Cut the primers
cutadapt <- "/usr/bin/cutadapt"
# system2(cutadapt, args = "--version")


path.cut <- file.path(wd, "cutPrimers")
if(!dir.exists(path.cut)) {
  dir.create(path.cut)
}
rF.cut <- file.path(path.cut, basename(rF.fN))


## Trim FWD and the reverse-complement of REV off of R1 &  REV and the
## reverse-complement of FWD off of R2
FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

R1.flags <- paste("-g", FWD, "-a", REV.RC)
R2.flags <- paste("-g", REV, "-a", FWD.RC)
# R2.flags <- paste("-G", REV, "-A", FWD.RC)


## Run Cutadapt
for(i in seq_along(rF)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2,
                             "-o", rF.cut[i], # output
                             rF.fN[i])) # input
}


## Compare reads before and after
(reads <- ceiling(runif(1, 1, length(rF.cut))))
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = rF.fN[[reads]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = rF.fN[[reads]]))
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = rF.cut[[reads]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = rF.cut[[reads]]))


### Filter and trim. Place filtered files in filtered/ subdirectory
rF.cut.f <- file.path(paste0(wd, "/filtered"),
                      paste0(sample.names, "_filt.fastq.gz"))

names(rF.cut.f) <- sample.names


## Filtering, min and max values as in Boynton et al. 2019 #####################
filt <- filterAndTrim(rF.cut, rF.cut.f, minLen = 200, maxLen = 1000,
                      maxN = 0, rm.phix = TRUE, truncQ = 2,
                      compress = TRUE, multithread = ncore, verbose = TRUE)


### Read quality profiles exemplarily
set.seed(74398)
startPlot <- ceiling(runif(1, 0, length(rF) - 3))

plot.rF <- plotQualityProfile(rF[startPlot:(startPlot + 2)])
plot.rF.cut.f <- plotQualityProfile(rF.cut.f[startPlot:(startPlot + 2)])

grid.arrange(plot.rF, plot.rF.cut.f, nrow = 2)


## Estimate and plot the error rates ###########################################
errF <- learnErrors(rF.cut.f,
                    multithread = ncore,
                    verbose = TRUE)

plotErrors(errF, nominalQ = TRUE)


## Dereplicate identical reads #################################################
derep.rF <- derepFastq(rF.cut.f,
                       verbose = TRUE)
names(derep.rF) <- sample.names


## Core sample inference algorithm #############################################
dadaF <- dada(derep.rF,
              errF,
              HOMOPOLYMER_GAP_PENALTY = -1,
              BAND_SIZE = 32,
              multithread = ncore)
# dadaF[[1]]


saveRDS(dadaF, paste0(wd, "/dadaF_PB.rds"))


## Construct ASV table #########################################################
seqtab <- makeSequenceTable(dadaF)
dim(seqtab)


### Distribution of sequence lengths
table(nchar(getSequences(seqtab)))
plot(table(nchar(getSequences(seqtab))))


## Chimera detection ###########################################################
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                    multithread = ncore, verbose = TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim) / sum(seqtab)


saveRDS(seqtab.nochim, paste0(wd, "/seqtab.nochim_PB.rds"))


## Track reads through the pipeline ############################################
getN <- function(x) {
  sum(getUniques(x))
}

track <- cbind(filt, sapply(dadaF, getN), rowSums(seqtab.nochim))

colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track) <- sample.names
head(track)


## Assign taxonomy #############################################################
### RDP classifier #############################################################
# seqtab.nochim <- readRDS(paste0(wd, "/seqtab.nochim_PB.rds"))

taxa <- assignTaxonomy(seqtab.nochim,
                       "refs/sh_general_release_dynamic_25.07.2023.fasta",
                       multithread = TRUE, verbose = TRUE)


saveRDS(taxa, paste0(wd, "/taxa_PB.rds"))


### IdTaxa and SILVA ###########################################################
dna <- DNAStringSet(getSequences(seqtab.nochim))

load("refs/UNITE_v2023_July2023.RData")


ids <- IdTaxa(dna, trainingSet, strand = "both", processors = ncore,
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

saveRDS(taxa.id, paste0(wd, "/taxa.id_PB.rds"))
