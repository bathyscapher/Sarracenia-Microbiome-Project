################################################################################
################################################################################
################################################################################
################################################################################
### SMP Read taxa (ASV), add metadata and much more
### Authors: Rachel Korn
### korn@cumulonimbus.at University of Fribourg 2019/2020
################################################################################


library("phyloseq")
rm(list = ls())


################################################################################
### Arguments

smpDir <- "~/Sarracenia-Microbiome-Project/Thesis/"
asvDir <- "~/Seafile/SMP_results"


source("~/Sarracenia-Microbiome-Project/Thesis/SMP_readMetadataAndASV.R")


################################################################################
### Read pitcher taxa
### Prokaryotes
prok.l <- readTaxa(primer = "16S", mosses = "without", glom = TRUE)
prok.l


list2env(prok.l, globalenv())

prok <- smp
prok.a <- smp.a
prok.s <- smp.s
prok.p <- smp.p
rm(smp, smp.a, smp.s, smp.p, prok.l)


### Eukaryotes
euk.l <- readTaxa(primer = "18S", mosses = "without", glom = TRUE)
euk.l

list2env(euk.l, globalenv())


euk <- smp
euk.a <- smp.a
euk.s <- smp.s
euk.p <- smp.p
rm(smp, smp.a, smp.s, smp.p, euk.l)


## Save
saveRDS(prok.a, "rds/SMP_prok.a.RDS")
saveRDS(euk.a, "rds/SMP_euk.a.RDS")


### Merge both
## Renumber taxa
taxa_names(prok.a) <- paste0("ASV", seq(ntaxa(prok.a)))
taxa_names(prok.a)


taxa_names(euk.a) <- paste0("ASV", seq(ntaxa(prok.a) + 1,
                                       ntaxa(prok.a) + ntaxa(euk.a)))
taxa_names(euk.a)


smp <- merge_phyloseq(euk.a, prok.a)
tax_table(smp)
smp


### Save
saveRDS(smp, "rds/SMP_smp.RDS")


################################################################################
### Read moss taxa summarized by sector replicate
### Prokaryotes
moss.prok.l <- readTaxa(primer = "16S", mosses = "only", merge.mosses = TRUE,
                        glom = TRUE)
moss.prok.l


list2env(moss.prok.l, globalenv())


moss.prok <- moss
moss.prok.a <- moss.a
moss.prok.s <- moss.s
rm(moss.prok.l, moss, moss.a, moss.s)


### Eukaryotes
moss.euk.l <- readTaxa(primer = "18S", mosses = "only", merge.mosses = TRUE,
                       glom = TRUE)
moss.euk.l


list2env(moss.euk.l, globalenv())


moss.euk <- moss
moss.euk.a <- moss.a
moss.euk.s <- moss.s
rm(moss.euk.l, moss, moss.a, moss.s)


## Save
saveRDS(moss.prok.a, "rds/SMP_moss20.prok.a.RDS")
saveRDS(moss.euk.a, "rds/SMP_moss20.euk.a.RDS")


### Merge both
## Renumber taxa
taxa_names(moss.prok.a) <- paste0("ASV", seq(ntaxa(moss.prok.a)))
taxa_names(moss.prok.a)


taxa_names(moss.euk.a) <- paste0("ASV", seq(ntaxa(moss.prok.a) + 1,
                                            ntaxa(moss.prok.a) +
                                                    ntaxa(moss.euk.a)))
taxa_names(moss.euk.a)


moss <- merge_phyloseq(moss.prok.a, moss.euk.a)
tax_table(moss)
moss


## Save
saveRDS(moss, "rds/SMP_moss20.RDS")


################################################################################
### Read moss taxa separately by sector replicate
### Prokaryotes
moss.prok.l <- readTaxa(primer = "16S", mosses = "only", merge.mosses = FALSE,
                        glom = TRUE)
moss.prok.l


list2env(moss.prok.l, globalenv())


moss.prok <- moss
moss.prok.a <- moss.a
moss.prok.s <- moss.s
rm(moss.prok.l, moss, moss.a, moss.s)


### Eukaryotes
moss.euk.l <- readTaxa(primer = "18S", mosses = "only", merge.mosses = FALSE,
                       glom = TRUE)
moss.euk.l


list2env(moss.euk.l, globalenv())


moss.euk <- moss
moss.euk.a <- moss.a
moss.euk.s <- moss.s
rm(moss.euk.l, moss, moss.a, moss.s)

## Save
saveRDS(moss.prok.a, "rds/SMP_moss40.prok.a.RDS")
saveRDS(moss.euk.a, "rds/SMP_moss40.euk.a.RDS")


### Merge both
## Renumber taxa
taxa_names(moss.prok.a) <- paste0("ASV", seq(ntaxa(moss.prok.a)))
taxa_names(moss.prok.a)


taxa_names(moss.euk.a) <- paste0("ASV", seq(ntaxa(moss.prok.a) + 1,
                                            ntaxa(moss.prok.a) +
                                                    ntaxa(moss.euk.a)))
taxa_names(moss.euk.a)


moss <- merge_phyloseq(moss.prok.a, moss.euk.a)
tax_table(moss)
moss


## Save
saveRDS(moss, "rds/SMP_moss40.RDS")


################################################################################
################################################################################
################################################################################
################################################################################
