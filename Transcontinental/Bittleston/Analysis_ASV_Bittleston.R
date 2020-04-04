################################################################################
################################################################################
################################################################################
################################################################################
### SMP ASV Bittleston data
### Author: korn@cumulonimbus.at University of Fribourg 2020
################################################################################

rm(list = ls())


library("phyloseq")
library("ggplot2")
theme_set(theme_bw(base_size = 20))


setwd("/scratch/pSMP/16S_raw_data_2015/FASTQ/")
# setwd("/scratch/pSMP/18S_raw_data_2015/FASTQ/")


################################################################################
### ASV
seqtab.nochim <- readRDS("seqtab.nochim_LB.rds")
tax.id <- readRDS("taxa.id_LB.rds")


eb <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),
               tax_table(tax.id))


dna <- Biostrings::DNAStringSet(taxa_names(eb))
names(dna) <- taxa_names(eb)
eb <- merge_phyloseq(eb, dna)
taxa_names(eb) <- paste0("ASV", seq(ntaxa(eb)))

head(tax_table(eb))

colnames(tax_table(eb)) <- c("Domain", "Phylum", "Class", "Order", "Family",
                              "Genus", "Species")
eb

ebMeta <- read.table("~/Sarracenia-Microbiome-Project/Transcontinental/Bittleston/SMP_Pool_eb.csv",
                     sep = "\t", header = TRUE)

ebMeta <- sample_data(ebMeta)
row.names(ebMeta) <- ebMeta$Sample

eb <- merge_phyloseq(eb, ebMeta)
# sample_data(eb)
sample_variables(eb)


################################################################################
## Only Sarracenia and related samples
eb <- prune_samples(eb@sam_data$Country !=  "Malaysia" &
                      eb@sam_data$Country !=  "Singapore"
                       , eb)
eb@sam_data$Country

## Remove empty taxa
any(taxa_sums(eb) == 0)
sum(taxa_sums(eb) == 0)
eb <- prune_taxa(taxa_sums(eb) > 0, eb)
eb


################################################################################
### Plot overview

plot_bar(eb, x = "Phylum", fill = "Dataset") +
  # facet_wrap( ~ Class, scales = "free", nrow = 1) +
  geom_bar(aes(color = Dataset, fill = Dataset),
           stat = "identity", position = "stack") +
  theme(legend.position = "top") +
  ylab("Abundance") +
  coord_flip()
# ggsave("LB_16S_Phyla.pdf", width = 20, height = 20)
ggsave("LB_18S_Phyla.pdf", width = 20, height = 20)


################################################################################
################################################################################
################################################################################
################################################################################
