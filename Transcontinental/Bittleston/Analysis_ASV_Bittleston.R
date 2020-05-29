################################################################################
################################################################################
### SMP ASV Bittleston data
################################################################################


library("phyloseq")
library("ggplot2")
theme_set(theme_bw(base_size = 20))


rm(list = ls())


## Choose either pro- or eukaryotes
setwd("/scratch/pSMP/16S_raw_data_2015/FASTQ/")
# setwd("/scratch/pSMP/18S_raw_data_2015/FASTQ/")


################################################################################
### ASV
seqtab.nochim <- readRDS("seqtab.nochim_LB.rds")
taxa.id <- readRDS("taxa_LB.rds")


lb <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),
               tax_table(taxa.id))


dna <- Biostrings::DNAStringSet(taxa_names(lb))
names(dna) <- taxa_names(lb)
lb <- merge_phyloseq(lb, dna)
taxa_names(lb) <- paste0("ASV", seq(ntaxa(lb)))

head(tax_table(lb))

colnames(tax_table(lb)) <- c("Domain", "Phylum", "Class", "Order", "Family",
                              "Genus")
lb

lbMeta <- read.table("~/Sarracenia-Microbiome-Project/Transcontinental/Bittleston/SMP_Pool_lb.csv",
                     sep = "\t", header = TRUE)

lbMeta <- sample_data(lbMeta)
row.names(lbMeta) <- lbMeta$ID

lb <- merge_phyloseq(lb, lbMeta)
# sample_data(lb)
sample_variables(lb)


################################################################################
## Only Sarracenia and related samples
lb <- prune_samples(lb@sam_data$Country !=  "Malaysia" &
                      lb@sam_data$Country !=  "Singapore", lb)
lb@sam_data$Country


## Remove empty taxa
any(taxa_sums(lb) == 0)
sum(taxa_sums(lb) == 0)
lb <- prune_taxa(taxa_sums(lb) > 0, lb)
lb


################################################################################
### Plot overview

plot_bar(lb, x = "Phylum", fill = "Habitat") +
  # facet_wrap( ~ Class, scales = "free", nrow = 1) +
  geom_bar(aes(color = Habitat, fill = Habitat),
           stat = "identity", position = "stack") +
  theme(legend.position = "top") +
  ylab("Abundance") +
  coord_flip()
ggsave("LB_16S_Phyla.pdf", width = 20, height = 20)
# ggsave("LB_18S_Phyla.pdf", width = 20, height = 20)


################################################################################
################################################################################
