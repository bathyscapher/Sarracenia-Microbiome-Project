################################################################################
################################################################################
### SMP ASV Young greenhouse data
################################################################################


library("phyloseq")
library("ggplot2")
theme_set(theme_bw(base_size = 20))


rm(list = ls())


setwd("/scratch/pSMP/Young_Greenhouse/")


################################################################################
### ASV
seqtab.nochim <- readRDS("seqtab.nochim_EY_GH.rds")
tax.id <- readRDS("taxa_EY_GH.rds")


eb <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),
               tax_table(tax.id))


dna <- Biostrings::DNAStringSet(taxa_names(eb))
names(dna) <- taxa_names(eb)
eb <- merge_phyloseq(eb, dna)
taxa_names(eb) <- paste0("ASV", seq(ntaxa(eb)))

head(tax_table(eb))

colnames(tax_table(eb)) <- c("Domain", "Phylum", "Class", "Order", "Family",
                              "Genus")
eb


ebMeta <- read.table("~/Sarracenia-Microbiome-Project/Transcontinental/Young_GH/SMP_Pool_ey_gh.csv",
                     sep = "\t", header = TRUE)

ebMeta <- sample_data(ebMeta)
row.names(ebMeta) <- ebMeta$ID

eb <- merge_phyloseq(eb, ebMeta)
# sample_data(eb)
sample_variables(eb)


################################################################################
### Plot overview

plot_bar(eb, x = "Phylum", fill = "Habitat") +
  # facet_wrap( ~ Class, scales = "free", nrow = 1) +
  geom_bar(aes(color = Habitat, fill = Habitat),
           stat = "identity", position = "stack") +
  theme(legend.position = "top") +
  ylab("Abundance [%]") +
  coord_flip()
ggsave("EY_GH_Phyla.pdf", width = 20, height = 20, bg = "transparent")


################################################################################
################################################################################
