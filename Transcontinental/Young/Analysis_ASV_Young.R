################################################################################
################################################################################
### SMP ASV Young data
################################################################################


library("phyloseq")
library("ggplot2")
theme_set(theme_bw(base_size = 20))


rm(list = ls())

setwd("/scratch/pSMP/Young/prok/")
# setwd("/scratch/pSMP/Young/euk/")


################################################################################
### ASV
seqtab.nochim <- readRDS("seqtab.nochim_EY.rds")
tax.id <- readRDS("taxa_EY.rds")


ey <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),
               tax_table(tax.id))


dna <- Biostrings::DNAStringSet(taxa_names(ey))
names(dna) <- taxa_names(ey)
ey <- merge_phyloseq(ey, dna)
taxa_names(ey) <- paste0("ASV", seq(ntaxa(ey)))

head(tax_table(ey))


colnames(tax_table(ey)) <- c("Domain", "Phylum", "Class", "Order", "Family",
                              "Genus")
ey

eyMeta <- read.table("~/Sarracenia-Microbiome-Project/Transcontinental/Young/SMP_Pool_ey.csv",
                     sep = "\t", header = TRUE)

eyMeta <- sample_data(eyMeta)
row.names(eyMeta) <- eyMeta$ID


ey <- merge_phyloseq(ey, eyMeta)
# sample_data(ey)
sample_variables(ey)


################################################################################
### Plot overview

plot_bar(ey, x = "Phylum", fill = "Habitat") +
  # facet_wrap( ~ Class, scales = "free", nrow = 1) +
  geom_bar(aes(color = Habitat, fill = Habitat),
           stat = "identity", position = "stack") +
  theme(legend.position = "top") +
  ylab("Abundance [%]") +
  coord_flip()
# ggsave("EY_16S_Phyla.pdf", width = 20, height = 20, bg = "transparent")
# ggsave("EY_18S_Phyla.pdf", width = 20, height = 20, bg = "transparent")


################################################################################
################################################################################
