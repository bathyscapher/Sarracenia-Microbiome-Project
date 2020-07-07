################################################################################
################################################################################
### SMP ASV Young data
################################################################################


library("phyloseq")
library("ggplot2")
theme_set(theme_bw(base_size = 20))


rm(list = ls())


# setwd("/scratch/pSMP/Freedman/Samp55/")
setwd("/scratch/pSMP/Freedman/Samp1/")


################################################################################
### ASV
seqtab.nochim <- readRDS("seqtab.nochim_ZF.rds")
tax.id <- readRDS("taxa_ZF.rds")


zf <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),
               tax_table(tax.id))


dna <- Biostrings::DNAStringSet(taxa_names(zf))
names(dna) <- taxa_names(zf)
zf <- merge_phyloseq(zf, dna)
taxa_names(zf) <- paste0("ASV", seq(ntaxa(zf)))

head(tax_table(zf))


colnames(tax_table(zf)) <- c("Domain", "Phylum", "Class", "Order", "Family",
                              "Genus")
zf

zfMeta <- read.table("~/Sarracenia-Microbiome-Project/Transcontinental/Young/SMP_Pool_zf.csv",
                     sep = "\t", header = TRUE)

zfMeta <- sample_data(zfMeta)
row.names(zfMeta) <- zfMeta$ID


zf <- merge_phyloseq(zf, zfMeta)
# sample_data(zf)
sample_variables(zf)


################################################################################
### Plot overview

plot_bar(zf, x = "Phylum", fill = "Habitat") +
  # facet_wrap( ~ Class, scales = "free", nrow = 1) +
  geom_bar(aes(color = Habitat, fill = Habitat),
           stat = "identity", position = "stack") +
  theme(legend.position = "top") +
  ylab("Abundance [%]") +
  coord_flip()
# ggsave("ZF_16S_Phyla.pdf", width = 20, height = 20, bg = "transparent")
# ggsave("ZF_18S_Phyla.pdf", width = 20, height = 20, bg = "transparent")


################################################################################
################################################################################
