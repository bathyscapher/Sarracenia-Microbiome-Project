################################################################################
################################################################################
### SMP ASV Korn data
################################################################################

rm(list = ls())


library("phyloseq")
library("ggplot2")
theme_set(theme_bw(base_size = 20) +
            theme(rect = element_rect(fill = "transparent")))


setwd("/scratch/pSMP/Korn/")


################################################################################
### ASV
seqtab.nochim <- readRDS("seqtab.nochim_RK.rds")
tax.id <- readRDS("taxa.id_RK.rds")


rk <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),
               tax_table(tax.id))


dna <- Biostrings::DNAStringSet(taxa_names(rk))
names(dna) <- taxa_names(rk)
rk <- merge_phyloseq(rk, dna)
taxa_names(rk) <- paste0("ASV", seq(ntaxa(rk)))

head(tax_table(rk))

colnames(tax_table(rk)) <- c("Domain", "Phylum", "Class", "Order", "Family",
                              "Genus", "Species")
rk

rkMeta <- read.table("~/Sarracenia-Microbiome-Project/Transcontinental/Korn/SMP_Pool_rk.csv",
                     sep = "\t", header = TRUE)

rkMeta <- sample_data(rkMeta)
row.names(rkMeta) <- rkMeta$Sample

rk <- merge_phyloseq(rk, rkMeta)
# sample_data(rk)
sample_variables(rk)


################################################################################
### Plot overview

plot_bar(rk, x = "Phylum", fill = "Dataset") +
  # facet_wrap( ~ Class, scales = "free", nrow = 1) +
  geom_bar(aes(color = Dataset, fill = Dataset),
           stat = "identity", position = "stack") +
  theme(legend.position = "top") +
  ylab("Abundance [%]") +
  coord_flip()
ggsave("RK_18S_Phyla.pdf", width = 20, height = 20, bg = "transparent")


################################################################################
################################################################################
