################################################################################
################################################################################
### SMP ASV Korn data
################################################################################


library("phyloseq")
library("ggplot2")
theme_set(theme_bw(base_size = 20) +
            theme(rect = element_rect(fill = "transparent")))


rm(list = ls())
gc()


# primer <- "16S"
primer <- "18S"

if (primer == "16S") {
  wd <- "Korn/16S"
} else if (primer == "18S") {
  wd <- "Korn/18S"
}


################################################################################
### ASV
seqtab.nochim <- readRDS(paste0(wd, "/seqtab.nochim_RK_", primer, ".rds"))
taxa.id <- readRDS(paste0(wd, "/taxa_RK_", primer, ".rds"))


rk <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),
               tax_table(taxa.id))


dna <- Biostrings::DNAStringSet(taxa_names(rk))
names(dna) <- taxa_names(rk)
rk <- merge_phyloseq(rk, dna)
taxa_names(rk) <- paste0("ASV", seq(ntaxa(rk)))

head(tax_table(rk))

colnames(tax_table(rk)) <- c("Domain", "Phylum", "Class", "Order", "Family",
                              "Genus")
rk

sample_names(rk)


rkMeta <- read.table("Korn/SMP_Pool_rk.csv",
                     sep = "\t", header = TRUE)

rkMeta <- sample_data(rkMeta)
row.names(rkMeta) <- rkMeta$ID


rk <- merge_phyloseq(rk, rkMeta)
sample_data(rk)
sample_variables(rk)
rk


################################################################################
### Plot overview

plot_bar(rk, x = "Phylum", fill = "Habitat") +
  # facet_wrap( ~ Class, scales = "free", nrow = 1) +
  geom_bar(aes(color = Habitat, fill = Habitat),
           stat = "identity", position = "stack") +
  theme(legend.position = "top") +
  ylab("Abundance [%]") +
  coord_flip()
ggsave("RK_16S_Phyla.pdf", width = 20, height = 20, bg = "transparent")
# ggsave("RK_18S_Phyla.pdf", width = 20, height = 20, bg = "transparent")


################################################################################
################################################################################
