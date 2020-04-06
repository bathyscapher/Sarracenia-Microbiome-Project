################################################################################
################################################################################
### SMP ASV Boynton data
################################################################################

rm(list = ls())


library("phyloseq")
library("ggplot2")
theme_set(theme_bw(base_size = 20))


# setwd("/scratch/pSMP/Boynton/geo/")
setwd("/scratch/pSMP/Boynton/temp/")


################################################################################
### ASV
seqtab.nochim <- readRDS("seqtab.nochim_PB.rds")
taxa.id <- readRDS("taxa.id_PB.rds")


pb <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),
               tax_table(taxa.id))


dna <- Biostrings::DNAStringSet(taxa_names(pb))
names(dna) <- taxa_names(pb)
pb <- merge_phyloseq(pb, dna)
taxa_names(pb) <- paste0("ASV", seq(ntaxa(pb)))

head(tax_table(pb))

colnames(tax_table(pb)) <- c("Domain", "Phylum", "Class", "Order", "Family",
                              "Genus", "Species")
pb

pbMeta <- read.table("~/Sarracenia-Microbiome-Project/Transcontinental/Boynton/SMP_Pool_pb.csv",
                     sep = "\t", header = TRUE)

pbMeta <- sample_data(pbMeta)
row.names(pbMeta) <- pbMeta$Sample

pb <- merge_phyloseq(pb, pbMeta)
sample_data(pb)
sample_variables(pb)


################################################################################
### Plot overview

plot_bar(pb, x = "Domain", fill = "Country") +
  # facet_wrap( ~ Class, scales = "free", nrow = 1) +
  geom_bar(aes(color = Dataset, fill = Dataset),
           stat = "identity", position = "stack") +
  theme(legend.position = "top") +
  ylab("Abundance [%]") +
  coord_flip()


################################################################################
################################################################################
