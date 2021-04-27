################################################################################
################################################################################
################################################################################
################################################################################
### SMP 16S MC OTU & ASV
### Author: korn@cumulonimbus.at University of Fribourg 2020
################################################################################


library("phyloseq")
library("ggplot2")
theme_set(theme_bw(base_size = 12) +
            theme(rect = element_rect(fill = "transparent")))
library("dada2")
library("DECIPHER")


################################################################################
rm(list = ls())
setwd("~/Sarracenia-Microbiome-Project/Thesis")


set.seed(34706)

source("~/Sarracenia-Microbiome-Project/Thesis/SMP_readMetadataAndASV.R")


################################################################################
### Arguments
mcDir <- "~/Sarracenia-Microbiome-Project/Thesis/"
asvDir <- "~/Seafile/SMP_results"

## Chose "16S" or "18S"
# primer <- "16S"
primer <- "18S"
glom <- TRUE


##############################################################################
### Read taxa (ASV)
setwd(asvDir)

if (primer == "16S")
  {
  seqtab.nochim <- readRDS(paste(asvDir,
                                 "SMP_16S_reseq/seqtab.nochim_SMP16S.rds",
                                 sep = "/"))
  taxaid <- readRDS(paste(asvDir,
                          "SMP_16S_reseq/taxa_SMP16S.rds",
                          sep = "/"))
  # taxaid <- readRDS("taxid_SMP16S.rds")
  }
if (primer == "18S")
  {
  seqtab.nochim <- readRDS(paste(asvDir,
                                 "SMP_18S_reseq/seqtab.nochim_SMP18S.rds",
                                 sep = "/"))
  taxaid <- readRDS("SMP_18S_reseq/taxa_SMP18S.rds")
  # taxaid <- readRDS("taxid_SMP18S.rds")
  }

mc <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),
              tax_table(taxaid))


dna <- Biostrings::DNAStringSet(taxa_names(mc))
names(dna) <- taxa_names(mc)
mc <- merge_phyloseq(mc, dna)
taxa_names(mc) <- paste0("ASV", seq(ntaxa(mc)))


### Rename taxonomic ranks
colnames(tax_table(mc)) <- c("Domain", "Phylum", "Class", "Order",
                            "Family", "Genus")


##############################################################################
### Transpose (sometimes the OTU table is transposed... :-|)
if (taxa_are_rows(mc))
  {otu_table(mc) <- t(otu_table(mc))}


##############################################################################
### Read metadata
mcMeta <- readMetadata(ps = mc, primer = primer)
mc <- merge_phyloseq(mc, mcMeta)


### Merge resequenced samples
sample_data(mc)$FullID <- gsub("r", "", sample_data(mc)$FullID)
mc <- merge_samples(mc, "FullID")


sample_data(mc) <- mcMeta


############################################################################
### Aggregate taxa by genus
if(glom == TRUE){
  mc.g  <- tax_glom(mc, taxrank = rank_names(mc)[6],
                     NArm = FALSE, bad_empty = c("", " ", "\t"))
  } else {
  mc.g <- mc
  }


##############################################################################
### Remove unwanted samples
if (primer == "16S")
  {
mc.p <- prune_samples(mc.g@sam_data$Succession == "NC Pilot" |
                        mc.g@sam_data$Succession == "NC SMP" |
                        mc.g@sam_data$Succession == "16S MC Pilot" |
                        mc.g@sam_data$Succession == "16S MC SMP", mc.g)
  }
if (primer == "18S")
  {
  mc.p <- prune_samples(mc.g@sam_data$Succession ==  "NC Pilot" |
                          mc.g@sam_data$Succession ==  "NC SMP" |
                          mc.g@sam_data$Succession ==  "18S MC Pilot" |
                          mc.g@sam_data$Succession == "18S MC SMP", mc.g)
  }


#############################################################################
### Remove empty taxa
any(taxa_sums(mc.p) == 0)
sum(taxa_sums(mc.p) == 0)
mc.s <- prune_taxa(taxa_sums(mc.p) > 0, mc.p)


## Rename levels
mc.s@sam_data$Leaf <- as.factor(mc.s@sam_data$Leaf)
levels(mc.s@sam_data$Leaf) <- c("100 %", "20 %", "NC")

names(mc.s@sam_data)[8] <- "Control"


## log-transform
mc.l <- transform_sample_counts(mc.s, function(otu) {log1p(otu)})


plot_bar(mc.l, x = "Genus", fill = "Control") +
  facet_grid( ~ Succession) +
  geom_bar(aes(color = Control, fill = Control), stat = "identity",
           position = "stack") +
  scale_fill_manual(values = c("#D55E00", "#E69F00", "#56B4E9")) +
  scale_color_manual(values = c("#D55E00", "#E69F00", "#56B4E9")) +
  theme(legend.position = "top") +
  xlab("") +
  ylab(expression(paste(ln[e], "(Abundance + 1)"))) +
  coord_flip()
# ggsave("SMP_16SMCNC.pdf", width = 11.69, height = 6)
# ggsave("SMP_18SMCNC.pdf", width = 11.69, height = 4)


############################################################################
### Convert abundance to relative abundance
mc.r <- transform_sample_counts(mc.s, function(otu) {otu / sum(otu)})


############################################################################
### Abundance filtering
## 16S
if(primer == "16S")
  {
  mc.a <- filter_taxa(mc.r, function(otu) {mean(otu) > 0.01}, prune = TRUE)
  } else {
  mc.a <- mc.r
  }

mc.a
plot(rowSums(otu_table(mc.a)))


## Remove negative control
mc.a <- prune_samples(mc.a@sam_data$Control != "NC", mc.a)


plot_bar(mc.a, x = "Genus", fill = "Control") +
  facet_grid(Control ~ Succession) +
  geom_bar(aes(color = Control, fill = Control), stat = "identity",
           position = "stack") +
  scale_fill_manual(values = c("#D55E00", "#E69F00")) +
  scale_color_manual(values = c("#D55E00", "#E69F00")) +
  theme(legend.position = "none") +
  xlab("") +
  ylab("Abundance [%]") +
  coord_flip()
# ggsave("SMP_16SMC.pdf", width = 11.69, height = 6)
# ggsave("SMP_18SMC.pdf", width = 11.69, height = 4)


################################################################################
################################################################################
################################################################################
################################################################################
