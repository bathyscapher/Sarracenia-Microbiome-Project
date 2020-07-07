################################################################################
################################################################################
### SMP ASV data
################################################################################


library("phyloseq")
library("ggplot2")
theme_set(theme_bw(base_size = 12))
library("reshape2")


rm(list = ls())


################################################################################
### Read metadata
psMeta.df <- read.table("~/Sarracenia-Microbiome-Project/Transcontinental/SMP_Pool.csv",
                        sep = "\t", header = TRUE)
psMeta <- sample_data(psMeta.df)
row.names(psMeta) <- psMeta$ID


### List all working directories and abbreviations
wd <- c("/scratch/pSMP/16S_raw_data_2015/FASTQ/",
        "/scratch/pSMP/18S_raw_data_2015/FASTQ/",
        "/scratch/pSMP/Boynton/geo/",
        "/scratch/pSMP/Boynton/temp/",
        "/scratch/pSMP/Korn/SMP_16S_reseq/",
        "/scratch/pSMP/Korn/SMP_18S_reseq/",
        "/scratch/pSMP/Young/prok/",
        "/scratch/pSMP/Young/euk/",
        "/scratch/pSMP/Young_Greenhouse/")

ps.l <- c("eb.p", "eb.e", "pb.g", "pb.t", "rk.p", "rk.e", "ey.p", "ey.e",
          "ey.g")


# i <- 5
### Read the ASV datasets
for (i in 1:length(wd)) {

  ## Set workingdirectory and grab the sequences and taxonomy file
  setwd(wd[i])
  rds <- list.files(pattern = "seqtab|taxa_")
  seqtab.nochim <- readRDS(rds[1])
  taxa.id <- readRDS(rds[2])

  ## Create phyloseq object
  ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),
                 tax_table(taxa.id))

  dna <- Biostrings::DNAStringSet(taxa_names(ps))
  names(dna) <- taxa_names(ps)
  ps <- merge_phyloseq(ps, dna)

  colnames(tax_table(ps)) <- c("Domain", "Phylum", "Class", "Order", "Family",
                               "Genus")

  ## Transfer data to specific naming and attach metadata
  assign(ps.l[i], merge_phyloseq(ps, psMeta))
  sample_data(eval(as.name(ps.l[i])))
  sample_variables(eval(as.name(ps.l[i])))

  ## Clean workspace
  rm(rds, seqtab.nochim, taxa.id, dna, ps)
}


setwd("~/Sarracenia-Microbiome-Project/Transcontinental/")


################################################################################
### Merge
## Prokaryotes
prok <- merge_phyloseq(eb.p, ey.p, ey.g, rk.p)
taxa_names(prok) <- paste0("ASV", seq(ntaxa(prok)))
head(tax_table(prok))
prok


## Eukaryotes
euk <- merge_phyloseq(eb.e, ey.e, rk.e, pb.g, pb.t)
taxa_names(euk) <- paste0("ASV", seq(ntaxa(euk)))
head(tax_table(euk))
euk


## Eukaryotes 18S (without ITS)
# euk <-  merge_phyloseq(eb.e, ey.e, rk.e)
# euk


## Clean workspace
rm(eb.p, ey.p, ey.g, rk.p, eb.e, ey.e, rk.e, pb.g, pb.t, wd, i, ps.l, psMeta)


################################################################################
### Only Sarracenia and related samples
## Prokaryotes
prok <- prune_samples(prok@sam_data$Country != "Malaysia" &
                        prok@sam_data$Country != "Singapore" #&
                        # prok@sam_data$Habitat != "Bog" &
                        # prok@sam_data$Habitat != "Soil"
                      , prok)
prok@sam_data$Country
levels(prok@sam_data$Habitat)


## Remove empty taxa
any(taxa_sums(prok) == 0)
sum(taxa_sums(prok) == 0)
prok <- prune_taxa(taxa_sums(prok) > 0, prok)
prok


## Eukaryotes
euk <- prune_samples(euk@sam_data$Country != "Malaysia" &
                       euk@sam_data$Country != "Singapore" #&
                       # euk@sam_data$Habitat != "Bog" &
                       # euk@sam_data$Habitat != "Soil"
                     , euk)
euk@sam_data$Country
levels(euk@sam_data$Habitat)


## Remove empty taxa
any(taxa_sums(euk) == 0)
sum(taxa_sums(euk) == 0)
euk <- prune_taxa(taxa_sums(euk) > 0, euk)
euk


################################################################################
### Taxonomic filtering
colnames(tax_table(prok))


## Prokaryotes
prok.taxa <- get_taxa_unique(prok, taxonomic.rank = rank_names(prok)[4])
prok.taxa
taxa.g <- grep("chloro", prok.taxa, ignore.case = TRUE)
prok.taxa[taxa.g]


prok.s <- subset_taxa(prok, !(Domain %in% c("unclassified_Root") |
                                Order %in% c("Chloroplast") |
                                Family %in% c("Mitochondria")))
prok.s


## Remove empty taxa
any(taxa_sums(prok.s) == 0)
# sum(taxa_sums(prok.s) == 0)
# prok.s <- prune_taxa(taxa_sums(prok.s) > 0, prok.s)
# prok.s

## Remove empty samples (1)
prok.s <- prune_samples(sample_sums(prok.s) > 0, prok.s)


### Percentual abundance
prok.s <- transform_sample_counts(prok.s, function(otu) {otu / sum(otu)})
plot(rowSums(otu_table(prok.s)))


### Abundance filtering
prok.a <- filter_taxa(prok.s, function(otu) {mean(otu) > 0.00001}, prune = TRUE)
prok.a
plot(rowSums(otu_table(prok.a)))


## Plot by country
prok.t <- data.frame(otu_table(prok.a))
prok.t$Sum <- rowSums(prok.t)
prok.t$ID <- rownames(prok.t)
prok.t <- merge(prok.t, psMeta.df[, 1:5])


ggplot(prok.t) +
  geom_point(aes(x = ID, y = Sum, color = Country), pch = 21) +
  theme(legend.position = "top",
                axis.text.x = element_blank()) +
  xlab("Sample") +
  ylab("Abundance [%]")
# ggsave("SMP_prok_AbundanceAfterFiltering.pdf", width = 11.69, height = 8.27)


any(taxa_sums(prok.s) == 0)


################################################################################
## Eukaryotes
euk.taxa <- get_taxa_unique(euk, taxonomic.rank = rank_names(euk)[1])
euk.taxa
taxa.g <- grep("Tremellales", euk.taxa, ignore.case = TRUE)
euk.taxa[taxa.g]


euk.s <- subset_taxa(euk, !(Domain %in% c(NA, "unclassified_Root", "Bacteria") |
                                Phylum %in% c(NA, "unclassified_Ecdysozoa",
                                              "unclassified_Dikarya",
                                              "unclassified_Fungi",
                                              # "Myxogastria",
                                              "Vertebrata", "Mollusca",
                                              "Microsporidia", "Zoopagomycota",
                                              "Mucoromycota") |
                                Class %in% c("Embryophyta",
                                             "Agaricomycetes",
                                             "unclassified_Agaricomycotina",
                                             "Pucciniomycetes",
                                             "unclassified_Pucciniomycotina",
                                             "Agaricostilbomycetes",
                                             "Exobasidiomycetes",
                                             "Spiculogloeomycetes",
                                             "Dothideomycetes", "Leotiomycetes",
                                             "Ustilaginomycetes",
                                             "unclassified_Pezizomycotina",
                                             "Eurotiomycetes",
                                             "Entomophthoromycetes",
                                             "Microbotryomycetes",
                                             "Sordariomycetes",
                                             "Malasseziomycetes", "Diplopoda",
                                             "Chilopoda", "Diplura", "Ellipura",
                                             "unclassified_Hexapoda", "Cestoda",
                                             "unclassified_Chelicerata",
                                             "Arachnida", "Trematoda",
                                             "Insecta") |
                                Order %in% c("Gregarinasina",
                                             "Trichosporonales",
                                             "Trypanosomatida",
                                             "Filobasidiales", "Tremellales") |
                                Family %in% c("Adeleorina", "Trichosporonaceae",
                                              "Mrakiaceae",
                                              "Legeriomycetaceae")))
euk.s
# get_taxa_unique(euk.s, taxonomic.rank = rank_names(euk)[2])


## Remove empty taxa
any(taxa_sums(euk.s) == 0)
# sum(taxa_sums(euk.s) == 0)
# euk.s <- prune_taxa(taxa_sums(euk.s) > 0, euk.s)
# euk.s


## Remove empty samples (> 200)
sample_names(prune_samples(sample_sums(euk.s) == 0, euk.s))
euk.s <- prune_samples(sample_sums(euk.s) > 0, euk.s)


### Percentual abundance
euk.s <- transform_sample_counts(euk.s, function(otu) {otu / sum(otu)})
plot(rowSums(otu_table(euk.s)))


### Abundance filtering
euk.a <- filter_taxa(euk.s, function(otu) {mean(otu) > 0.00001}, prune = TRUE)
euk.a
plot(rowSums(otu_table(euk.a)))


## Plot by country
euk.t <- data.frame(otu_table(euk.a))
euk.t$Sum <- rowSums(euk.t)
euk.t$ID <- rownames(euk.t)
euk.t <- merge(euk.t, psMeta.df[, 1:5])


ggplot(euk.t) +
  geom_point(aes(x = ID, y = Sum, color = Country), pch = 21) +
  theme(legend.position = "top",
        axis.text.x = element_blank()) +
  xlab("Sample") +
  ylab("Abundance [%]")
# ggsave("SMP_euk_AbundanceAfterFiltering.pdf", width = 11.69, height = 8.27)

rm(prok.t, euk.t)

any(taxa_sums(euk.a) == 0)


################################################################################
### Plot overview
plot_bar(prok.a, x = "Phylum", fill = "Habitat") +
  facet_wrap( ~ Country) +
  geom_bar(aes(color = Habitat, fill = Habitat),
           stat = "identity", position = "stack") +
  theme(legend.position = "top") +
  xlab("") +
  ylab("Abundance [%]") +
  coord_flip()
ggsave("SMP_prok_Phyla.pdf", width = 11.69, height = 8.27)


plot_bar(euk.a, x = "Phylum", fill = "Habitat") +
  facet_wrap( ~ Country) +
  geom_bar(aes(color = Habitat, fill = Habitat),
           stat = "identity", position = "stack") +
  theme(legend.position = "top") +
  xlab("") +
  ylab("Abundance [%]") +
  coord_flip()
ggsave("SMP_euk_Phyla.pdf", width = 11.69, height = 8.27)


################################################################################
### Ordination
## Prokaryotes
prok.log <- transform_sample_counts(prok.a, function(otu) {log1p(otu)})
prok.nmds <- ordinate(prok.log, method = "PCoA", distance = "bray")
# prok.nmds <- ordinate(prok.log, method = "NMDS", distance = "bray", k = 4)
prok.nmds


plot_ordination(prok.a, prok.nmds, shape = "Country",
                # label = "Site",
                color = "Site", title = NULL) +
  # stat_ellipse(aes(group = Habitat), type = "norm", linetype = 2, size = 0.2) +
  # scale_colour_manual(values = gradCol) +
  coord_fixed(ratio = 1) +
  theme(legend.position = "top", legend.direction = "horizontal",
        legend.box = "vertical") +
  # geom_text(mapping = aes(label = ID), size = 1.5, vjust = 1.5) +
  guides(fill = guide_legend(nrow = 3, byrow = TRUE))
  # xlab("PCoA/MDS1 (95 % CI)") +
  # ylab("PCoA/MDS2 (95 % CI)")
# ggsave("SMP_prok_nMDS.pdf", width = 12, height = 8.27, bg = "transparent")
ggsave("SMP_prok_PCoA.pdf", width = 12, height = 8.27, bg = "transparent")
# ggsave("SMP_prok_PCoA_Country.pdf", width = 10, height = 8.27, bg = "transparent")
# ggsave("SMP_prok_PCoA_WithSurroundings.pdf", width = 10, height = 8.27, bg = "transparent")


### Eukaryotes
euk.log <- transform_sample_counts(euk.a, function(otu) {log1p(otu)})
euk.nmds <- ordinate(euk.log, method = "PCoA", distance = "bray")
# euk.nmds <- ordinate(euk.log, method = "NMDS", distance = "bray", k = 4)
euk.nmds


plot_ordination(euk.a, euk.nmds, shape = "Country",
                color = "Site", title = NULL) +
  stat_ellipse(aes(group = Country), type = "norm", linetype = 2, size = 0.2) +
  # scale_colour_manual(values = gradCol) +
  coord_fixed(ratio = 1) +
  theme(legend.position = "top", legend.direction = "horizontal",
        legend.box = "vertical") +
  geom_text(mapping = aes(label = ID), size = 1.5, vjust = 1.5) +
  guides(fill = guide_legend(nrow = 3, byrow = TRUE))
# xlab("PCoA/MDS1 (95 % CI)") +
# ylab("PCoA/MDS2 (95 % CI)")
ggsave("SMP_euk_PCoA.pdf", width = 10, height = 8.27, bg = "transparent")
# ggsave("SMP_euk_PCoA_Country.pdf", width = 10, height = 8.27, bg = "transparent")
# ggsave("SMP_euk_PCoA_18SOnly.pdf", width = 10, height = 8.27, bg = "transparent")
# ggsave("SMP_euk_PCoA_WithSurroundings.pdf", width = 10, height = 8.27, bg = "transparent")


################################################################################
################################################################################
