################################################################################
################################################################################
################################################################################
################################################################################
### SMP Read taxa (OTU or ASV), add metadata and much more
### Authors: Rachel Korn
### korn@cumulonimbus.at University of Fribourg 2019/2020
################################################################################


library("phyloseq")
library("ggplot2")
library("GGally")
theme_set(theme_bw(base_size = 20) +
            theme(rect = element_rect(fill = "transparent")))
library("reshape2")


rm(list = ls())
setwd("~/Sarracenia-Microbiome-Project/Thesis")


set.seed(34706)


################################################################################
## Colorblind temperature scale for the 5 sites
gradCol <- c("#D55E00", "#E69F00", "#F0E442", "#009E73", "#56B4E9")


################################################################################
### Pitcher data
prok.a <- readRDS("rds/SMP_prok.a.RDS")
euk.a <- readRDS("rds/SMP_euk.a.RDS")


################################################################################
### Moss taxa
## Merged
moss.prok.a <- readRDS("rds/SMP_moss20.prok.a.RDS")
moss.euk.a <- readRDS("rds/SMP_moss20.euk.a.RDS")


################################################################################
### nMDS pitchers
## Choose pro- or eukaryotes
# primer <- "16S"
primer <- "18S"


if(primer == "16S")
  {smp <- prok.a}
if(primer == "18S")
  {smp <- euk.a}


## log-transform
smp.log <- transform_sample_counts(smp, function(otu) {log1p(otu)})


## nMDS
smp.nmds <- ordinate(smp.log, method = "NMDS", distance = "bray", k = 4,
                     autotransform = FALSE, trymax = 100)
smp.nmds
stressplot(smp.nmds)


## Plot nMDS
plot_ordination(smp, smp.nmds, shape = "Succession", color = "Site",
                title = NULL, axes = 1:2) +
  stat_ellipse(aes(group = Site), type = "t", linetype = 2, size = 0.2) +
  geom_point(size = 3) +
  scale_colour_manual(values = gradCol) +
  coord_fixed(ratio = 1) +
  theme(legend.position = "top", legend.direction = "horizontal",
        legend.box = "vertical") +
  xlab("nMDS1") +
  ylab("nMDS2")
# ggsave("SMP_16S_nMDS_ASV.pdf", width = 8.27, height = 8.27)
# ggsave("SMP_18S_nMDS_ASV.pdf", width = 8.27, height = 8.27)


### nMDS pitchers and mosses
## Choose pro- or eukaryotes
# primer <- "16S"
primer <- "18S"

if(primer == "16S")
  {smp.moss <- merge_phyloseq(prok.a, moss.prok.a)}
if(primer == "18S")
  {smp.moss <- merge_phyloseq(euk.a, moss.euk.a)}


## log-transform
smp.moss.log <- transform_sample_counts(smp.moss, function(otu) {log1p(otu)})


## nMDS
smp.moss.nmds <- ordinate(smp.moss.log, method = "NMDS", distance = "bray",
                          k = 4, autotransform = FALSE, trymax = 100,
                          weakties = FALSE)
smp.moss.nmds
stressplot(smp.moss.nmds)


## Plot nMDS
plot_ordination(smp.moss, smp.moss.nmds, shape = "Site", color = "Succession",
                title = NULL, axes = 1:2) +
  stat_ellipse(aes(group = Succession), type = "t", linetype = 2, size = 0.2) +
  geom_point(size = 3) +
  scale_colour_manual(values = gradCol[c(2, 1, 4)]) +
  coord_fixed(ratio = 1) +
  theme(legend.position = "top", legend.direction = "horizontal",
        legend.box = "vertical") +
  xlab("nMDS1") +
  ylab("nMDS2")
# ggsave("SMP_16S_nMDS_withMoss.pdf", width = 8.27, height = 8.27)
# ggsave("SMP_18S_nMDS_withMoss.pdf", width = 8.27, height = 8.27)


################################################################################
### Plot taxa
## Load pro- and eukaryotic taxa together
smp <- readRDS("rds/SMP_smp.RDS")


### Dataframe with taxa and abundances
## Pitchers
tax.l <- as.data.frame(smp@tax_table@.Data)
otu.l <- t(as.data.frame(otu_table(smp)))

tax.otu.l <- merge(tax.l, otu.l, by = 0, all = TRUE)
rownames(tax.otu.l) <- tax.otu.l$Row.names
colnames(tax.otu.l)[1] <- "Taxa"

rm(tax.l, otu.l)


### log-transform abundances
tax.otu.l[, -c(1:7)] <- log1p(tax.otu.l[, -c(1:7)])


tax.otu.lm <- reshape2::melt(tax.otu.l, id.vars = c("Taxa", "Domain",
                                                      "Phylum", "Class",
                                                      "Order", "Family",
                                                      "Genus"),
                               variable.name = "FullID", value.name = "Count")
tax.otu.lm$Site <- substr(tax.otu.lm$FullID, 1, 2)
tax.otu.lm$Site <- ordered(tax.otu.lm$Site,
                           levels = c("CB", "LT", "LE", "LV", "LM"))
tax.otu.lm$Succession <- as.factor(substr(tax.otu.lm$FullID, 7, 7))
levels(tax.otu.lm$Succession) <- c("Early", "Late")
tax.otu.lm$Domain[tax.otu.lm$Domain == "Archaea"] <- "A."


## Reorder by taxonomy
tax.otu.lm[, c(2:7)] <- lapply(tax.otu.lm[, c(2:7)], as.factor)


tax <- c("Crenarchaeota", ## Archaea
         ## Proteobacteria
         "Bdellovibrionota", "Campilobacterota", "Desulfobacterota",
         "Proteobacteria",
         "Myxococcota", ## Deltaproteobacteria
         "Planctomycetota", "Verrucomicrobiota", ## PVC group
         "Gemmatimonadota", ## FCB group
         ## Singles
         "Acidobacteriota", "Bacteroidota", "Chloroflexi", "Dependentiae",
         ## Terrabacteria
         "Actinobacteriota", "Armatimonadota", "Cyanobacteria", "Firmicutes",
         "WPS-2", ## ?
         "Amoebozoa_ph", "Schizoplasmodiida", ## Amoeba
         ## SAR:Heterokonta
         "Diatomea", "Bicosoecida", "Ochrophyta_ph", "Peronosporomycetes",
         "Cercozoa", "Retaria", ## SAR: Rhizaria
         ## Harosa:Alveolata
         "Ciliophora", "Dinoflagellata", "Protalveolata",
         "Euglenozoa", "Heterolobosea", ## Excavata
         "Chlorophyta_ph", "Cryptophyceae_ph", ## Algae
         "Nucleariidae_and_Fonticula_group", ## Holomycota
         ## Fungi
         "Ascomycota", "Basidiomycota", "Blastocladiomycota", "Chytridiomycota",
         "Cryptomycota", "Zoopagomycota",
         "Holozoa_ph", ## Choanozoa
         ## Metazoa
         "Annelida", "Platyhelminthes", "Nematozoa", "Rotifera", "Tardigrada",
         "Arthropoda")


## Pronounce small values even more
tax.otu.lm$Count[tax.otu.lm$Count < 0.1 & tax.otu.lm$Count > 0] <- 0.1

summary(tax.otu.lm$Count)


## Plot
ggplot(tax.otu.lm, aes(x = factor(Phylum, level = rev(tax)), y = Count,
                       fill = Site)) +
  geom_bar(stat = "identity") +
  facet_grid(Domain ~ Succession, scales = "free_y", space = "free_y") +
  scale_fill_manual(values = gradCol) +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab(expression(paste(ln[e], "(Relative abundance + 1)"))) +
  xlab("") +
  coord_flip() +
  theme_bw(base_size = 15)
# ggsave("SMP_PhylaOverview.pdf", width = 8.27, height = 11.69)


# aggregate(Count ~ Class, data = tax.otu.lm, FUN = sum)


### Inspect taxa and frequencies
df <- tax.otu.l[tax.otu.l$Phylum == "Basidiomycota", ]
summary(colSums(df[, -c(1:7)]) == 0)
plot(colSums(df[, -c(1:7)]))

100 / 160 * 5


### Plot incidence
## Presence/absence
tax.otu.pa <- tax.otu.l
tax.otu.pa[, -c(1:7)][tax.otu.pa[, -c(1:7)] > 0] <- 1



## Sum by Phylum and replace with 1
inc <- aggregate(tax.otu.pa[, -c(1:7)],
                 by = list(tax.otu.pa$Domain, tax.otu.pa$Phylum), FUN = sum)
names(inc)[1:2] <- c("Domain", "Phylum")

inc[, -c(1:2)][inc[, -c(1:2)] > 0] <- 1
inc$Incidence <- rowSums(inc[, -c(1:2)])


## Sort
inc <- inc[order(inc$Incidence), ]
inc$Phylum <- factor(inc$Phylum, levels = unique(as.character(inc$Phylum)))
inc$IncidencePercent <- 100 / 160 * inc$Incidence


## Abbreviate archaea
inc$Domain[inc$Domain == "Archaea"] <- "A."


## Plot
ggplot(inc, aes(x = Phylum, y = IncidencePercent)) +
  geom_point() +
  facet_grid(Domain ~ ., scales = "free", space = "free_y") +
  theme(legend.position = "none") +
  ylab("Incidence [%]") +
  xlab("") +
  coord_flip() +
  theme_bw(base_size = 10)
# ggsave("SMP_Incidence.pdf", width = 8.27, height = 11.69 / 2)


################################################################################
################################################################################
################################################################################
################################################################################
