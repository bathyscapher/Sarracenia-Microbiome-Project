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
library("vegan")
library("ape")
library("gplots")
library("reshape2")
# library("scales")
# library("cluster")
# library("apcluster")


# setwd("~/Sarracenia-Microbiome-Project/Thesis")
setwd("~/Seafile/_FNS_2016-NicheDynamCoexist/SarraceniaMicrobiomeProject")

# rm(list = ls())


################################################################################
## Colorblind temperature scale for the 5 sites
gradCol <- c("#D55E00", "#E69F00", "#F0E442", "#009E73", "#56B4E9")


################################################################################
### Pitcher data
tst.prok.a <- readRDS("SMP_prok.a.RDS")
tst.euk.a <- readRDS("SMP_euk.a.RDS")
tst.smp <- readRDS("SMP_smp.RDS")


################################################################################
### Moss taxa
## Merged
# moss.prok.a <- readRDS("SMP_moss.prok.a.RDS")
# moss.euk.a <- readRDS("SMP_moss.euk.a.RDS")
# moss <- readRDS("SMP_moss.RDS")


## Separately
moss.prok.a <- readRDS("Mosses/SMP_moss40.prok.a.RDS")
moss.euk.a <- readRDS("Mosses/SMP_moss40.euk.a.RDS")
moss <- readRDS("Mosses/SMP_moss40.RDS")


################################################################################
### Export dataframe with taxa and abundances
## Pitchers
tax.l <- as.data.frame(smp@tax_table@.Data)
otu.l <- t(as.data.frame(otu_table(smp)))

tax.otu.l <- merge(tax.l, otu.l, by = 0, all = TRUE)
rownames(tax.otu.l) <- tax.otu.l$Row.names
colnames(tax.otu.l)[1] <- "Taxa"

rm(tax.l, otu.l)

# write.table(tax.otu.l, "SMP_OTU_Leaf.csv", sep = "\t", row.names = FALSE)
# write.table(tax.otu.l, "SMP_ASV_Leaf.csv", sep = "\t", row.names = FALSE)


## Mosses
tax.m <- as.data.frame(moss@tax_table@.Data)
otu.m <- as.data.frame(t(otu_table(moss)))

tax.otu.m <- merge(tax.m, otu.m, by = 0, all = TRUE)
rownames(tax.otu.m) <- tax.otu.m$Row.names
tax.otu.m$Row.names <- NULL

rm(tax.m, otu.m)


# write.table(tax.otu.m, "SMP_OTU_16S_Moss.csv", sep = "\t", row.names = FALSE)
# write.table(tax.otu.m, "SMP_OTU_18S_Moss.csv", sep = "\t", row.names = FALSE)


################################################################################
### nMDS pitchers
smp.log <- transform_sample_counts(euk.a, function(otu) {log1p(otu)})


## Eukaryotes only: remove outlier with only one Saccharomyces (LM3a08L)
smp.log <- prune_samples(smp.log@sam_data$FullID != "LM3a08L", smp.log)
any(taxa_sums(smp.log) == 0)

smp.nmds <- ordinate(smp.log, method = "NMDS", distance = "bray", k = 4,
                     autotransform = FALSE, trymax = 100)
# smp.nmds <- ordinate(smp.log, method = "PCoA", distance = "bray")
smp.nmds


stressplot(smp.nmds)
# barplot(smp.nmds$values$Relative_eig)
# biplot(smp.nmds, data.frame(otu_table(smp.a)))


plot_ordination(smp, smp.nmds, shape = "Succession",
                color = "Site", title = NULL, #label = "FullID",
                axes = 1:2) +
  stat_ellipse(aes(group = Site), type = "t", linetype = 2, size = 0.2) +
  geom_point(size = 3) +
  scale_colour_manual(values = gradCol) +
  coord_fixed(ratio = 1) +
  theme(legend.position = "top", legend.direction = "horizontal",
        legend.box="vertical") +
  # guides(fill = guide_legend(nrow = 2, byrow = TRUE)) + # grrr
  xlab("nMDS1") +
  ylab("nMDS2")
# ggsave("SMP_16S_nMDS_ASV.pdf", width = 8.27, height = 8.27)
# ggsave("SMP_18S_nMDS_ASV.pdf", width = 8.27, height = 8.27)
# ggsave("SMP_nMDS_ASV.pdf", width = 8.27, height = 8.27)


### nMDS pitchers and mosses
smp.moss <- merge_phyloseq(prok.a, moss.prok.a)
smp.moss <- merge_phyloseq(euk.a, moss.euk.a)


smp.moss.log <- transform_sample_counts(smp.moss, function(otu) {log1p(otu)})


smp.moss.nmds <- ordinate(smp.moss.log, method = "NMDS", distance = "bray",
                          k = 4, autotransform = FALSE, trymax = 100,
                          weakties = FALSE)
smp.moss.nmds


stressplot(smp.moss.nmds)
# barplot(smp.moss.nmds$values$Relative_eig)
# biplot(smp.moss.nmds, data.frame(otu_table(smp.a)))


plot_ordination(smp.moss, smp.moss.nmds, shape = "Succession",
                color = "Site", title = NULL, #label = "FullID",
                axes = 1:2) +
  stat_ellipse(aes(group = Succession), type = "t", linetype = 2, size = 0.2) +
  geom_point(size = 3) +
  scale_colour_manual(values = gradCol) +
  coord_fixed(ratio = 1) +
  theme(legend.position = "top", legend.direction = "horizontal") +
  guides(linetype = guide_legend(nrow = 2)) +
  xlab("nMDS1") +
  ylab("nMDS2")
# ggsave("SMP_16S_nMDS_withMoss.pdf", width = 8.27, height = 8.27)
# ggsave("SMP_18S_nMDS_withMoss.pdf", width = 8.27, height = 8.27)


################################################################################
### PERMANOVA
## Prokaryotes
tax.l <- as.data.frame(prok.a@tax_table@.Data)
otu.l <- t(as.data.frame(otu_table(prok.a)))

tax.otu.prok <- merge(tax.l, otu.l, by = 0, all = TRUE)
rownames(tax.otu.prok) <- tax.otu.prok$Row.names
colnames(tax.otu.prok)[1] <- "Taxa"

tax.otu.prok.t <- data.frame(t(tax.otu.prok[c(8:167)]))
tax.otu.prok.t$Site <- gsub(".{5}$", "", rownames(tax.otu.prok.t))

rm(tax.l, otu.l, tax.otu.prok)


prok.adonis <- adonis(tax.otu.prok.t[, -346] ~ tax.otu.prok.t$Site,
                      method = "bray", perm = 999)
summary(prok.adonis)
prok.adonis


## Eukaryotes
tax.l <- as.data.frame(euk.a@tax_table@.Data)
otu.l <- t(as.data.frame(otu_table(euk.a)))

tax.otu.euk <- merge(tax.l, otu.l, by = 0, all = TRUE)
rownames(tax.otu.euk) <- tax.otu.euk$Row.names
colnames(tax.otu.euk)[1] <- "Taxa"

tax.otu.euk.t <- data.frame(t(tax.otu.euk[c(8:166)]))
tax.otu.euk.t$Site <- gsub(".{5}$", "", rownames(tax.otu.euk.t))

rm(tax.l, otu.l, tax.otu.euk)


euk.adonis <- adonis(tax.otu.euk.t[, -115] ~ tax.otu.euk.t$Site, method = "bray",
                     perm = 999)
euk.adonis


## Combined
tax.otu.t <- data.frame(t(tax.otu.l[c(8:167)]))
tax.otu.t$Site <- gsub(".{5}$", "", rownames(tax.otu.t))


smp.adonis <- adonis(tax.otu.t[, -460] ~ tax.otu.t$Site, method = "bray",
                     perm = 999)
smp.adonis


################################################################################
### Clustering
## UPGMA
smp.bc <- vegdist(t(tax.otu.prok[, -c(1:7)]), method = "bray")

smp.bc.UPGMA <- hclust(smp.bc, method = "average")
plot(smp.bc.UPGMA, main = "", cex = 0.5)


smp.bc.UPGMA.coph <- cophenetic(smp.bc.UPGMA)
cor(smp.bc, smp.bc.UPGMA.coph)


plot(smp.bc, smp.bc.UPGMA.coph, xlab = "Bray-Curtis dissimilarity",
     ylab="Cophenetic distance", asp=1, xlim=c(0,1), ylim=c(0,1),
     main=c(paste("UPGMA \n", "Cophenetic correlation =",
                  round(cor(smp.bc, smp.bc.UPGMA.coph),3))),
     cex=0.8, cex.main=0.9)
abline(0,1)
lines(lowess(smp.bc, smp.bc.UPGMA.coph), col="red")


asw <- numeric(nrow(t(tax.otu.prok[, -c(1:7)])))

for (k in 2:(nrow(t(tax.otu.prok[, -c(1:7)])) - 1)) {
  sil <- silhouette(cutree(smp.bc.UPGMA, k = k), smp.bc)
  asw[k] <- summary(sil)$avg.width
}
k.best <- which.max(asw)

plot(1:nrow(t(tax.otu.prok[, -c(1:7)])), asw, type = "h",
     xlab = "k", ylab = "Average silhouette width")
axis(1, k.best, paste(k.best, sep = "\n"), col = "red", font = 2,
     col.axis = "red")
points(k.best, max(asw), pch = 16, col = "red", cex = 1.5)
cat("", "Silhouette-optimal number of clusters k =", k.best, "\n",
    "with an average silhouette width of", max(asw), "\n")


## AP Clustering
d.apclus <- apcluster(negDistMat(r = 1), as.matrix(smp.bc))

cat("affinity propogation optimal number of clusters:",
    length(d.apclus@clusters), "\n")


# pdf("SMP_APCluster_prok.pdf", height = 8.27, width = 8.27)
# pdf("SMP_APCluster_euk.pdf", height = 8.27, width = 8.27)
# pdf("SMP_APCluster_all.pdf", height = 8.27, width = 8.27)
heatmap(d.apclus, as.matrix(smp.bc), margin = c(3.5, 3.5),
        cexRow = 0.3, cexCol = 0.3)
dev.off()


plot(d.apclus, as.matrix(smp.bc))


################################################################################
### Plot taxa
smp.log <- transform_sample_counts(smp, function(otu) {log1p(otu)})

plot_bar(smp.log, x = "Class", fill = "Site") +
  geom_bar(aes(color = Site, fill = Site), stat = "identity",
           position = "stack") +
  facet_grid( ~ Domain, scales = "free", space = "free_y") +
  theme(legend.position = "top") +
  scale_color_manual(values = gradCol) +
  scale_fill_manual(values = gradCol) +
  xlab("") +
  ylab("ln+1(Abundance) [%]") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_flip()
# ggsave("SMP_Phyla_prok.pdf", width = 8.27, height = 11.69)
# ggsave("SMP_Phyla_euk.pdf", width = 8.27, height = 11.69)


### Without S4-class
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


## Reorder taxonomically
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
         "Bicosoecida", "Ochrophyta_ph", "Peronosporomycetes", ## SAR:Heterokonta
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
tax.otu.lm$Count[tax.otu.lm$Count < 0.01 & tax.otu.lm$Count > 0] <- 0.1

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


aggregate(Count ~ Class, data = tax.otu.lm, FUN = sum)


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
inc$Phylum <- factor(inc$Phylum,
                                    levels = unique(as.character(
                                      inc$Phylum)))
inc$IncidencePercent <- 100 / 160 * inc$Incidence


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
ggsave("SMP_Incidence.pdf", width = 8.27, height = 11.69 / 2)



## Control
# tax.log <- as.data.frame(smp.log@tax_table@.Data)
# otu.log <- t(as.data.frame(otu_table(smp.log)))
#
# tax.otu.log <- merge(tax.log, otu.log, by = 0, all = TRUE)
# rownames(tax.otu.log) <- tax.otu.l$Row.names
# colnames(tax.otu.log)[1] <- "Taxa"
#
# tax.otu.logm <- reshape2::melt(tax.otu.log, id.vars = c("Taxa", "Domain",
#                                                         "Phylum", "Class",
#                                                         "Order", "Family",
#                                                         "Genus"),
#                                variable.name = "FullID", value.name = "Count")
#
# aggregate(Count ~ Class, data = tax.otu.logm, FUN = sum)


### Check for (spurious) taxa
rank_names(smp)
smp.taxa <- get_taxa_unique(smp, taxonomic.rank = rank_names(smp))
# moss.taxa <- get_taxa_unique(moss.a, taxonomic.rank = rank_names(moss.a)[5])

taxa.g <- grep("unknown", smp.taxa, ignore.case = TRUE)
smp.taxa[taxa.g]


## Explore taxa
smp.sub <- subset_taxa(smp, (Phylum %in% c("Rotifera")))
tax.df <- data.frame(tax_table(smp.sub))

summary(is.na(tax.df$Genus))
dim(tax.df)

otu.df <- data.frame(otu_table(smp.sub))
otu.df[otu.df > 1] <- 1
rowSums(otu.df)
summary(rowSums(otu.df) == 0)

get_taxa_unique(smp.sub, taxonomic.rank = "Genus")


################################################################################
################################################################################
################################################################################
################################################################################
