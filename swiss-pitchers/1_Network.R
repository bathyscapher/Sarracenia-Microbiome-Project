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
theme_set(theme_bw(base_size = 20) +
            theme(rect = element_rect(fill = "transparent")))
library("reshape2")
library("igraph")
library("Hmisc")


setwd("~/Sarracenia-Microbiome-Project/Thesis/")


rm(list = ls())
set.seed(30008)


################################################################################
## Colorblind temperature scale for the 5 sites
gradCol <- c("#D55E00", "#E69F00", "#F0E442", "#009E73", "#56B4E9")


################################################################################
### Pitcher data
smp <- readRDS("rds/SMP_smp.RDS")


################################################################################
### Approach with rank correlation
otus <- data.frame(t(otu_table(smp)))


## Presence/absence and taxa that occur at least in 10 leaves and then, subset
otus.pa <- data.frame(t(otu_table(smp)))
otus.pa[otus.pa > 0] <- 1
otus.pa <- otus.pa[rowSums(otus.pa) > 9, ]
summary(rowSums(otus.pa))
subset.inc <- rownames(otus.pa)

otus <- otus[rownames(otus) %in% subset.inc, ]


## Correlation analysis based on spearman's co-efficient
otus.dist <- rcorr(t(otus), type = "spearman")
otus.cor <- otus.dist$r
otus.cor.p <- otus.dist$P


## Multiple testing correction using Benjamini-Hochberg standard FDR correction
otus.cor.p <- p.adjust(otus.cor.p, method = "BH")


## Positive and netagive cooccurence at given coefficient and p-value cutoff
cor.cutoff <- 0.4
p.cutoff <- 0.01


otus.cor[which(otus.cor >= (- cor.cutoff) & otus.cor <= cor.cutoff)] <- 0
otus.cor[which(otus.cor.p > p.cutoff)] <- 0


## Delete rows and columns with sum = 0
otus.cor <- otus.cor[which(rowSums(otus.cor) != 1), ]
otus.cor <- otus.cor[, which(colSums(otus.cor) != 0)]


## Create a graph
g3 <- graph.adjacency(otus.cor, weight = TRUE, mode = "undirected")
g3 <- simplify(g3) # remove duplicate and loop edges


### Add taxonomy
## Extract taxonomy and merge
tax <- data.frame(tax_table(smp))
otus.tax <- merge(tax, otus, by = 0)
names(otus.tax)[1] <- "ASV"


## Fill unannotated genera with higher-level taxa
otus.tax$Genus <- ifelse(is.na(otus.tax$Genus), otus.tax$Family, otus.tax$Genus)
otus.tax$Genus <- ifelse(is.na(otus.tax$Genus), otus.tax$Order, otus.tax$Genus)
otus.tax$Genus <- ifelse(is.na(otus.tax$Genus), otus.tax$Class, otus.tax$Genus)
otus.tax$Genus <- ifelse(is.na(otus.tax$Genus), otus.tax$Phylum, otus.tax$Genus)


## Total abundance
otus.tax$Abundance <- rowSums(otus.tax[, 8:167])
# otus.tax$Abundance[otus.tax$Genus == "Rubrivivax"]
# min(otus.tax$Abundance)


## Get taxa in network and merge with taxonomy
target <- data.frame(ASV = V(g3)$name,
                     to_sort = seq(1:length(V(g3)$name)))


otu.target <- merge(target, otus.tax, by = "ASV")
otu.target <- otu.target[order(otu.target$to_sort), ]


E(g3)
V(g3)
V(g3)$name
V(g3)$degree <- degree(g3) # connectivity


## Get pairwise correlations
otus.cor[lower.tri(otus.cor, diag = TRUE)] <- NA
otus.cor.m <- melt(otus.cor, value.name = "rho")
otus.cor.m$rho[otus.cor.m$rho == 0] <- NA
otus.cor.m <- otus.cor.m[complete.cases(otus.cor.m), ]


### Add edge attributes
## Rho as weight
E(g3)$width <- otus.cor.m$rho


## Get sign
E(g3)$sign <- ifelse(E(g3)$width < 0, "Negative", "Positive")


## Convert to positive numbers for scaled display
E(g3)$scale <- ifelse(E(g3)$width < 0, E(g3)$width * -1, E(g3)$width)


### Add vertex attributes
g3 <- set_vertex_attr(g3, "ASV", value = otu.target$ASV)
V(g3)$ASV

g3 <- set_vertex_attr(g3, "Genus", value = otu.target$Genus)
V(g3)$Genus
V(g3)$name <- V(g3)$Genus

g3 <- set_vertex_attr(g3, "Domain", value = otu.target$Domain)
V(g3)$Domain

g3 <- set_vertex_attr(g3, "Abundance", value = otu.target$Abundance)
V(g3)$Abundance


### Louvain community detection
## Set negative relations to zero
E(g3)$width <- ifelse(E(g3)$width < 0, 0, E(g3)$width)
E(g3)$width

g3.louv <- cluster_louvain(g3, weights = E(g3)$width)
sizes(g3.louv)


V(g3)$louvain.member <- g3.louv$membership


### Export graph
write_graph(g3, "SMP_Co-occurrence.graphml", format = "graphml")


## Network statistics
vcount(g3)
ecount(g3)

table(V(g3)$Domain)
table(degree(g3))
mean(degree(g3))
table(degree(g3))
summary(degree(g3))


## Negative and positive edges
summary(otus.cor.m$rho > 0)


## Subnetworks
components(g3)


## Clusters
sizes(g3.louv)
table(sizes(g3.louv))


################################################################################
################################################################################
################################################################################
################################################################################
