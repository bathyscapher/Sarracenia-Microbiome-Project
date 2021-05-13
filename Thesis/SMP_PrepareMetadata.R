################################################################################
################################################################################
################################################################################
################################################################################
### SMP Prepare metadata
### Authors: Rachel Korn
### korn@cumulonimbus.at University of Fribourg 2020
################################################################################


library("phyloseq")
library("ggplot2")
library("GGally")
theme_set(theme_bw(base_size = 20) +
            theme(rect = element_rect(fill = "transparent")))
library("vegan")
library("gplots")
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
### Prepare metadata for SEM
smpMeta.df <- data.frame(sample_data(prok.a))


## Merge resequenced samples and double-spiked mosses
smpMeta.df <- smpMeta.df[!grepl("r$", smpMeta.df$FullID), ]
mossMeta.df <- smpMeta.df[smpMeta.df$Succession == "Moss", ]
mossMeta.df <- mossMeta.df[!grepl("Q$", mossMeta.df$FullID), ]


smpMeta.df <- smpMeta.df[smpMeta.df$Site != "NC" &
                           smpMeta.df$Site != "MC" &
                           smpMeta.df$Succession != "Moss", ]


### Prey composition
prey <- read.table("csv/SMP_Prey_clean.csv", sep = "\t", header = TRUE)


rownames(prey) <- prey$FullID
prey$FullID <- NULL
prey$Site <- NULL


## Remove empty columns
prey <- log1p(prey)


## Reduce to 2 dimensions
prey.nmds <- vegan:::metaMDS(prey, distance = "euclidean", k = 2, trymax = 100,
                             autotransform = FALSE)
prey.nmds
stressplot(prey.nmds)
plot(prey.nmds)
orditorp(prey.nmds, display = "species", col = "red", air = 0.01)
orditorp(prey.nmds, display = "sites", cex = 0.75, air = 0.01)


prey.nmds.sc <- as.data.frame(scores(prey.nmds)[, 1])
names(prey.nmds.sc) <- "prey"
prey.nmds.sc$FullID <- rownames(prey.nmds.sc)

smpMeta.df <- merge(smpMeta.df, prey.nmds.sc, by = "FullID", all = TRUE)


### Leaf microbiome
## Prokaryotes
leaf.prok.mb <- as.data.frame(otu_table(prok.a))
leaf.prok.mb <- log1p(leaf.prok.mb)


leaf.prok.nmds <- vegan:::metaMDS(leaf.prok.mb, distance = "bray", k = 4,
                                  trymax = 100, autotransform = FALSE)
leaf.prok.nmds
stressplot(leaf.prok.nmds)
# plot(goodness(leaf.prok.nmds))


leaf.prok.nmds.sc <- as.data.frame(scores(leaf.prok.nmds)[, 1])
names(leaf.prok.nmds.sc) <- "mb.leaf.prok"
leaf.prok.nmds.sc$FullID <- rownames(leaf.prok.nmds.sc)


smpMeta.df <- merge(smpMeta.df, leaf.prok.nmds.sc, by = "FullID", all = TRUE)


## Eukaryotes
leaf.euk.mb <- as.data.frame(otu_table(euk.a))
leaf.euk.mb <- log1p(leaf.euk.mb)


leaf.euk.nmds <- vegan:::metaMDS(leaf.euk.mb, distance = "bray", k = 4,
                                 trymax = 100, autotransform = FALSE)
leaf.euk.nmds
stressplot(leaf.euk.nmds)
# plot(goodness(leaf.euk.nmds))


leaf.euk.nmds.sc <- as.data.frame(scores(leaf.euk.nmds)[, 1])
names(leaf.euk.nmds.sc) <- "mb.leaf.euk"
leaf.euk.nmds.sc$FullID <- rownames(leaf.euk.nmds.sc)


smpMeta.df <- merge(smpMeta.df, leaf.euk.nmds.sc, by = "FullID", all = TRUE)


################################################################################
### Moss microbiome
## Prokaryotes
moss.prok.mb <- as.data.frame(otu_table(moss.prok.a))
moss.prok.mb <- log1p(moss.prok.mb)


moss.prok.nmds <- vegan:::metaMDS(moss.prok.mb, distance = "bray",
                                  trymax = 50, k = 2, autotransform = FALSE)
moss.prok.nmds
stressplot(moss.prok.nmds)


moss.prok.nmds.sc <- as.data.frame(scores(moss.prok.nmds)[, 1])
names(moss.prok.nmds.sc) <- "mb.moss.prok"
moss.prok.nmds.sc$FullID <- rownames(moss.prok.nmds.sc)
sample_data(moss.prok.a)$mb.moss.prok <- moss.prok.nmds.sc$mb.moss.prok


smpMeta.df <- merge(smpMeta.df, sample_data(moss.prok.a)[, c(2, 6, 24)],
                    by = c("Site", "Sector"), all = TRUE)


## Eukaryotes
moss.euk.mb <- as.data.frame(otu_table(moss.euk.a))
moss.euk.mb <- log1p(moss.euk.mb)


moss.euk.nmds <- vegan:::metaMDS(moss.euk.mb, distance = "bray",
                                 trymax = 50, k = 2, autotransform = FALSE)
moss.euk.nmds
stressplot(moss.euk.nmds)


moss.euk.nmds.sc <- as.data.frame(scores(moss.euk.nmds)[, 1])
names(moss.euk.nmds.sc) <- "mb.moss.euk"
moss.euk.nmds.sc$FullID <- rownames(moss.euk.nmds.sc)
sample_data(moss.euk.a)$mb.moss.euk <- moss.euk.nmds.sc$mb.moss.euk


smpMeta.df <- merge(smpMeta.df, sample_data(moss.euk.a)[, c(2, 6, 24)],
                    by = c("Site", "Sector"), all = TRUE)


## Convert factor to dummy
smpMeta.df$SampleColor <- as.factor(smpMeta.df$SampleColor)
smpMeta.df$SampleColor <- ordered(smpMeta.df$SampleColor,
                                  levels = c("clear", "cloudy", "green",
                                             "brown"))
smpMeta.df$SampleColorDummy <- as.integer(unclass(smpMeta.df$SampleColor))


### Write unscaled data
# write.table(smpMeta.df, "csv/SMP_Metadata_SEM_raw.csv", sep = "\t")


### Scale numeric variables
str(smpMeta.df[, -c(1:8, 10, 18)])

par(mfrow = c(2, 1))
boxplot(smpMeta.df[, -c(1:8, 10, 18)])
abline(h = 0, lty = 2, col = "gray")


smpMeta.df[, -c(1:8, 10, 18)] <- apply(smpMeta.df[, -c(1:8, 10, 18)], 2, scale)


boxplot(smpMeta.df[, -c(1:8, 10, 18)])
abline(h = 0, lty = 2, col = "gray")
dev.off()


### Write and read the data
# write.table(smpMeta.df, "csv/SMP_Metadata_SEM_scaled.csv", sep = "\t")
# smpMeta.df <- read.table("csv/SMP_Metadata_SEM_scaled.csv", sep = "\t")


################################################################################
## Overview plots
### Plot the pairs
smpMeta.df$prey.items <- log1p(smpMeta.df$prey.items)


p <- ggpairs(smpMeta.df[, - c(2:8, 10, 13, 15:18)],
             aes(colour = Site, alpha = 0.4),
             upper = list(continuous = wrap("cor", size = 2)),
             lower = list(continuous = wrap("points", alpha = 0.3, size = 0.7),
                          combo = wrap("dot", alpha = 0.4, size = 0.2))) +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_bw(base_size = 10)


for(i in 1:p$nrow) {
  for(j in 1:p$ncol){
    p[i,j] <- p[i,j] +
      scale_fill_manual(values = gradCol) +
      scale_color_manual(values = gradCol)
  }
}

p
# ggsave("SMP_MetadataPairs.pdf", width = 11.69, height = 11.69)


################################################################################
### Leaf metrics
p <- ggpairs(smpMeta.df[, c(1, 9, 11, 13, 15:17)],
             aes(colour = Site, alpha = 0.4),
             upper = list(continuous = wrap("cor"))) +
  theme_bw(base_size = 10) +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1))


for(i in 1:p$nrow) {
  for(j in 1:p$ncol){
    p[i,j] <- p[i,j] +
      scale_fill_manual(values = gradCol) +
      scale_color_manual(values = gradCol)
  }
}

p
# ggsave("SMP_MetadataPairs_Morphology.pdf", width = 8.27, height = 8.27)


## PCA
leafPCA <- prcomp(~ SampleVolume_mL + PotentialVolume_mL + PitcherLength +
                    KeelWidth + MouthWidth + PitcherWidth,
                  data = smpMeta.df, scale = TRUE)
leafPCA$rotation # loadings
print(summary(leafPCA), digits = 5)


pdf("SMP_Leafmorphometrics.pdf", height = 4, width = 8.27)
par(mfrow = c(1, 2))
barplot(summary(leafPCA)$importance[3, ], las = 2,
        ylab = "Cumulative variance explained")
biplot(leafPCA, col = c("gray", "black"), cex = 0.7, #xlim = c(-0.2, 0.4),
       asp = 1, las = 1, xlim = c(-0.45, 0.2),
       xlab = "PC1 (total variance explained)",
       ylab = "PC2 (total variance explained)")
abline(h = 0, lty = 2, col = "black")
abline(v = 0, lty = 2, col = "black")
dev.off()


################################################################################
## Unscaled data
smpMeta.df <- read.table("csv/SMP_Metadata_SEM_scaled.csv", sep = "\t")


smpMeta.df$Sector <- as.factor(smpMeta.df$Sector)
smpMeta.df$Site <- ordered(smpMeta.df$Site,
                           levels = c("CB", "LT", "LE", "LV", "LM" ))


### Correlation and covariance
my_palette <- colorRampPalette(c("#D55E00", "white",  "#56B4E9"))(n = 19)


## Get covariance and sort data.frame
covs <- data.frame(cov(smpMeta.df[, -c(1:8, 10:11, 13, 15:18)],
                       use = "complete"))
covs <- covs[c("mb.leaf.prok", "mb.leaf.euk", "mb.moss.prok", "mb.moss.euk",
               "prey", "prey.items", "Altitude", "MeanTemperature", "Distance",
               "CanopyCover",  "Age_d", "SampleVolume_mL", "pH",
               "SampleColorDummy")]
covs <- covs[colnames(covs), ]


pdf("SMP_Metadata_Covariance.pdf", height = 8.27, width = 8.27)
heatmap.2(as.matrix(covs), dendrogram = "none", key = TRUE, key.title = "",
          cellnote = round(covs, digits = 2), notecol = "black", notecex = 0.65,
          trace = "none", distfun = function(x) {x}, cexCol = 1.4, cexRow = 1.4,
          symm = TRUE, margins = c(16, 16), col = my_palette,
          tracecol = "black", Rowv = FALSE, Colv = FALSE)
dev.off()


################################################################################
### Plot pH
ggplot(smpMeta.df) +
  geom_hline(yintercept = 7, linetype = "dashed", color = "gray") +
  geom_boxplot(aes(x = Site, y = pH), color = "black",
               size = 0.2, outlier.shape = NA, na.rm = TRUE) +
  geom_jitter(aes(x = Site, y = pH, fill = Site), width = 0.3, height = 0,
              pch = 21, size = 2, na.rm = TRUE) +
  facet_grid( ~ Site, scales = "free") +
  scale_fill_manual(values = gradCol) +
  theme(legend.position = "none", legend.direction = "horizontal",
        axis.text.x = element_blank()) +
  xlab("") +
  ylab("pH")
ggsave("SMP_pH.pdf", width = 11.69, height = 5)


################################################################################
### Plot sample color
smpMeta.df$SampleColor <- ordered(smpMeta.df$SampleColor,
                                  levels = c("clear", "cloudy", "green",
                                             "brown"))

## Color scale
sampleCol <- c("#EBE7E0", "#56B4E9", "#009E73", "#C59434")


ggplot(smpMeta.df) +
  geom_bar(aes(x = SampleColor, fill = SampleColor)) +
  facet_grid(Sector ~ Site) +
  scale_fill_manual(values = sampleCol) +
  theme(legend.position = "top", legend.direction = "horizontal",
        axis.text.x = element_blank()) +
  xlab("") +
  ylab("Frequency")
# ggsave("SMP_SampleColors.pdf", width = 11.69, height = 8.27)


################################################################################
### Plot canopy openness
ggplot(smpMeta.df) +
  geom_boxplot(aes(x = Site, y = CanopyCover), color = "black",
               size = 0.2, outlier.shape = NA) +
  geom_jitter(aes(x = Site, y = CanopyCover, fill = Site), na.rm = TRUE,
              height = 0, pch = 21, size = 2) +
  facet_grid( ~ Site, scales = "free") +
  scale_fill_manual(values = gradCol) +
  theme(legend.position = "none", legend.direction = "horizontal",
        axis.text.x = element_blank()) +
  xlab("") +
  ylab("Canopy cover [%]")
# ggsave("SMP_CanopyCover.pdf", width = 11.69, height = 5)


################################################################################
################################################################################
################################################################################
################################################################################
