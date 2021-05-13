################################################################################
################################################################################
################################################################################
################################################################################
### SMP Distance decay
### Authors: Rachel Korn
### korn@cumulonimbus.at University of Fribourg 2019/2020
################################################################################


library("phyloseq")
library("ggplot2")
theme_set(theme_bw(base_size = 20) +
            theme(rect = element_rect(fill = "transparent")))
library("vegan")
library("reshape2")
library("scales")
library("lme4")


rm(list = ls())
setwd("~/Sarracenia-Microbiome-Project/Thesis")


################################################################################
## Colorblind temperature scale for the 5 sites
gradCol <- c("#D55E00", "#E69F00", "#F0E442", "#009E73", "#56B4E9")



################################################################################
### Pitcher data
prok.a <- readRDS("rds/SMP_prok.a.RDS")
euk.a <- readRDS("rds/SMP_euk.a.RDS")


### As data.frames
## Prokaryotes
tax.p <- as.data.frame(prok.a@tax_table@.Data)
otu.p <- t(as.data.frame(otu_table(prok.a)))

tax.otu.p <- merge(tax.p, otu.p, by = 0, all = TRUE)
rownames(tax.otu.p) <- tax.otu.p$Row.names
tax.otu.p$Row.names <- NULL

rm(tax.p, otu.p)


## Eukaryotes
tax.e <- as.data.frame(euk.a@tax_table@.Data)
otu.e <- t(as.data.frame(otu_table(euk.a)))

tax.otu.e <- merge(tax.e, otu.e, by = 0, all = TRUE)
rownames(tax.otu.e) <- tax.otu.e$Row.names
tax.otu.e$Row.names <- NULL
# tax.otu.e$LT4a06E <- 0

rm(tax.e, otu.e)


## Merge
# tax.otu.df <- rbind(tax.otu.p, tax.otu.e)


################################################################################
## Bray Curtis distance
tax.otu.bc.p <- vegdist(t(as.matrix(tax.otu.p[, -c(1:6)])), method = "bray",
                        upper = TRUE)
tax.otu.bc.e <- vegdist(t(as.matrix(tax.otu.e[, -c(1:6)])), method = "bray",
                        upper = TRUE)


################################################################################
### Distance decay with Bray-Curtis distances
### Geographic distance
smpMeta.df <- as(sample_data(prok.a), "data.frame")
distance <- read.table("csv/SMP_2018_DistanceBetweenPlants.csv",
                       header = TRUE, sep = ",", check.names = FALSE)


distance <- merge(distance, smpMeta.df[, c(1, 4)], by.x = c("ID"),
                  by.y = "PlantID")
rownames(distance) <- distance$FullID
distance <- distance[with(distance, order(FullID)), ]
distance$ID <- NULL
distance$FullID <- NULL


distance <- as.data.frame(t(distance))
distance <- cbind("ID" = row.names(distance), distance)


distance <- merge(distance, smpMeta.df[, c(1, 4)], by.x = "ID",
                  by.y = "PlantID")
rownames(distance) <- distance$FullID
distance <- distance[with(distance, order(FullID)), ]
distance$ID <- NULL
distance$FullID <- NULL

distance[lower.tri(distance)] <- NA
diag(distance) <- NA

distance.m <- melt(as.matrix(distance), value.name = "Geo_distance")
distance.m <- distance.m[complete.cases(distance.m), ]
distance.m$Pair <- as.factor(paste(distance.m$Var1, distance.m$Var2, sep = "-"))


### Community distances
## Prokaryotes
bc.p <- as.matrix(tax.otu.bc.p)
bc.p[lower.tri(bc.p)] <- NA
diag(bc.p) <- NA


bc.pm <- melt(bc.p, value.name = "BC_distance")
bc.pm <- bc.pm[complete.cases(bc.pm), ]
bc.pm$Pair <- as.factor(paste(bc.pm$Var1, bc.pm$Var2, sep = "-"))
bc.pm$Domain <- "Prokaryotes"


## Eukaryotes
bc.e <- as.matrix(tax.otu.bc.e)
bc.e[lower.tri(bc.e)] <- NA
diag(bc.e) <- NA


bc.em <- melt(bc.e, value.name = "BC_distance")
bc.em <- bc.em[complete.cases(bc.em), ]
bc.em$Pair <- as.factor(paste(bc.em$Var1, bc.em$Var2, sep = "-"))
bc.em$Domain <- "Eukaryotes"


## Merge domains
bc.m <- rbind(bc.pm, bc.em)
bc.m$Domain <- ordered(bc.m$Domain, levels = c("Prokaryotes", "Eukaryotes"))


### Merge with distance
bc.distance <- merge(bc.m, distance.m, by = "Pair")
rm(bc.em, bc.pm, bc.m, distance, distance.m)

bc.distance$Pairs <- as.factor(paste(gsub(".{5}$", "", bc.distance$Var1.x),
                                     gsub(".{5}$", "", bc.distance$Var2.x),
                                     sep = "-"))


## Add succession
bc.distance$Succession <- as.factor(paste(gsub("^.{6}", "", bc.distance$Var1.x),
                                          gsub("^.{6}", "", bc.distance$Var2.x),
                                          sep = "-"))
levels(bc.distance$Succession) <- c("Early", "Mix", "Mix", "Late")


## Only within sites
bc.distance.s <- subset(bc.distance, bc.distance$Pairs == "CB-CB" |
                          bc.distance$Pairs == "LE-LE" |
                          bc.distance$Pairs == "LM-LM" |
                          bc.distance$Pairs == "LV-LV" |
                          bc.distance$Pairs == "LT-LT")
bc.distance.s$Pairs <- droplevels(bc.distance.s$Pairs)
bc.distance.s$Pairs <- ordered(bc.distance.s$Pairs,
                               levels = c("CB-CB", "LT-LT", "LE-LE", "LV-LV",
                                          "LM-LM"))
bc.distance.s$Site <- gsub(".{3}", "", bc.distance.s$Pairs)
bc.distance.s$Site <- ordered(bc.distance.s$Site,
                              levels = c("CB", "LT", "LE", "LV", "LM"))


## Colorblind temperature scale for the 5 sites
levels(bc.distance$Pairs)
ramp <- hue_pal()(15)
ramp[1] <- "#D55E00" # CB
ramp[13] <- "#E69F00" # LT
ramp[6] <- "#F0E442" # LE
ramp[15] <- "#009E73" # LV
ramp[10] <- "#56B4E9" # LM
ramp[-c(1, 13, 6, 15, 10)] <- "gray"


## Binned on geo distance
ggplot(bc.distance, aes(x = round(log1p(Geo_distance), digits = 1),
                        y = BC_distance,
                        group = round(log1p(Geo_distance), digits = 1))) +
  geom_boxplot() +
  geom_smooth(aes(x = log1p(Geo_distance), y = BC_distance, group = 1),
              method = "lm", formula = y ~ x, se = FALSE, color = "gray") +
  geom_smooth(data = bc.distance.s,
              aes(x = log1p(Geo_distance), y = BC_distance,
                  group = Site, color = Site),
              method = 'lm', formula = y ~ x, se = FALSE) +
  facet_grid(~ Domain) +
  scale_colour_manual(values = gradCol) +
  ylim(0, 1) +
  theme(legend.position = "top") +
  xlab(expression(paste(italic(ln), "(Geographical distance + 1) [m]")))  +
  ylab("Bray-Curtis dissimilarity")
# ggsave("SMP_DistanceDecay_bin.pdf", width = 11.69, height = 6.5)


## Regression coefficients
bc.distance.d <- subset(bc.distance.s, bc.distance$Domain == "Eukaryotes")


dd.fit <- lmList(BC_distance ~ log1p(Geo_distance) | Pairs,
                 data = bc.distance.d, na.action = na.exclude)
summary(dd.fit)


### BC distance sorted by median
## Subset by domain (change for Prok- and Eukaryotes, respectively)
bc.distance.d <- bc.distance[bc.distance$Domain == "Prokaryotes", ]
# bc.distance.d <- bc.distance[bc.distance$Domain == "Eukaryotes", ]


## Sort by median
fac <- with(bc.distance.d, reorder(Pairs, BC_distance, median, order = TRUE))
bc.distance.d$Pairs <- factor(bc.distance.d$Pairs, levels = levels(fac))


## Color by site
bc.distance.d$col.box <- "#ffffff"
bc.distance.d$col.box[bc.distance.d$Pairs == "CB-CB"] <- "#D55E00"
bc.distance.d$col.box[bc.distance.d$Pairs == "LT-LT"] <- "#E69F00"
bc.distance.d$col.box[bc.distance.d$Pairs == "LE-LE"] <- "#F0E442"
bc.distance.d$col.box[bc.distance.d$Pairs == "LV-LV"] <- "#009E73"
bc.distance.d$col.box[bc.distance.d$Pairs == "LM-LM"] <- "#56B4E9"


## Only early or late succession pairs
bc.distance.d <- bc.distance.d[!bc.distance.d$Succession == "Mix", ]


ggplot(bc.distance.d, aes(x = Pairs, y = BC_distance, fill = col.box)) +
  geom_boxplot(outlier.shape = 21) +
  facet_wrap(~ Domain) +
  scale_colour_manual(values = ramp) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        # axis.text.y = element_blank(), # outcomment for eukaryotes
        legend.position = "top") +
  ylim(0, 1) +
  xlab("") +
  ylab("Bray-Curtis dissimilarity") +
  scale_fill_identity()
# ggsave("SMP_16S_BC-dissimilarity.pdf", width = 6.5, height = 8.27)
# ggsave("SMP_18S_BC-dissimilarity.pdf", width = 6.5, height = 8.27)


################################################################################
################################################################################
################################################################################
################################################################################
