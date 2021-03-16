################################################################################
################################################################################
################################################################################
################################################################################
### SMP Pitchers and mosses
### Authors: Rachel Korn
### korn@cumulonimbus.at University of Fribourg 2020
################################################################################


library("phyloseq")
library("ggplot2")
theme_set(theme_bw(base_size = 20) +
            theme(rect = element_rect(fill = "transparent")))
library("VennDiagram")
library("gridExtra")
library("ggalluvial")


rm(list = ls())
setwd("~/Sarracenia-Microbiome-Project/Thesis")


################################################################################
## Colorblind temperature scale for the 5 sites
gradCol <- c("#D55E00", "#E69F00", "#F0E442", "#009E73", "#56B4E9")


################################################################################
### Pitcher data
prok.a <- readRDS("rds/SMP_prok.a.RDS")
euk.a <- readRDS("rds/SMP_euk.a.RDS")


################################################################################
### Read moss taxa
moss.prok.a <- readRDS("rds/SMP_moss20.prok.a.RDS")
moss.euk.a <- readRDS("rds/SMP_moss20.euk.a.RDS")


################################################################################
### Pitcher overlap with moss in percent pitcher-wise
## Prokaryotes
tax.l <- data.frame(t(otu_table(prok.a)))
tax.m <- data.frame(t(otu_table(moss.prok.a)))

pitchers <- colnames(tax.l)
mosses <- colnames(tax.m)


## Empty dataframe to populate
overlaps.p <- data.frame(Pitcher = rep(NA_character_, length(pitchers)),
                         Moss = rep(NA_character_, length(pitchers)),
                         Overlap = rep(NA_real_, length(pitchers)))


k <- 1
# i <- 1; j <- 1
for (i in 1:length(mosses)){
  site.sector <- gsub(".{4}$", "", mosses[i])
  print(site.sector)

  site.sector.l <- grep(site.sector, pitchers, value = TRUE)

  # print(site.sector.l)
  for (j in 1:length(site.sector.l)){
    # print(site.sector.l[j])
    pitch <- row.names(tax.l)[which(tax.l[[site.sector.l[j]]] > 0)]
    moss <- row.names(tax.m)[which(tax.m[[mosses[i]]] > 0)]

    overlap <- 100 / length(pitch) * length(intersect(pitch, moss))
    print(sprintf("%s overlaps with %s %f", site.sector.l[j], mosses[i],
                  overlap))

    ## Fill dataframe
    overlaps.p[k, ] <- c(site.sector.l[j], mosses[i], as.numeric(overlap))
    print(k)
    k <- k + 1
  }
}

overlaps.p$Domain <- "Prokaryotes"


## Eukaryote
tax.l <- data.frame(t(otu_table(euk.a)))
tax.m <- data.frame(t(otu_table(moss.euk.a)))

pitchers <- colnames(tax.l)
mosses <- colnames(tax.m)


## Empty dataframe to populate
overlaps.e <- data.frame(Pitcher = rep(NA_character_, length(pitchers)),
                         Moss = rep(NA_character_, length(pitchers)),
                         Overlap = rep(NA_real_, length(pitchers)))


k <- 1
# i <- 1; j <- 1
for (i in 1:length(mosses)){
  site.sector <- gsub(".{4}$", "", mosses[i])
  print(site.sector)

  site.sector.l <- grep(site.sector, pitchers, value = TRUE)

  # print(site.sector.l)
  for (j in 1:length(site.sector.l)){
    # print(site.sector.l[j])
    pitch <- row.names(tax.l)[which(tax.l[[site.sector.l[j]]] > 0)]
    moss <- row.names(tax.m)[which(tax.m[[mosses[i]]] > 0)]

    overlap <- 100 / length(pitch) * length(intersect(pitch, moss))
    print(sprintf("%s overlaps with %s %f", site.sector.l[j], mosses[i],
                  overlap))

    ## Fill dataframe
    overlaps.e[k, ] <- c(site.sector.l[j], mosses[i], as.numeric(overlap))
    print(k)
    k <- k + 1
  }
}


overlaps.e$Domain <- "Eukaryotes"


## Merge
overlaps <- rbind(overlaps.p, overlaps.e)


overlaps$Overlap <- as.numeric(overlaps$Overlap)
overlaps$Site <- gsub(".{5}$", "", overlaps$Pitcher)
overlaps$Site <- ordered(overlaps$Site,
                         levels = c("CB", "LT", "LE", "LV", "LM"))
overlaps$Sector <- substr(overlaps$Pitcher, 3, 3)
overlaps$Succession <- as.factor(substr(overlaps$Pitcher, 7, 7))
levels(overlaps$Succession) <- c("Early", "Late")
overlaps$Domain <- ordered(overlaps$Domain,
                           levels = c("Prokaryotes", "Eukaryotes"))


ggplot(overlaps) +
  geom_boxplot(aes(x = Succession, y = Overlap), size = 0.2,
               outlier.shape = NA, color = "black") +
  geom_jitter(aes(x = Succession, y = Overlap, group = Sector, fill = Site),
              size = 2, height = 0, width = 0.3, pch = 21) +
  facet_grid(Domain ~ Site, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = gradCol) +
  theme(legend.position = "none", legend.direction = "horizontal") +
  ylim(0, 100) +
  xlab("") +
  ylab("Overlap pitcher- with moss taxa [%]")
# ggsave("SMP_Overlap.pdf", width = 11.69, height = 6)


## Get statistics
aggregate(Overlap ~ Site + Domain, data = overlaps, FUN = range)
aggregate(Overlap ~ Site + Domain, data = overlaps, FUN = mean)


aggregate(Overlap ~ Domain, data = overlaps, FUN = range)
aggregate(Overlap ~ Domain, data = overlaps, FUN = mean)
aggregate(Overlap ~ Domain, data = overlaps, FUN = sd)


################################################################################
### Alluvial diagrams = taxonomic overlap
## Pitcher and moss
# smp.moss <- merge_phyloseq(prok.a, moss.prok.a) # Prokaryotes
smp.moss <- merge_phyloseq(euk.a, moss.euk.a) # Eukaryotes


tax.sm <- as.data.frame(smp.moss@tax_table@.Data)
otu.sm <- as.data.frame(t(otu_table(smp.moss)))

tax.otu.sm <- base:::merge(tax.sm, otu.sm, by = 0, all = TRUE)
rownames(tax.otu.sm) <- tax.otu.sm$Row.names
tax.otu.sm$Row.names <- NULL

tax.otu.sm <- droplevels(tax.otu.sm)
tax.otu.sm$Genus <- as.character(tax.otu.sm$Genus)
str(tax.otu.sm)


rm(tax.sm, otu.sm)


## Pitcher and moss taxa
# mb.pitcher <- colnames(otu_table(prok.a))
# mb.moss <- colnames(otu_table(moss.prok.a))
mb.pitcher <- colnames(otu_table(euk.a))
mb.moss <- colnames(otu_table(moss.euk.a))

shared <- intersect(mb.pitcher, mb.moss)
only.p <- setdiff(mb.pitcher, mb.moss)
only.m <- setdiff(mb.moss, mb.pitcher)


## Rewrite NAs
tax.otu.sm$Class[is.na(tax.otu.sm$Class)]  <- "Unknown class"
tax.otu.sm$Genus[is.na(tax.otu.sm$Genus)]  <- "Unknown genus"

## Sort by class
tax.otu.sm <- tax.otu.sm[order(tax.otu.sm$Class), ]


tax.otu.sm$Freq <- 1
tax.otu.sm$OTU <- rownames(tax.otu.sm)


tax.otu.sm$Habitat <- NA
tax.otu.sm[tax.otu.sm$OTU %in% shared, ]$Habitat <- "Pitcher & moss"
tax.otu.sm[tax.otu.sm$OTU %in% only.p, ]$Habitat <- "Pitcher"
tax.otu.sm[tax.otu.sm$OTU %in% only.m, ]$Habitat <- "Moss"


## Exclude mosses
tax.otu.sm <- tax.otu.sm[!(tax.otu.sm$Habitat == "Moss"), ]


## Rename as they match with Order levels otherwise
tax.otu.sm$Genus <- gsub("unclassified", "unclassified_ge", tax.otu.sm$Genus)
# intersect(unique(tax.otu.sm$Genus), unique(tax.otu.sm$Order))

tax.otu.sm$Genus <- factor(tax.otu.sm$Genus,
                           levels = as.character(unique(tax.otu.sm$Genus)))


ggplot(data = tax.otu.sm,
       aes(axis1 = Class, axis2 = Genus, axis3 = Habitat, y = Freq)) +
  stat_alluvium(aes(fill = Habitat)) +
  geom_stratum(linetype = 1, lwd = 0.01) +
  ## prok
  # geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 0.85) +
  ## euk
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 1.8) +
  scale_fill_manual(values = c("#D55E00", "#009E73")) +
  theme_void() +
  theme(legend.position = "none")
# ggsave("SMP_16S_Alluvial_ASV.pdf", width = 8.27, height = 11.69)
ggsave("SMP_18S_Alluvial_ASV.pdf", width = 8.27, height = 11.69)


################################################################################
### Venn diagrams for site-level overlap comparison
## Pitcher and moss taxa
mb.pitcher <- colnames(otu_table(prok.a))
mb.moss <- colnames(otu_table(moss.prok.a))
# mb.pitcher <- colnames(otu_table(euk.a))
# mb.moss <- colnames(otu_table(moss.euk.a))


## Venn global and site-wise
grid.newpage()
vennAll <-
  draw.pairwise.venn(area1 = length(mb.pitcher),
                     area2 = length(mb.moss),
                     cross.area = length(intersect(mb.pitcher, mb.moss)),
                     category = c("",
                                  sprintf("%s (%s) %%", "All sites",
                                          format(round(100 /
                                                         length(mb.pitcher) *
                                                         length(intersect(
                                                           mb.pitcher,
                                                           mb.moss)),
                                                       digits = 1),
                                                 nsmall = 1))),
                     lty = rep("blank", 2), fill = c("#D55E00", "#009E73"),
                     alpha = rep(0.4, 2), cat.pos = c(180, 180),
                     cex = rep(1.8, 3), cat.cex = 2.2, cat.dist = rep(0.025, 2),
                     scaled = TRUE)
dev.off()


################################################################################
## Site-wise
smp.a <- prok.a
moss.a <- moss.prok.a
# smp.a <- euk.a
# moss.a <- moss.euk.a


sites <- levels(smp.a@sam_data$Site)

for (i in 1:length(sites)) {
  ## Get pitcher mb for site i
  smp.site <- prune_samples(smp.a@sam_data$Site == sites[i], smp.a)
  smp.site <- prune_taxa(taxa_sums(smp.site) > 0, smp.site)

  ## Get moss mb for site i
  moss.site <- prune_samples(moss.a@sam_data$Site == sites[i], moss.a)
  moss.site <- prune_taxa(taxa_sums(moss.site) > 0, moss.site)

  ## Get OTUs
  mb.pitcher <- colnames(otu_table(smp.site))
  mb.moss <- colnames(otu_table(moss.site))

  ## Draw the venns
  assign(paste0("venn", sites[i]),
         draw.pairwise.venn(area1 = length(mb.pitcher),
                            area2 = length(mb.moss),
                            cross.area = length(intersect(mb.pitcher, mb.moss)),
                            category = c("",
                                         sprintf("%s (%s %%)", sites[i],
                                                 format(round(100 /
                                                                length(
                                                                  mb.pitcher) *
                                                                length(
                                                                  intersect(
                                                                    mb.pitcher,
                                                                    mb.moss)),
                                                              digits = 1),
                                                        nsmall = 1))),
                            lty = rep("blank", 2),
                            fill = c("#D55E00", "#009E73"),
                            alpha = rep(0.4, 2), cat.pos = c(180, 180),
                            cex = rep(1.8, 3), cat.cex = 2.2,
                            cat.dist = rep(0.025, 2), scaled = TRUE))
}


pdf("SMP_venn_16S_Sitewise.pdf", height = 7, width = 8.27)
# pdf("SMP_venn_18S_Sitewise.pdf", height = 7, width = 8.27)
grid.arrange(grobTree(vennAll), grobTree(vennCB), grobTree(vennLT),
             grobTree(vennLE), grobTree(vennLV), grobTree(vennLM), ncol = 2)
dev.off()


################################################################################
################################################################################
################################################################################
################################################################################
