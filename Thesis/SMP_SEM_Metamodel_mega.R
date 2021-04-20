################################################################################
################################################################################
################################################################################
################################################################################
### SMP SEM Metamodel v3: full model as the bestest model
### Author: korn@cumulonimbus.at University of Fribourg 2020
################################################################################


library("lavaan")
library("semPlot")
# library("blavaan")
# library("rjags")
# library("AICcmodavg")
# library("piecewiseSEM")

# sessionInfo()


rm(list = ls())
setwd("~/Sarracenia-Microbiome-Project/Thesis")


set.seed(34706)


## Read data
smpMeta.df <- read.table("csv/SMP_Metadata_SEM_scaled.csv", sep = "\t")


dim(smpMeta.df) # should be a data.frame with 160 rows and 29 columns
colnames(smpMeta.df)
# smpMeta.df$FullID


## ln+1-transform prey.items to remove the curvy relationship with prey
smpMeta.df$prey.items <- log1p(smpMeta.df$prey.items)


################################################################################
### Microhabitat
## Variables:
## * PitcherLength = proxy for habitat size
## * Age_d = age of the pitcher in days
## * prey.items = number of prey items
## * prey = 1st nMDS axis of Euclidean distance of prey on order level
## * SampleColorDummy = color scale as proxy for nutrient content
## * pH = pH
## Response:
## * mb.leaf.prok = 1st nMDS axis of Bray-Curtis distance of prokaryotes from
##   the pitchers
## * mb.moss.prok = 1st nMDS axis of Bray-Curtis distance of prokaryotes from
##   the mosses

### Macrohabitat
## Variables:
## * Distance = geographical distance between plants
## * MeanTemperature = mean temperature measured during field season
## * Altitude
## * CanopyCover = canopy cover in %
## Response:
## * prey = 1st nMDS axis of Euclidean distance of prey on order level
## * mb.leaf.prok = 1st nMDS axis of Bray-Curtis distance of prokaryotes from
##   the pitchers
## * mb.moss.prok = 1st nMDS axis of Bray-Curtis distance of prokaryotes from
##   the moss


################################################################################
################################################################################
### Run mega SEM
mega <- '
prey.items ~ Age_d + prey
prey ~ Age_d + prey.items + SampleVolume_mL
pH ~ mb.leaf.prok + mb.leaf.euk + Age_d + SampleVolume_mL + prey.items + prey
SampleColorDummy ~ mb.leaf.prok + mb.leaf.euk + prey + Age_d

prey ~ Distance + Altitude + MeanTemperature + CanopyCover
mb.moss.prok ~ Distance + Altitude + MeanTemperature + CanopyCover
mb.moss.euk ~ Distance + Altitude + MeanTemperature + CanopyCover

mb.leaf.prok ~ Distance + Altitude + MeanTemperature + CanopyCover +
  SampleVolume_mL + Age_d + SampleColorDummy + prey + prey.items + pH +
  mb.moss.prok #+ mb.leaf.euk

mb.leaf.euk ~ Distance + Altitude + MeanTemperature + CanopyCover +
  SampleVolume_mL + Age_d + SampleColorDummy + prey + prey.items + pH +
  mb.moss.euk #+ mb.leaf.prok
'


fit.mega <- sem(mega, data = smpMeta.df, estimator = "MLM")
summary(fit.mega, rsq = TRUE, fit.measures = TRUE)


modindices(fit.mega, minimum.value = 6)


### Plot bestest model
pdf(file = "SMP_mega.pdf", height = 8.27, width = 11.69)
semPaths(fit.mega, what = "std", whatLabels = "std", residuals = TRUE,
         sizeMan = 12, sizeMan2 = 8, edge.label.cex = 0.8,
         fade = FALSE, layout = "tree2", style = "mx", nCharNodes = 0,
         rotation = 3, layoutSplit = TRUE, curve = 2, curvature = 2)
dev.off()


################################################################################
################################################################################
################################################################################
################################################################################
