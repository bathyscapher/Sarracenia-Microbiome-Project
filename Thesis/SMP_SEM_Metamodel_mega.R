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
### Full model
### Composite "microhabitat"
full.micro <-
'mb.leaf.prok ~ PotentialVolume_mL + Age_d + prey.items + SampleColorDummy + pH
mb.leaf.euk ~ PotentialVolume_mL + Age_d + prey.items + SampleColorDummy + pH'


fit.full.micro <- sem(full.micro, data = smpMeta.df)
summary(fit.full.micro, rsq = TRUE)


## Extract coefficients
PotentialVolume.p <- coef(fit.full.micro)[[1]]
Age_d.p  <- coef(fit.full.micro)[[2]]
prey.items.p  <- coef(fit.full.micro)[[3]]
SampleColorDummy.p  <- coef(fit.full.micro)[[4]]
pH.p <- coef(fit.full.micro)[[5]]


## Extract coefficients
PotentialVolume.e <- coef(fit.full.micro)[[6]]
Age_d.e <- coef(fit.full.micro)[[7]]
prey.items.e <- coef(fit.full.micro)[[8]]
SampleColorDummy.e  <- coef(fit.full.micro)[[9]]
pH.e <- coef(fit.full.micro)[[10]]

coef(fit.full.micro); PotentialVolume.p; Age_d.p; prey.items.p; SampleColorDummy.p; pH.p; PotentialVolume.e; Age_d.e; prey.items.e; SampleColorDummy.e; pH.e


## Compute composite "microhabitat"
smpMeta.df$full.comp.micro.p <- Age_d.p * smpMeta.df$Age_d +
  PotentialVolume.p * smpMeta.df$PotentialVolume_mL +
  prey.items.p * smpMeta.df$prey.items +
  SampleColorDummy.p * smpMeta.df$SampleColorDummy +
  pH.p * smpMeta.df$pH

smpMeta.df$full.comp.micro.e <- Age_d.e * smpMeta.df$Age_d +
  PotentialVolume.e * smpMeta.df$PotentialVolume_mL +
  prey.items.e * smpMeta.df$prey.items +
  SampleColorDummy.e * smpMeta.df$SampleColorDummy +
  pH.e * smpMeta.df$pH


rm(PotentialVolume.p, Age_d.p, prey.items.p, SampleColorDummy.p, pH.p,
   PotentialVolume.e, Age_d.e, prey.items.e, SampleColorDummy.e, pH.e)


### Run SEM with the manually computed composite
full.comp.micro <-
'mb.leaf.prok ~ full.comp.micro.p
mb.leaf.euk ~ full.comp.micro.e'
fit.full.micro.comp <- sem(full.comp.micro, data = smpMeta.df)


## Control if the standard coefficient and z-value (standardised coefficient)
## are similar
summary(fit.full.micro.comp, rsq = TRUE, standardized = TRUE)
summary(fit.full.micro, rsq = TRUE, standardized = TRUE)


rm(fit.full.micro.comp, full.micro, full.comp.micro)


### Composite "macrohabitat": estimate separate composites for the responses
full.macro <-
'mb.leaf.prok ~ Distance + MeanTemperature + Altitude + CanopyCover
mb.moss.prok ~ Distance + MeanTemperature + Altitude + CanopyCover

mb.leaf.euk ~ Distance + MeanTemperature + Altitude + CanopyCover
mb.moss.euk ~ Distance + MeanTemperature + Altitude + CanopyCover

prey ~ Distance + MeanTemperature + Altitude + CanopyCover'

fit.full.macro <- sem(full.macro, data = smpMeta.df)
summary(fit.full.macro, rsq = TRUE, fit.measures = TRUE)


## Extract coefficients
## For the pitchers/leaves
Distance.lp <- coef(fit.full.macro)[[1]]
meanTemperature.lp  <- coef(fit.full.macro)[[2]]
Altitude.lp  <- coef(fit.full.macro)[[3]]
CanopyCover.lp  <- coef(fit.full.macro)[[4]]

coef(fit.full.macro)[1:4]; Distance.lp; meanTemperature.lp; Altitude.lp; CanopyCover.lp


## For the mosses
Distance.mp <- coef(fit.full.macro)[[5]]
meanTemperature.mp  <- coef(fit.full.macro)[[6]]
Altitude.mp <- coef(fit.full.macro)[[7]]
CanopyCover.mp <- coef(fit.full.macro)[[8]]

coef(fit.full.macro)[5:8]; Distance.mp; meanTemperature.mp; Altitude.mp; CanopyCover.mp


## For the pitchers/leaves
Distance.le <- coef(fit.full.macro)[[9]]
meanTemperature.le  <- coef(fit.full.macro)[[10]]
Altitude.le  <- coef(fit.full.macro)[[11]]
CanopyCover.le  <- coef(fit.full.macro)[[12]]

coef(fit.full.macro)[9:12]; Distance.le; meanTemperature.le; Altitude.le; CanopyCover.le


## For the mosses
Distance.me <- coef(fit.full.macro)[[13]]
meanTemperature.me  <- coef(fit.full.macro)[[14]]
Altitude.me <- coef(fit.full.macro)[[15]]
CanopyCover.me <- coef(fit.full.macro)[[16]]

coef(fit.full.macro)[13:16]; Distance.me; meanTemperature.me; Altitude.me; CanopyCover.me


## For the prey
Distance.p <- coef(fit.full.macro)[[17]]
meanTemperature.p  <- coef(fit.full.macro)[[18]]
Altitude.p  <- coef(fit.full.macro)[[19]]
CanopyCover.p  <- coef(fit.full.macro)[[20]]

coef(fit.full.macro)[17:20]; Distance.p; meanTemperature.p; Altitude.p; CanopyCover.p


### Compute composites "macrohabitat" for the leaves, the prey and the moss
## For the pitchers/leaves
smpMeta.df$full.comp.macro.lp <- Distance.lp * smpMeta.df$Distance +
  meanTemperature.lp * smpMeta.df$MeanTemperature +
  Altitude.lp * smpMeta.df$Altitude +
  CanopyCover.lp * smpMeta.df$CanopyCover

smpMeta.df$full.comp.macro.le <- Distance.le * smpMeta.df$Distance +
  meanTemperature.le * smpMeta.df$MeanTemperature +
  Altitude.le * smpMeta.df$Altitude +
  CanopyCover.le * smpMeta.df$CanopyCover

## For the moss
smpMeta.df$full.comp.macro.mp <- Distance.mp * smpMeta.df$Distance +
  meanTemperature.mp * smpMeta.df$MeanTemperature +
  Altitude.mp * smpMeta.df$Altitude +
  CanopyCover.mp * smpMeta.df$CanopyCover

smpMeta.df$full.comp.macro.me <- Distance.me * smpMeta.df$Distance +
  meanTemperature.me * smpMeta.df$MeanTemperature +
  Altitude.me * smpMeta.df$Altitude +
  CanopyCover.me * smpMeta.df$CanopyCover

## For the prey
smpMeta.df$full.comp.macro.prey <- Distance.p * smpMeta.df$Distance +
  meanTemperature.p * smpMeta.df$MeanTemperature +
  Altitude.p * smpMeta.df$Altitude +
  CanopyCover.p * smpMeta.df$CanopyCover


rm(Distance.lp, meanTemperature.lp, Altitude.lp, CanopyCover.lp, Distance.mp,
   meanTemperature.mp, Altitude.mp, CanopyCover.mp, Distance.le,
   meanTemperature.le, Altitude.le, CanopyCover.le, Distance.me,
   meanTemperature.me, Altitude.me, CanopyCover.me, Distance.p,
   meanTemperature.p, Altitude.p,  CanopyCover.p)


### SEM with composite macrohabitat
full.comp.macro <-
'mb.leaf.prok ~ full.comp.macro.lp
# mb.moss.prok ~ full.comp.macro.mp
mb.leaf.euk ~ full.comp.macro.le
mb.moss.euk ~ full.comp.macro.me
prey ~ full.comp.macro.prey'

fit.full.macro.comp <- sem(full.comp.macro, data = smpMeta.df)


## Control if the standard coefficient and z-value (standardised coefficient)
## are similar
summary(fit.full.macro.comp, rsq = TRUE, fit.measures = TRUE)
summary(fit.full.macro, rsq = TRUE, fit.measures = TRUE)


rm(fit.full.macro.comp, full.macro, full.comp.macro)


### Run mega SEM
# mega <-
# 'full.comp.micro.prok ~ full.comp.macro.lp
# full.comp.micro.euk ~ full.comp.macro.le
#
# mb.moss.prok ~ full.comp.macro.mp
# mb.moss.euk ~ full.comp.macro.me
# #
# # prey ~ full.comp.macro.prey
# #
# mb.leaf.prok ~ full.comp.micro.prok + mb.moss.prok + full.comp.macro.lp + prey
# mb.leaf.euk ~ full.comp.micro.euk + mb.moss.euk + full.comp.macro.le + prey
#
# mb.leaf.prok ~~ mb.leaf.euk
# full.comp.micro.prok ~ full.comp.micro.euk
# full.comp.micro.euk ~ full.comp.micro.prok
# '


### Simple model
mega <- '
prey.items ~ Age_d + prey
prey ~ Age_d + SampleVolume_mL + prey.items
pH ~ mb.leaf.prok + mb.leaf.euk + Age_d + SampleVolume_mL + prey.items
SampleColorDummy ~ mb.leaf.prok + mb.leaf.euk + prey + Age_d

prey ~ Distance + Altitude + MeanTemperature + CanopyCover
mb.moss.prok ~ Distance + Altitude + MeanTemperature + CanopyCover
mb.moss.euk ~ Distance + Altitude + MeanTemperature + CanopyCover

mb.leaf.prok ~ Distance + Altitude + MeanTemperature + CanopyCover +
  SampleVolume_mL + Age_d + SampleColorDummy + prey + prey.items + pH +
  mb.moss.euk #+ mb.leaf.euk

mb.leaf.euk ~ Distance + Altitude + MeanTemperature + CanopyCover +
  SampleVolume_mL + Age_d + SampleColorDummy + prey + prey.items + pH +
  mb.moss.euk #+ mb.leaf.prok

# prey ~~ prey.items
# mb.leaf.prok ~~ mb.leaf.euk
'


fit.mega <- sem(mega, data = smpMeta.df)
summary(fit.mega, rsq = TRUE, fit.measures = TRUE)


modindices(fit.mega, minimum.value = 3)


## Update model
fit.mega.up <- update(fit.mega, add = "prey ~ prey.items")
summary(fit.mega.up, fit.measures = TRUE, rsq = TRUE)

modificationindices(fit.mega.up, minimum.value = 3)


## Update model 2
fit.mega.up2 <- update(fit.mega.up, add = "prey ~ Age_d")
summary(fit.mega.up2, fit.measures = TRUE, rsq = TRUE)

modificationindices(fit.mega.up2, minimum.value = 3)


### Plot bestest model
pdf(file = "SMP_micro_prok.pdf", height = 8.27, width = 11.69)
semPaths(fit.mega, what = "std", whatLabels = "std", residuals = TRUE,
         sizeMan = 12, sizeMan2 = 8, edge.label.cex = 0.8,
         fade = FALSE, layout = "spring", style = "mx", nCharNodes = 0,
         rotation = 2, layoutSplit = TRUE, curve = 2, curvature = 3)
dev.off()


pdf(file = "SMP_macro_prok.pdf", height = 8, width = 11)
semPaths(fit.full.macro, what = "std", whatLabels = "std", residuals = TRUE,
         sizeMan = 12, sizeMan2 = 8, edge.label.cex = 0.8,
         fade = FALSE, layout = "tree", style = "mx", nCharNodes = 0,
         rotation = 2, layoutSplit = TRUE, curve = 2, curvature = 2)
dev.off()


pdf(file = "SMP_full_prok.pdf", height = 8.27, width = 11.69)
semPaths(fit.mega.up2, what = "std", whatLabels = "std", residuals = TRUE,
         sizeMan = 12, sizeMan2 = 8, edge.label.cex = 0.8,
         fade = FALSE, layout = "tree2", style = "mx", nCharNodes = 0,
         rotation = 3, layoutSplit = TRUE, curve = 2, curvature = 2)
dev.off()


################################################################################
################################################################################
################################################################################
################################################################################
