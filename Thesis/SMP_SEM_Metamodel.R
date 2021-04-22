################################################################################
################################################################################
################################################################################
################################################################################
### SMP Structural equation modeling
### Author: korn@cumulonimbus.at University of Fribourg 2020
################################################################################


library("GGally")
library("lavaan")
library("semPlot")
library("blavaan")


rm(list = ls())
setwd("~/Sarracenia-Microbiome-Project/Thesis")


set.seed(34706)


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


smpMeta.df <- read.table("csv/SMP_Metadata_SEM_scaled.csv", sep = "\t")


dim(smpMeta.df) # should be a data.frame with 160 rows and 29 columns
colnames(smpMeta.df)
# smpMeta.df$FullID # 160 sample names


## ln+1-transform prey.items to remove the curvy relationship with prey
smpMeta.df$prey.items <- log1p(smpMeta.df$prey.items)


################################################################################
################################################################################
### Model with composite variables
### Composite "microhabitat"
micro <-
'mb.leaf.prok ~ SampleVolume_mL + Age_d + SampleColorDummy + pH + prey.items
mb.leaf.euk ~ SampleVolume_mL + Age_d + SampleColorDummy + pH + prey.items'


fit.micro <- sem(micro, data = smpMeta.df)
summary(fit.micro, rsq = TRUE)


## Extract coefficients
SampleVolume.lp <- coef(fit.micro)[[1]]
Age_d.lp <- coef(fit.micro)[[2]]
SampleColorDummy.lp <- coef(fit.micro)[[3]]
pH.lp <- coef(fit.micro)[[4]]
prey.items.lp <- coef(fit.micro)[[5]]


coef(fit.micro)[1:5]; SampleVolume.lp; Age_d.lp; SampleColorDummy.lp; pH.lp; prey.items.lp


SampleVolume.le <- coef(fit.micro)[[6]]
Age_d.le <- coef(fit.micro)[[7]]
SampleColorDummy.le <- coef(fit.micro)[[8]]
pH.le <- coef(fit.micro)[[9]]
prey.items.le <- coef(fit.micro)[[10]]

coef(fit.micro)[6:10]; SampleVolume.le; Age_d.le; SampleColorDummy.le; pH.le; prey.items.le


## Compute composite "microhabitat"
smpMeta.df$comp.micro.euk <- SampleVolume.le * smpMeta.df$SampleVolume_mL +
  Age_d.le * smpMeta.df$Age_d +
  prey.items.le * smpMeta.df$prey.items +
  SampleColorDummy.le * smpMeta.df$SampleColorDummy +
  pH.le * smpMeta.df$pH


smpMeta.df$comp.micro.prok <- SampleVolume.lp * smpMeta.df$SampleVolume_mL +
  Age_d.lp * smpMeta.df$Age_d +
  prey.items.lp * smpMeta.df$prey.items +
  SampleColorDummy.lp * smpMeta.df$SampleColorDummy +
  pH.lp * smpMeta.df$pH


rm(SampleVolume.le, Age_d.le, prey.items.le, SampleColorDummy.le, pH.le,
   SampleVolume.lp, Age_d.lp, prey.items.lp, SampleColorDummy.lp, pH.lp)


### Run SEM with the manually computed composite
comp.micro <-
'mb.leaf.prok ~ comp.micro.prok
mb.leaf.euk ~ comp.micro.euk'


fit.micro.comp <- sem(comp.micro, data = smpMeta.df)


## Control if the standard coefficient and z-value (standardised coefficient)
## are similar
summary(fit.micro.comp, rsq = TRUE, standardized = TRUE)
summary(fit.micro, rsq = TRUE, standardized = TRUE)


rm(fit.micro.comp, micro, comp.micro)


### Composite "macrohabitat": estimate separate composites for the responses
macro <-
'mb.leaf.prok ~ Distance + MeanTemperature + Altitude + CanopyCover
mb.leaf.euk ~ Distance + MeanTemperature + Altitude + CanopyCover'


fit.macro <- sem(macro, data = smpMeta.df)
summary(fit.macro, rsq = TRUE, fit.measures = TRUE)


## Extract coefficients
## For the pitchers/leaves
Distance.lp <- coef(fit.macro)[[1]]
meanTemperature.lp <- coef(fit.macro)[[2]]
Altitude.lp <- coef(fit.macro)[[3]]
CanopyCover.lp <- coef(fit.macro)[[4]]

coef(fit.macro)[1:4]; Distance.lp; meanTemperature.lp; Altitude.lp; CanopyCover.lp


Distance.le <- coef(fit.macro)[[5]]
meanTemperature.le <- coef(fit.macro)[[6]]
Altitude.le <- coef(fit.macro)[[7]]
CanopyCover.le <- coef(fit.macro)[[8]]

coef(fit.macro)[5:8]; Distance.le; meanTemperature.le; Altitude.le; CanopyCover.le


### Compute composites "macrohabitat" for the leaves, the prey and the moss
## For the pitchers/leaves
smpMeta.df$comp.macro.le <- Distance.le * smpMeta.df$Distance +
  meanTemperature.le * smpMeta.df$MeanTemperature +
  Altitude.le * smpMeta.df$Altitude +
  CanopyCover.le * smpMeta.df$CanopyCover

smpMeta.df$comp.macro.lp <- Distance.lp * smpMeta.df$Distance +
  meanTemperature.lp * smpMeta.df$MeanTemperature +
  Altitude.lp * smpMeta.df$Altitude +
  CanopyCover.lp * smpMeta.df$CanopyCover


rm(Distance.le, meanTemperature.le, Altitude.le, CanopyCover.le, Distance.lp,
   meanTemperature.lp, Altitude.lp, CanopyCover.lp)


### SEM with composite macrohabitat
comp.macro <-
'mb.leaf.prok ~ comp.macro.lp
mb.leaf.euk ~ comp.macro.le'


fit.macro.comp <- sem(comp.macro, data = smpMeta.df)


## Control if the standard coefficient and z-value (standardised coefficient)
## are similar
summary(fit.macro.comp, rsq = TRUE, fit.measures = TRUE)
summary(fit.macro, rsq = TRUE, fit.measures = TRUE)


rm(macro, comp.macro)


################################################################################
check.vars <- c("Site",
                "mb.leaf.prok", "comp.macro.lp", "comp.micro.prok",
                "mb.leaf.euk", "comp.macro.le", "comp.micro.euk")
smpMeta.df.sub <- smpMeta.df[check.vars]
print(round(cor(smpMeta.df.sub[, -1], use = "complete"), digits = 1))


pairsWithAbline <- function(x, y){
  points(x, y, pch = 1, col = "black")
  abline(lm(y ~ x), lty = 2, col = "red")
}


pairs(smpMeta.df.sub[, -1], panel = pairsWithAbline)
pairs(smpMeta.df.sub[, -1], panel = panel.smooth)


## With GGally
## Colorblind temperature scale for the 5 sites
gradCol <- c("#D55E00", "#E69F00", "#F0E442", "#009E73", "#56B4E9")


smpMeta.df.sub$Site <- ordered(smpMeta.df.sub$Site,
                               levels = c("CB", "LT", "LE", "LV", "LM"))


p <- ggpairs(smpMeta.df.sub, aes(colour = Site, alpha = 0.5),
             upper = list(continuous = wrap("cor", size = 1.7)),
             lower = list(continuous = wrap("points", alpha = 0.3,
                                            size = 0.4))) +
  theme_bw(base_size = 5) +
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
# ggsave("SMP_SEM_MetadataPairs.pdf", width = 8.27, height = 8.27)

rm(check.vars, pairsWithAbline, smpMeta.df.sub)


################################################################################
### Run SEM with the two composites
comp <- '
mb.leaf.prok ~ comp.micro.prok + mb.moss.prok + comp.macro.lp + prey
mb.leaf.euk ~ comp.micro.euk + mb.moss.euk + comp.macro.le + prey

mb.leaf.euk ~ mb.leaf.prok
mb.leaf.prok ~ mb.leaf.euk'


fit.comp <- sem(comp, data = smpMeta.df, estimator = "MLM")
summary(fit.comp, rsq = TRUE, fit.measures = TRUE)


resid(fit.comp, "cor") # misfit of the bivariate associations
modindices(fit.comp, minimum.value = 3, op = "~")


fit.comp.b <- bsem(comp, data = smpMeta.df,
                        n.chains = 4, bcontrol = list(cores = 6))
summary(fit.comp.b, rsq = TRUE, fit.measures = TRUE)


################################################################################
### Plot model
pdf(file = "SMP_SEM_comp.pdf", height = 8.27, width = 6)
semPaths(fit.comp.b, what = "est", whatLabels = "est", residuals = FALSE,
         intercepts = FALSE, sizeMan = 5, sizeMan2 = 3, edge.label.cex = 0.4,
         fade = FALSE, layout = "tree", style = "mx", nCharNodes = 0,
         posCol = "#009e73ff", negCol = "#d55e00ff", edge.label.color = "black",
         layoutSplit = TRUE, curve = 1, curvature = 1, fixedStyle = 1,
         exoCov = FALSE, rotation = 1)
dev.off()


### Plot composites
pdf(file = "SMP_micro.pdf", height = 8.27, width = 6)
semPaths(fit.micro, what = "est", whatLabels = "est", residuals = FALSE,
         intercepts = FALSE, sizeMan = 5, sizeMan2 = 3, edge.label.cex = 0.4,
         fade = FALSE, layout = "tree", style = "mx", nCharNodes = 0,
         posCol = "#009e73ff", negCol = "#d55e00ff", edge.label.color = "black",
         layoutSplit = TRUE, curve = 1, curvature = 1, fixedStyle = 1,
         exoCov = FALSE, rotation = 1)
dev.off()


pdf(file = "SMP_macro.pdf", height = 8.27, width = 6)
semPaths(fit.macro, what = "est", whatLabels = "est", residuals = FALSE,
         intercepts = FALSE, sizeMan = 5, sizeMan2 = 3, edge.label.cex = 0.4,
         fade = FALSE, layout = "tree", style = "mx", nCharNodes = 0,
         posCol = "#009e73ff", negCol = "#d55e00ff", edge.label.color = "black",
         layoutSplit = TRUE, curve = 1, curvature = 1, fixedStyle = 1,
         exoCov = FALSE, rotation = 1)
dev.off()


################################################################################
################################################################################
################################################################################
### Without composite variables, but linking the individual variables with each
### other if justified

# mega <- '
# prey.items ~ Age_d + prey
# prey ~ Age_d + prey.items + SampleVolume_mL
# pH ~ mb.leaf.prok + mb.leaf.euk + Age_d + SampleVolume_mL + prey.items + prey
# SampleColorDummy ~ mb.leaf.prok + mb.leaf.euk + prey + Age_d
#
# prey ~ Distance + Altitude + MeanTemperature + CanopyCover
# mb.moss.prok ~ Distance + Altitude + MeanTemperature + CanopyCover
# mb.moss.euk ~ Distance + Altitude + MeanTemperature + CanopyCover
#
# mb.leaf.prok ~ Distance + Altitude + MeanTemperature + CanopyCover +
#   SampleVolume_mL + Age_d + SampleColorDummy + prey + prey.items + pH +
#   mb.moss.prok #+ mb.leaf.euk
#
# mb.leaf.euk ~ Distance + Altitude + MeanTemperature + CanopyCover +
#   SampleVolume_mL + Age_d + SampleColorDummy + prey + prey.items + pH +
#   mb.moss.euk #+ mb.leaf.prok
# '
#
#
# fit.mega <- sem(mega, data = smpMeta.df, estimator = "MLM")
# summary(fit.mega, rsq = TRUE, fit.measures = TRUE)
#
#
# modindices(fit.mega, minimum.value = 6)
#
#
# ### Plot bestest model
# pdf(file = "SMP_mega.pdf", height = 8.27, width = 11.69)
# semPaths(fit.mega, what = "est", whatLabels = "est", residuals = FALSE,
#          intercepts = FALSE, sizeMan = 5, sizeMan2 = 3, edge.label.cex = 0.4,
#          fade = FALSE, layout = "tree2", style = "mx", nCharNodes = 0,
#          posCol = "#009e73ff", negCol = "#d55e00ff", edge.label.color = "black",
#          layoutSplit = TRUE, curve = 1, curvature = 1, fixedStyle = 1,
#          exoCov = FALSE, rotation = 1)
# dev.off()



################################################################################
################################################################################
################################################################################
################################################################################
