################################################################################
################################################################################
################################################################################
################################################################################
### SMP SEM Metamodel v3 - lavaan + blavaan
### Author: korn@cumulonimbus.at University of Fribourg 2020
################################################################################


library("GGally")
library("lavaan")
library("semPlot")
library("blavaan")


rm(list = ls())
setwd("~/Sarracenia-Microbiome-Project/Thesis")


set.seed(34706)


smpMeta.df <- read.table("csv/SMP_Metadata_SEM_scaled.csv", sep = "\t")


dim(smpMeta.df) # should be a data.frame with 160 rows and 29 columns
colnames(smpMeta.df)
# smpMeta.df$FullID


## ln+1-transform prey.items to remove the curvy relationship with prey
smpMeta.df$prey.items <- log1p(smpMeta.df$prey.items)


################################################################################
################################################################################
### Full model
### Composite "microhabitat"
micro <-
'mb.leaf.prok ~ SampleVolume_mL + Age_d + SampleColorDummy + pH + prey.items
mb.leaf.euk ~ SampleVolume_mL + Age_d + SampleColorDummy + pH + prey.items
prey ~ SampleVolume_mL + Age_d + SampleColorDummy + pH + prey.items'


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


SampleVolume.p <- coef(fit.micro)[[11]]
Age_d.p <- coef(fit.micro)[[12]]
SampleColorDummy.p <- coef(fit.micro)[[13]]
pH.p <- coef(fit.micro)[[14]]
prey.items.p <- coef(fit.micro)[[15]]

coef(fit.micro)[11:15]; SampleVolume.p; Age_d.p; SampleColorDummy.p; pH.p; prey.items.p


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


smpMeta.df$comp.micro.prey <- SampleVolume.p * smpMeta.df$SampleVolume_mL +
  Age_d.p * smpMeta.df$Age_d +
  prey.items.p * smpMeta.df$prey.items +
  SampleColorDummy.p * smpMeta.df$SampleColorDummy +
  pH.p * smpMeta.df$pH


rm(SampleVolume.le, Age_d.le, prey.items.le, SampleColorDummy.le, pH.le,
   SampleVolume.lp, Age_d.lp, prey.items.lp, SampleColorDummy.lp, pH.lp,
   SampleVolume.p, Age_d.p, prey.items.p, SampleColorDummy.p, pH.p)


### Run SEM with the manually computed composite
comp.micro <-
'mb.leaf.prok ~ comp.micro.prok
mb.leaf.euk ~ comp.micro.euk
prey ~ comp.micro.prey'


fit.micro.comp <- sem(comp.micro, data = smpMeta.df)


## Control if the standard coefficient and z-value (standardised coefficient)
## are similar
summary(fit.micro.comp, rsq = TRUE, standardized = TRUE)
summary(fit.micro, rsq = TRUE, standardized = TRUE)


rm(fit.micro.comp, micro, comp.micro)


### Composite "macrohabitat": estimate separate composites for the responses
macro <-
'mb.leaf.prok ~ Distance + MeanTemperature + Altitude + CanopyCover
mb.leaf.euk ~ Distance + MeanTemperature + Altitude + CanopyCover
prey ~ Distance + MeanTemperature + Altitude + CanopyCover
mb.moss.prok ~ Distance + MeanTemperature + Altitude + CanopyCover
mb.moss.euk ~ Distance + MeanTemperature + Altitude + CanopyCover'


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


## For the prey
Distance.p <- coef(fit.macro)[[9]]
meanTemperature.p <- coef(fit.macro)[[10]]
Altitude.p <- coef(fit.macro)[[11]]
CanopyCover.p <- coef(fit.macro)[[12]]

coef(fit.macro)[9:12]; Distance.p; meanTemperature.p; Altitude.p; CanopyCover.p


## For the mosses
Distance.mp <- coef(fit.macro)[[13]]
meanTemperature.mp  <- coef(fit.macro)[[14]]
Altitude.mp <- coef(fit.macro)[[15]]
CanopyCover.mp <- coef(fit.macro)[[16]]

coef(fit.macro)[13:16]; Distance.mp; meanTemperature.mp; Altitude.mp; CanopyCover.mp


Distance.me <- coef(fit.macro)[[17]]
meanTemperature.me  <- coef(fit.macro)[[18]]
Altitude.me <- coef(fit.macro)[[19]]
CanopyCover.me <- coef(fit.macro)[[20]]

coef(fit.macro)[17:20]; Distance.me; meanTemperature.me; Altitude.me; CanopyCover.me


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

## For the prey
smpMeta.df$comp.macro.prey <- Distance.p * smpMeta.df$Distance +
  meanTemperature.p * smpMeta.df$MeanTemperature +
  Altitude.p * smpMeta.df$Altitude +
  CanopyCover.p * smpMeta.df$CanopyCover

## For the moss
smpMeta.df$comp.macro.mp <- Distance.mp * smpMeta.df$Distance +
  meanTemperature.mp * smpMeta.df$MeanTemperature +
  Altitude.mp * smpMeta.df$Altitude +
  CanopyCover.mp * smpMeta.df$CanopyCover

smpMeta.df$comp.macro.me <- Distance.me * smpMeta.df$Distance +
  meanTemperature.me * smpMeta.df$MeanTemperature +
  Altitude.me * smpMeta.df$Altitude +
  CanopyCover.me * smpMeta.df$CanopyCover


rm(Distance.le, meanTemperature.le, Altitude.le, CanopyCover.le, Distance.lp,
   meanTemperature.lp, Altitude.lp, CanopyCover.lp, Distance.p,
   meanTemperature.p, Altitude.p, Distance.mp, meanTemperature.mp, Altitude.mp,
   CanopyCover.mp, Distance.me, meanTemperature.me, Altitude.me,
   CanopyCover.me, CanopyCover.p)


### SEM with composite macrohabitat
comp.macro1 <-
'mb.leaf.prok ~ comp.macro.lp
mb.leaf.euk ~ comp.macro.le
prey ~ comp.macro.prey'

comp.macro2 <-
'mb.moss.prok ~ comp.macro.mp
mb.moss.euk ~ comp.macro.me'


fit.macro.comp1 <- sem(comp.macro1, data = smpMeta.df)
fit.macro.comp2 <- sem(comp.macro2, data = smpMeta.df)


## Control if the standard coefficient and z-value (standardised coefficient)
## are similar
summary(fit.macro.comp1, rsq = TRUE, fit.measures = TRUE)
summary(fit.macro.comp2, rsq = TRUE, fit.measures = TRUE)
summary(fit.macro, rsq = TRUE, fit.measures = TRUE)


rm(macro, comp.macro1, comp.macro2)


################################################################################
check1.vars <- c("Site",
                 "mb.leaf.prok", "comp.macro.lp", "comp.micro.prok",
                 "mb.leaf.euk", "comp.macro.le", "comp.micro.euk",
                 "mb.moss.prok", "comp.macro.mp",
                 "mb.moss.euk", "comp.macro.me",
                 "prey", "comp.macro.prey", "comp.micro.prey")
newdata2 <- smpMeta.df[check1.vars]
print(round(cor(newdata2[, -1], use = "complete"), digits = 1))


# pairsWithAbline <- function(x, y){
#   points(x, y, pch = 1, col = "black")
#   abline(lm(y ~ x), lty = 2, col = "red")
# }
#
# pairs(newdata2, panel = pairsWithAbline)
# pairs(newdata2, panel = panel.smooth)
#
# rm(check1.vars, pairsWithAbline, newdata2)


## With GGally
## Colorblind temperature scale for the 5 sites
gradCol <- c("#D55E00", "#E69F00", "#F0E442", "#009E73", "#56B4E9")


newdata2$Site <- ordered(newdata2$Site,
                         levels = c("CB", "LT", "LE", "LV", "LM"))


# p <- ggpairs(newdata2, aes(colour = Site, alpha = 0.5),
#              upper = list(continuous = wrap("cor", size = 1.7)),
#              lower = list(continuous = wrap("points", alpha = 0.3,
#                                             size = 0.4))) +
#   theme_bw(base_size = 5) +
#   theme(legend.position = "top",
#         axis.text.x = element_text(angle = 45, hjust = 1))
#
#
# for(i in 1:p$nrow) {
#   for(j in 1:p$ncol){
#     p[i,j] <- p[i,j] +
#       scale_fill_manual(values = gradCol) +
#       scale_color_manual(values = gradCol)
#   }
# }
#
# p
# ggsave("SMP_SEM_MetadataPairs.pdf", width = 8.27, height = 8.27)


################################################################################
### Run SEM with the two composites
comp <- '
# mb.moss.prok ~ comp.macro.mp
# mb.moss.euk ~ comp.macro.me
# prey ~ comp.macro.prey + comp.micro.prey

mb.leaf.prok ~ comp.micro.prok + mb.moss.prok + comp.macro.lp + prey
mb.leaf.euk ~ comp.micro.euk + mb.moss.euk + comp.macro.le + prey

mb.leaf.euk ~ mb.leaf.prok
mb.leaf.prok ~ mb.leaf.euk'


fit.comp <- sem(comp, data = smpMeta.df, estimator = "MLM")
summary(fit.comp, rsq = TRUE, fit.measures = TRUE)


resid(fit.comp, "cor") # misfit of the bivariate associations
modindices(fit.comp, minimum.value = 3, op = "~")


# ## Update model
# fit.comp.up <- update(fit.comp,
#                            add = "prey ~ comp.micro.prok")
# summary(fit.comp.up, fit.measures = TRUE, rsq = TRUE)
#
#
# resid(fit.comp.up, "cor") # misfit of the bivariate associations
# modificationindices(fit.comp.up, minimum.value = 3, op = "~")


fit.comp.b <- bsem(comp, data = smpMeta.df,
                        n.chains = 4, #burnin = 2000, sample = 2000,
                        bcontrol = list(cores = 6)#,
                        # control = list(adapt_delta = 0.9)
)
summary(fit.comp.b, rsq = TRUE, fit.measures = TRUE)


################################################################################
### Plot model
pdf(file = "SMP_SEM.pdf", height = 8.27, width = 6)
semPaths(fit.comp.b, what = "est", whatLabels = "est", residuals = FALSE,
         intercepts = FALSE, sizeMan = 5, sizeMan2 = 3, edge.label.cex = 0.35,
         fade = FALSE, layout = "tree", style = "mx", nCharNodes = 0,
         posCol = "#009e73ff", negCol = "#d55e00ff", edge.label.color = "black",
         layoutSplit = TRUE, curve = 1, curvature = 1, fixedStyle = 1,
         exoCov = FALSE, rotation = 1)
dev.off()


### Plot composites
pdf(file = "SMP_micro_leaf_prok.pdf", height = 8.27, width = 6)
semPaths(fit.micro, what = "est", whatLabels = "est", residuals = FALSE,
         intercepts = FALSE, sizeMan = 5, sizeMan2 = 3, edge.label.cex = 0.35,
         fade = FALSE, layout = "tree", style = "mx", nCharNodes = 0,
         posCol = "#009e73ff", negCol = "#d55e00ff", edge.label.color = "black",
         layoutSplit = TRUE, curve = 1, curvature = 1, fixedStyle = 1,
         exoCov = FALSE, rotation = 1)
dev.off()


pdf(file = "SMP_macro_leaf_prok.pdf", height = 8.27, width = 6)
semPaths(fit.macro, what = "est", whatLabels = "est", residuals = FALSE,
         intercepts = FALSE, sizeMan = 5, sizeMan2 = 3, edge.label.cex = 0.35,
         fade = FALSE, layout = "tree", style = "mx", nCharNodes = 0,
         posCol = "#009e73ff", negCol = "#d55e00ff", edge.label.color = "black",
         layoutSplit = TRUE, curve = 1, curvature = 1, fixedStyle = 1,
         exoCov = FALSE, rotation = 1)
dev.off()


################################################################################
################################################################################
################################################################################
################################################################################
