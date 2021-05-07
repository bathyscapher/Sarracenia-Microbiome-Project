################################################################################
################################################################################
################################################################################
################################################################################
### Mosses SEM
### Author: korn@cumulonimbus.at University of Fribourg 2020
################################################################################


library("lavaan")
library("semPlot")
library("blavaan")
library("gplots")


rm(list = ls())
setwd("~/Sarracenia-Microbiome-Project/Thesis/")


mossMeta <- read.table("csv/Moss_Metadata_scaled.csv", sep = "\t",
                       header = TRUE)
dim(mossMeta)
colnames(mossMeta)
mossMeta$FullID


################################################################################
### Control for highly correlated variables
check.vars <- c("Site", "Altitude", "Light", "Temperature", "Humidity",
               "Soil.reaction", "Nitrogen", "CanopyCover", "Moss_g",
               "Distance", "MossHabitat", "mb.moss.prok", "mb.moss.euk")
mossMeta.sub <- mossMeta[check.vars]
print(round(cor(mossMeta.sub[, -1], use = "complete"), digits = 1))


mossMeta.sub$Site <- ordered(mossMeta.sub$Site,
                             levels = c("CB", "LE", "LV", "LT", "LM"))


### Correlation and covariance
my_palette <- colorRampPalette(c("#D55E00", "white",  "#56B4E9"))(n = 19)


## Sort (out) correlation data.frame
covs <- data.frame(cov(mossMeta[, -c(1:9, 11:15, 18, 22)], use = "complete"))
covs <- covs[c("Altitude", "Distance", "CanopyCover", "Light", "Moss_g",
               "MossHabitat", "Temperature", "Humidity", "Soil.reaction",
               "Nitrogen", "mb.moss.prok", "mb.moss.euk")]
covs <- covs[colnames(covs), ]


# pdf("Moss_Metadata_Covariance.pdf", height = 8.27, width = 8.27)
# heatmap.2(as.matrix(covs), dendrogram = "none", key = TRUE, key.title = "",
#           cellnote = round(covs, digits = 2), notecol = "black", notecex = 0.7,
#           trace = "none", distfun = function(x) {x}, cexCol = 1.4, cexRow = 1.4,
#           symm = TRUE, margins = c(16, 16), col = my_palette,
#           tracecol = "black", Rowv = FALSE, Colv = FALSE)
# dev.off()


################################################################################
### SEM
## EV composite
ev <-
'mb.moss.prok ~ Light + Temperature + Humidity + Soil.reaction
mb.moss.euk ~ Light + Temperature + Humidity + Soil.reaction'


fit.ev <- sem(ev, data = mossMeta)
summary(fit.ev, rsq = TRUE)


## Extract coefficients
Light.p <- coef(fit.ev)[[1]]
Temperature.p <- coef(fit.ev)[[2]]
Humidity.p <- coef(fit.ev)[[3]]
Soil.reaction.p <- coef(fit.ev)[[4]]


coef(fit.ev)[1:4]; Light.p; Temperature.p; Humidity.p; Soil.reaction.p


Light.e <- coef(fit.ev)[[5]]
Temperature.e <- coef(fit.ev)[[6]]
Humidity.e <- coef(fit.ev)[[7]]
Soil.reaction.e <- coef(fit.ev)[[8]]


coef(fit.ev)[5:8]; Light.e; Temperature.e; Humidity.e; Soil.reaction.e


## Compute composite "ev"
mossMeta$ev.p <- Light.p * mossMeta$Light +
  Temperature.p * mossMeta$Temperature +
  Humidity.p * mossMeta$Humidity +
  Soil.reaction.p * mossMeta$Soil.reaction

mossMeta$ev.e <- Light.e * mossMeta$Light +
  Temperature.e * mossMeta$Temperature +
  Humidity.e  * mossMeta$Humidity +
  Soil.reaction.e * mossMeta$Soil.reaction


### Run SEM with the manually computed composite
ev.comp <-
'mb.moss.prok ~ ev.p
mb.moss.euk ~ ev.e'

fit.ev.comp <- sem(ev.comp, data = mossMeta)


## Control if the standard coefficient and z-value (standardised coefficient)
## are similar
summary(fit.ev, rsq = TRUE, standardized = TRUE)
summary(fit.ev.comp, rsq = TRUE, standardized = TRUE)


rm(fit.ev.comp, ev, ev.comp, Light.p, Temperature.p, Humidity.p,
   Soil.reaction.p, Light.e, Temperature.e, Humidity.e, Soil.reaction.e)


################################################################################
### SEM
set.seed(12614)


moss <-
'MossHabitat ~ Altitude + CanopyCover + Distance
ev.p ~ Altitude + CanopyCover + Distance
ev.e ~ Altitude + CanopyCover + Distance

mb.moss.prok ~ Altitude + CanopyCover + Distance + ev.p + MossHabitat + Moss_g
mb.moss.euk ~ Altitude + CanopyCover + Distance + ev.e + MossHabitat + Moss_g

mb.moss.prok ~ mb.moss.euk
mb.moss.euk ~ mb.moss.prok

ev.p ~~ MossHabitat
ev.e ~~ MossHabitat
ev.e ~~ ev.p'


fit.moss <- sem(moss, data = mossMeta, estimator = "MLM")
summary(fit.moss, rsq = TRUE, fit.measures = TRUE)


modificationindices(fit.moss, minimum.value = 3)


### Plot
pdf(file = "Moss_SEM.pdf", height = 8.27, width = 8.27)
semPaths(fit.moss, what = "est", whatLabels = "est", residuals = FALSE,
         intercepts = FALSE, sizeMan2 = 3, edge.label.cex = 0.35,
         fade = FALSE, layout = "tree", style = "mx", nCharNodes = 0,
         posCol = "#009e73ff", negCol = "#d55e00ff", edge.label.color = "black",
         layoutSplit = TRUE, curve = 1, curvature = 1, fixedStyle = 1,
         exoCov = FALSE, rotation = 1)
dev.off()


pdf(file = "Moss_SEM_composites.pdf", height = 8.27, width = 8.27)
semPaths(fit.ev, what = "est", whatLabels = "est", residuals = FALSE,
         intercepts = FALSE, sizeMan = 5, sizeMan2 = 3, edge.label.cex = 0.5,
         fade = FALSE, layout = "tree", style = "mx", nCharNodes = 0,
         posCol = "#009e73ff", negCol = "#d55e00ff", edge.label.color = "black",
         layoutSplit = TRUE, curve = 2, curvature = 2, fixedStyle = 1,
         exoCov = FALSE, rotation = 4)
dev.off()


################################################################################
################################################################################
################################################################################
################################################################################
