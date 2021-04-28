################################################################################
################################################################################
################################################################################
################################################################################
### SMP Biodiversity metrics
### Authors: Rachel Korn
### korn@cumulonimbus.at University of Fribourg 2020
################################################################################


library("phyloseq")
library("ggplot2")
theme_set(theme_bw(base_size = 20) +
            theme(rect = element_rect(fill = "transparent")))
library("vegan")
library("reshape2")
library("nlme")


rm(list = ls())
setwd("~/Sarracenia-Microbiome-Project/Thesis/")


################################################################################
### Read data
prok.a <- readRDS("rds/SMP_prok.a.RDS")
euk.a <- readRDS("rds/SMP_euk.a.RDS")
smp <- readRDS("rds/SMP_smp.RDS")


################################################################################
## Colorblind temperature scale for the 5 sites
gradCol <- c("#D55E00", "#E69F00", "#F0E442", "#009E73", "#56B4E9")


################################################################################
### Export dataframe with taxa and abundances
## Pitchers
tax.l <- as.data.frame(smp@tax_table@.Data)
otu.l <- as.data.frame(t(otu_table(smp)))

tax.otu.l <- merge(tax.l, otu.l, by = 0, all = TRUE)
rownames(tax.otu.l) <- tax.otu.l$Row.names
colnames(tax.otu.l)[1] <- "Taxa"


rm(tax.l, otu.l)


################################################################################
### Alpha-diversity
## Prokaryotes
alpha.p <- data.frame(otu_table(transform_sample_counts(prok.a,
                                                        function(abund) {
                                                          1 * (abund > 0)
                                                          })))

alpha.p$alpha <- rowSums(alpha.p)
alpha.p <- alpha.p[c("alpha")]
alpha.p$FullID <- rownames(alpha.p)
alpha.p$Domain <- "Prokaryotes"


## Eukaryotes
alpha.e <- data.frame(otu_table(transform_sample_counts(euk.a,
                                                        function(abund) {
                                                          1 * (abund > 0)
                                                          })))


## Add empty sample
alpha.e <- rbind(alpha.e, "LT4a06E" = rep(0, dim(alpha.e)[2]))


alpha.e$alpha <- rowSums(alpha.e)
alpha.e <- alpha.e[c("alpha")]
alpha.e$FullID <- rownames(alpha.e)
alpha.e$Domain <- "Eukaryotes"


## Merge
alpha <- rbind(alpha.p, alpha.e)
rm(alpha.p, alpha.e)


## Merge with metadata
meta.l <- data.frame(sample_data(prok.a))

alpha <- merge(alpha, meta.l, by = "FullID", all = TRUE)
# rownames(alpha) == rownames(meta.l)
alpha$Domain <- ordered(alpha$Domain, levels = c("Prokaryotes", "Eukaryotes"))


## Alpha-diversity
ggplot(alpha, aes(x = Succession, y = alpha, color = Site,
                  shape = Succession)) +
  geom_boxplot(color = "black", size = 0.2, outlier.shape = NA) +
  geom_jitter(height = 0, width = 0.3) +
  facet_grid(Domain ~ Site, scales = "free") +
  scale_colour_manual(values = gradCol) +
  xlab("") +
  ylab(expression(alpha-diversity)) +
  theme(legend.position = "top", legend.direction = "horizontal",
        axis.text.x = element_blank())
# ggsave("SMP_Alpha.pdf", width = 11.69, height = 6)


## Alpha by temperature
ggplot(alpha, aes(x = MeanTemperature, y = alpha, color = Site)) +
  geom_smooth(aes(group = 1), color = "gray",
              method = "lm", formula = y ~ poly(x, 2), se = FALSE) +
  geom_point(pch = 21) +
  facet_grid(Domain ~ Succession, scales = "free_y") +
  scale_colour_manual(values = gradCol) +
  theme(legend.position = "top", legend.direction = "horizontal") +
  xlab("Mean temperature [°C]") +
  ylab(expression(paste(alpha, "-diversity")))
# ggsave("SMP_AlphaByTemperature.pdf", width = 11.69, height = 6)


## Calculate fit
alpha.d <- alpha[alpha$Domain == "Eukaryotes", ]


a.fit <- lmList(alpha ~ poly(MeanTemperature, 2, raw = TRUE) | Succession,
                data = alpha.d, na.action = na.exclude)
summary(a.fit)
coef(a.fit)


################################################################################
### Gamma-diversity
## Prokaryotes
otus <- as.data.frame(t(otu_table(prok.a)))


## Get site names: fill with .* for regex matching
sites <- unique(paste(substr(names(otus), 1, 2), substr(names(otus), 7, 7),
                      sep = ".*"))


## Get rowSums by site and convert to presence/absence
# grep(sites[10], names(otus), value = TRUE) # Control match
gamma.site.p <- sapply(sites, function(xx) {rowSums(otus[, grep(xx,
                                                                   names(otus)),
                                                       drop = FALSE])})
gamma.site.p[gamma.site.p > 0] <- 1
gamma.site.p <- data.frame(t(gamma.site.p))
gamma.site.p$gamma <- rowSums(gamma.site.p)

gamma.site.p <- gamma.site.p["gamma"]
gamma.site.p$Domain <- "Prokaryotes"


### Eukaryotes
otus <- as.data.frame(t(otu_table(euk.a)))


## Get rowSums by site and convert to presence/absence
gamma.site.e <- sapply(sites, function(xx) {rowSums(otus[, grep(xx,
                                                                   names(otus)),
                                                       drop = FALSE])})
gamma.site.e[gamma.site.e > 0] <- 1
gamma.site.e <- data.frame(t(gamma.site.e))
gamma.site.e$gamma <- rowSums(gamma.site.e)


gamma.site.e <- gamma.site.e["gamma"]
gamma.site.e$Domain <- "Eukaryotes"


## Merge both
gamma.site <- rbind(gamma.site.p, gamma.site.e)
gamma.site$Site <- substr(rownames(gamma.site), 1, 2)
gamma.site$Sector <- "Site-wise"
gamma.site$Succession <- as.factor(substr(rownames(gamma.site), 5, 5))
levels(gamma.site$Succession) <- c("Early", "Late")
row.names(gamma.site) <- NULL


rm(gamma.site.p, gamma.site.e, sites, otus)


### Gamma-diversity site-, sector- and succession-wise
## Prokaryotes
otus <- as.data.frame(t(otu_table(prok.a)))


## Get site+sector names
site.sec <- unique(paste(substr(names(otus), 1, 3), substr(names(otus), 7, 7),
                         sep = ".*"))


## Get rowSums by site and convert to presence/absence
# grep(site.sec[7], names(otus), value = TRUE) # Control match
gamma.sec.p <- sapply(site.sec, function(xx) {rowSums(otus[, grep(xx,
                                                                  names(otus)),
                                                           drop = FALSE])})
gamma.sec.p[gamma.sec.p > 0] <- 1

gamma.sec.p <- data.frame(t(gamma.sec.p))
gamma.sec.p$gamma <- rowSums(gamma.sec.p)


gamma.sec.p <- gamma.sec.p["gamma"]
gamma.sec.p$Domain <- "Prokaryotes"


### Eukaryotes
otus <- as.data.frame(t(otu_table(euk.a)))


## Get rowSums by site and convert to presence/absence
gamma.sec.e <- sapply(site.sec, function(xx) {rowSums(otus[, grep(xx,
                                                                  names(otus)),
                                                       drop = FALSE])})
gamma.sec.e[gamma.sec.e > 0] <- 1

gamma.sec.e <- data.frame(t(gamma.sec.e))
gamma.sec.e$gamma <- rowSums(gamma.sec.e)


gamma.sec.e <- gamma.sec.e["gamma"]
gamma.sec.e$Domain <- "Eukaryotes"


## Merge both
gamma.sec <- rbind(gamma.sec.p, gamma.sec.e)
gamma.sec$Site <- substr(rownames(gamma.sec), 1, 2)
gamma.sec$Sector <- substr(rownames(gamma.sec), 3, 3)
gamma.sec$Succession <- as.factor(substr(rownames(gamma.sec), 6, 6))
levels(gamma.sec$Succession) <- c("Early", "Late")
row.names(gamma.sec) <- NULL


rm(gamma.sec.p, gamma.sec.e, site.sec, otus)


### Altogether
gamma <- rbind(gamma.site, gamma.sec)


gamma.sec$Domain <- factor(gamma.sec$Domain,
                           levels = c("Prokaryotes", "Eukaryotes"))
gamma.sec$Site <- factor(gamma.sec$Site,
                         levels = c("CB", "LT", "LE", "LV", "LM"))


ggplot(gamma.sec, aes(x = Succession, y = gamma)) +
  geom_boxplot(aes(fill = Site),
               color = "black", size = 0.2, outlier.shape = NA) +
  geom_point(aes(group = Site), fill = "white",
             pch = 21, position = position_dodge(0.75)) +
  stat_summary(aes(group = Site), fun = mean, position = position_dodge(0.75),
               colour = "black", geom = "point",
               shape = 8, size = 2, show.legend = FALSE) +
  facet_wrap( ~ Domain, scales = "free") +
  scale_fill_manual(values = gradCol) +
  xlab("Successional stage") +
  ylab(expression(gamma-diversity)) +
  theme(legend.position = "top", legend.direction = "horizontal")
# ggsave("SMP_Gamma.pdf", width = 11.69, height = 5)


################################################################################
### Alpha vs gamma
## Get average alphas
alpha.mean.site <- aggregate(alpha$alpha, by = list(alpha$Site, alpha$Domain,
                                                    alpha$Succession),
                             FUN = mean, simplify = TRUE)
names(alpha.mean.site) <- c("Site", "Domain", "Succession", "mean.alpha")
alpha.mean.site$Sector <- "Site-wise"


alpha.mean.sec <- aggregate(alpha$alpha,
                            by = list(alpha$Site, alpha$Sector, alpha$Domain,
                                      alpha$Succession),
                            FUN = mean, simplify = TRUE)
names(alpha.mean.sec) <- c("Site", "Sector", "Domain", "Succession",
                           "mean.alpha")


alpha.mean <- rbind(alpha.mean.site, alpha.mean.sec)


## Merge mean alpha with gamma
ag <- merge(alpha.mean, gamma, by = c("Site", "Sector", "Domain", "Succession"))


## Remove site-level means
am.g.s <- ag[ag$Sector != "Site-wise", ]


## Subset data
am.g.pE <- am.g.s[am.g.s$Domain == "Prokaryotes" &
                    am.g.s$Succession == "Early", ]
am.g.pL <- am.g.s[am.g.s$Domain == "Prokaryotes" &
                    am.g.s$Succession == "Late", ]
am.g.eE <- am.g.s[am.g.s$Domain == "Eukaryotes" &
                    am.g.s$Succession == "Early", ]
am.g.eL <- am.g.s[am.g.s$Domain == "Eukaryotes" & am.g.s$Succession == "Late", ]


am.g.list <- list(am.g.pE, am.g.pL, am.g.eE, am.g.eL)
names(am.g.list) <- c("am.g.pE", "am.g.pL", "am.g.eE", "am.g.eL")
rm(am.g.pE, am.g.pL, am.g.eE, am.g.eL)


### Fit Michaelis-Menten with one parameter
## Global functions and values
fit.alpha <- function(gamma, a){
  fit.mean.alpha <- a * gamma / (a + gamma)
  return(fit.mean.alpha)
}


## Empty list for pro- and eukaryotes and both successional stages
am.g.list.e <- vector("list", length(names(am.g.list)))
names(am.g.list.e) <- names(am.g.list)


get.a <- function(alpha.gamma){
  ## Michaelis-Menten
  fit.res <- function(a){
    ss <- sum((alpha.gamma$mean.alpha -
                 fit.alpha(alpha.gamma$gamma, a))^2)
    return(ss)}

  ## Minimize a
  out <- optim(par = 20, fn = fit.res, method = "Brent", lower = 0,
               upper = 1000)
  a <- out$par

  ## Write to data.frame and to list
  x <- 0:(max(alpha.gamma$gamma) + 5)
  fit <- data.frame(x = x,
                    y = a * x / (a + x),
                    Domain = rep(alpha.gamma$Domain[1], length(x)),
                    Succession = rep(alpha.gamma$Succession[1], length(x)))

  am.g.list.e <- fit
  return(am.g.list.e)
}


am.g.fit <- lapply(am.g.list, FUN = get.a)
am.g.fit <- do.call("rbind", am.g.fit)


## Plot
ggplot(am.g.s) +
  geom_line(data = am.g.fit, aes(x, y), color = "black", lty = 1) +
  geom_smooth(aes(x = gamma, y = mean.alpha), color = "gray", lty = 2,
              size = 0.5, method = "lm", formula = y ~ x, fullrange = TRUE,
              se = FALSE) +
  geom_point(aes(x = gamma, y = mean.alpha, color = Site, shape = Sector),
             size = 2.5) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  facet_grid(Succession ~ Domain, scales = "free") +
  scale_colour_manual(values = gradCol) +
  ylab(expression(bar(alpha)-diversity)) +
  xlab(expression(gamma-diversity)) +
  ylim(0, 100) +
  theme(legend.position = "top", legend.direction = "horizontal")
# ggsave("SMP_GammaVsAlpha.pdf", width = 11.69, height = 6)


################################################################################
### Beta diversity
beta <- merge(alpha, gamma.sec[, 1:5],
              by = c("Site", "Sector", "Succession", "Domain"))
beta$beta <- beta$gamma / beta$alpha


## Plot
ggplot(beta, aes(x = Succession, y = log(beta), color = Site,
                 shape = Succession)) +
  geom_boxplot(color = "black", size = 0.2, outlier.shape = NA) +
  geom_jitter(height = 0, width = 0.3) +
  facet_grid(Domain ~ Site) +
  scale_colour_manual(values = gradCol) +
  xlab("") +
  ylim(0, 3) +
  ylab(expression(log[e](beta-diversity))) +
  theme(legend.position = "top", legend.direction = "horizontal",
        axis.text.x = element_blank())
# ggsave("SMP_Beta.pdf", width = 11.69, height = 6)


################################################################################
### Pilou's Evenness
## Prokaryotes
even.p <- as.data.frame(otu_table(prok.a))


even.p$H <- diversity(even.p, index = "shannon", MARGIN = 1, base = exp(1))
even.p$H_max <- log(specnumber(even.p[, -346]))
even.p$J <- even.p$H / even.p$H_max
even.p <- even.p[c("H", "H_max", "J")]
even.p$Domain <- "Prokaryotes"


## Eukaryotes
even.e <- as.data.frame(otu_table(euk.a))


## Add empty sample
even.e <- rbind(even.e, "LT4a06E" = rep(0, dim(even.e)[2]))


even.e$H <- diversity(even.e, index = "shannon", MARGIN = 1, base = exp(1))
even.e$H_max <- log(specnumber(even.e[, -108]))
even.e$J <- even.e$H / even.e$H_max
even.e <- even.e[c("H", "H_max", "J")]
even.e$Domain <- "Eukaryotes"


## Merge
even.p$FullID <- rownames(even.p)
even.e$FullID <- rownames(even.e)

even <- rbind(even.p, even.e)


even$Site <- substr(rownames(even), 1, 2)
even$Sector <- substr(rownames(even), 3, 3)
even$Succession <- as.factor(substr(rownames(even), 7, 7))
levels(even$Succession) <- c("Early", "Late")


even$Domain <- ordered(even$Domain, levels = c("Prokaryotes", "Eukaryotes"))
even$Site <- ordered(even$Site, levels = c("CB", "LT", "LE", "LV", "LM"))


## Plot
ggplot(even, aes(x = Succession, y = J, color = Site,
                 shape = Succession)) +
  geom_boxplot(color = "black", size = 0.2, outlier.shape = NA) +
  geom_jitter(height = 0, width = 0.3) +
  facet_grid(Domain ~ Site) +
  scale_colour_manual(values = gradCol) +
  xlab("") +
  ylab(expression(paste("Pielou's evenness ", italic(J)))) +
  theme(legend.position = "top", legend.direction = "horizontal",
        axis.text.x = element_blank())
# ggsave("SMP_Evenness.pdf", width = 11.69, height = 6)


## Plot evenness by temperature
even <- merge(even, alpha[, c(1, 3, 21)], by = c("Domain", "FullID"), all = TRUE)


ggplot(even, aes(x = MeanTemperature, y = J, color = Site)) +
  geom_smooth(aes(group = 1), color = "gray",
              method = "lm", formula = y ~ x, se = FALSE) +
  geom_point(pch = 21) +
  facet_grid(Domain ~ Succession) +
  scale_colour_manual(values = gradCol) +
  xlab("Mean temperature [°C]") +
  ylab(expression(paste("Pielou's evenness ", italic(J)))) +
  ylim(0, 1) +
  theme(legend.position = "top", legend.direction = "horizontal")
# ggsave("SMP_EvennessByTemperature.pdf", width = 11.69, height = 6)


## Calculate fit
even.d <- even[even$Domain == "Prokaryotes", ]


e.fit <- lmList(J ~ MeanTemperature | Succession,
                data = even.d, na.action = na.exclude)
summary(e.fit)
coef(e.fit)


################################################################################
### Renyi profiles
## Prokaryotes
otus <- data.frame(otu_table(prok.a))

renyi.p <- renyi(otus)
renyi.p$Domain <- "Prokaryotes"
renyi.p$Site <- substr(rownames(renyi.p), 1, 2)
renyi.p$Succession <- as.factor(substr(rownames(renyi.p), 7, 7))
levels(renyi.p$Succession) <- c("Early", "Late")


renyi.p.m <- melt(renyi.p, id.vars = c("Domain", "Site", "Succession"),
                  variable.name = "alpha", value.name = "RenyiEntropy")


## Eukaryotes
otus <- data.frame(otu_table(euk.a))

renyi.e <- renyi(otus)
renyi.e$Domain <- "Eukaryotes"
renyi.e$Site <- substr(rownames(renyi.e), 1, 2)
renyi.e$Succession <- as.factor(substr(rownames(renyi.e), 7, 7))
levels(renyi.e$Succession) <- c("Early", "Late")



renyi.e.m <- melt(renyi.e, id.vars = c("Domain", "Site", "Succession"),
                  variable.name = "alpha", value.name = "RenyiEntropy")


## Merge both
renyi.m <- rbind(renyi.p.m, renyi.e.m)
renyi.m$alpha <- factor(renyi.m$alpha)


renyi.m$Domain <- factor(renyi.m$Domain,
                         levels = c("Prokaryotes", "Eukaryotes"))
renyi.m$Site <- factor(renyi.m$Site, levels = c("CB", "LT", "LE", "LV", "LM"))


## Mean, max and min by domain, site and succession
renyi.mean <- aggregate(renyi.m$RenyiEntropy,
                        by = list(renyi.m$Domain, renyi.m$alpha,
                                  renyi.m$Succession, renyi.m$Site), FUN = mean)
names(renyi.mean) <- c("Domain", "alpha", "Succession", "Site", "Mean")

renyi.min <- aggregate(renyi.m$RenyiEntropy,
                       by = list(renyi.m$Domain, renyi.m$alpha,
                                 renyi.m$Succession, renyi.m$Site), FUN = min)
names(renyi.min) <- c("Domain", "alpha", "Succession", "Site", "Min")

renyi.max <- aggregate(renyi.m$RenyiEntropy,
                       by = list(renyi.m$Domain, renyi.m$alpha,
                                 renyi.m$Succession, renyi.m$Site), FUN = max)
names(renyi.max) <- c("Domain", "alpha", "Succession", "Site", "Max")
# all(renyi.mean[, 1:3] == renyi.max[, 1:3])
# all(renyi.max[, 1:3] == renyi.min[, 1:3])


renyi.stats <- cbind(renyi.mean,
                     Min = renyi.min$Min,
                     Max =  renyi.max$Max)
rm(renyi.mean, renyi.min, renyi.max)
renyi.stats$alpha <- factor(renyi.stats$alpha)


renyi.stats$Domain <- factor(renyi.stats$Domain,
                             levels = c("Prokaryotes", "Eukaryotes"))


## Plot
ggplot() +
  geom_line(data = renyi.stats, aes(x = alpha, y = Mean, group = Site,
                                    color = Site)) +
  geom_boxplot(data = renyi.m, aes(x = alpha, y = RenyiEntropy, fill = Site),
               outlier.shape = 21, outlier.size = 0.7) +
  facet_grid(Succession ~ Domain) +
  scale_fill_manual(values = gradCol) +
  scale_colour_manual(values = gradCol) +
  xlab(expression(alpha)) +
  ylab("Rényi entropy") +
  theme(legend.position = "top", legend.direction = "horizontal",
        axis.text.x = element_text(angle = 45, hjust = 1))
# ggsave("SMP_Renyi.pdf", width = 11.69, height = 6)


################################################################################
################################################################################
################################################################################
################################################################################
