################################################################################
################################################################################
################################################################################
################################################################################
### SMP Biodiversity metrics linear models prokaryotes
################################################################################


library("phyloseq")
library("ggplot2")
theme_set(theme_bw(base_size = 20) +
            theme(rect = element_rect(fill = "transparent")))
library("vegan")
library("reshape2")
library("brms")


rm(list = ls())
setwd("~/Sarracenia-Microbiome-Project/Thesis/")


seed <- 34306


################################################################################
## Colorblind temperature scale for the 5 sites
gradCol <- c("#D55E00", "#E69F00", "#F0E442", "#009E73", "#56B4E9")


################################################################################
### Read data
prok.a <- readRDS("rds/SMP_prok.a.RDS")
euk.a <- readRDS("rds/SMP_euk.a.RDS")


################################################################################
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


## Sort factors and scale numeric variables
alpha$Domain <- ordered(alpha$Domain, levels = c("Prokaryotes", "Eukaryotes"))
alpha$SampleColor <- factor(alpha$SampleColor,
                            levels = c("clear", "cloudy", "green", "brown"))


## Prokaryotes and Eukaryotes
alpha.p <- alpha[alpha$Domain == "Prokaryotes", ]
alpha.e <- alpha[alpha$Domain == "Eukaryotes", ]


## Add alpha of other domain
alpha.p$alpha.OtherDomain <- alpha.e$alpha
# alpha.e$alpha.OtherDomain <- alpha.p$alpha
# summary(alpha.p[, 1:2] == alpha.e[, c(1, 26)])
# summary(alpha.e[, 1:2] == alpha.p[, c(1, 26)])

rm(alpha)


################################################################################
### Gamma-diversity site-, sector- and succession-wise
## Prokaryotes
otus <- as.data.frame(t(otu_table(prok.a)))


## Get site+sector names
site.sec <- unique(paste(substr(names(otus), 1, 3), substr(names(otus), 7, 7),
                         sep = ".*"))


## Get rowSums by site and convert to presence/absence
gamma.sec.p <- sapply(site.sec, function(xx) {rowSums(otus[, grep(xx,
                                                                  names(otus)),
                                                           drop = FALSE])})
gamma.sec.p[gamma.sec.p > 0] <- 1

gamma.sec.p <- data.frame(t(gamma.sec.p))
gamma.sec.p$gamma <- rowSums(gamma.sec.p)


gamma.sec.p <- gamma.sec.p["gamma"]
gamma.sec.p$Domain <- "Prokaryotes"
gamma.sec.p$Site <- substr(rownames(gamma.sec.p), 1, 2)
gamma.sec.p$Sector <- substr(rownames(gamma.sec.p), 3, 3)
gamma.sec.p$Succession <- as.factor(substr(rownames(gamma.sec.p), 6, 6))
levels(gamma.sec.p$Succession) <- c("Early", "Late")
row.names(gamma.sec.p) <- NULL


## Eukaryotes
otus <- as.data.frame(t(otu_table(euk.a)))


## Get site+sector names
site.sec <- unique(paste(substr(names(otus), 1, 3), substr(names(otus), 7, 7),
                         sep = ".*"))


## Get rowSums by site and convert to presence/absence
gamma.sec.e <- sapply(site.sec, function(xx) {rowSums(otus[, grep(xx,
                                                                  names(otus)),
                                                           drop = FALSE])})
gamma.sec.e[gamma.sec.e > 0] <- 1

gamma.sec.e <- data.frame(t(gamma.sec.e))
gamma.sec.e$gamma <- rowSums(gamma.sec.e)


gamma.sec.e <- gamma.sec.e["gamma"]
gamma.sec.e$Domain <- "Eukaryotes"
gamma.sec.e$Site <- substr(rownames(gamma.sec.e), 1, 2)
gamma.sec.e$Sector <- substr(rownames(gamma.sec.e), 3, 3)
gamma.sec.e$Succession <- as.factor(substr(rownames(gamma.sec.e), 6, 6))
levels(gamma.sec.e$Succession) <- c("Early", "Late")
row.names(gamma.sec.e) <- NULL


## Add gamma of other domain
gamma.sec.p$gamma.OtherDomain <- gamma.sec.e$gamma
# summary(gamma.sec.p[, c(3:6)] == gamma.sec.e[, c(3:5, 1)])


## Add metadata
meta.l.sec <- aggregate(meta.l[, c(9, 11:17, 19:23)],
                        by = list(meta.l$Site, meta.l$Sector,
                                  meta.l$Succession),
                        FUN = mean, na.rm = TRUE)
names(meta.l.sec)[1:3] <- c("Site", "Sector", "Succession")


gamma.sec.p <- merge(gamma.sec.p, meta.l.sec,
                     by = c("Site", "Sector", "Succession"), all = TRUE)


rm(site.sec, otus, meta.l.sec)


################################################################################
### Beta diversity
## Prokaryotes
## Average alpha on site-sector level
alpha.sec.p <- aggregate(alpha.p[, -c(1, 3:10, 12, 20)],
                         by = list(meta.l$Site, meta.l$Sector,
                                   meta.l$Succession),
                         FUN = mean, na.rm = TRUE)
names(alpha.sec.p)[1:4] <- c("Site", "Sector", "Succession", "MeanAlpha")
alpha.sec.p$Domain <- "Prokaryotes"


## Merge and calculate beta
beta.p <- merge(alpha.sec.p, gamma.sec.p[, 1:5],
                by = c("Site", "Sector", "Succession", "Domain"))
beta.p$beta <- beta.p$gamma / beta.p$MeanAlpha


## Eukaryotes
## Average alpha on site-sector level
alpha.sec.e <- aggregate(alpha.e[, -c(1, 3:10, 12, 20)],
                         by = list(meta.l$Site, meta.l$Sector,
                                   meta.l$Succession),
                         FUN = mean, na.rm = TRUE)
names(alpha.sec.e)[1:4] <- c("Site", "Sector", "Succession", "MeanAlpha")
alpha.sec.e$Domain <- "Eukaryotes"


## Merge and calculate beta
beta.e <- merge(alpha.sec.e, gamma.sec.e,
                by = c("Site", "Sector", "Succession", "Domain"))
beta.e$beta <- beta.e$gamma / beta.e$MeanAlpha


## Add beta of other domain
beta.p$beta.OtherDomain <- beta.e$beta
# summary(beta.p[, c(1:3, 22)] == beta.e[, c(1:3, 20)])


################################################################################
### Pilou's Evenness
## Prokaryotes
even.p <- as.data.frame(otu_table(prok.a))


even.p$H <- diversity(even.p, index = "shannon", MARGIN = 1, base = exp(1))
even.p$H_max <- log(specnumber(even.p[, -346]))
even.p$J <- even.p$H / even.p$H_max
even.p <- even.p[c("H", "H_max", "J")]


even.p$FullID <- rownames(even.p)
even.p$Site <- substr(rownames(even.p), 1, 2)
even.p$Sector <- substr(rownames(even.p), 3, 3)
even.p$Succession <- as.factor(substr(rownames(even.p), 7, 7))
levels(even.p$Succession) <- c("Early", "Late")


## Eukaryotes
even.e <- as.data.frame(otu_table(euk.a))


## Add empty sample
even.e <- rbind(even.e, "LT4a06E" = rep(0, dim(even.e)[2]))


even.e$H <- diversity(even.e, index = "shannon", MARGIN = 1, base = exp(1))
even.e$H_max <- log(specnumber(even.e[, -108]))
even.e$J <- even.e$H / even.e$H_max
even.e <- even.e[c("H", "H_max", "J")]


## Add evenness of other domain
even.e <- even.e[match(rownames(even.p), rownames(even.e)), ]
# rownames(even.e) == rownames(even.p)

even.p$J.OtherDomain <- even.e$J
# all(rbind(rownames(even.p), even.p$J.OtherDomain) ==
#       rbind(rownames(even.e), even.e$J))


## Merge with metadata
even.p <- merge(even.p, meta.l,
                by = c("FullID", "Site", "Sector", "Succession"), all = TRUE)
even.p$SampleColor <- factor(even.p$SampleColor,
                             levels = c("clear", "cloudy", "green", "brown"))


rm(meta.l, alpha.e, alpha.sec.e, beta.e, gamma.sec.e, even.e, prok.a, euk.a)


################################################################################
#### Linear models
### Alpha-diversity
## Scale
alpha.p[, c(2, 11, 13:19, 21:26)] <- apply(alpha.p[, c(2, 11, 13:19, 21:26)],
                                           2, scale)


### Full model with hierarchical sampling design
priors.nest <- get_prior(alpha ~ pH + prey.items + SampleColor + Age_d +
                           SampleVolume_mL + CanopyCover + alpha.OtherDomain +
                           (1 | Site / Sector),
                         family = gaussian(), data = alpha.p)
alpha.nest <- brm(alpha ~ pH + prey.items + SampleColor + Age_d +
                  SampleVolume_mL + CanopyCover + alpha.OtherDomain +
                  (1 | Site / Sector),
                  family = gaussian(), prior = priors.nest, cores = 8,
                  iter = 8000, save_all_pars = TRUE, seed = seed,
                  control = list(adapt_delta = 0.99999999999,
                                 max_treedepth = 15),
                  data = alpha.p)
summary(alpha.nest)
alpha.nest$fit

alpha.nest <- add_criterion(alpha.nest, criterion = "waic")
alpha.nest$criteria$waic
alpha.nest.loo <- loo(alpha.nest, reloo = TRUE)


## Save it
save(alpha.nest, file = "StanLM/SMP_LM_alpha.nest_prok.gz", compress = "gzip")


### With distance instead of site/sector as random effect
priors.dist <- get_prior(alpha ~ pH + prey.items + SampleColor + Age_d +
                           SampleVolume_mL + CanopyCover + alpha.OtherDomain +
                           Distance,
                         family = gaussian(), data = alpha.p)
alpha.dist <- brm(alpha ~ pH + prey.items + SampleColor + Age_d +
                  SampleVolume_mL + CanopyCover + alpha.OtherDomain +
                  Distance,
                  family = gaussian(), prior = priors.dist, cores = 8,
                  iter = 8000, save_all_pars = TRUE, seed = seed,
                  control = list(adapt_delta = 0.9, max_treedepth = 15),
                  data = alpha.p)
summary(alpha.dist)
alpha.dist$fit


alpha.dist <- add_criterion(alpha.dist, criterion = "waic")
alpha.dist$criteria$waic
alpha.dist.loo <- loo(alpha.dist, reloo = TRUE)


## Save it
save(alpha.dist, file = "StanLM/SMP_LM_alpha.dist_prok.gz", compress = "gzip")


### With interactions
priors.int <- get_prior(alpha ~ pH * prey.items * SampleColor * Age_d +
                           SampleVolume_mL + CanopyCover + alpha.OtherDomain +
                           (1 | Site / Sector),
                         family = gaussian(), data = alpha.p)
alpha.int <- brm(alpha ~ pH * prey.items * SampleColor * Age_d +
                 SampleVolume_mL + CanopyCover + alpha.OtherDomain +
                 (1 | Site / Sector),
                 save_all_pars = TRUE, seed = seed, family = gaussian(),
                 prior = priors.int, cores = 8, iter = 8000,
                 control = list(adapt_delta = 0.999, max_treedepth = 15),
                 data = alpha.p)
summary(alpha.int)
alpha.int$fit


alpha.int <- add_criterion(alpha.int, criterion = "waic")
alpha.int$criteria$waic
# alpha.int.loo <- loo(alpha.int, reloo = TRUE) # takes too long


## Save it
save(alpha.int, file = "StanLM/SMP_LM_alpha.int_prok.gz", compress = "gzip")


### Model comparison
# load("StanLM/SMP_LM_alpha.nest_prok.gz")
# load("StanLM/SMP_LM_alpha.dist_prok.gz")
# load("StanLM/SMP_LM_alpha.int_prok.gz")


loo_compare(alpha.nest, alpha.dist, alpha.int, criterion = "waic")
loo_compare(alpha.nest.loo, alpha.dist.loo)


## Plot random effects
alpha.nest.ran <- data.frame(ranef(alpha.nest)[2])
names(alpha.nest.ran) <- c("Estimate", "SDEstimate", "CrI2.5", "CrI97.5")
alpha.nest.ran$Site <- substr(rownames(alpha.nest.ran), 1, 2)
alpha.nest.ran$Sector <- substr(rownames(alpha.nest.ran), 4, 4)
alpha.nest.ran$SiteSector <- paste(alpha.nest.ran$Site, alpha.nest.ran$Sector)


alpha.nest.ran <- alpha.nest.ran[order(alpha.nest.ran$Estimate), ]
alpha.nest.ran$SiteSector <- factor(alpha.nest.ran$SiteSector,
                                    levels = unique(as.character(
                                      alpha.nest.ran$SiteSector)))
alpha.nest.ran$Site <- factor(alpha.nest.ran$Site,
                              levels = c("CB", "LT", "LE", "LV", "LM"))


ggplot(alpha.nest.ran, aes(x = SiteSector, y = Estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(ymin = CrI2.5,
                    ymax = CrI97.5), width = 0, color = "black") +
  geom_errorbar(aes(ymin = Estimate - SDEstimate,
                    ymax = Estimate + SDEstimate),
                width = 0, color = "black", size = 1.5) +
  geom_point(aes(fill = Site), pch = 21, size = 4) +
  scale_fill_manual(values = gradCol) +
  theme(legend.position = "none", legend.direction = "horizontal") +
  xlab("") +
  ylab(expression(paste("Estimates ", alpha-diversity))) +
  coord_flip()
# ggsave("SMP_Alpha_ranef_prok.pdf", width = 11.69 / 2, height = 6)


## Plot fixed effects
alpha.nest.fix <- data.frame(fixef(alpha.nest))
names(alpha.nest.fix) <- c("Estimate", "SDEstimate", "CrI2.5", "CrI97.5")
alpha.nest.fix$Variable <- rownames(alpha.nest.fix)


## Rename intercept by SampleColorclear
alpha.nest.fix$Variable[alpha.nest.fix$Variable == "Intercept"] <- "SampleColorclear"


alpha.nest.fix <- alpha.nest.fix[order(alpha.nest.fix$Estimate), ]
alpha.nest.fix$Variable <- factor(alpha.nest.fix$Variable,
                                    levels = unique(as.character(
                                      alpha.nest.fix$Variable)))


ggplot(alpha.nest.fix, aes(x = Variable, y = Estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(ymin = CrI2.5,
                    ymax = CrI97.5), width = 0, color = "black") +
  geom_errorbar(aes(ymin = Estimate - SDEstimate,
                    ymax = Estimate + SDEstimate),
                width = 0, color = "black", size = 1.5) +
  geom_point(pch = 21, fill = "white", size = 4) +
  xlab("") +
  ylab(expression(paste("Estimate ", alpha-diversity))) +
  coord_flip()
# ggsave("SMP_Alpha_fixef_prok.pdf", width = 11.69 / 2, height = 8.27 / 2)


################################################################################
### Beta diversity
## Scale numeric variables
beta.p[, -c(1:4)] <- apply(beta.p[, -c(1:4)], 2, scale)


### Full model with hierarchical sampling design
priors.nest <- get_prior(beta ~ pH + prey.items + Age_d + SampleVolume_mL +
                           CanopyCover + beta.OtherDomain + (1 | Site),
                         family = gaussian(), data = beta.p)
beta.nest <- brm(beta ~ pH + prey.items + Age_d + SampleVolume_mL +
                   CanopyCover + beta.OtherDomain + (1 | Site),
                 family = gaussian(), prior = priors.nest, cores = 8,
                 iter = 8000, save_all_pars = TRUE, seed = seed,
                 control = list(adapt_delta = 0.999999, max_treedepth = 15),
                 data = beta.p)
summary(beta.nest)
beta.nest$fit


beta.nest <- add_criterion(beta.nest, criterion = "waic")
beta.nest$criteria$waic
beta.nest.loo <- loo(beta.nest, reloo = TRUE)


## Save it
save(beta.nest, file = "StanLM/SMP_LM_beta.nest_prok.gz", compress = "gzip")


### With distance instead of site/sector as random effect
priors.dist <- get_prior(beta ~ pH + prey.items + Age_d + SampleVolume_mL +
                           CanopyCover + beta.OtherDomain + Distance,
                         family = gaussian(), data = beta.p)
beta.dist <- brm(beta ~ pH + prey.items + Age_d + SampleVolume_mL +
                   CanopyCover+ beta.OtherDomain + Distance,
                 family = gaussian(), prior = priors.dist, cores = 8,
                 iter = 8000, save_all_pars = TRUE, seed = seed,
                 control = list(adapt_delta = 0.9, max_treedepth = 15),
                 data = beta.p)
summary(beta.dist)
beta.dist$fit


beta.dist <- add_criterion(beta.dist, criterion = "waic")
beta.dist$criteria$waic
beta.dist.loo <- loo(beta.dist, reloo = TRUE)


## Save it
save(beta.dist, file = "StanLM/SMP_LM_beta.dist_prok.gz", compress = "gzip")


### With interactions
priors.int <- get_prior(beta ~ pH * prey.items * Age_d + SampleVolume_mL +
                          CanopyCover + beta.OtherDomain + (1 | Site),
                         family = gaussian(), data = beta.p)
beta.int <- brm(beta ~ pH * prey.items * Age_d + SampleVolume_mL + CanopyCover +
                  beta.OtherDomain + (1 | Site),
                save_all_pars = TRUE, seed = seed, family = gaussian(),
                prior = priors.int, cores = 8, iter = 8000,
                control = list(adapt_delta = 0.99, max_treedepth = 15),
                data = beta.p)
summary(beta.int)
beta.int$fit


beta.int <- add_criterion(beta.int, criterion = "waic")
beta.int$criteria$waic
beta.int.loo <- loo(beta.int, reloo = TRUE)


## Save it
save(beta.int, file = "StanLM/SMP_LM_beta.int_prok.gz", compress = "gzip")


### Model comparison
# load("StanLM/SMP_LM_beta.nest_prok.gz")
# load("StanLM/SMP_LM_beta.dist_prok,gz")
# load("StanLM/SMP_LM_beta.int_prok.gz")


loo_compare(beta.nest, beta.dist, beta.int, criterion = "waic")
loo_compare(beta.nest.loo, beta.dist.loo, beta.int.loo)


## Plot random effects
beta.nest.ran <- data.frame(ranef(beta.nest))
names(beta.nest.ran) <- c("Estimate", "SDEstimate", "CrI2.5", "CrI97.5")
beta.nest.ran$Site <- substr(rownames(beta.nest.ran), 1, 2)
beta.nest.ran$Sector <- substr(rownames(beta.nest.ran), 4, 4)
beta.nest.ran$SiteSector <- paste(beta.nest.ran$Site, beta.nest.ran$Sector)


beta.nest.ran <- beta.nest.ran[order(beta.nest.ran$Estimate), ]
beta.nest.ran$SiteSector <- factor(beta.nest.ran$SiteSector,
                                    levels = unique(as.character(
                                      beta.nest.ran$SiteSector)))


ggplot(beta.nest.ran, aes(x = SiteSector, y = Estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(ymin = CrI2.5,
                    ymax = CrI97.5), width = 0, color = "black") +
  geom_errorbar(aes(ymin = Estimate - SDEstimate,
                    ymax = Estimate + SDEstimate),
                width = 0, color = "black", size = 1.5) +
  geom_point(aes(fill = Site), pch = 21, size = 4) +
  scale_fill_manual(values = c("#D55E00", "#F0E442", "#56B4E9", "#E69F00",
                    "#009E73")) +
  theme(legend.position = "none", legend.direction = "horizontal") +
  xlab("") +
  ylab(expression(paste("Estimate ", beta-diversity))) +
  coord_flip()
# ggsave("SMP_Beta_ranef_prok.pdf", width = 11.69 / 2, height = 8.27 / 3)


## Plot fixed effects
beta.nest.fix <- data.frame(fixef(beta.nest))
names(beta.nest.fix) <- c("Estimate", "SDEstimate", "CrI2.5", "CrI97.5")
beta.nest.fix$Variable <- rownames(beta.nest.fix)

beta.nest.fix <- beta.nest.fix[!rownames(beta.nest.fix) %in% "Intercept", ]


beta.nest.fix <- beta.nest.fix[order(beta.nest.fix$Estimate), ]
beta.nest.fix$Variable <- factor(beta.nest.fix$Variable,
                                 levels = unique(as.character(
                                   beta.nest.fix$Variable)))


ggplot(beta.nest.fix, aes(x = Variable, y = Estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(ymin = CrI2.5,
                    ymax = CrI97.5), width = 0, color = "black") +
  geom_errorbar(aes(ymin = Estimate - SDEstimate,
                    ymax = Estimate + SDEstimate),
                width = 0, color = "black", size = 1.5) +
  geom_point(pch = 21, fill = "white", size = 4) +
  xlab("") +
  ylab(expression(paste("Estimate ", beta-diversity))) +
  coord_flip()
# ggsave("SMP_Beta_fixef_prok.pdf", width = 11.69 / 2, height = 8.27 / 3)


################################################################################
### Gamma diversity
## Scale numeric variables
gamma.sec.p[, c(4, 6:19)] <- apply(gamma.sec.p[, c(4, 6:19)], 2, scale)


### Full model with hierarchical sampling design
priors.nest <- get_prior(gamma ~ pH + prey.items + Age_d + SampleVolume_mL +
                           CanopyCover + gamma.OtherDomain + (1 | Site),
                         family = gaussian(), data = gamma.sec.p)
gamma.nest <- brm(gamma ~ pH + prey.items + Age_d + SampleVolume_mL +
                    CanopyCover + gamma.OtherDomain + (1 | Site),
                  family = gaussian(), prior = priors.nest, cores = 8,
                  iter = 8000, save_all_pars = TRUE, seed = seed,
                  control = list(adapt_delta = 0.999, max_treedepth = 15),
                  data = gamma.sec.p)
summary(gamma.nest)
gamma.nest$fit


gamma.nest <- add_criterion(gamma.nest, criterion = "waic")
gamma.nest$criteria$waic
gamma.nest.loo <- loo(gamma.nest, reloo = TRUE)


## Save it
save(gamma.nest, file = "StanLM/SMP_LM_gamma.nest_prok.gz", compress = "gzip")


### dist model with distance
priors.dist <- get_prior(gamma ~ pH + prey.items + Age_d + SampleVolume_mL +
                           CanopyCover + gamma.OtherDomain + Distance,
                         family = gaussian(), data = gamma.sec.p)
gamma.dist <- brm(gamma ~ pH + prey.items + Age_d + SampleVolume_mL +
                    CanopyCover + gamma.OtherDomain + Distance,
                  family = gaussian(), prior = priors.dist, cores = 8,
                  iter = 8000, save_all_pars = TRUE, seed = seed,
                  control = list(adapt_delta = 0.999, max_treedepth = 15),
                  data = gamma.sec.p)
summary(gamma.dist)
gamma.dist$fit


gamma.dist <- add_criterion(gamma.dist, criterion = "waic")
gamma.dist$criteria$waic
gamma.dist.loo <- loo(gamma.dist, reloo = TRUE)


## Save it
save(gamma.dist, file = "StanLM/SMP_LM_gamma.dist_prok.gz", compress = "gzip")


### With interactions
priors.int <- get_prior(gamma ~ pH * prey.items * Age_d + SampleVolume_mL +
                          CanopyCover + gamma.OtherDomain + (1 | Site),
                        family = gaussian(), data = gamma.sec.p)
gamma.int <- brm(gamma ~ pH * prey.items * Age_d + SampleVolume_mL +
                   CanopyCover + gamma.OtherDomain + (1 | Site),
                 family = gaussian(), prior = priors.int, cores = 8,
                 iter = 8000, save_all_pars = TRUE, seed = seed,
                 control = list(adapt_delta = 0.999, max_treedepth = 20),
                 data = gamma.sec.p)
summary(gamma.int)
gamma.int$fit


gamma.int <- add_criterion(gamma.int, criterion = "waic")
gamma.int$criteria$waic
gamma.int.loo <- loo(gamma.int, reloo = TRUE)


## Save it
save(gamma.int, file = "StanLM/SMP_LM_gamma.int_prok.gz", compress = "gzip")


### Model comparison
# load("StanLM/SMP_LM_gamma.nest_prok.gz")
# load("StanLM/SMP_LM_gamma.dist_prok.gz")
# load("StanLM/SMP_LM_gamma.int_prok.gz")


loo_compare(gamma.nest, gamma.dist, gamma.int, criterion = "waic")
loo_compare(gamma.nest.loo, gamma.dist.loo, gamma.int.loo)


## Plot random effects
gamma.nest.ran <- data.frame(ranef(gamma.nest)[1])
names(gamma.nest.ran) <- c("Estimate", "SDEstimate", "CrI2.5", "CrI97.5")
gamma.nest.ran$Site <- substr(rownames(gamma.nest.ran), 1, 2)


gamma.nest.ran <- gamma.nest.ran[order(gamma.nest.ran$Estimate), ]
gamma.nest.ran$Site <- factor(gamma.nest.ran$Site,
                                   levels = unique(as.character(
                                     gamma.nest.ran$Site)))


ggplot(gamma.nest.ran, aes(x = Site, y = Estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(ymin = CrI2.5,
                    ymax = CrI97.5), width = 0, color = "black") +
  geom_errorbar(aes(ymin = Estimate - SDEstimate,
                    ymax = Estimate + SDEstimate),
                width = 0, color = "black", size = 1.5) +
  geom_point(aes(fill = Site), pch = 21, size = 4) +
  scale_fill_manual(values = c("#E69F00", "#009E73", "#D55E00", "#F0E442",
                               "#56B4E9")) +
  theme(legend.position = "none", legend.direction = "horizontal") +
  xlab("") +
  ylab(expression(paste("Estimate ", gamma-diversity))) +
  coord_flip()
# ggsave("SMP_Gamma_ranef_prok.pdf", width = 11.69 / 2, height = 8.27 / 3)


## Plot fixed effects
gamma.nest.fix <- data.frame(fixef(gamma.nest))
names(gamma.nest.fix) <- c("Estimate", "SDEstimate", "CrI2.5", "CrI97.5")
gamma.nest.fix$Variable <- rownames(gamma.nest.fix)


gamma.nest.fix <- gamma.nest.fix[!rownames(gamma.nest.fix) %in% "Intercept", ]


gamma.nest.fix <- gamma.nest.fix[order(gamma.nest.fix$Estimate), ]
gamma.nest.fix$Variable <- factor(gamma.nest.fix$Variable,
                                 levels = unique(as.character(
                                   gamma.nest.fix$Variable)))


ggplot(gamma.nest.fix, aes(x = Variable, y = Estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(ymin = CrI2.5,
                    ymax = CrI97.5), width = 0, color = "black") +
  geom_errorbar(aes(ymin = Estimate - SDEstimate,
                    ymax = Estimate + SDEstimate),
                width = 0, color = "black", size = 1.5) +
  geom_point(pch = 21, fill = "white", size = 4) +
  xlab("") +
  ylab(expression(paste("Estimate ", gamma-diversity))) +
  coord_flip()
# ggsave("SMP_Gamma_fixef_prok.pdf", width = 11.69 / 2, height = 8.27 / 3)


################################################################################
### Evenness
## Scale: note: Pielou's evenness J not scaled as it returns NaN
even.p[, -c(1:4, 9:12, 14, 22)] <- apply(even.p[, -c(1:4, 9:12, 14, 22)],
                                         2, scale)


### Full model with hierarchical sampling design
priors.nest <- get_prior(J ~ pH + prey.items + SampleColor + Age_d +
                           SampleVolume_mL + CanopyCover + J.OtherDomain +
                           (1 | Site / Sector),
                         family = gaussian(), data = even.p)
even.nest <- brm(J ~ pH + prey.items + SampleColor + Age_d +
                   SampleVolume_mL + CanopyCover + J.OtherDomain +
                   (1 | Site / Sector),
                 family = gaussian(), prior = priors.nest, cores = 8,
                 iter = 8000, save_all_pars = TRUE, seed = seed,
                 control = list(adapt_delta = 0.999999, max_treedepth = 15),
                 data = even.p)
summary(even.nest)
even.nest$fit

even.nest <- add_criterion(even.nest, criterion = "waic")
even.nest$criteria$waic
even.nest.loo <- loo(even.nest, reloo = TRUE)


## Save it
save(even.nest, file = "StanLM/SMP_LM_even.nest_prok.gz", compress = "gzip")


### With distance instead of site/sector as random effect
priors.dist <- get_prior(J ~ pH + prey.items + SampleColor + Age_d +
                           SampleVolume_mL + CanopyCover + J.OtherDomain +
                           Distance,
                         family = gaussian(), data = even.p)
even.dist <- brm(J ~ pH + prey.items + SampleColor + Age_d +
                   SampleVolume_mL + CanopyCover + J.OtherDomain + Distance,
                 family = gaussian(), prior = priors.dist, cores = 8,
                 iter = 8000, save_all_pars = TRUE, seed = seed,
                 control = list(adapt_delta = 0.9, max_treedepth = 15),
                 data = even.p)
summary(even.dist)
even.dist$fit


even.dist <- add_criterion(even.dist, criterion = "waic")
even.dist$criteria$waic
even.dist.loo <- loo(even.dist, reloo = TRUE)


## Save it
save(even.dist, file = "StanLM/SMP_LM_even.dist_prok.gz", compress = "gzip")


### With interactions
priors.int <- get_prior(J ~ pH * prey.items * SampleColor * Age_d +
                          SampleVolume_mL + CanopyCover + J.OtherDomain +
                          (1 | Site / Sector),
                        family = gaussian(), data = even.p)
even.int <- brm(J ~ pH * prey.items * SampleColor * Age_d + SampleVolume_mL +
                  CanopyCover + J.OtherDomain + (1 | Site / Sector),
                save_all_pars = TRUE, seed = seed, family = gaussian(),
                prior = priors.int, cores = 8, iter = 8000,
                control = list(adapt_delta = 0.99, max_treedepth = 15),
                data = even.p)
summary(even.int)
even.int$fit


even.int <- add_criterion(even.int, criterion = "waic")
even.int$criteria$waic
# even.int.loo <- loo(even.int, reloo = TRUE) # takes too long...


## Save it
save(even.int, file = "StanLM/SMP_LM_even.int_prok.gz", compress = "gzip")


### Model comparison
# load("StanLM/SMP_LM_even.nest_prok.gz")
# load("StanLM/SMP_LM_even.dist_prok.gz")
# load("StanLM/SMP_LM_even.int_prok.gz")


loo_compare(even.nest, even.dist, even.int, criterion = "waic")
loo_compare(even.nest.loo, even.dist.loo)


## Plot random effects
even.nest.ran <- data.frame(ranef(even.nest)[2])
names(even.nest.ran) <- c("Estimate", "SDEstimate", "CrI2.5", "CrI97.5")
even.nest.ran$Site <- substr(rownames(even.nest.ran), 1, 2)
even.nest.ran$Sector <- substr(rownames(even.nest.ran), 4, 4)
even.nest.ran$SiteSector <- paste(even.nest.ran$Site, even.nest.ran$Sector)


even.nest.ran <- even.nest.ran[order(even.nest.ran$Estimate), ]
even.nest.ran$SiteSector <- factor(even.nest.ran$SiteSector,
                                   levels = unique(as.character(
                                     even.nest.ran$SiteSector)))
even.nest.ran$Site <- factor(even.nest.ran$Site,
                             levels = c("CB", "LT", "LE", "LV", "LM"))


ggplot(even.nest.ran, aes(x = SiteSector, y = Estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(ymin = CrI2.5,
                    ymax = CrI97.5), width = 0, color = "black") +
  geom_errorbar(aes(ymin = Estimate - SDEstimate,
                    ymax = Estimate + SDEstimate),
                width = 0, color = "black", size = 1.5) +
  geom_point(aes(fill = Site), pch = 21, size = 4) +
  scale_fill_manual(values = gradCol) +
  theme(legend.position = "none", legend.direction = "horizontal") +
  xlab("") +
  ylab("Estimate Evenness") +
  coord_flip()
# ggsave("SMP_Even_ranef_prok.pdf", width = 11.69 / 2, height = 6)


## Plot fixed effects
even.nest.fix <- data.frame(fixef(even.nest))
names(even.nest.fix) <- c("Estimate", "SDEstimate", "CrI2.5", "CrI97.5")
even.nest.fix$Variable <- rownames(even.nest.fix)


## Rename intercept by SampleColorclear
even.nest.fix$Variable[even.nest.fix$Variable == "Intercept"] <- "SampleColorclear"


even.nest.fix <- even.nest.fix[order(even.nest.fix$Estimate), ]
even.nest.fix$Variable <- factor(even.nest.fix$Variable,
                                 levels = unique(as.character(
                                   even.nest.fix$Variable)))


ggplot(even.nest.fix, aes(x = Variable, y = Estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(ymin = CrI2.5,
                    ymax = CrI97.5), width = 0, color = "black") +
  geom_errorbar(aes(ymin = Estimate - SDEstimate,
                    ymax = Estimate + SDEstimate),
                width = 0, color = "black", size = 1.5) +
  geom_point(pch = 21, fill = "white", size = 4) +
  xlab("") +
  ylab("Estimate evenness") +
  coord_flip()
# ggsave("SMP_Even_fixef_prok.pdf", width = 11.69 / 2, height = 8.27 / 2)


################################################################################
################################################################################
################################################################################
################################################################################
