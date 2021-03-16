################################################################################
################################################################################
################################################################################
################################################################################
### SMP samples prey: format and export
### Author: korn@cumulonimbus.at University of Fribourg 2020
################################################################################


library("ggplot2")
theme_set(theme_bw(base_size = 20) +
            theme(rect = element_rect(fill = "transparent")))
library("plyr")
library("data.table")
library("reshape2")


rm(list = ls())
setwd("~/Sarracenia-Microbiome-Project/Thesis/csv")


set.seed(34706)


################################################################################
## Colorblind temperature scale for the 5 sites
gradCol <- c("#D55E00", "#E69F00", "#F0E442", "#009E73", "#56B4E9")


## Read csv files
prey.l <- lapply(list.files(pattern = glob2rx("SMP_2020*.csv")), read.delim)


## Reformat and bring into a single df
prey.l <- lapply(prey.l, function(df){
  names(df) <- df[1, ]

  df <- as.data.frame(t(df[-c(1:3), ]))

  ## Remove empty columns, rename columns and remove first three rows
  df <- Filter(function(x) !(all(x == "")), df)
  names(df) <- df[1, ]
  df <- df[-c(1:3), ]

  ## Remove comments and whitespace, replace NA with 0 & convert into integer
  df[] <- lapply(df, gsub, pattern = '[[:alpha:]]|_|,|\\s', replacement = '')

  ## Type conversion
  df[] <- lapply(df, as.integer)
  # str(df)

  ## Sum by order
  df <- as.data.frame(t(rowsum(t(df), group = colnames(df), na.rm = TRUE)))

  ## Add ID
  df$FullID <- rownames(df)
  df$FullID <- gsub("_|-", "", df$FullID)

  return(df)
  })


lapply(prey.l, names)
lapply(prey.l, rownames)
class(prey.l[[1]])


## Merge into one df
prey <- rbind.fill(prey.l[[1]], prey.l[[2]], prey.l[[3]], prey.l[[4]],
                   prey.l[[5]])
prey[is.na(prey)] <- 0
names(prey)
str(prey)


## Add site
prey$Site <- ordered(substr(prey$FullID, 1, 2),
                     levels = c("CB", "LT", "LE", "LV", "LM"))


## Show counts per order
colSums(prey[, !names(prey) %in% c("FullID", "Site")], na.rm = TRUE)


## Remove empty columns (Orthoptera)
prey <- Filter(function(x) !(all(x == 0)), prey)


## Export
# write.table(prey, "SMP_Prey_clean.csv", sep = "\t", row.names = FALSE)


## Convert into long form
prey.m <- reshape2::melt(prey, id.vars = c("Site", "FullID"))
names(prey.m)[3:4] <- c("Order", "Count")
names(prey.m)


str(prey.m)
levels(prey.m$Order)


### Prey per site and order
ggplot(data = prey.m) +
  geom_bar(aes(x = reorder(Order, -Count, sum), y = Count, fill = Site),
           stat = "identity") +
  theme(legend.position = "top", legend.direction = "horizontal",
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_trans(x = "identity", y = "log1p") +
  scale_fill_manual(values = gradCol) +
  scale_y_continuous(breaks = c(0, 1, 10, 100, 1000, 2000)) +
  ylab(expression(log[e](Counts))) +
  xlab("")
# ggsave("../SMP_Prey.pdf", width = 11.69, height = 6.27)


## Number of prey per sample
prey.sum <- aggregate(prey.m$Count, by = list(FullID = prey.m$FullID),
                      FUN = sum, na.rm = TRUE)
prey.sum$Site <- gsub(".{5}$", "", prey.sum$FullID)
prey.sum$Site <- ordered(prey.sum$Site,
                         levels = c("CB", "LT", "LE", "LV", "LM"))
prey.sum$Succession <- as.factor(gsub("^.{6}", "", prey.sum$FullID))
levels(prey.sum$Succession) <- c("Early", "Late")


ggplot(prey.sum, aes(x = Succession, y = x)) +
  geom_boxplot(size = 0.2, outlier.shape = NA) +
  geom_jitter(aes(col = Site, shape = Succession), height = 0, width = 0.3) +
  facet_grid( ~ Site, scales = "free_x") +
  scale_color_manual(values = gradCol) +
  theme(legend.position = "top", legend.direction = "horizontal",
        axis.text.x = element_blank()) +
  scale_y_continuous(breaks = seq(0, 90, 20)) +
  xlab("") +
  ylab("Count")
# ggsave("../SMP_Prey_Counts.pdf", width = 11.69, height = 5)


################################################################################
################################################################################
################################################################################
################################################################################
