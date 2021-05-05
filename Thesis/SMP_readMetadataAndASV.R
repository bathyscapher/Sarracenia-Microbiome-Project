################################################################################
################################################################################
################################################################################
################################################################################
### SMP read OTUs & ASV
### Authors: Rachel Korn
### korn@cumulonimbus.at University of Fribourg 2019/2020
################################################################################
### Returned items are (processed in the following order):
## * smp = full dataset
## * smp.g = taxa merged on genus level
## * smp.p = subset to needed samples
## * smp.s = pruned from noise taxa
## * smp.f = remove empty samples (for 18S only)
## * smp.a = abundance filtered
### NOTE: mosses == "with" is useless, as the filter criteria are different


################################################################################
## Colorblind temperature scale for the 5 sites
gradCol <- c("#D55E00", "#E69F00", "#F0E442", "#009E73", "#56B4E9")


################################################################################
readTaxa <- function (primer = c("16S", "18S"),
                      mosses = c("with", "without", "only"),
                      merge.mosses = c(TRUE, FALSE), glom = c(TRUE, FALSE)) {

  ##############################################################################
  ### Metadata: ps = phyloseq object
  readMetadata <- function (ps, primer = c("16S", "18S")) {

    ## Leaf,  moss and DNA extraction metadata
    setwd(smpDir)

    ## Leaf samples
    leaves <- read.table("csv/SMP_LeafSamples_2018.csv",
                         header = TRUE, sep = ",", fill = FALSE)

    leaves$SampleColor <- ordered(leaves$SampleColor,
                                  levels = c("clear", "cloudy", "green", "red",
                                             "brown"))


    ## Canopy cover
    cover <- read.table("csv/SMP_CanopyCover.csv",
                        header = TRUE, sep = "\t", fill = FALSE,
                        check.names = TRUE)

    leaves <- merge(leaves, cover[, c(1, 2, 8)],
                    by.x = c("Site", 'PlantOld'), by.y = c("Site", 'PlantOld'))


    ## Leaf morphometry
    morph <- read.table("csv/SMP_LeafMorphometrics_2018.csv",
                        header = TRUE, sep = "\t", fill = FALSE)


    ## Moss samples
    mosses <- read.table("csv/SMP_MossSamples_2018.csv",
                         header = TRUE, sep = "\t", fill = FALSE)


    ## Positive and negative controls
    controls <- read.table("csv/SMP_2018_MC_NC.csv",
                           header = TRUE, sep = "\t", fill = FALSE)


    ## Geographical position: reduce to 1D with MDS (and plot if needed)
    xy <- read.table("csv/SMP_Plants_XY_m_CH1903.csv", header = TRUE,
                     sep = "\t")
    rownames(xy) <- xy$IDPlant

    xy.mds <- cmdscale(dist(xy[, c(1:2)]), eig = TRUE, k = 1)
    # xy.mds$GOF # should be > 0.8
    distance <- as.data.frame(xy.mds$points)
    colnames(distance) <- "Distance"
    distance$PlantID <- row.names(distance)
    ##
    # xy.mds.df <- as.data.frame(xy.mds$points)
    # colnames(xy.mds.df) <- "Points"
    # xy.mds.df$Plant <- rownames(xy.mds.df)
    # ggplot(xy.mds.df) +
    #   geom_text(aes(x = 0, y = Points, label = Plant), na.rm = TRUE) +
    #   xlim(-0.1, 0.1) +
    #   scale_x_continuous(breaks = seq(-1, 1, 1)) +
    #   xlab("") +
    #   ylab("")
    # ggsave("SMP_Distance_1D.pdf", width = 3, height = 6)


    ## Temperature from data logger
    logger <- read.table("csv/SMP_allLoggersDateTrimmed.csv",
                         header = TRUE, sep = ";", check.names = TRUE)

    options(digits.secs = 0)
    temp.mu <- aggregate(logger$Temperature, by = list(logger$Site), FUN = mean)
    colnames(temp.mu) <- c("Site", "MeanTemperature")


    ## Complete leaf morphology
    leaves.up <- merge(leaves, morph[, c(3, 9:13)],
                       by.x = 'FullIDOld', by.y = 'FullIDOld')


    ## Complete moss data
    colnames(mosses)[11] <- "SampleVolume_mL"
    mosses$Leaf <- paste0(0, mosses$Replicate)
    mosses$SampleColor <- NA
    mosses$PotentialVolume_mL <- NA_real_
    mosses$pH <- NA_real_
    mosses$PitcherLength <- NA_real_
    mosses$CanopyCover <- NA_real_
    mosses$KeelWidth <- NA_real_
    mosses$PitcherWidth <- NA_real_
    mosses$MouthWidth <- NA_real_
    mosses$Comments <- NA


    mosses.leaves <- rbind(mosses[, c(1:4, 6:7, 11, 13:22)],
                           leaves.up[, c(1, 2, 6:13, 16, 44:49)])
    mosses.leaves$Leaf <- as.factor(mosses.leaves$Leaf)


    ## Complete controls
    controls$SampleVolume_mL <- NA_real_
    controls$SampleColor <- NA
    controls$PotentialVolume_mL <- NA_real_
    controls$pH <- NA_real_
    controls$PitcherLength <- NA_real_
    controls$CanopyCover <- NA_real_
    controls$KeelWidth <- NA_real_
    controls$PitcherWidth <- NA_real_
    controls$MouthWidth <- NA_real_
    controls$Comments <- NA


    ## Merge all
    mosses.leaves.ctr <- rbind(controls[, c(1:2, 5:7, 9:10, 27:36)],
                               mosses.leaves)


    smpMeta <- merge(mosses.leaves.ctr, temp.mu, by = "Site", all.x = TRUE)
    smpMeta$PlantID <- paste(smpMeta$Site, smpMeta$SectorPlant, sep = "-")
    smpMeta <- merge(smpMeta, distance, by = "PlantID", all.x = TRUE)


    ## Remove 16S (= P and R) or 18S (= U and V) MC
    if (primer == "16S")
      {smpMeta <- smpMeta[!grepl("U|V", smpMeta$Succession), ]}
    if (primer == "18S")
      {smpMeta <- smpMeta[!grepl("P|R", smpMeta$Succession), ]}

    smpMeta <- droplevels(smpMeta)


    smpMeta$Site <- ordered(smpMeta$Site, levels = c("CB", "LT", "LE", "LV",
                                                     "LM", "MC", "NC"))
    smpMeta$Sector <- as.factor(smpMeta$Sector)


    ## Add altitude
    sites <- read.table("csv/SMP_SarraceniaField.csv", header = TRUE, sep = ";")
    smpMeta <- merge(smpMeta, sites[, c(1, 8)], by = "Site", all.x = TRUE)


    ## Replace Succession with habitat age [d] (= days since first rain event)
    age <- read.table("csv/SMP_SampleDates.csv", header = TRUE,
                      sep = "\t", check.names = TRUE)
    smpMeta <- merge(smpMeta, age, by = c("Site", "Succession"), all.x = TRUE)
    smpMeta$Succession <- as.factor(smpMeta$Succession)


    ## Remove spacers in sample names
    smpMeta$FullID <- gsub("[-|_]", "", smpMeta$FullID)


    ## Prey items
    prey <- read.table("csv/SMP_Prey_clean.csv", sep = "\t", header = TRUE)


    prey$prey.items <- rowSums(prey[, -c(12, 20)])
    prey


    smpMeta <- merge(smpMeta, prey[, c(12, 21)], by = "FullID", all = TRUE)


    ## Duplicate first all moss samples, then all samples and add a resequenced
    ##  tag ("r") to the end of the sample names
    moss <- smpMeta[smpMeta$Succession == "M", ]
    moss$FullID <- gsub("M$", "Q", moss$FullID)

    smpMeta <- rbind(smpMeta, moss)


    smpMeta.r <- smpMeta
    smpMeta.r$FullID <- gsub("$", "r", smpMeta.r$FullID)

    smpMeta <- rbind(smpMeta, smpMeta.r)


    ## Rename succession
    if (primer == "16S")
      {
      levels(smpMeta$Succession) <- c("Early", "Late", "Moss",
                                      "NC Pilot", "NC SMP",
                                      "16S MC Pilot", "16S MC SMP")
      }
    if (primer == "18S")
      {
      levels(smpMeta$Succession) <- c("Early", "Late", "Moss",
                                      "NC Pilot", "NC SMP",
                                      "18S MC Pilot", "18S MC SMP")
      }


    ## Match with sample names from dada2
    smpMeta <- smpMeta[match(sample_names(smp), smpMeta$FullID), ]


    ## Convert data.frame to sample_data and add row.names
    smpMeta <- sample_data(smpMeta)
    row.names(smpMeta) <- smpMeta$FullID

    return(smpMeta)
  }


  ##############################################################################
  ### Read taxa (ASV)
  setwd(asvDir)

  if (primer == "16S")
  {
    seqtab.nochim <- readRDS(paste(asvDir,
                                   "SMP_16S_reseq/seqtab.nochim_SMP16S.rds",
                                   sep = "/"))
    taxaid <- readRDS(paste(asvDir,
                            "SMP_16S_reseq/taxa_SMP16S.rds",
                            sep = "/"))
    # taxaid <- readRDS("taxid_SMP16S.rds")
  }
  if (primer == "18S")
  {
    seqtab.nochim <- readRDS(paste(asvDir,
                                   "SMP_18S_reseq/seqtab.nochim_SMP18S.rds",
                                   sep = "/"))
    taxaid <- readRDS("SMP_18S_reseq/taxa_SMP18S.rds")
    # taxaid <- readRDS("taxid_SMP18S.rds")
  }

  smp <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),
                  tax_table(taxaid))


  dna <- Biostrings::DNAStringSet(taxa_names(smp))
  names(dna) <- taxa_names(smp)
  smp <- merge_phyloseq(smp, dna)
  taxa_names(smp) <- paste0("ASV", seq(ntaxa(smp)))


  ### Rename taxonomic ranks
  colnames(tax_table(smp)) <- c("Domain", "Phylum", "Class", "Order",
                                "Family", "Genus")


  ##############################################################################
  ### Transpose (sometimes the OTU table is transposed... :-|)
  if (taxa_are_rows(smp))
    {otu_table(smp) <- t(otu_table(smp))}


  ##############################################################################
  ### Read metadata
  smpMeta <- readMetadata(ps = smp, primer = primer)

  smp <- merge_phyloseq(smp, smpMeta)


  ### Merge resequenced samples (for ASV only)
  sample_data(smp)$FullID <- gsub("r", "", sample_data(smp)$FullID)
  smp <- merge_samples(smp, "FullID")


  sample_data(smp) <- smpMeta


  ############################################################################
  ### Aggregate taxa by genus
  if(glom == TRUE){
    smp.g  <- tax_glom(smp, taxrank = rank_names(smp)[6],
                       NArm = FALSE, bad_empty = c("", " ", "\t"))
  } else {
    smp.g <- smp
  }

  ##############################################################################
  ### Remove unwanted samples
  if (mosses == "with")
  {
    ### Merge resequenced moss samples
    sample_data(smp.g)$FullID <- gsub("Q", "M", sample_data(smp.g)$FullID)
    smp.g <- merge_samples(smp.g, "FullID")


    sample_data(smp.g) <- smpMeta


    if (primer == "16S")
    {
      smp.p <- prune_samples(smp.g@sam_data$Succession != "NC Pilot" &
                               smp.g@sam_data$Succession != "NC SMP" &
                               smp.g@sam_data$Succession != "16S MC Pilot" &
                               smp.g@sam_data$Succession != "16S MC SMP", smp.g)
    }
    if (primer == "18S")
    {
      smp.p <- prune_samples(smp.g@sam_data$Succession != "NC Pilot" &
                               smp.g@sam_data$Succession != "NC SMP" &
                               smp.g@sam_data$Succession != "18S MC Pilot" &
                               smp.g@sam_data$Succession != "18S MC SMP", smp.g)
    }
  }

  if (mosses == "without")
  {
    if (primer == "16S")
    {
      smp.p <- prune_samples(smp.g@sam_data$Succession != "NC Pilot" &
                               smp.g@sam_data$Succession != "NC SMP" &
                               smp.g@sam_data$Succession != "16S MC Pilot" &
                               smp.g@sam_data$Succession != "16S MC SMP" &
                               smp.g@sam_data$Succession != "Moss", smp.g)
    }
    if (primer == "18S")
    {
      smp.p <- prune_samples(smp.g@sam_data$Succession !=  "NC Pilot" &
                               smp.g@sam_data$Succession !=  "NC SMP" &
                               smp.g@sam_data$Succession !=  "18S MC Pilot" &
                               smp.g@sam_data$Succession != "18S MC SMP" &
                               smp.g@sam_data$Succession != "Moss", smp.g)
    }
  }


  ##############################################################################
  ## Mosses

  if(mosses == "with" |  mosses == "without"){

    ### Remove spurious taxa
    if(primer == "16S")
    {
      smp.s <- subset_taxa(smp.p, !(Domain %in% c("unknown", "Eukaryota", NA) |
                                      Phylum %in% c("Eukaryota_unclassified",
                                                    NA) |
                                      Order %in% c("Chloroplast") |
                                      Family %in% c("Mitochondria")))
    }
    if(primer == "18S")
    {
      smp.s <- subset_taxa(smp.p, !(Domain %in% c("Bacteria", "unknown") |
                                      Phylum %in% c("Eukaryota_unclassified",
                                                    "Myxogastria",
                                                    "Apicomplexa",
                                                    "Mollusca", "Vertebrata",
                                                    "Microsporidia",
                                                    "Mucoromycota",
                                                    "Archaeorhizomycetes", NA) |
                                      Class %in% c("Insecta", "Ellipura",
                                                   "Embryophyta", "Arachnida",
                                                   "Heterophyidae",
                                                   "Ichthyophonae",
                                                   "Arthropoda_unclassified",
                                                   "unclassified_Hexapoda",
                                                   "Ascomycota_unclassified",
                                                   "Agaricomycetes",
                                                   "Basidiomycota_unclassified",
                                                   "Exobasidiomycetes",
                                                   "Microbotryomycetes",
                                                   "Dothideomycetes",
                                                   "Pucciniomycetes",
                                                   "Spiculogloeomycetes",
                                                   "Ustilaginomycetes",
                                                   "Taphrinomycetes",
                                                   "Lecanoromycetes",
                                                   "Basidiomycota_unclassified",
                                                   "Sordariomycetes",
                                                   "Pezizomycetes",
                                                   "Entorrhizomycetes",
                                                   "Agaricostilbomycetes",
                                                   "Orbiliomycetes",
                                                   "Eurotiomycetes",
                                                   "Leotiomycetes",
                                                   "Entomophthoromycetes",
                                                   "Dacrymycetes") |
                                      Order %in% c("Filobasidiales",
                                                   "Spizellomycetales",
                                                   "Rhytismatales",
                                                   "Dothideomycetes",
                                                   "Mytilinidiales",
                                                   "Naohideales",
                                                   "Pleosporales",
                                                   "Kickxellales",
                                                   "Archaeorhizomycetes",
                                                   "Atractiellomycetes",
                                                   "Trichosporonales") |
                                      Family %in% c("Pleosporaceae",
                                                    "Mrakiaceae",
                                                    "Teratosphaeriaceae",
                                                    "Leotiaceae",
                                                    "Pleosporales_fa",
                                                    "Didymellaceae",
                                                    "Phaeosphaeriaceae",
                                                    "Rhynchogastremataceae",
                                                    "Phaeotremellaceae") |
                                      Genus %in% c("Lecophagus", "Genolevuria",
                                                   "Trimorphomyces") |
                                      Phylum == "Arthropoda" & Class %in% NA |
                                      Phylum == "Ascomycota" & Class %in% NA |
                                      Phylum == "Basidiomycota" & Class %in% NA))


      ## Remove the samples without target taxa
      smp.s <- prune_samples(sample_sums(smp.s) > 0, smp.s)
    }


    ############################################################################
    ### Remove empty taxa
    any(taxa_sums(smp.s) == 0)
    sum(taxa_sums(smp.s) == 0)
    smp.s <- prune_taxa(taxa_sums(smp.s) > 0, smp.s)


    ############################################################################
    ### Convert abundance to relative abundance
    smp.r <- transform_sample_counts(smp.s, function(otu) {otu / sum(otu)})
  }


  ############################################################################
  ### Abundance filtering
  if(mosses == "without"){
  # ifelse(method == "OTU", threshold <- 0.0001, threshold <- 0.00001)

  smp.a <- filter_taxa(smp.r, function(otu) {mean(otu) > 0.00001},
                       prune = TRUE)
  }
  # smp.a
  # plot(rowSums(otu_table(smp.a)))


  if(mosses == "with"){
    ## Split by moss and pitcher samples
    smp.r.moss <- prune_samples(smp.r@sam_data$Succession ==  "Moss", smp.r)
    smp.r.pitch <- prune_samples(smp.r@sam_data$Succession !=  "Moss", smp.r)


    ## Filter pitcher samples
    # ifelse(method == "OTU", threshold <- 0.0001, threshold <- 0.00001)

    smp.a.pitch <- filter_taxa(smp.r.pitch,
                               function(otu) {mean(otu) > 0.00001},
                               prune = TRUE)


    ## Filter moss samples
    # ifelse(method == "OTU", threshold <- 0.00001, threshold <- 0.000001)

    smp.a.moss <- filter_taxa(smp.r.moss, function(otu) {mean(otu) > 0.000001},
                              prune = TRUE)

    ## Merge them
    smp.a <- merge_phyloseq(smp.a.moss, smp.a.pitch)

  }

  ##############################################################################
  ### Same for mosses
  if(mosses == "only")
  {
    moss <- prune_samples(smp@sam_data$Succession == "Moss", smp.g)


    ### Merge resequenced moss samples (samples *M and *Q) and restore metadata
    sample_data(moss)$FullID <- gsub("Q", "M", sample_data(moss)$FullID)
    moss <- merge_samples(moss, "FullID")


    sample_data(moss) <- smpMeta


    if (merge.mosses == TRUE)
    {
      ### Merge sector-level replicates of moss samples: create new ID to merge
      ### replicates, merge and restore metadata
      sample_data(moss)$FullIDMerge <- paste0(sample_data(moss)$Site,
                                              sample_data(moss)$SectorPlant,
                                              "00M")


      moss.m <- merge_samples(moss, "FullIDMerge")

      sample_names(moss.m) <- gsub("0M$", "1M", sample_names(moss.m))
      sample_data(moss.m) <- smpMeta
    }
    if (merge.mosses == "FALSE")
    {
      moss.m <- moss
    }


    ### Remove spurious taxa
    if(primer == "16S")
    {
      moss.s <- subset_taxa(moss.m, !(Domain %in% c("unknown", "Eukaryota",
                                                    NA) |
                                        Phylum %in% c("Eukaryota_unclassified",
                                                      NA) |
                                        Order %in% c("Chloroplast") |
                                        Family %in% c("Mitochondria")))
    }

    if(primer == "18S")
    {
      moss.s <- subset_taxa(moss.m, !(Domain %in% c("Bacteria", "unknown") |
                                        Phylum %in% c("Eukaryota_unclassified",
                                                      "Myxogastria",
                                                      "Apicomplexa",
                                                      "Mollusca", "Vertebrata",
                                                      "Microsporidia",
                                                      "Mucoromycota",
                                                      "Archaeorhizomycetes", NA) |
                                        Class %in% c("Insecta", "Ellipura",
                                                     "Embryophyta", "Arachnida",
                                                     "Heterophyidae",
                                                     "Ichthyophonae",
                                                     "Arthropoda_unclassified",
                                                     "unclassified_Hexapoda",
                                                     "Ascomycota_unclassified",
                                                     "Agaricomycetes",
                                                     "Basidiomycota_unclassified",
                                                     "Exobasidiomycetes",
                                                     "Microbotryomycetes",
                                                     "Dothideomycetes",
                                                     "Pucciniomycetes",
                                                     "Spiculogloeomycetes",
                                                     "Ustilaginomycetes",
                                                     "Taphrinomycetes",
                                                     "Lecanoromycetes",
                                                     "Basidiomycota_unclassified",
                                                     "Sordariomycetes",
                                                     "Pezizomycetes",
                                                     "Entorrhizomycetes",
                                                     "Agaricostilbomycetes",
                                                     "Orbiliomycetes",
                                                     "Eurotiomycetes",
                                                     "Leotiomycetes",
                                                     "Entomophthoromycetes",
                                                     "Dacrymycetes") |
                                        Order %in% c("Filobasidiales",
                                                     "Spizellomycetales",
                                                     "Rhytismatales",
                                                     "Dothideomycetes",
                                                     "Mytilinidiales",
                                                     "Naohideales",
                                                     "Pleosporales",
                                                     "Kickxellales",
                                                     "Archaeorhizomycetes",
                                                     "Atractiellomycetes",
                                                     "Trichosporonales") |
                                        Family %in% c("Pleosporaceae",
                                                      "Mrakiaceae",
                                                      "Teratosphaeriaceae",
                                                      "Leotiaceae",
                                                      "Pleosporales_fa",
                                                      "Didymellaceae",
                                                      "Phaeosphaeriaceae",
                                                      "Rhynchogastremataceae",
                                                      "Phaeotremellaceae") |
                                        Genus %in% c("Lecophagus",
                                                     "Genolevuria",
                                                     "Trimorphomyces") |
                                        Phylum == "Arthropoda" & Class %in% NA |
                                        Phylum == "Ascomycota" & Class %in% NA |
                                        Phylum == "Basidiomycota" &
                                          Class %in% NA))
    }


    ## Remove empty taxa
    any(taxa_sums(moss.s) == 0)
    sum(taxa_sums(moss.s) == 0)
    moss.s <- prune_taxa(taxa_sums(moss.s) > 0, moss.s)


    ## Relative abundance
    moss.r <- transform_sample_counts(moss.s, function(otu) {otu / sum(otu)})


    ## Abundance filtering
    # ifelse(method == "OTU", threshold <- 0.00001, threshold <- 0.000001)

    moss.a <- filter_taxa(moss.r, function(x) {mean(x) > 0.000001}, TRUE)


    ## Remove empty taxa
    moss.a <- prune_taxa(taxa_sums(moss.a) > 0, moss.a)
    }


  ##############################################################################
  ### Return variables
  if(mosses == "only")
  {
    moss.l <- list(smpMeta = smpMeta, moss = moss,
                   moss.m = moss.m, moss.s = moss.s,  moss.r = moss.r,
                   moss.a = moss.a)
  } else {
    smp.l <- list(smpMeta = smpMeta, smp = smp, smp.g = smp.g, smp.p = smp.p,
                  smp.s = smp.s, smp.r = smp.r, smp.a = smp.a)
  }

  ifelse(mosses == "only", return(moss.l), return(smp.l))
  }


################################################################################
### Read metadata
readMetadata <- function (ps, primer = c("16S", "18S")) {
  ## Leaf,  moss and DNA extraction metadata
  setwd(mcDir)

  ## Leaf samples
  leaves <- read.table("csv/SMP_LeafSamples_2018.csv",
                       header = TRUE, sep = ",", fill = FALSE)

  leaves$SampleColor <- ordered(leaves$SampleColor,
                                levels = c("clear", "cloudy", "green", "red",
                                           "brown"))

  ## Canopy cover
  cover <- read.table("csv/SMP_CanopyCover.csv",
                      header = TRUE, sep = "\t", fill = FALSE,
                      check.names = TRUE)

  leaves <- merge(leaves, cover[, c(1, 2, 8)],
                  by.x = c("Site", 'PlantOld'), by.y = c("Site", 'PlantOld'))

  ## Leaf morphometry
  morph <- read.table("csv/SMP_LeafMorphometrics_2018.csv",
                      header = TRUE, sep = "\t", fill = FALSE)

  ## Moss samples
  mosses <- read.table("csv/SMP_MossSamples_2018.csv",
                       header = TRUE, sep = "\t", fill = FALSE)

  ## Positive and negative controls
  controls <- read.table("csv/SMP_2018_MC_NC.csv",
                         header = TRUE, sep = "\t", fill = FALSE)

  ## Geographical position: reduce to 1D with MDS (and plot if needed)
  xy <- read.table("csv/SMP_Plants_XY_m_CH1903.csv", header = TRUE,
                   sep = "\t")
  rownames(xy) <- xy$IDPlant

  xy.mds <- cmdscale(dist(xy[, c(1:2)]), eig = TRUE, k = 1)
  # xy.mds$GOF # should be > 0.8
  distance <- as.data.frame(xy.mds$points)
  colnames(distance) <- "Distance"
  distance$PlantID <- row.names(distance)

  ## Temperature from data logger
  logger <- read.table("csv/SMP_allLoggersDateTrimmed.csv",
                       header = TRUE, sep = ";", check.names = TRUE)

  options(digits.secs = 0)
  temp.mu <- aggregate(logger$Temperature, by = list(logger$Site), FUN = mean)
  colnames(temp.mu) <- c("Site", "MeanTemperature")

  ## Complete leaf morphology
  leaves.up <- merge(leaves, morph[, c(3, 9:13)],
                     by.x = 'FullIDOld', by.y = 'FullIDOld')

  ## Complete moss data
  colnames(mosses)[11] <- "SampleVolume_mL"
  mosses$Leaf <- paste0(0, mosses$Replicate)
  mosses$SampleColor <- NA
  mosses$PotentialVolume_mL <- NA_real_
  mosses$pH <- NA_real_
  mosses$PitcherLength <- NA_real_
  mosses$CanopyCover <- NA_real_
  mosses$KeelWidth <- NA_real_
  mosses$PitcherWidth <- NA_real_
  mosses$MouthWidth <- NA_real_
  mosses$Comments <- NA

  mosses.leaves <- rbind(mosses[, c(1:4, 6:7, 11, 13:22)],
                         leaves.up[, c(1, 2, 6:13, 16, 44:49)])
  mosses.leaves$Leaf <- as.factor(mosses.leaves$Leaf)

  ## Complete controls
  controls$SampleVolume_mL <- NA_real_
  controls$SampleColor <- NA
  controls$PotentialVolume_mL <- NA_real_
  controls$pH <- NA_real_
  controls$PitcherLength <- NA_real_
  controls$CanopyCover <- NA_real_
  controls$KeelWidth <- NA_real_
  controls$PitcherWidth <- NA_real_
  controls$MouthWidth <- NA_real_
  controls$Comments <- NA

  ## Merge all
  mosses.leaves.ctr <- rbind(controls[, c(1:2, 5:7, 9:10, 27:36)],
                             mosses.leaves)

  mcMeta <- merge(mosses.leaves.ctr, temp.mu, by = "Site", all.x = TRUE)
  mcMeta$PlantID <- paste(mcMeta$Site, mcMeta$SectorPlant, sep = "-")
  mcMeta <- merge(mcMeta, distance, by = "PlantID", all.x = TRUE)

  ## Remove 16S (= P and R) or 18S (= U and V) MC
  if (primer == "16S")
  {mcMeta <- mcMeta[!grepl("U|V", mcMeta$Succession), ]}
  if (primer == "18S")
  {mcMeta <- mcMeta[!grepl("P|R", mcMeta$Succession), ]}

  mcMeta <- droplevels(mcMeta)


  mcMeta$Site <- ordered(mcMeta$Site, levels = c("CB", "LT", "LE", "LV",
                                                 "LM", "MC", "NC"))
  mcMeta$Sector <- as.factor(mcMeta$Sector)

  ## Add altitude
  sites <- read.table("csv/SMP_SarraceniaField.csv", header = TRUE, sep = ";")
  mcMeta <- merge(mcMeta, sites[, c(1, 8)], by = "Site", all.x = TRUE)


  ## Replace Succession with habitat age [d] (= days since first rain event)
  age <- read.table("csv/SMP_SampleDates.csv", header = TRUE,
                    sep = "\t", check.names = TRUE)
  mcMeta <- merge(mcMeta, age, by = c("Site", "Succession"), all.x = TRUE)
  mcMeta$Succession <- as.factor(mcMeta$Succession)

  ## Remove spacers in sample names
  mcMeta$FullID <- gsub("[-|_]", "", mcMeta$FullID)

  ## Prey items
  prey <- read.table("csv/SMP_Prey_clean.csv", sep = "\t", header = TRUE)

  prey$prey.items <- rowSums(prey[, -c(12, 20)])
  prey

  mcMeta <- merge(mcMeta, prey[, c(12, 21)], by = "FullID", all = TRUE)

  ## Duplicate first all moss samples, then all samples and add a resequenced
  ##  tag ("r") to the end of the sample names
  moss <- mcMeta[mcMeta$Succession == "M", ]
  moss$FullID <- gsub("M$", "Q", moss$FullID)

  mcMeta <- rbind(mcMeta, moss)

  mcMeta.r <- mcMeta
  mcMeta.r$FullID <- gsub("$", "r", mcMeta.r$FullID)

  mcMeta <- rbind(mcMeta, mcMeta.r)


  ## Rename succession
  if (primer == "16S")
  {
    levels(mcMeta$Succession) <- c("Early", "Late", "Moss",
                                   "NC Pilot", "NC SMP",
                                   "16S MC Pilot", "16S MC SMP")
  }
  if (primer == "18S")
  {
    levels(mcMeta$Succession) <- c("Early", "Late", "Moss",
                                   "NC Pilot", "NC SMP",
                                   "18S MC Pilot", "18S MC SMP")
  }


  ## Match with sample names from dada2
  mcMeta <- mcMeta[match(sample_names(mc), mcMeta$FullID), ]

  ## Convert data.frame to sample_data and add row.names
  mcMeta <- sample_data(mcMeta)
  row.names(mcMeta) <- mcMeta$FullID

  return(mcMeta)
}


################################################################################
################################################################################
################################################################################
################################################################################
