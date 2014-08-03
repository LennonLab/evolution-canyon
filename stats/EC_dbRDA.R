################################################################################
#                                                                              #
# Evolution Canyon Project: Microbial Community RDA                            #
#                                                                              #
################################################################################
#                                                                              #
# Written by: Jay Lennon and Mario Muscarella                                  #
#                                                                              #
# Last update: 2014/08/01                                                      #
#                                                                              #
################################################################################
#                                                                              #
# Notes: This code performs the RDA function on the Evolution Canyon        #
#        emperical data                                                        #
#                                                                              #
# Recent Changes:                                                              #
#         1. Calculate Variance Explained                                      #
#         2. Confirm Normalizations                                            #
#         3. Removed Singletons                                                #
#                                                                              #
# Future Changes (To-Do List):                                                 #
#         1. ID influential Taxa                                               #
#         2. Load Taxonomy Data                                                #
#         4. Confirm Problamatic Taxa (EC-2A-D, EC-2A-R, EC-2C-R, EC-2D-R)     #
#         5. Compare with PERMANOVA (See other code)                           #
#         6. Paired PERMANOVA       (See other code)                           #
#                                                                              #
################################################################################

# Setup Work Environment
rm(list=ls())
setwd("~/GitHub/evolution-canyon")
source("./bin/DiversityFunctions.r")
require("vegan")

shared <- "./mothur/EC.bac.final.shared"
design <- "./data/design.txt"
level  <-  "0.03"



  # Import Site by OTU Matrix
  ec_data <- read.otu(shared, "0.03")
  design <- read.delim(design, header=T, row.names=1)

  # Note: owing to amplification issues, we only sequenced 76 of 80 samples.
  # The four samples not included are C-1E-R, EC-2G-R, EC-2J-R, EC-6I-D

  # Remove problematic samples (EC_2A_D, EC_2A_R, EC_2C_R, EC_2D_R)
  ec_data.tmp <- ec_data[-c(9,20:21,24:27,32,37,74),]

  # Remove OTUs with less than ten observations
  ec_data_red <- ec_data.tmp[,colSums(ec_data.tmp) >= 100]

  # design matrix w/o problem samples & pairs
  design_red <- design[rownames(ec_data_red),]
  design_red$paired <- c(rep(seq(1:14), each=2), rep(seq(1:19), each=2))

  # Create factors for model
  slope <- design_red$slope # factor 1
  molecule <- design_red$molecule # factor 2
  paired <- design_red$paired
  site <- design_red$site
  station <- design_red$station
  slope.molecule <- data.frame(cbind(as.character(slope),
    as.character(molecule))) # Y matrix with factor 1 and 2
  slope.molecule.concat <- do.call(paste, c(slope.molecule[c("X1", "X2")],
    sep = "")) # create unique treat ID vector
  pair.station <- c(rep(seq(1:9), each=2), rep(seq(1:5), each=2),
    rep(seq(1:10), each=2), rep(seq(1:9), each=2))
  # To include "questionable" data uncomment the following
  # pair.station <- c(rep(seq(1:9), each=2), rep(seq(1:8), each=2),
  #    rep(seq(1:10), each=10), rep(seq(1:9), each=2))
  # Create a vector of molecules by station
  station.molecule.concat <- paste(station, molecule, sep = "")

  # Calculate Presense Absence
  dataPA <- (ec_data_red > 0)*1

  # Calculating Relative Abundance
  dataREL <- ec_data_red
  for(i in 1:nrow(ec_data_red)){
    dataREL[i,] = ec_data_red[i,]/sum(ec_data_red[i,])
    }

  # Log Transform Relative Abundance
  dataREL.l <- decostand(dataREL,method="log")

  # Chord Transformation
  dataREL.c <- decostand(dataREL, method="normalize")

  # Hellinger Transformation
  dataREL.h <- decostand(dataREL, method="hellinger")

  # Create Distance Matrix with bray (deafault), manhattan, euclidean, canberra, bray, kulczynski, jaccard, gower, altGower, morisita, horn, mountford, raup, binomial, or chao. Most should be part of vegan, but possilbly 'labdsv' or 'BiodiversityR' packages
  samplePA.dist <- vegdist(dataHell,method="bray")
  sampleREL.dist <- vegdist(dataHell,method="bray")

  # Distance Based Redundancy Analysis
  dbRDA <- capscale(dataPA ~ slope + molecule, distance="bray")

  head(summary(dbRDA))
  RsquareAdj(dbRDA)

  plot(dbRDA, display = c("wa", "cn"))
  ordiellipse(dbRDA, groups=slope.molecule.concat, label=T, kind="sd", conf=0.95)

  anova(dbRDA, by="terms", permu=500)
  
  varpart(dataREL.l,  ~ slope, ~ molecule)