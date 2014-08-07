################################################################################
#                                                                              #
#	Evolution Canyon Project: Microbial Community Annalysis by                   #
#                           Distance Based Redundancy Analysis                 #
#                                                                              #
# This code was used to develop the simulation results figures for ESA         #
#                                                                              #
################################################################################
#                                                                              #
#	Written by: Mario Muscarella                                                 #
#                                                                              #
#	Last update: 2014/08/04                                                      #
#                                                                              #
################################################################################

rm(list=ls())
setwd("~/GitHub/evolution-canyon")
source("./bin/DiversityFunctions.r")
require("vegan")

# Plot 1 - Blank Ordination for descriptive purposes ###########################
png(file="./plots/blank_ordination.png",
    width=600, height=600, antialias = "cleartype", res=96)

  # Plot Parameters
  par(mfrow=c(1,1), mar=c(5,6,1,1))

  # Initiate Plot
  plot(NA,
    xlim = c(-0.1, 0.1),ylim= c(-0.1, 0.1), pch=16, cex=2.0,
    type="n",ann=FALSE, xaxt="n", yaxt = "n", cex.lab=1.5, cex.axis=1.2)
  axis(side=1, las=1)
  axis(side=2, las=1)
  abline(h=0, lty="dotted")
  abline(v=0, lty="dotted")
  mtext(side = 1, text = "Axis 1 (% Variation Explained)", line = 3, cex=1.5)
  mtext(side = 2, text = "Axis 2 (% Variation Explained)", line = 3.5, cex=1.5)
  box(lwd=2)

dev.off()

# Plot 2 - Condition 1 #########################################################
png(file="./plots/Condition1_ESA.png",
    width=600, height=600, antialias = "cleartype", res=96)

  # Inputs
  shared     = "./microbide/SbyS/Condition1.txt"
  cutoff     = "0.03"
  design     = "./data/simmy.design.txt"

  # Import Site by OTU Matrix
  ec_data.sim <- read.otu(shared, "0.03")
  design <- read.delim(design, header=T, row.names=1)

  # Remove OTUs with less than ten observations
  ec_data_red <- ec_data.sim[,colSums(ec_data.sim) >= 10]

  # Create factors for model
  slope <- design$slope # factor 1
  molecule <- design$molecule # factor 2
  paired <- design$paired
  slope.molecule <- data.frame(cbind(as.character(slope),
    as.character(molecule))) # Y matrix with factor 1 and 2
  slope.molecule.concat <- do.call(paste, c(slope.molecule[c("X1", "X2")],
    sep = "")) # create unique treat ID vector

  # Calculate Presense Absence
  dataPA <- (ec_data_red > 0)*1

  # Calculating Relative Abundance
  dataREL <- ec_data_red
  for(i in 1:nrow(ec_data_red)){
    dataREL[i,] = ec_data_red[i,]/sum(ec_data_red[i,])
    }

  # Log Transform Relative Abundance
  dataREL.l <- decostand(dataREL,method="log")

  # Distance Based Redundancy Analysis
  dbRDA <- capscale(dataREL ~ slope + molecule+ Condition(paired), distance="bray")

  # Percent Variance Explained Using PCoA (Axis 1,2,3)
  explainvar1 <- round(dbRDA$CCA$eig[1]/sum(dbRDA$CCA$eig, dbRDA$CA$eig)*100,2)
  explainvar2 <- round(dbRDA$CCA$eig[2]/sum(dbRDA$CCA$eig, dbRDA$CA$eig)*100,2)

  RDA <- as.data.frame(dbRDA$CCA$wa.eig)
  RDA$molecule <- design$molecule
  RDA$slope <- design$slope
  RDA$labs <- slope.molecule.concat


  # Plot Parameters
  par(mfrow=c(1,1), mar=c(5,6,1,1))

  # Initiate Plot
  plot(RDA$RDA1, RDA$RDA2,
    xlim = c(-0.1, 0.1),ylim= c(-0.1, 0.1), pch=16, cex=2.0,
    type="n",ann=FALSE, xaxt="n", yaxt = "n", cex.lab=1.5, cex.axis=1.2)
  axis(side=1, las=1)
  axis(side=2, las=1)
  abline(h=0, lty="dotted")
  abline(v=0, lty="dotted")
  mtext(side = 1, text = paste("Axis 1 (",explainvar1, "%)", sep=""), line = 3, cex=1.5)
  mtext(side = 2, text = paste("Axis 2 (",explainvar2, "%)", sep=""), line = 3.5, cex=1.5)
  box(lwd=2)
  mol.shape <- rep(NA, dim(RDA)[1])
    for (i in 1:length(mol.shape)){
      if (RDA$molecule[i] == "DNA"){mol.shape[i] = 21}
      else {mol.shape[i] = 22}
    }
  slope.color <- rep(NA, dim(RDA)[1])
    for (i in 1:length(slope.color)){
      if (RDA$slope[i] == "North") {slope.color[i] = "brown"}
      else {slope.color[i] = "green3"}
    }
  points(RDA$RDA1, RDA$RDA2, pch=mol.shape, cex=2.0, col="black", bg=slope.color, lwd=2)
  ordiellipse(cbind(RDA$RDA1, RDA$RDA2), RDA$labs, kind="sd", conf=0.95,
    lwd=2, lty=3, draw = "lines", col = "black", label=FALSE)

dev.off()


# Plot 3 - Condition 2 #########################################################
png(file="./plots/Condition2_ESA.png",
    width=600, height=600, antialias = "cleartype", res=96)

  # Inputs
  shared     = "./microbide/SbyS/Condition2.txt"
  cutoff     = "0.03"
  design     = "./data/simmy.design.txt"

  # Import Site by OTU Matrix
  ec_data.sim <- read.otu(shared, "0.03")
  design <- read.delim(design, header=T, row.names=1)

  # Remove OTUs with less than ten observations
  ec_data_red <- ec_data.sim[,colSums(ec_data.sim) >= 10]

  # Create factors for model
  slope <- design$slope # factor 1
  molecule <- design$molecule # factor 2
  paired <- design$paired
  slope.molecule <- data.frame(cbind(as.character(slope),
    as.character(molecule))) # Y matrix with factor 1 and 2
  slope.molecule.concat <- do.call(paste, c(slope.molecule[c("X1", "X2")],
    sep = "")) # create unique treat ID vector

  # Calculate Presense Absence
  dataPA <- (ec_data_red > 0)*1

  # Calculating Relative Abundance
  dataREL <- ec_data_red
  for(i in 1:nrow(ec_data_red)){
    dataREL[i,] = ec_data_red[i,]/sum(ec_data_red[i,])
    }

  # Log Transform Relative Abundance
  dataREL.l <- decostand(dataREL,method="log")

  # Distance Based Redundancy Analysis
  dbRDA <- capscale(dataREL ~ slope + molecule+ Condition(paired), distance="bray")

  # Percent Variance Explained Using PCoA (Axis 1,2,3)
  explainvar1 <- round(dbRDA$CCA$eig[1]/sum(dbRDA$CCA$eig, dbRDA$CA$eig)*100,2)
  explainvar2 <- round(dbRDA$CCA$eig[2]/sum(dbRDA$CCA$eig, dbRDA$CA$eig)*100,2)

  RDA <- as.data.frame(dbRDA$CCA$wa.eig)
  RDA$molecule <- design$molecule
  RDA$slope <- design$slope
  RDA$labs <- slope.molecule.concat


  # Plot Parameters
  par(mfrow=c(1,1), mar=c(5,6,1,1))

  # Initiate Plot
  plot(RDA$RDA1, RDA$RDA2,
    xlim = c(-0.1, 0.1),ylim= c(-0.1, 0.1), pch=16, cex=2.0,
    type="n",ann=FALSE, xaxt="n", yaxt = "n", cex.lab=1.5, cex.axis=1.2)
  axis(side=1, las=1)
  axis(side=2, las=1)
  abline(h=0, lty="dotted")
  abline(v=0, lty="dotted")
  mtext(side = 1, text = paste("Axis 1 (",explainvar1, "%)", sep=""), line = 3, cex=1.5)
  mtext(side = 2, text = paste("Axis 2 (",explainvar2, "%)", sep=""), line = 3.5, cex=1.5)
  box(lwd=2)
  mol.shape <- rep(NA, dim(RDA)[1])
    for (i in 1:length(mol.shape)){
      if (RDA$molecule[i] == "DNA"){mol.shape[i] = 21}
      else {mol.shape[i] = 22}
    }
  slope.color <- rep(NA, dim(RDA)[1])
    for (i in 1:length(slope.color)){
      if (RDA$slope[i] == "North") {slope.color[i] = "brown"}
      else {slope.color[i] = "green3"}
    }
  points(RDA$RDA1, RDA$RDA2, pch=mol.shape, cex=2.0, col="black", bg=slope.color, lwd=2)
  ordiellipse(cbind(RDA$RDA1, RDA$RDA2), RDA$labs, kind="sd", conf=0.95,
    lwd=2, lty=3, draw = "lines", col = "black", label=FALSE)

dev.off()



# Plot 4 - Condition 3 #########################################################
png(file="./plots/Condition3_ESA.png",
    width=600, height=600, antialias = "cleartype", res=96)

  # Inputs
  shared     = "./microbide/SbyS/Condition3.txt"
  cutoff     = "0.03"
  design     = "./data/simmy.design.txt"

  # Import Site by OTU Matrix
  ec_data.sim <- read.otu(shared, "0.03")
  design <- read.delim(design, header=T, row.names=1)

  # Remove OTUs with less than ten observations
  ec_data_red <- ec_data.sim[,colSums(ec_data.sim) >= 10]

  # Create factors for model
  slope <- design$slope # factor 1
  molecule <- design$molecule # factor 2
  paired <- design$paired
  slope.molecule <- data.frame(cbind(as.character(slope),
    as.character(molecule))) # Y matrix with factor 1 and 2
  slope.molecule.concat <- do.call(paste, c(slope.molecule[c("X1", "X2")],
    sep = "")) # create unique treat ID vector

  # Calculate Presense Absence
  dataPA <- (ec_data_red > 0)*1

  # Calculating Relative Abundance
  dataREL <- ec_data_red
  for(i in 1:nrow(ec_data_red)){
    dataREL[i,] = ec_data_red[i,]/sum(ec_data_red[i,])
    }

  # Log Transform Relative Abundance
  dataREL.l <- decostand(dataREL,method="log")

  # Distance Based Redundancy Analysis
  dbRDA <- capscale(dataREL ~ slope + molecule+ Condition(paired), distance="bray")

  # Percent Variance Explained Using PCoA (Axis 1,2,3)
  explainvar1 <- round(dbRDA$CCA$eig[1]/sum(dbRDA$CCA$eig, dbRDA$CA$eig)*100,2)
  explainvar2 <- round(dbRDA$CCA$eig[2]/sum(dbRDA$CCA$eig, dbRDA$CA$eig)*100,2)

  RDA <- as.data.frame(dbRDA$CCA$wa.eig)
  RDA$molecule <- design$molecule
  RDA$slope <- design$slope
  RDA$labs <- slope.molecule.concat


  # Plot Parameters
  par(mfrow=c(1,1), mar=c(5,6,1,1))

  # Initiate Plot
  plot(RDA$RDA2, RDA$RDA1,
    xlim = c(-0.1, 0.1),ylim= c(-0.1, 0.1), pch=16, cex=2.0,
    type="n",ann=FALSE, xaxt="n", yaxt = "n", cex.lab=1.5, cex.axis=1.2)
  axis(side=1, las=1)
  axis(side=2, las=1)
  abline(h=0, lty="dotted")
  abline(v=0, lty="dotted")
  mtext(side = 1, text = paste("Axis 1 (",explainvar2, "%)", sep=""), line = 3, cex=1.5)
  mtext(side = 2, text = paste("Axis 2 (",explainvar1, "%)", sep=""), line = 3.5, cex=1.5)
  box(lwd=2)
  mol.shape <- rep(NA, dim(RDA)[1])
    for (i in 1:length(mol.shape)){
      if (RDA$molecule[i] == "DNA"){mol.shape[i] = 21}
      else {mol.shape[i] = 22}
    }
  slope.color <- rep(NA, dim(RDA)[1])
    for (i in 1:length(slope.color)){
      if (RDA$slope[i] == "North") {slope.color[i] = "brown"}
      else {slope.color[i] = "green3"}
    }
  points(RDA$RDA2, RDA$RDA1, pch=mol.shape, cex=2.0, col="black", bg=slope.color, lwd=2)
  ordiellipse(cbind(RDA$RDA2, RDA$RDA1), RDA$labs, kind="sd", conf=0.95,
    lwd=2, lty=3, draw = "lines", col = "black", label=FALSE)

dev.off()







points(RDA$RDA1, RDA$RDA2, pch=mol.shape, cex=2.0, col="black", bg=slope.color, lwd=2)
ordiellipse(cbind(RDA$RDA1, RDA$RDA2), RDA$labs, kind="sd", conf=0.95,
  lwd=2, lty=3, draw = "lines", col = "black", label=FALSE)


par(mar=c(0, 3, 0, 0))
plot.new()
legend("center", c(paste("All; ",levels(RDA$slope)[1]," Slope", sep=""),
  paste("All; ",levels(RDA$slope)[2]," Slope", sep=""),
  paste("Active; ",levels(RDA$slope)[1]," Slope", sep=""),
  paste("Active; ",levels(RDA$slope)[2]," Slope", sep="")),
  pt.lwd=2, col="black", pt.bg=c("brown", "green3", "brown",
  "green3"), pch=c(21,21,22,22), bty='n', ncol=2, cex=1.5, pt.cex=2)





  shared     = "./microbide/SbyS/Condition6.txt"
  cutoff     = "0.03"
  design     = "./data/simmy.design.txt"



  source("bin/DiversityFunctions.r")
  require(vegan)

  # Import Site by OTU Matrix
  ec_data.sim <- read.otu(shared, "0.03")
  design <- read.delim(design, header=T, row.names=1)

  # Remove OTUs with less than ten observations
  ec_data_red <- ec_data.sim[,colSums(ec_data.sim) >= 10]

  # Create factors for model
  slope <- design$slope # factor 1
  molecule <- design$molecule # factor 2
  paired <- design$paired
  site <- design$site
  station <- design$station
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
  samplePA.dist <- vegdist(dataREL.h,method="bray")
  sampleREL.dist <- vegdist(dataREL.h,method="bray")

  # Distance Based Redundancy Analysis
  dbRDA <- capscale(dataREL.l ~ slope + molecule, distance="bray")
# Distance Based Redundancy Analysis
dbRDA <- capscale(dataREL ~ slope + molecule+ Condition(paired), distance="bray")

head(summary(dbRDA))
RsquareAdj(dbRDA)
anova(dbRDA, by="terms", permu=999)

varpart(dataREL,  ~ slope, ~ molecule)


  # Percent Variance Explained Using PCoA (Axis 1,2,3)
  explainvar1 <- round(dbRDA$CCA$eig[1]/sum(dbRDA$CCA$eig, dbRDA$CA$eig)*100,2)
  explainvar2 <- round(dbRDA$CCA$eig[2]/sum(dbRDA$CCA$eig, dbRDA$CA$eig)*100,2)



RDA <- as.data.frame(dbRDA$CCA$wa.eig)
RDA$molecule <- design$molecule
RDA$slope <- design$slope
RDA$labs <- slope.molecule.concat


# Plot Parameters
par(mfrow=c(1,1), mar=c(5,5,1,1))
layout(rbind(1, 2), height=c(7, 1))
x.dim <- c(min(RDA$RDA1)+min(RDA$RDA1)*0.2,max(RDA$RDA1)+max(RDA$RDA1)*0.2)
y.dim <- c(min(RDA$RDA2)+min(RDA$RDA2)*0.2,max(RDA$RDA2)+max(RDA$RDA2)*0.2)

# Initiate Plot
plot(RDA$RDA1, RDA$RDA2,
  xlab = paste("RDA Axis 1 (",explainvar1, "%)", sep=""),
  ylab = paste("RDA Axis 2 (",explainvar2, "%)", sep=""),
  xlim = c(-0.2, 0.2),ylim= c(-0.1, 0.1), pch=16, cex=2.0, type="n",xaxt="n",
  yaxt = "n", cex.lab=1.5, cex.axis=1.2)
axis(side=1, las=1)
axis(side=2, las=1)
abline(h=0, lty="dotted")
abline(v=0, lty="dotted")
mol.shape <- rep(NA, dim(RDA)[1])
  for (i in 1:length(mol.shape)){
    if (RDA$molecule[i] == "DNA"){mol.shape[i] = 21}
    else {mol.shape[i] = 22}
  }
slope.color <- rep(NA, dim(RDA)[1])
  for (i in 1:length(slope.color)){
    if (RDA$slope[i] == "North") {slope.color[i] = "brown"}
    else {slope.color[i] = "green3"}
  }
points(RDA$RDA1, RDA$RDA2, pch=mol.shape, cex=2.0, col="black", bg=slope.color, lwd=2)
ordiellipse(cbind(RDA$RDA1, RDA$RDA2), RDA$labs, kind="sd", conf=0.95,
  lwd=2, lty=3, draw = "lines", col = "black", label=FALSE)

box(lwd=2)
par(mar=c(0, 3, 0, 0))
plot.new()
legend("center", c(paste("All; ",levels(RDA$slope)[1]," Slope", sep=""),
  paste("All; ",levels(RDA$slope)[2]," Slope", sep=""),
  paste("Active; ",levels(RDA$slope)[1]," Slope", sep=""),
  paste("Active; ",levels(RDA$slope)[2]," Slope", sep="")),
  pt.lwd=2, col="black", pt.bg=c("brown", "green3", "brown",
  "green3"), pch=c(21,21,22,22), bty='n', ncol=2, cex=1.5, pt.cex=2)


anova(dbRDA, by="terms", permu=999)
adonis(dataREL.c ~ slope + molecule, strata=paired, method="bray", permutations=1000)

varpart(dataREL,  ~ slope, ~ molecule)