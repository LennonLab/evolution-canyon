################################################################################
#                                                                              #
# Evolution Canyon Project: Microbial Community RDA                            #
#                                                                              #
################################################################################
#                                                                              #
# Written by: Jay Lennon and Mario Muscarella                                  #
#                                                                              #
# Last update: 2014/08/04                                                      #
#                                                                              #
################################################################################
#                                                                              #
# Notes: This code performs the RDA function on the Evolution Canyon           #
#        emperical data                                                        #
#                                                                              #
# Recent Changes:                                                              #
#         1. Calculate Variance Explained                                      #
#         2. Confirm Normalizations                                            #
#         3. Removed Singletons                                                #
#         4. Clean up plots                                                    #
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

# Create Distance Matrix with bray (deafault), manhattan, euclidean, canberra, 
#  bray, kulczynski, jaccard, gower, altGower, morisita, horn, mountford, raup, 
#  binomial, or chao. Most should be part of vegan, but possilbly 'labdsv' or 
#  'BiodiversityR' packages
samplePA.dist <- vegdist(dataREL.h,method="bray")
sampleREL.dist <- vegdist(dataREL.h,method="bray")

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
RDA$molecule <- design_red$molecule
RDA$slope <- design_red$slope
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
  xlim = x.dim,ylim= y.dim, pch=16, cex=2.0, type="n",xaxt="n",
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
    if (RDA$slope[i] == "AF") {slope.color[i] = "brown"}
    else {slope.color[i] = "green3"}
  } 
points(RDA$RDA1, RDA$RDA2, pch=mol.shape, cex=2.0, col="black", bg=slope.color, lwd=2)   
ordiellipse(cbind(RDA$RDA1, RDA$RDA2), RDA$labs, kind="sd", conf=0.95,
  lwd=2, lty=3, draw = "lines", col = "black", label=TRUE)

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



# Overlay Plot
# Plot Parameters
par(mar=c(5,5,1,1)) 
 
x.dim <- c(min(RDA$RDA1)+min(RDA$RDA1)*0.2,max(RDA$RDA1)+max(RDA$RDA1)*0.2)
y.dim <- c(min(RDA$RDA2)+min(RDA$RDA2)*0.2,max(RDA$RDA2)+max(RDA$RDA2)*0.2)
                       
# Initiate Plot
plot(RDA$RDA1, RDA$RDA2, 
  xlab = "",
  ylab = "", 
  xlim = x.dim,ylim= y.dim, pch=16, cex=2.0, type="n",xaxt="n",
  yaxt = "n", cex.lab=1.5, cex.axis=1.2)  

mol.shape <- rep(NA, dim(RDA)[1])
  for (i in 1:length(mol.shape)){
    if (RDA$molecule[i] == "DNA"){mol.shape[i] = 21}
    else {mol.shape[i] = 22}
  }
slope.color <- rep(NA, dim(RDA)[1])
  for (i in 1:length(slope.color)){
    if (RDA$slope[i] == "AF") {slope.color[i] = "brown"}
    else {slope.color[i] = "green3"}
  } 
points(RDA$RDA1, RDA$RDA2, pch=mol.shape, cex=2.0, col="black", bg=slope.color, lwd=2)   
ordiellipse(cbind(RDA$RDA1, RDA$RDA2), RDA$labs, kind="sd", conf=0.95,
  lwd=2, lty=3, draw = "lines", col = "black", label=FALSE)

box(lwd=2)
