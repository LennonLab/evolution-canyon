################################################################################
#                                                                              #
# Evolution Canyon Project: Microbial Community PL-SDA                         #
#                                                                              #
################################################################################
#                                                                              #
# Written by: Jay Lennon and Mario Muscarella                                  #
#                                                                              #
# Last update: 2014/07/25                                                      #
#                                                                              #
################################################################################
#                                                                              #
# Notes: This code performs the PLS-DA function on the Evolution Canyon        #
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
require("vegan")||install.packages("vegan");require("vegan")
require("mixOmics")||install.packages("mixOmics");require("mixOmics")

shared <- "./mothur/EC.bac.final.shared" 
design <- "./data/design.txt"
 
#1 -- Import data (site by OTU matrix) and deal w/ "problem" samples
ec_data <- t(read.otu(shared, "0.03")) # takes a long time (10 mins)
design <- read.delim(design, header=T, row.names=1) 
 
  # Note: owing to amplification issues, we only submitted 76 of 80 samples
  # The four samples NOT submitted were C-1E-R, EC-2G-R, EC-2J-R, EC-6I-D
  # So, removed their pairs for PLS-DA: EC-1E-D, EC-2G-D, EC-2J-D,EC-6I-R
  # 4 samples had crappy quality in past: EC-2A-D, EC-2A-R, EC-2C-R, EC-2D-R
  # Remove pairs EC-2C-D, EC-2D-R, for a total of 6 "questionable" samples
  # Results in a total of 66 samples (14 problematics)
  # Use following code to identify columns to remove
  
# Remove problamatic samples and pairs
samps <- (colnames(ec_data))
which(samps == "EC-2A-D")
ec_data_red <- ec_data[,-c(9,20:21,24:27,32,37,74)]

# To keeps in the a priori "questionable" samples that were run:
# ec_data_red <- ec_data[,-c(9,32,37,74)]
  
# Remove any OTUs with <2 OTUs
ec_data_red <- ec_data_red[which(rowSums(ec_data_red) >= 2),]

# Remove problamatic samples and pairs from design
samps_red <- colnames(ec_data_red)
design_red <- design[samps_red,] 

# Design matrix with renumbered pairs after removing problem pairs
design_red$paired <- c(rep(seq(1:14), each=2), rep(seq(1:19), each=2))
# To keep questionable samples, uncomment below
# design_red$paired <- c(rep(seq(1:17), each=2), rep(seq(1:19), each=2))

  
#2 -- Options for transformation, relativation, and normalization
# "Transformation 1": transpose and log10 transform, also remove zero taxa
Xt <- t(ec_data_red)[,which(colSums(ec_data_red) !=0)] 
Xlogt <- decostand(Xt,method="log")[,which(colSums(Xt) !=0)]
  
# "Transformation 2": relativize, transpose, log10-transform
Xrel<-ec_data_red
  for(i in 1:ncol(ec_data_red)){
    Xrel[,i]=ec_data_red[,i]/sum(ec_data_red[,i])
  } 
Xrelt<-t(Xrel)
Xlogrelt <- decostand(Xrelt,method="log")[,which(colSums(Xrelt) !=0)]
  
# "Transformation 3": relativize, transpose, Hellinger-transform 
Xrel<-ec_data_red
  for(i in 1:ncol(ec_data_red)){
    Xrel[,i]=ec_data_red[,i]/sum(ec_data_red[,i])
  } 
Xrelt<-t(Xrel)
Xhellit <- decostand(Xrelt,method="hellinger")[,which(colSums(Xrelt) !=0)]

#3 -- Create factors for multilevel model
slope    <- design_red$slope    # PLS-DA factor 1
molecule <- design_red$molecule # PLS-DA factor 2 
site     <- design_red$site
station  <- design_red$station

slope.molecule <- data.frame(cbind(as.character(slope),
    as.character(molecule))) # Y matrix with factor 1 and 2
slope.molecule.concat <- do.call(paste, c(slope.molecule[c("X1", "X2")],
    sep = "")) # create unique treat ID vector
pair.station <- c(rep(seq(1:9), each=2), rep(seq(1:5), each=2),
    rep(seq(1:10), each=2), rep(seq(1:9), each=2))
    
  # To include "questionable" data uncomment the following
  # pair.station <- c(rep(seq(1:9), each=2), rep(seq(1:8), each=2),
  #    rep(seq(1:10), each=10), rep(seq(1:9), each=2))

station.molecule.concat <- paste(station, molecule, sep = "")

#4 -- Run multilevel models
one.way        <- multilevel(Xt, cond = slope.molecule.concat,
                    sample = pair.station, ncomp = 2, method = 'splsda')
two.way.log    <- multilevel(Xlogt, cond = slope.molecule,
                    sample = pair.station, ncomp = 2, method = 'splsda')
two.way.logrel <- multilevel(Xlogrelt, cond = slope.molecule,
                    sample = pair.station, ncomp = 2, method = 'splsda')
two.way.hel    <- multilevel(Xhellit, cond = slope.molecule,
                    sample = pair.station, ncomp = 2, method = 'splsda')
  
# Indicate which model you want to graph
EC_multilevel <- two.way.logrel

#5 -- Variation Explained by axis
# Calculate distance between samples (Bray Curtis or Euclidean?)
X.dist  <- vegdist(Xlogrelt,method="euclidean")

# Calculate distance between samples in reduced (ordination) space
plsda.1 <- dist(EC_multilevel$variates$X[,1],method="euclidean")
plsda.2 <- dist(EC_multilevel$variates$X[,2],method="euclidean")
# plsda.3 <- dist(EC_multilevel$variates$X[,3],method="euclidean")

# Calculate variation explained
var1 <- round(cor(X.dist, plsda.1), 2)*100
var2 <- round(cor(X.dist, plsda.2), 2)*100
# var3 <- cor(X.dist, plsda.3)

#6 -- Plot PLSDA ordination
# Plot Parameters
par(mar=c(4.5,4.5,0,0), oma=c(1,1,1,1)+0.1 )

points <- EC_multilevel$variates$X
mol.shape <- rep(NA, dim(points)[1])
  for (i in 1:length(mol.shape)){
    if (molecule[i] == "DNA"){mol.shape[i] = 21}
    else {mol.shape[i] = 22}
  } # identifies symbol shape based on molecule
slope.color <- rep(NA, dim(points)[1])
  for (i in 1:length(slope.color)){
    if (slope[i] == levels(slope)[1]) {slope.color[i] = "brown"}
    else {slope.color[i] = "green3"}
  } 
  
# Initiate Plot
plot(points[,1], points[,2], 
     xlab=paste("PLS-DA Axis 1 (",var1, "%)", sep=""), 
     ylab=paste("PLS-DA Axis 2 (",var2, "%)", sep=""), 
     xlim=c(min(points[,1])+min(points[,1])*0.1,
            max(points[,1])+max(points[,1])*0.1),
     ylim=c(min(points[,2])+min(points[,2])*0.1,
            max(points[,2])+max(points[,2])*0.1),
     pch=16, cex=2.0, type="n",xaxt="n", yaxt="n", 
     cex.lab=1.5, cex.axis=1.2) 
axis(side=1, las=1, at=c(seq(-0.2,2,by=0.1))) # add x-axis ticks and labels  
axis(side=2, las=1, at=c(seq(-0.2,2,by=0.1))) # add y-axis ticks and labels   
abline(h=0, lty="dotted") # add horizontal dashed line at 0 
abline(v=0, lty="dotted") # add vertical dashed line at 0       
points(points[,1], points[,2], 
       pch=mol.shape, cex=2.0, col="black", bg=slope.color, lwd=2)
ordiellipse(cbind(points[,1], points[,2]), 
            slope.molecule.concat, kind="sd", conf=0.95, lwd=2, lty=0, 
            draw = "lines", col = "black", label=TRUE, pos=3, offset=2.7)

#7 -- Contribution (%) variance of factors (Y) to PLS-DA axes
#  http://perso.math.univ-toulouse.fr/mixomics/faq/numerical-outputs/
# not sure what "Rd" means. U are the "variates"...

# Slope
Y <- slope 
# turn categorical vars into quantitative
Rd.YvsU = cor(as.numeric(as.factor(Y)),EC_multilevel$variates$X)
Rd.YvsU = apply(Rd.YvsU^2, 2, sum)
Rd.Y = cbind(Rd.YvsU, cumsum(Rd.YvsU))
colnames(Rd.Y) = c("Proportion", "Cumulative")
Rd.Y # percent of variance explained by each component

# Molecule
Y <- molecule 
Rd.YvsU = cor(as.numeric(as.factor(Y)),EC_multilevel$variates$X) 
Rd.YvsU = apply(Rd.YvsU^2, 2, sum)
Rd.Y = cbind(Rd.YvsU, cumsum(Rd.YvsU))
colnames(Rd.Y) = c("Proportion", "Cumulative")
Rd.Y # percent of variance explained by each component
    
#8 -- Other stuff: some note and tries based on following website:
# http://perso.math.univ-toulouse.fr/mixomics/methods/spls-da/   
# calculate the coefficients of the linear combinations
two.way.logrel <- multilevel(Xlogrelt, cond = slope.molecule,
  sample = pair.station, ncomp = 2, method = 'splsda')
  
ec_splsda <- spls(Xlogrelt, model.matrix(~ X1 + X2 - 1, slope.molecule), 
                  ncomp = 2)
  
pred <- predict(EC_multilevel, EC_multilevel$X[1:2, ])
pred$B.hat
# calculate R2 and Q2 values for sPLS-DA?
# Q2 = "Q2 is the square of the correlation between the actual 
#       and predicted response"
# warning = this is memory intenstive, but didn't crash
val <- valid(EC_multilevel, criterion = c("R2", "Q2"))
val <- valid(EC_multilevel, validation = "loo")
  

  
### Stuff below is leftover from PCoA. May want to consider 
 
# From P-Meso Project
# Combine with Taxonomy
hetero.tax.raw <- (read.delim("../mothur/output/tanks.bac.final.0.03.taxonomy"))
hetero.tax <- transform(hetero.tax.raw, Taxonomy=colsplit(hetero.tax.raw[,3],
 split="\\;", names=c("Domain","Phylum","Class","Order","Family","Genus")))
rownames(hetero.tax) <- hetero.tax[,1]
rownames(heteros.bio) <- gsub("Otu0", "Otu", rownames(heteros.bio)) 
hetero.data <- merge(heteros.bio, hetero.tax, by = "row.names", all.x=T )
hetero.data <- hetero.data[order(hetero.data$group, -hetero.data$indval), ]
hetero.data <- hetero.data[hetero.data$Taxonomy$Phylum != "unclassified(100)",]

