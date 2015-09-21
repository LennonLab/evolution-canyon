################################################################################
#                                                                              #
#	Evolution Canyon Project: Microbial Community Characterization               #
#                                                                              #
################################################################################
#                                                                              #
#	Written by: Jay Lennon (2014/02/03)                                          #
# Modified by: M. Muscarella                                                   #
#                                                                              #
#	Last update: 2014/03/11                                                      #
#                                                                              #
################################################################################


# Setup Work Environment
rm(list=ls())
setwd("~/GitHub/evolution-canyon/community")
se <- function(x){sd(x)/sqrt(length(x))}

# Code Dependencies
source("functions/DiversityFunctions.r")  
require(vegan)                                   
require(reshape)
require("BiodiversityR")
require("ecodist")

## Data Import and Manipulations ###############################################
################################################################################

# Site by OTU Matrix
EC_data <- read.otu("data/EC.bad.shared", "0.03") 
  # I manually edited to remove the unique and millions of blank otus
EC_data <- t(EC_data)
EC_data_red <- EC_data[,c(2:21,24:26,28,30:81)] 
  # select this way instead of with cbind - fixes loss of column name issue
	# EC_data_red removes the following samples that had independently 
  # documented problems with nucleic acid extraction and/or amplification

# Taxonomy Data
EC_tax_raw <- (read.delim("data/EC.bad.taxonomy")) # load genus-level data
EC_tax <- transform(EC_tax_raw,Taxonomy=colsplit(EC_tax_raw[,3],
  split="\\;",names=c("Domain","Phylum","Class","Order","Family","Genus")))

# Merge Taxonomy and Recovery; not really needed
mix_data <- cbind(EC_tax,EC_data_red)

# Calculate Presense Absence
dataPA <- (EC_data_red > 0)*1 

# Calculating Relative Abundance
dataREL <- EC_data_red
 for(i in 1:ncol(EC_data_red)){
  dataREL[,i] = EC_data_red[,i]/sum(EC_data_red[,i])
  }  

## Coverage and Diversity Measures #############################################
################################################################################

# These are only here because they identify the issues with the current data set

## Measure of Sample Coverage
#EC_coverage <- count.groups(read.otu("data/EC.shared", "0.03"))
#EC_coverage
#
## Taxon Richness with Resampling
#EC_Richness1 <- richness.iter("data/EC.shared", "0.03", 10000, 100)
#EC_Richness2 <- richness.iter("data/EC.shared", "0.03", 1000, 100)
#EC_Richness3 <- richness.iter("data/EC.shared", "0.03", 400, 100)
#EC_Richness4 <- richness.iter("data/EC.shared", "0.03", 200, 100)
#
## Taxon Diversity with Resampling
#EC_Shannon1 <- diversity.iter("data/EC.shared", "shannon", "0.03", 10000, 100)
#EC_Shannon2 <- diversity.iter("data/EC.shared", "shannon", "0.03", 1000, 100)
#EC_Shannon3 <- diversity.iter("data/EC.shared", "shannon", "0.03", 400, 100)
#EC_Shannon4 <- diversity.iter("data/EC.shared", "shannon", "0.03", 200, 100)
#
## Taxon Evenness with Resampling
#EC_Evenness1 <- evenness.iter("data/EC.shared", "0.03", 10000, 100)
#EC_Evenness2 <- evenness.iter("data/EC.shared", "0.03", 1000, 100)
#EC_Evenness3 <- evenness.iter("data/EC.shared", "0.03", 400, 100)
#EC_Evenness4 <- evenness.iter("data/EC.shared", "0.03", 200, 100)

## PCoA Analysis ###############################################################
################################################################################

# Create Distance Matrix
samplePA.dist <- vegdist(t(dataPA),method="bray")
sampleREL.dist <- vegdist(t(dataREL),method="bray")

# Principal Coordinates Analysis
EC_pcoa <- cmdscale(sampleREL.dist,k=3,eig=TRUE,add=FALSE) 
  # Classical (Metric) Multidimensional Scaling; returns PCoA coordinates
  # eig=TRUE returns eigenvalues; k = # of dimensions to calculate

# Responder Analysis Based on PCoA 
pcoaS <- add.spec.scores(EC_pcoa,t(dataREL),method="cor.scores",Rscale=TRUE,
  scaling=1,multi=1) 
  # retrieves correlation coefficient for each taxon's relative 
  # abudnace with respect to PCoA coordinates (k = 3)
  
# PCoA Axis 1 Responders
cor_spp_a1 <- cbind(EC_tax[,3],pcoaS$cproj[,1]) 
  # creates matrix of taxonomy and correlations (r)
colnames(cor_spp_a1)[7] <- "Corr" 
  # renames correlation column name
cor_spp_a1 <- cor_spp_a1[order(cor_spp_a1[,7],decreasing=TRUE),] 
  # sorts based on r
responders_a1 <- cor_spp_a1[a1 <- abs(cor_spp_a1[,7])>0.7,] 
  # subset based on r > |0.7| 

# PCoA Axis 2 Responders
cor_spp_a2 <- cbind(EC_tax[,3],pcoaS$cproj[,2]) 
  # creates matrix of taxonomy and correlations(r) axis 2
colnames(cor_spp_a2)[7] <- "Corr" 
  # renames correlation column name
cor_spp_a2 < -cor_spp_a2[order(cor_spp_a1[,7],decreasing=TRUE),] 
  # sorts based on r
responders_a2 <- cor_spp_a2[a2 <- abs(cor_spp_a2[,7])>0.7,] 
  # subset based on r > |0.7| 

# Percent Variance Explained Using PCoA 
Eig_sum <- sum(pcoaS$eig) 
  # sum of k eigenvalues
eig1 <- pcoaS$eig[1] 
  # eigenvalue for PCoA 1
explainvar1 <- pcoaS$eig[1]/Eig_sum 
  # eigenvalue for PCoA 1 divided by sum of k eigenvalues
explainvar1 
  # percent variance explained by PCoA 1
eig2 <- pcoaS$eig[2] 
  # eigenvalue for PCoA 2
explainvar2 <- pcoaS$eig[2]/Eig_sum 
  # eigenvalue for PCoA 2 divided by sum of k eigenvalues
explainvar2 
  # percent variance explained by PCoA 2
eig3 <- pcoaS$eig[3] 
  # eigenvalue for PCoA 3
explainvar3 <- pcoaS$eig[3]/Eig_sum 
  # eigenvalue for PCoA 3 divided by sum of k eigenvalues
explainvar3 
  # percent variance explained by PCoA 3

## PCoA Plots ##################################################################
################################################################################

# Organize Samples for Plotting
# Plot Codes:
# "St" = station; 1 and 2 = African (dry) slope; 5 and 6 = European (wet) slope 
# "Loc" = locations; 10 replicate plots on each station, identified (1, 2, 5, 6)

# Stripping PCoA Loadings out of cmdscale Output
pcoap <- as.data.frame(pcoaS$points)
site <- rep(NA,76)
site[c(1,3,5,7,9,11,13,15,17,19)] <- "AF1_DNA"
site[c(21,23,24,25,27,29,31,33,35)] <- "AF2_DNA"
site[c(37,39,41,43,45,47,49,51,53,55)] <- "EU1_DNA"
site[c(57,59,61,63,65,67,69,71,73,75)] <- "EU2_DNA"
site[c(2,4,6,8,10,12,14,16,18,20)] <- "AF1_RNA"
site[c(22,26,28,30,32,34,36)] <- "AF2_RNA"
site[c(38,40,42,44,46,48,50,52,54,56)] <- "EU1_RNA"
site[c(58,60,62,64,66,68,70,72,74,76)] <- "EU2_RNA"
pcoap <- cbind(pcoap, site)

#--> AF: Station 1 - DNA
AF1_DNA_1 <- pcoap[c(1,3,5,7,9,11,13,15,17,19),1]             #PCoA 1 scores
AF1_DNA_2 <- pcoap[c(1,3,5,7,9,11,13,15,17,19),2]             #PCoA 2 scores 
#--> AF: Station 2 - DNA
AF2_DNA_1 <- pcoap[c(21,23,24,25,27,29,31,33,35),1]           #PCoA 1 scores  
AF2_DNA_2 <- pcoap[c(21,23,24,25,27,29,31,33,35),2]           #PCoA 2 scores 
#--> EU: Station 5 - DNA
EU1_DNA_1 <- pcoap[c(37,39,41,43,45,47,49,51,53,55),1]        #PCoA 1 scores
EU1_DNA_2 <- pcoap[c(37,39,41,43,45,47,49,51,53,55),2]        #PCoA 2 scores 
#--> EU: Station 6 - DNA
EU2_DNA_1 <- pcoap[c(57,59,61,63,65,67,69,71,73,75),1]        #PCoA 1 scores
EU2_DNA_2 <- pcoap[c(57,59,61,63,65,67,69,71,73,75),2]        #PCoA 2 scores 
#--> AF: Station 1 - RNA
AF1_RNA_1 <- pcoap[c(2,4,6,8,10,12,14,16,18,20),1]            #PCoA 1 scores
AF1_RNA_2 <- pcoap[c(2,4,6,8,10,12,14,16,18,20),2]            #PCoA 2 scores 
#--> AF: Station 2 - RNA
AF2_RNA_1 <- pcoap[c(22,26,28,30,32,34,36),1]                 #PCoA 1 scores  
AF2_RNA_2 <- pcoap[c(22,26,28,30,32,34,36),2]                 #PCoA 2 scores  
#--> EU: Station 5 - RNA
EU1_RNA_1 <- pcoap[c(38,40,42,44,46,48,50,52,54,56),1]        #PCoA 1 scores 
EU1_RNA_2 <- pcoap[c(38,40,42,44,46,48,50,52,54,56),2]        #PCoA 2 scores   
#--> EU: Station 6 - RNA
EU2_RNA_1 <- pcoap[c(58,60,62,64,66,68,70,72,74,76),1]        #PCoA 1 scores 
EU2_RNA_2 <- pcoap[c(58,60,62,64,66,68,70,72,74,76),2]        #PCoA 2 scores  

# Plot Parameters
par(mfrow=c(1,1), mar=c(5,5,1,1)) 
x <- c(min(pcoap[,2])+min(pcoap[,2])*0.25,max(pcoap[,2])+max(pcoap[,2])*0.1)
y <- c(min(pcoap[,1])+min(pcoap[,1])*0.1,max(pcoap[,1])+max(pcoap[,1])*0.1)

# Initiate Plot
plot(pcoap[,2], pcoap[,1], xlab="PCoA Axis 2", ylab="PCoA Axis 1", 
  xlim=rev(range(x)),ylim= y, pch=16, cex=2.0, type='n',xaxt="n",yaxt="n", 
  cex.lab=1.5, cex.axis=1.2) 
    # Creates PCoA byplot for PCoA axis 1 and PCoA axis 2; adds axes labels, 
    # Adjusts axes lengths and size
    
  axis(side=1,at=c(-0.4,-0.2,0,0.2), las=1)    # adds x-axis ticks
  axis(side=2,at=c(-0.2,0,0.2,0.4), las=1)     # adds y-axis ticks
  segments(-1500, -0, 1500, 0, lty="dotted")   # adds horizontal reference line
  segments(0, -1500, 0, 1500, lty="dotted")    # adds vertical reference line
  
  points(AF1_DNA_2,AF1_DNA_1,pch=21,cex=2.0,col="black",bg="brown3",lwd=2)
  points(AF2_DNA_2,AF2_DNA_1,pch=21,cex=2.0,col="black",bg="brown3",lwd=2)
  points(EU1_DNA_2,EU1_DNA_1,pch=21,cex=2.0,col="black",bg="green3",lwd=2)
  points(EU2_DNA_2,EU2_DNA_1,pch=21,cex=2.0,col="black",bg="green3",lwd=2)
  points(AF1_RNA_2,AF1_RNA_1,pch=22,cex=2.0,col="black",bg="brown3",lwd=2)
  points(AF2_RNA_2,AF2_RNA_1,pch=22,cex=2.0,col="black",bg="brown3",lwd=2)
  points(EU1_RNA_2,EU1_RNA_1,pch=22,cex=2.0,col="black",bg="green3",lwd=2)
  points(EU2_RNA_2,EU2_RNA_1,pch=22,cex=2.0,col="black",bg="green3",lwd=2)
    # adds PCoA scores for each to existing graph
  
  ordiellipse(cbind(pcoap[,2], pcoap[,1]), pcoap$site, kind="sd", conf=0.95,
    lwd=2, lty=3, draw = "lines", col = "black", label=TRUE)  
  #legend(-0.25, -0.08,"AF-1-DNA", cex=2.0,col="red",pch=16,bty="n")  
  #legend(-0.25, -0.12,"AF-2-DNA", cex=2.0,col="pink",pch=16,bty="n")
  #legend(-0.25, -0.12,"EU-1-DNA", cex=1.5,col="darkgreen",pch=16,bty="n")
  #legend(-0.25, -0.12,"EU-2-DNA", cex=1,col="chartreuse4",pch=16,bty="n")
  #legend(-0.25, -0.12,"AF-1-RNA", cex=1,col="red",pch=16,bty="n")
  #legend(-0.25, -0.12,"AF-1-RNA", cex=1,col="pink",pch=16,bty="n")
  #legend(-0.25, -0.12,"AF-1-RNA", cex=1,col="red",pch=16,bty="n")
  #legend(-0.25, -0.12,"AF-1-RNA", cex=1,col="pink",pch=16,bty="n")
    # adds legends to existing graph
  box(lwd=2)

## Hypothesis Testing ##########################################################
################################################################################

site <- rep(NA,76)
site[c(1,3,5,7,9,11,13,15,17,19)] <- "AF1_DNA"
site[c(21,23,24,25,27,29,31,33,35)] <- "AF2_DNA"
site[c(37,39,41,43,45,47,49,51,53,55)] <- "EU1_DNA"
site[c(57,59,61,63,65,67,69,71,73,75)] <- "EU2_DNA"
site[c(2,4,6,8,10,12,14,16,18,20)] <- "AF1_RNA"
site[c(22,26,28,30,32,34,36)] <- "AF2_RNA"
site[c(38,40,42,44,46,48,50,52,54,56)] <- "EU1_RNA"
site[c(58,60,62,64,66,68,70,72,74,76)] <- "EU2_RNA"

station <- rep(NA, 76)
station[1:20] <- "AF1"
station[21:36] <- "AF2"
station[37:56] <- "EU1"
station[57:76] <- "EU2"

molecule <- rep(NA, 76)
molecule[c(1,3,5,7,9,11,13,15,17,19,21,23,24,25,27,29,31,33,35,
  37,39,41,43,45,47,49,51,53,55,57,59,61,63,65,67,69,71,73,75)] <- "DNA"
molecule[c(2,4,6,8,10,12,14,16,18,20,22,26,28,30,32,34,36,
  38,40,42,44,46,48,50,52,54,56,58,60,62,64,66,68,70,72,74,76)] <- "RNA"
  
slope <- rep(NA, 76)
slope[1:36] <- "AF"
slope[37:76] <- "EU"

# Adonis (PERMANOVA)
  # adonis runs a PERMANOVA (Created by Marti J. Anderson) this is very similar 
  # to ANOVA but for multivariate data. You can make very complex
  # experimental designs with it.
  # The default distance measure is bray-curtis, but other measures 
  # (Chao, Jaccard, Euclidean) can be used when specified  
Adonis <- adonis(sampleREL.dist ~ molecule*slope, method="bray", 
  permutations=1000)
Adonis