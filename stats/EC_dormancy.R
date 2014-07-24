################################################################################
#                                                                              #
# Evolution Canyon Project: Microbial Community Dormancy                       #
#                                                                              #
################################################################################
#                                                                              #
# Written by: Mario Muscarella                                                 #
#                                                                              #
# Last update: 2014/07/23                                                      #
#                                                                              #
################################################################################
#                                                                              #
# Notes: This code estimates the %Taxa Dormant                                 #
#                                                                              #
# Issues: None detected yet                                                    #
#                                                                              #
# Recent Changes:                                                              #
#         1. Just Started                                                      #
#                                                                              #
# Future Changes (To-Do List):                                                 #
#         1.                                                                   #
#         2.                                                                   #
#                                                                              #
################################################################################

# Setup Work Environment
rm(list=ls())
setwd("~/GitHub/evolution-canyon")
source("./bin/DiversityFunctions.R")

se <- function(x, ...){sd(x, ...)/sqrt(length(x))}

# Define Inputs and Parameters
design <- "./data/design.txt"
shared <- "./mothur/EC.bac.final.shared"
level  <- "0.03"
size   <- 50000

# Import Design File
design <- read.delim(design, header=T, row.names=1)
slope.molecule <- data.frame(cbind(as.character(design$slope),
  as.character(design$molecule))) # Y matrix with factor 1 and 2
design$slope.molecule <- do.call(paste, c(slope.molecule[c("X1", "X2")],
  sep = "")) # create unique treat ID vector
station.site <- data.frame(cbind(as.character(design$station), as.character(design$site)))
design$station.site <- do.call(paste, c(station.site[c("X1", "X2")], sep = ""))

groups <- data.frame(cbind(rownames(design), design$station.site, as.character(design$slope)))
groups$X1 <- as.character(groups$X1)


system.time({
# Import OTU Table
ec_data <- read.otu(shared, level)

# Calculate Sampling Depth & Coverage
counts <- count.groups(ec_data)
})
# Calculate Presense Absence
dataPA <- (ec_data > 0)*1

# Calculating Relative Abundance
dataREL <- ec_data
  for(i in 1:ncol(ec_data)){
    dataREL[,i] = ec_data[,i]/sum(ec_data[,i])
    }

# Calculate Dormancy

  sites <- levels(groups$X2)
  per.d <- rep(NA, length(sites))
  total <- rep(NA, length(sites))
  active <- rep(NA, length(sites))
  slope <- as.character(design[design$molecule == "DNA",]$slope)
   
  for (i in 1:length(sites)){
  
  site <- groups$X1[groups$X2 == sites[i]]
  
  if (length(site) == 1) {
  
    print(paste(sites[i], "is not paired", sep=" "))   
    } else {
    
      if (length(site) == 2) {
      

  X.tmp1 <- t(subset(dataPA, row.names(dataPA) %in% site)) 
  X.tmp2 <- X.tmp[which(rowSums(X.tmp1) != 0),]  
  tot <- dim(X.tmp2)[1]
  total[i] <- tot
  act <- sum(X.tmp2[,2])
  active[i] <- act
  per.d[i] <- 1 - (act/tot)

  
  } else {
  print(paste(sites[i], "contains more than two samples", sep=" "))
  }}}

out <- data.frame(sites, slope, as.numeric(total), as.numeric(active), per.d = as.numeric(round(per.d, 2)))

plot(out$slope, out$per.d, ylab = "Dormancy")







dormancy <- function(X = "", groupings = ""){

  dataPA <- (X > 0)*1
  sites <- levels(groupings$X2)
   
  
  
  tmp <- groups$X1[groups$X2 == i]
  
  if (length(tmp) == 1) {
  
    print(%i "is not paired")   
    } else {
    
      if (length(tmp) == 2 {
      
testX <- matrix(c(1,0,1,0,1,1), 3, 2)
for (1:dim(testX)[1]){}

sites <- levels(X$site)
per.d <- rep(NA, length(sites))
for i in (1:length(sites)){
  X.tmp <- subset(X, X$site == i)
  X.tmp <- X.tmp[,which(rowSums(X.tmp) != 0)]  
  tot <- dim(X.tmp)[1]
  act <- sum(X.tmp[,1])
  per.d[i] <- 1 - (active/total)
  }
  
  } else {
  print(%i "containts more than two samples")
  }}


