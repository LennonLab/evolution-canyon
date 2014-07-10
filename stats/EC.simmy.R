################################################################################
#                                                                              #
#	Evolution Canyon Project: Microbial Community PcoA and PERMANOVA             #
#   Analysis of MICROBIDE output 1                                             #
#                                                                              #
################################################################################
#                                                                              #
#	Written by: Mario Muscarella                                                 #
#                                                                              #
#	Last update: 2014/03/12                                                      #
#                                                                              #
################################################################################

# Setup Work Environment
rm(list=ls())
setwd("~/GitHub/evolution-canyon/community")
source("functions/ec.pcoa.r")

# Run analysis of MICROBIDE simulation data

input <- list.files("./data", "*optima.shared")

permanova.output <- list()

for (i in 1:length(input)){
  permanova.output[[input[i]]] = ec.pcoa(shared = paste("./data/", input[i], 
  sep=""), design = "./data/simmy.design.txt", plot.title = input[i])
  }
  
permanova.output # Prints PERMANOVA Results for All
