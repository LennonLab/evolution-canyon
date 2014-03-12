################################################################################
#                                                                              #
#	Evolution Canyon Project: Microbial Community PcoA and PERMANOVA             #
#   Analysis of simmulation data                                               #
#                                                                              #
################################################################################
#                                                                              #
#	Written by: Mario Muscarella                                                 #
#                                                                              #
#	Last update: 2014/03/11                                                      #
#                                                                              #
################################################################################

# Setup Work Environment
rm(list=ls())
setwd("~/GitHub/evolution-canyon/community")
source("functions/ec.pcoa.r")

# Run analysis of Bad Evolution Canyon data set
ec.simmy1 <- ec.pcoa(shared = "./data/simmy1.shared", design = "./data/simmy.design.txt",
  plot.title = "EC.simmy2")
  
ec.simmy1 # Prints PERMANOVA Results