################################################################################
#                                                                              #
#	Evolution Canyon Project: Microbial Community PcoA and PERMANOVA             #
#   Analysis of BAD sequncing run data                                         #
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
setwd("~/GitHub/evolution-canyon")
source("bin/ec.pcoa.r")

# Run analysis of Bad Evolution Canyon data set
ec.test <- ec.pcoa(shared = "./data/EC.bad.shared", design = "./data/design.txt",
  plot.title = "EC.test")
  
ec.test # Prints PERMANOVA Results

#testing for J with Git