################################################################################
#                                                                              #
#	Evolution Canyon Project: Microbial Community PLS-DA and PERMANOVA           #
#   Analysis of MICROBIDE model output                                         #
#                                                                              #
################################################################################
#                                                                              #
#	Written by: Mario Muscarella                                                 #
#                                                                              #
#	Last update: 2014/07/25                                                      #
#                                                                              #
################################################################################

# Setup Work Environment
rm(list=ls())
setwd("~/GitHub/evolution-canyon/")
source("bin/ec.plsda.fun.R")                              
source("bin/ec.pcoa.R")

# Run analysis of MICROBIDE simulation data

input <- list.files("./microbide/SbyS/", "*.txt")

test1 <- ec.plsda(shared     = "./microbide/SbyS/Condition1.txt",
                     cutoff     = "0.03",
                     design     = "./data/simmy.design.txt")
                     
test1plot <- ec.plsda.plot(plsda.in = "test1")

test1 <- ec.pcoa.fun(shared     = "./microbide/SbyS/Condition6.txt",
                     cutoff     = "0.03",
                     design     = "./data/simmy.design.txt")

test1plot <- ec.pcoa.plot(pcoa.in = test1)
                            


test1 <- ec.plsda(shared     = "./microbide/SbyS/Condition5.txt",
                     cutoff     = "0.03",
                     design     = "./data/simmy.design.txt")
                     
test1plot <- ec.plsda.plot(plsda.in = "test1")



permanova.output <- list()

for (i in 1:length(input)){
  permanova.output[[input[i]]] = ec.pcoa(shared = paste("./data/", input[i], 
  sep=""), design = "./data/simmy.design.txt", plot.title = input[i])
  }
  
permanova.output # Prints PERMANOVA Results for All
