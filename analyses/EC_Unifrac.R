################################################################################
#                                                                              #
# Evolution Canyon Project: Microbial Community UniFrac                        #
#                                                                              #
################################################################################
#                                                                              #
# Written by: Mario Muscarella                                                 #
#                                                                              #
# Last update: 2014/07/09                                                      #
#                                                                              #
################################################################################
#                                                                              #
# Notes: This code calculates UniFrac distances from the EC project            #
#        For more information on UniFrac see: Lozupone 2005                    #
#                                                                              #
# Issues: Currently, this code does not work; exploratory                      #
#                                                                              #
# Recent Changes:                                                              #
#         1. New Code                                                          #
#                                                                              #
# Future Changes (To-Do List):                                                 #
#         1. Pick R packages for UniFrac                                       #
#         2. Build Rooted Tree                                                 #
#         3. Compare with Mothur                                               #
#                                                                              #
################################################################################

# Setup Work Environment
rm(list=ls())
setwd("~/GitHub/evolution-canyon")
source("./bin/DiversityFunctions.r")

# Packages to pick from:
# 1. GUniFrac - GUniFrac package
# 2. unifrac - picante package
# 3. UniFrac - phyloseq package
