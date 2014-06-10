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

### Some reading on multilevel paried data analysis ###
	# PLSDA: http://link.springer.com/article/10.1007/s11306-009-0185-z/fulltext.html
	# mixOmics: http://cran.r-project.org/web/packages/mixOmics/mixOmics.pdf
	# mixOmics: http://perso.math.univ-toulouse.fr/mixomics/methods/spls-da/
	# mixOmics: http://perso.math.univ-toulouse.fr/mixomics/case-studies/multilevelliver-toxicity/3/
	# mixOmics: http://www.inside-r.org/packages/cran/mixOmics/docs/data.simu
	# mixOmics: http://perso.math.univ-toulouse.fr/mixomics/faq/pre-processing/
	# Standardizations: http://cc.oulu.fi/~jarioksa/softhelp/vegan/html/decostand.html

# Load packages
library(vegan)
library(mixOmics)
library(som)

### MULTILEVEL ANALYSIS ###

	# Creating matrices for mixOmics
	X<-t(data)
	X<-normalize(X,byrow=TRUE) # normalizing by row (mean = 0, variance = 1)
	tax<-t(taxonomy[,1])
	colnames(X)<-tax

	land<-as.factor(c("T1","T1","T1","T1","T1","T1","T7","T7","T7","T7","T7","T7","DF","DF","DF","DF","CF","CF","CF","CF","CF","CF")) # factor 1

	rewet<-as.factor(c("dry","wet","dry","wet","dry","wet","dry","wet","dry","wet","dry","wet","dry","wet","dry","wet","dry","wet","dry","wet","dry","wet")) # factor 2

	land.rewet<-data.frame(cbind(as.character(land),as.character(rewet))) # combine factor 1 and 2

	paired<-as.integer(c(1,1,2,2,3,3,1,1,2,2,3,3,1,1,2,2,1,1,2,2,3,3)) # ID for paired within land

SIP_multilevel <- multilevel(X, cond = land.rewet, sample = paired, ncomp = 3, method = 'splsda')
	# # sPLS-DA centers and scales the variables	
	# # centering, presumably subtracts mean value from each column
	# # scaling may be dividing centered by sd or root mean square; not 100% sure