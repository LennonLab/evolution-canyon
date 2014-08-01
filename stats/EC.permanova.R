################################################################################
#                                                                              #
#	Evolution Canyon Project: Microbial Community PERMANOVA                      #
#                                                                              #
################################################################################
#                                                                              #
#	Written by: Jay Lennon and Mario Muscarella                                  #
#                                                                              #
#	Last update: 2014/07/01                                                      #
#                                                                              #
################################################################################

#### Things left to do
	# use greenegenes or silva database
	# confirm factors are set up right (paired)
	# confirm normalization (relative recovery, log transformation)
	# figure out how to get explained total variance for each axis
	# figure out how to identify influential taxa using loadings
	# load taxonomy in (see Pmeso code on github)
	# consider removal of rare taxa
	# use Stuart's PERMANOVA that accounts for paired samples
	# confirm these are actually "bad": EC-2A-D, EC-2A-R, EC-2C-R, EC-2D-R

# ec.pcoa <- function(shared = " ", design = " ", plot.title = "test"){


# Setup Work Environment
rm(list=ls())
setwd("~/GitHub/evolution-canyon")
source("./bin/DiversityFunctions.r")
require("vegan")
require("som")

shared     = "./mothur/EC.bac.final.shared"
design     = "./data/design.txt"

#1 -- Import data (site by OTU matrix) and deal w/ "problem" samples
ec_data <- t(read.otu(shared, "0.03")) # takes a long time (10 mins)
design <- read.delim(design, header=T, row.names=1)


# Note: owing to amplification issues, we only submitted 76 of 80 samples.
# The four samples NOT submitted were C-1E-R, EC-2G-R, EC-2J-R, EC-6I-D
# So, removed their pairs for PLS-DA: EC-1E-D, EC-2G-D, EC-2J-D,EC-6I-R
# 4 samples have had crappy quality in past: EC-2A-D, EC-2A-R, EC-2C-R, EC-2D-R
# Automatically exlcuded, but may want to look again
# Remove their pairs EC-2C-D, EC-2D-R
# Results in a total of 66 samples (14 problematics)
# Use following code to identify columns to remove

samps <- (colnames(ec_data))
which(samps == "EC-2A-D")

# Remove problematic samples and their pairs
ec_data_red <- ec_data[,-c(9,20:21,24:27,32,37,74)]

# Remove any OTUs with <2 OTUs
ec_data_red <- ec_data_red[which(rowSums(ec_data_red) >= 2),]

samps_red <- colnames(ec_data_red) # recover samples from reduced dataset
design_red <- design[samps_red,] # design matrix w/o problem samples & pairs

design_red$paired <- c(rep(seq(1:14), each=2), rep(seq(1:19), each=2))

#2 -- Options for transformation, relativation, and normalization
Xt <- t(ec_data_red) # transpose sample x OTU matrix; necssary for 'multilelve'
Xlogt <- decostand(Xt,method="log")
Xpa <- (Xt > 0)*1 # Calculate Presense Absence
Xrel <- ec_data_red
 for(i in 1:ncol(ec_data_red)){
  Xrel[,i] = ec_data_red[,i]/sum(ec_data_red[,i])}
Xrelt <- t(Xrel)
Xlog <- decostand(ec_data_red,method="log")
Xrellog <- Xlog
 for(i in 1:ncol(Xlog)){
  Xrellog[,i] = Xlog[,i]/sum(Xlog[,i])}
Xrellogt<-t(Xrellog)
X <- normalize(Xrelt,byrow=TRUE) # normalizing by row (mean = 0, var = 1) with som package, not mixOmics

#3 -- Calculate PERMANOVA Stats, options for design below

X.in <- Xt    # Pick the transformation you want

# Molecule and Slope (no interactin, no grouping) No Transformation
Adonis.a <- adonis(t(ec_data_red) ~ design_red$molecule + design_red$slope)

Adonis.b <- adonis(t(ec_data_red) ~ design_red$molecule * design_red$slope)


# Molecule and Slope (no interactin, no grouping)  Transformation
Adonis.a <- adonis(Xrellogt ~ design_red$molecule + design_red$slope)

Adonis.b <- adonis(Xrellogt ~ design_red$molecule * design_red$slope)



# Using Strats
Adonis.c <- adonis(Xrellogt ~ design_red$molecule * design_red$slope, strata=design_red$paired)
Adonis.c <- adonis(Xrellogt ~ design_red$molecule * design_red$slope, strata=design_red$station)



# Manual permutation based on biobucket to test time effect --> results is similar to adonis
# Is this really what we want????


print(fit<-adonis(dist(X)~time,permutations=1))

# number of perms
B<-999

pop<-rep(NA,B+1)

pop[1]<-fit$aov.tab[1,5]

#### this is where I think the biobucket example goes wrong. The biobucket example constrained permutations with in "site", but if observations are paired, I think transformations should be constrained at "individual or replicate" level
#ctrl<-permControl(strata=land,within=Within(type="series",mirror=FALSE))
ctrl<-permControl(strata=landRep,within=Within(type="series",mirror=FALSE))

nobs=nrow(X)

shuffle(nobs,control=ctrl)

for(i in 2:(B+1)){
	idx<-shuffle(nobs,control=ctrl)
	fit.rand<-adonis(dist(X)~time[idx],permutations=1)
	pop[i]<-fit.rand$aov.tab[1,5]
}

print(pval<-sum(pop>=pop[1])/(B+1))

hist(pop,xlab="Population R2")
abline(v=pop[1],col=2,lty=3)
text(0.08,300,paste("true R2,\np=",pval,sep=""))=