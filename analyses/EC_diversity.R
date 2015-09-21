################################################################################
#                                                                              #
# Evolution Canyon Project: Microbial Community Alpha Diversity                #
#                                                                              #
################################################################################
#                                                                              #
# Written by: Mario Muscarella                                                 #
#                                                                              #
# Last update: 2014/07/11                                                      #
#                                                                              #
################################################################################
#                                                                              #
# Notes: This code estimates species richness in the EC Project                #
#                                                                              #
# Issues: None detected yet                                                    #
#                                                                              #
# Recent Changes:                                                              #
#         1. Just Started                                                      #
#                                                                              #
# Future Changes (To-Do List):                                                 #
#         1. Add Richness Estimation                                           #
#         2. Plot Richness                                                     #
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
level <- "0.03"
size   <- 50000

# Import Design File
design <- read.delim(design, header=T, row.names=1)
slope.molecule <- data.frame(cbind(as.character(design$slope),
  as.character(design$molecule))) # Y matrix with factor 1 and 2
design$slope.molecule.concat <- do.call(paste, c(slope.molecule[c("X1", "X2")],
  sep = "")) # create unique treat ID vector

# Calculate Sampling Depth & Coverage
counts <- count.groups(read.otu(shared, cutoff = level))
coverage <- coverage(input= shared, cutoff = level, size = size)
counts2 <- counts[counts > 50000]

design <- design[names(counts2), ]
  
# Alpha Diversity with Resampling
rich <- round(richness.iter(input = shared, cutoff = level, size = size,
  iters = 100), 3)
rich <- rich[names(counts2),]
even <- round(evenness.iter(input = shared, cutoff = level, size = size,
  iters = 100), 3)
even <- even[names(counts2),]
shan <- round(diversity.iter(input = shared, index = "shannon",
  cutoff = level, size = size, iters = 100), 3)
shan <- shan[names(counts2),]

# Richness Summary
rich_data <- merge(design[,1:4], rich, by="row.names")
rich_data$mean <- round(apply(rich, 1, mean, na.rm = TRUE),3)
rich_data$se <- round(apply(rich, 1, se, na.rm = TRUE), 3)
rich_data$Design <- as.factor(design$slope.molecule.concat)

# Evenness Summary
even_data <- merge(design[,1:4], even, by="row.names")
even_data$mean <- round(apply(even, 1, mean, na.rm = TRUE),3)
even_data$se <- round(apply(even, 1, se, na.rm = TRUE),3)
even_data$Design <- as.factor(design$slope.molecule.concat)

# Shannon Diversity Summary
shan_data <- merge(design[,1:4], shan, by="row.names")
shan_data$mean <- round(apply(shan, 1, mean, na.rm = TRUE),3)
shan_data$se <- round(apply(shan, 1, se, na.rm = TRUE),3)
shan_data$Design <- as.factor(design$slope.molecule.concat)

# Plot Parameters
windows.options(width=6, height=12)
par(mfrow=c(3,1), mar=c(0.25,5,0.25,1), oma=c(3,1,1,1)+0.1 )

# Richness Plot A
rich.means <- tapply(rich_data$mean, rich_data$Design, mean)
rich.se <- tapply(rich_data$mean, rich_data$Design, se)
dim.r <- round(c(min(rich.means)-(max(rich.se)*1.1),
  max(rich.means)+(max(rich.se)*1.1)),2)
plot(rich.means, ylim=dim.r, xlim=c(0.5,4.5), pch=15, cex=2, yaxt="n",
  cex.lab=2, cex.axis=1.25, las=1, xaxt="n", ylab="")
arrows(x0 = 1:4, y0 = rich.means, y1 = rich.means - rich.se, angle = 90,
  length=0.05, lwd = 2)
arrows(x0 = 1:4, y0 = rich.means, y1 = rich.means + rich.se, angle = 90,
  length=0.05, lwd = 2)
title(ylab = "Taxa Richness", cex.lab = 2, line = 3.5)
box(lwd=2)
axis(side = 1, labels=F, lwd.ticks=2)
axis(side = 2, labels=T, lwd.ticks=2, las=1, cex.axis=1.25)
axis(side = 1, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 2, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 3, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 4, tck=0.02, labels=F, lwd.ticks=2)

# Evenness Plot A
even.means <- tapply(even_data$mean, even_data$Design, mean)
even.se <- tapply(even_data$mean, even_data$Design, se)
dim.e <- round(c(min(even.means)-(max(even.se)*1.1),
  max(even.means)+(max(even.se)*1.1)),2)
plot(even.means, ylim=dim.e, xlim=c(0.5,4.5), pch=15, cex=2, yaxt="n",
  cex.lab=2, cex.axis=1.25, las=1, xaxt="n", ylab="")
arrows(x0 = 1:4, y0 = even.means, y1 = even.means - even.se, angle = 90,
  length=0.05, lwd = 2)
arrows(x0 = 1:4, y0 = even.means, y1 = even.means + even.se, angle = 90,
  length=0.05, lwd = 2)
title(ylab = "Taxa Evenness", cex.lab = 2, line = 3.5)
box(lwd=2)
axis(side = 1, labels=F, lwd.ticks=2)
axis(side = 2, labels=T, lwd.ticks=2, las=1, cex.axis=1.25)
axis(side = 1, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 2, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 3, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 4, tck=0.02, labels=F, lwd.ticks=2)

# Shannon Plot A
shan.means <- tapply(shan_data$mean, shan_data$Design, mean)
shan.se <- tapply(shan_data$mean, shan_data$Design, se)
dim.s <- round(c(min(shan.means)-(max(shan.se)*1.1),
  max(shan.means)+(max(shan.se)*1.1)),2)
plot(shan.means, ylim=dim.s, xlim=c(0.5,4.5), pch=15, cex=2, yaxt="n",
  cex.lab=2, cex.axis=1.25, las=1, xaxt="n", ylab="")
points(shan.means, pch=15, cex=2)
arrows(x0 = 1:4, y0 = shan.means, y1 = shan.means - shan.se, angle = 90,
  length=0.05, lwd = 2)
arrows(x0 = 1:4, y0 = shan.means, y1 = shan.means + shan.se, angle = 90,
  length=0.05, lwd = 2)
title(ylab = "Shannon Diversity", cex.lab = 2, line = 3.5)
box(lwd=2)
axis(side = 1, labels=c('AFDNA', 'AFRNA', 'EUDNA', 'EURNA'), at=1:4,
  lwd.ticks=2, cex.axis=1.25)
axis(side = 2, labels=T, lwd.ticks=2, las=1, cex.axis=1.25)
axis(side = 1, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 2, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 3, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 4, tck=0.02, labels=F, lwd.ticks=2)


# Plot Parameters
dev.off()
windows.options(reset=T)
windows.options(width=6, height=6)
par(mar=c(0.25,5.5,0.25,0.25), oma=c(3,1,1,1)+0.1 )

# Richness Plot A
rich.means <- tapply(rich_data$mean, rich_data$Design, mean)
rich.se <- tapply(rich_data$mean, rich_data$Design, se)
dim.r <- round(c(min(rich.means)-(max(rich.se)*1.1),
  max(rich.means)+(max(rich.se)*1.1)),2)
plot(rich.means, ylim=dim.r, xlim=c(0.5,4.5), pch=15, cex=2, yaxt="n",
  cex.lab=2, cex.axis=1.25, las=1, xaxt="n", ylab="")
arrows(x0 = 1:4, y0 = rich.means, y1 = rich.means - rich.se, angle = 90,
  length=0.05, lwd = 2)
arrows(x0 = 1:4, y0 = rich.means, y1 = rich.means + rich.se, angle = 90,
  length=0.05, lwd = 2)
title(ylab = "Taxa Richness", cex.lab = 2, line = 4)
axis(side = 1, labels=c('AFDNA', 'AFRNA', 'EUDNA', 'EURNA'), at=1:4,
  lwd.ticks=2, cex.axis=1.25)
box(lwd=2)
axis(side = 1, labels=F, lwd.ticks=2)
axis(side = 2, labels=T, lwd.ticks=2, las=1, cex.axis=1.25)
axis(side = 1, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 2, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 3, tck=0.02, labels=F, lwd.ticks=2)
axis(side = 4, tck=0.02, labels=F, lwd.ticks=2)
