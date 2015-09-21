################################################################################
#                                                                              #
# Evolution Canyon Project: Soil Gravametic Water                              #
#                                                                              #
################################################################################
#                                                                              #
# Written by: Mario Muscarella                                                 #
#                                                                              #
# Last update: 2014/08/05                                                      #
#                                                                              #
################################################################################
#                                                                              #
# Notes: This code plots soil moisture                                         #
#                                                                              #
################################################################################

# Setup Work Environment
rm(list=ls())
setwd("~/GitHub/evolution-canyon")
se <- function(x, ...){sd(x)/sqrt(length(x))}

ec.soil <- read.delim("./data/ec_moisture.txt", header=T)

slope <- ec.soil[,1]
GWC <- ec.soil[,7]

gwc.es <- GWC[1:3]
gwc.as <- GWC[4:6]

EC.slopes <- c("Mesic", "Xeric")
EC.means <- c(mean(gwc.es), mean(gwc.as))
EC.ses   <- c(se(gwc.es), se(gwc.as))

png(file="./plots/ec_gwc.png", width=800, height=600, antialias = "cleartype")
par(mar=c(0.5,3,4,0.5), oma=c(1,1,1,1)+0.1, lwd=2)
ec_plot <- barplot(EC.means, names.arg=EC.slopes, ylim=c(0,0.4), col="gray",lwd=3, yaxt="n", xaxt="n", width=0.5, space=c(2,1))
arrows(x0 = ec_plot, y0 = EC.means, y1 = EC.means - EC.ses, angle = 90, length=0.1, lwd = 2)
arrows(x0 = ec_plot, y0 = EC.means, y1 = EC.means + EC.ses, angle = 90, length=0.1, lwd = 2)
axis(side = 2, labels=T, lwd.ticks=2, las=2, lwd=3, cex=2, cex.axis=2)
mtext("Gravimetic Water Content", side = 3, cex = 3, line = 1)
mtext(expression(paste("(g H"["2"], "O * g Soil"^" -1",")")), side=3, cex=2, line = -1.5)
dev.off()

