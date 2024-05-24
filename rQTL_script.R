# Load the software
library(qtl2)
library(qtl)
library(tidyverse)

DU <- read.cross(format="csvs", dir="D://chpt4_barbneom_QTL/QTL/20240521_scripts_that_work/20240524_Q30_manually_ordered/", 
                 genfile="Allchromos.20240524_Q30.manually.ordered.genotypes.csv",
                 phefile="20240524_F2_trait_averages_lnnectvol.csv",
                 estimate.map = F)

summary(DU)
par(mfrow = c(1,1))
plot.map(DU, show.marker.names = FALSE)

#to drop the duplicate markers
dupmark <- findDupMarkers(DU)
DU_nodups <- drop.markers(DU, dupmark)
summary(DU_nodups)

#summaryMap() provides a summary of the average inter-marker distance and 
# the largest gap on each chromosome
summary.map(DU_nodups)

#plot the map
plot.map(DU_nodups)




#########
# rQTL2 #
#########

DU <- DU_nodups

# First convert cross into RQTL2 format:
DU_Q <- convert2cross2(DU)

# insert pseudomarkers
map_DU <- insert_pseudomarkers(DU_Q$gmap, step=1)

# calculate genotype probabilities
pr_DU <- calc_genoprob(DU_Q, map_DU, error_prob = 0.0002)
apr <- genoprob_to_alleleprob(pr_DU)

# Performing a genome scan
out <- scan1(pr_DU, DU_Q$pheno)
par(mar=c(5.1, 4.1, 1.1, 1.1))
ymx <- maxlod(out) # overall maximum LOD score

# Plot LOD Score vs Position
plot(out, map_DU, lodcolumn=1, col="darkviolet", ylim=c(0, ymx*1.02))
plot(out, map_DU, lodcolumn=2, col="darkturquoise", add=TRUE)
plot(out, map_DU, lodcolumn=3, col="yellow", add=TRUE)
plot(out, map_DU, lodcolumn=4, col="seagreen", add=TRUE)
plot(out, map_DU, lodcolumn=5, col="navy", add=TRUE)
plot(out, map_DU, lodcolumn=6, col="darkred", add=TRUE)
plot(out, map_DU, lodcolumn=7, col="orange", add=TRUE)
legend("topright", lwd=2, 
       col=c('darkviolet', 'darkturquoise', 'yellow', 'seagreen', 'navy', 'darkred', 'orange'), 
       colnames(out), bg="gray90")

# Find peaks
find_peaks(out, map_DU, threshold=4, drop=1.5)

# LOD support
lod_int(out, map_DU, lodcolumn=1, peakdrop=1.8, drop=1.5)

# Bayes support
bayes_int(out, map_DU, lodcolumn=1, prob=0.95)

# Permutation test
operm <- scan1perm(pr_DU, DU_Q$pheno, n_perm=1000)
summary(operm)
