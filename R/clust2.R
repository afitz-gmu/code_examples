################################################################################
##################################### INIT #####################################
################################################################################
rm(list=ls(all=TRUE))
# load required libraries
library(igraph)   # to cluster
library(raster)   # to cluster
library(spatstat) # for smoothing

setwd("~/Dropbox/Research/Frazier Research/Nigeria_Paper")
load("lga_matrix.RData")

################################################################################
################################### FUNCTION ###################################
################################################################################

# helper function to plot non-reversed images
image2 <- function(mat, ...) {
  mat <- apply(mat, 2, rev)
  image(t(mat), ...)
}

# inputs:
# M - a matrix with population numbers in cells.
# dat - 3 column data frame with x-coordinate, y-coordinate and population.
# cutoff - cutoff value (between 0 and 100) to cut the density and for clusters.
# doPlot - if TRUE then plots are produced. If FALSE then no plotting is done.
smoothCluster <- function(M, dat, cutoff=1000, sigma=1, plot1="pic1.png", plot2="pic2.png") {
  # estimate the density using kde function from package ks
  fit <- as.matrix(blur(as.im(M), sigma=sigma, bleed=FALSE))
  # cut the population at specified cutoff level.
  # this produces a matrix of TRUE/FALSE values. TRUE if the density is above
  # cutoff and false otherwise.
  A <- fit > cutoff
  # cluster the produced "islands" using the clump function from raster package.
  clusters <- as.matrix(clump(raster(A)))
  clusters[is.na(clusters) & !is.na(fit)] <- 0

  # produce plots
  png(plot1, width=400, height=400)
  # pick pretty colors
  colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F",
                               "yellow", "#FF7F00", "red", "#7F0000"
                               ))(100)
  levels <- seq(0, 20000, length.out=100)

  # parameters for plot (2 plots on one page + margin settings)
  savedPar <- par(mar=c(2,1,2,1))
  layout(matrix(c(1,2), nrow=1), widths=c(0.8,0.2))
  # plot estimated densities)
  image2(fit, col=colors, xaxt='n', yaxt='n', zlim=c(0, 20000))
  mtext("Average Smoothed Population", 1) # add text at the bottom margin
  # add legend for density levels
  par(mar=c(1,0,2,5))
  plot(NA, xlim=c(0,0), ylim=c(0,20000), frame=FALSE, axes=F, xaxs="i", yaxs="i", ylab="", xlab="")
  rect(0, levels[-length(levels)], 1, levels[-1L], col=colors, border=FALSE)
  axis(4, round(levels[seq(0,100,10)]), labels=TRUE, las=2)
  # add cut point
  abline(h=cutoff, lwd=2, col="purple")
  dev.off()

  # produce an image of clusters
  png(plot2, width=400, height=400)
  par(mar=c(2,1,2,1))
  img <- clusters
  img[is.na(img)] <- -1
  img[clusters!=0] <- 1
  image2(img, xaxt="n", yaxt="n", col=c("white", "grey", "orange"))
  mtext("Clustered Cities", 1) # add text at the bottom margin
  # restore the changed plot paramteres to saved dfaults
  par(savedPar)
  dev.off()

  # return the answer
  dat <- cbind(dat, t(clusters)[!is.na(t(clusters))])

  # check if merging went OK
  if(!all.equal(dat[,3], t(M)[!is.na(t(M))])) {
    stop("something is wrong")
  }
  list(clusters=clusters, data=dat)
}

################################################################################
##################################### RUN ######################################
################################################################################

############################### obtain clusters ################################

M <- as.matrix(lga_mask_05) # population matrix
res <- smoothCluster(M, lga_mask_05_matrix, cutoff=750, sigma=1, plot1="fig1.png", plot2="fig2.png")

df  <- as.data.frame(res$data)
names(df) <- c("X", "Y", "Population", "Cluster")
subsetted_df <- subset(df, Cluster > 0)
clust <- res$clusters

test_ordering <- tapply(df[,3], df[,4], sum)
test_ordering <- as.data.frame(test_ordering)
data_ordered <- as.data.frame(test_ordering[order(-test_ordering), ])
data_ordered$Rank <- 1:nrow(data_ordered)
names(data_ordered) <- c("Total_Cluster_Population", "Rank")

##################################### Cool Interactive Plot #####################################

library(rgl)
cols <- c("darkblue", "blue", "blue", "cornflowerblue", "cornflowerblue", "lightblue")
cols <- colorRampPalette(cols)(100)
colors <- cols[ceiling(M/max(M[clust==0],na.rm=TRUE)*100)]
colors[clust>0] <- "red"
persp3d(M, col=colors, lit=FALSE, specular="black")

