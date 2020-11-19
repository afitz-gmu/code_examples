################################################################################
##################################### INIT #####################################
################################################################################

# load required libraries
library(ks)      # to estimate the density
library(igraph)  # to cluster
library(raster)  # to cluster
library(MASS)    # just for the toy example below (geyser)

setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/Nigeria_Paper")
load("lga_matrix.RData")
load("LGAtotals.RData")

################################################################################
################################### FUNCTION ###################################
################################################################################

# helper function to plot non-reversed images
image2 <- function(mat, ...) {
  mat <- apply(mat, 2, rev)
  image(t(mat), ...)
}

# function to obtain density
getDensity <- function(x, bw= bw) {
  lim1 <- range(x[,1]) # range of x coordinates
  lim2 <- range(x[,2]) # range of y coordinates
  weights <- x[,3] * nrow(x) / sum(x[,3])  # rescale weights
  # estimate the density using kde function from package ks
  H <- cbind(c(bw,0),c(0,bw))
  fit <- kde(x=x[,1:2], binned=TRUE, compute.cont=TRUE, xmin=c(lim1[1], lim2[1]),
             xmax=c(lim1[2], lim2[2]), bgridsize=c(1441, 1155),
             w=weights, H=H
             )
  # normalize to a score between 0 and 100
  fit$estimate <- fit$estimate - min(fit$estimate)
  fit$estimate <- fit$estimate / max(fit$estimate) * 100
  fit
}

densityCluster <- function(x, fit, cutoff= cutoff) {
  # cut the density at specified cutoff level.
  # this produces a matrix of TRUE/FALSE values. TRUE if the density is above
  # cutoff and false otherwise.
  A <- fit$estimate > cutoff
  # cluster the produced "islands" using the clump function from raster package.
  clusters <- as.matrix(clump(raster(A)))
  clusters[is.na(clusters)] <- 0
  # traverse data rows and assign them to clusters
  clust <- numeric(nrow(x)) # clusters will be stroed in this vector
  for(i in 1:nrow(x)) {
    minx1 <- which.min(abs(x[i,1]-fit$eval.points[[1]])) # closest estimated density element by x coordinate
    minx2 <- which.min(abs(x[i,2]-fit$eval.points[[2]])) # closest estimated density element by y coordinate
    clust[i] <- clusters[minx1, minx2] # assigned cluster based on the value of density that is closest to the point
  }
  # return the answer
  cbind(x, clust)
}

plotDensities <- function(dens, mask, cutmark) {
  E <- apply(t(dens$estimate), 2, rev)
  M <- as.matrix(mask)
  E[is.na(M)] <- NA

  colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F",
                               "yellow", "#FF7F00", "red", "#7F0000"
                               ))(100)
  levels <- seq(0, max(pretty(E)), length.out=100)

  # parameters for plot (2 plots on one page + margin settings)
  savedPar <- par(mar=c(2,1,2,1))
  layout(matrix(c(1,2), nrow=1), widths=c(0.8,0.2))
  # plot estimated densities)
  image2(E, col=colors, xaxt='n', yaxt='n', zlim=c(0, max(levels)))

  par(mar=c(1,0,2,5))
  plot(NA, xlim=c(0,0), ylim=c(0,max(levels)), frame=FALSE, axes=F, xaxs="i", yaxs="i", ylab="", xlab="")
  rect(0, levels[-length(levels)], 1, levels[-1L], col=colors, border=FALSE)
  axis(4, round(levels[seq(0,100,10)]), labels=TRUE, las=2)
  # add cut point
  abline(h=cutmark, lwd=2, col="purple")
  par(savedPar)
}

plotClusters <- function(fit, mask, cutoff) {
  A <- apply(t(fit$estimate), 2, rev) > cutoff
  # cluster the produced "islands" using the clump function from raster package.
  clusters <- as.matrix(clump(raster(A)))
  # plot
  savedPar <- par(mar=c(2,1,2,1))
  clusters[clusters!=0] <- 1
  clusters[is.na(clusters)] <- 0
  clusters[is.na(as.matrix(mask))] <- -1
  image2(clusters, xaxt="n", yaxt="n", col=c("white", "grey", "orange"))
  mtext("formed clusters", 1) # add text at the bottom margin
  par(savedPar)
}

################################################################################
##################################### RUN ######################################
################################################################################

# obtain density
# obtain clusters
cutoff <- 6
bw <- 0.001

dens <- getDensity(lga_mask_05_matrix, bw=bw)
res <- densityCluster(lga_mask_05_matrix, dens, cutoff=cutoff)
totals_cluster <- tapply(res[,3], res[,4], sum)

#################################### plot 1 ####################################

png("fig1.png", width=400, height=400)
plotDensities(dens, lga_mask_05, cutmark=cutoff)
dev.off()

png("fig2.png", width=400, height=400)
plotClusters(dens, lga_mask_05, cutoff=cutoff)
dev.off()

################################### plot 3d ####################################

# plot data and color it by cluster
plot(x[,1], x[,2], col=res[,4])

library(rgl)
cols <- ifelse(res[,4]==0, "grey", "red")
plot3d(res[,1], res[,2], res[,3], col=cols, type="h")


