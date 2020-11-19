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

################################################################################
################################### FUNCTION ###################################
################################################################################

# inputs:
# x - a matrix with 3 columns - x coordinates, y coordinates and weights
# cutoff - cutoff value (between 0 and 100) to cut the density and for clusters
# doPlot - if TRUE then plots are produced. If FALSE then no plotting is done.
densityCluster <- function(x, cutoff=5, bw=0.005, plot1 = "picbw1.png", plot2 = "picbw2.png") {
  lim1 <- range(x[,1]) # range of x coordinates
  lim2 <- range(x[,2]) # range of y coordinates
  weights <- x[,3] * nrow(x) / sum(x[,3])  # rescale weights
  # estimate the density using kde function from package ks
  H <- cbind(c(bw,0),c(0,bw))
  fit <- kde(x=x[,1:2], binned=TRUE, compute.cont=TRUE, xmin=c(lim1[1], lim2[1]),
             xmax=c(lim1[2], lim2[2]), bgridsize=rep(ceiling(sqrt(nrow(x))), 2),
             w=weights, H=H
             )
  # normalize to a score between 0 and 100
  fit$estimate <- fit$estimate - min(fit$estimate)
  fit$estimate <- fit$estimate / max(fit$estimate) * 100
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
  # produce plots
    png(plotbw1, width = 400, height = 400)
    # pick pretty colors
    colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F",
                                "yellow", "#FF7F00", "red", "#7F0000"
                                ))(100)
    levels <- seq(min(fit$estimate), max(fit$estimate), length.out=100)

    # parameters for plot (2 plots on one page + margin settings)
    savedPar <- par(mar=c(2,2,2,2))
    layout(matrix(c(1,2,3), nrow=1), widths=c(0.4,0.1,0.5))
    # plot estimated densities)
    plot(NA, xlim=range(x[,1]), ylim=range(x[,2]), axes=FALSE, xlab="", ylab="")
    .filled.contour(fit$eval.points[[1]], fit$eval.points[[2]], fit$estimate,
                    levels, colors
                    )
    mtext("estimated density", 1) # add text at the bottom margin
    # add legend for density levels
    par(mar=c(1,0,2,3))
    plot(NA, xlim=c(0,0), ylim=c(0,100), frame=FALSE, axes=F, xaxs="i", yaxs="i", ylab="", xlab="")
    rect(0, levels[-length(levels)], 1, levels[-1L], col=colors, border=FALSE)
    axis(4, seq(0,100,10), labels=TRUE, las=2)
    # add cut point
    abline(h=cutoff, lwd=2, col="purple")

    # produce an image of clusters
    png(plotbw2, width = 400, heigh = 400)
    par(mar=c(2,1,2,1))
    img <- clusters
    img[clusters!=0] <- 1
    image(img, xaxt="n", yaxt="n")
    mtext("formed clusters", 1) # add text at the bottom margin
    # restore the changed plot paramteres to saved dfaults
    par(savedPar)
    dev.off()
  
  # return the answer
  cbind(x, clust)
}


################################################################################
##################################### RUN ######################################
################################################################################

# obtain clusters
res <- densityCluster(lga_mask_05_matrix, cutoff=5, bw=.005)
totals_cluster <- tapply(res[,3], res[,4], sum)

# plot data and color it by cluster
plot(x[,1], x[,2], col=res[,4])

library(rgl)
cols <- ifelse(res[,4]==0, "grey", "red")
plot3d(res[,1], res[,2], res[,3], col=cols, type="h")


