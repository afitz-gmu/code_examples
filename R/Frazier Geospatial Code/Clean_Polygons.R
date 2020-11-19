#Cleaned script from 2/5/17 - KM - Contour Clustering 

setwd("~/Dropbox/Research/Frazier Research/Nigeria_Cluster")

# install.packages("rgdal", dependencies = TRUE)
# install.packages("sp", dependencies = TRUE)
# install.packages("raster", dependencies = TRUE)
# install.packages("rgeos", dependencies = TRUE)
# install.packages("spatstat", dependencies = TRUE)
# install.packages("maptools", dependencies = TRUE)
# install.packages("ggmap", dependencies = TRUE)
# install.packages("spdep", dependencies = TRUE)
# install.packages("sf", dependencies = TRUE)
# install.packages("devtools")
# install_github("arcdiagram", username = "gastonstat")

library(rgdal)
library(sp)
library(raster)
library(rgeos)
library(spatstat)
library(maptools)
library(ggmap)
library(devtools)
library(arcdiagram)
library(igraph)
library(spdep)
library(sf)

############
## Import/Prep ##
############

wpop_lbr14 <- raster("LBR14adjv1/LBR14adjv1.tif")
clans <- readOGR(dsn = 'spatial_data', layer = "LBR_adm3", stringsAsFactors=FALSE, verbose=FALSE)

prj_lbr <- ("+proj=utm +zone=29 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
prj_wpop <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"


wpop_lbr14 <- projectRaster(wpop_lbr14, crs = prj_lbr)
clans <- spTransform(clans, CRS = prj_lbr)

bain <- subset(clans, NAME_3 == "Bain")

proj4string(wpop_lbr14)
proj4string(clans)
proj4string(bain)

wpop_bain14 <- crop(wpop_lbr14, bain)

wpop_bain14 <- mask(wpop_bain14, bain)

#################
## Description ##
#################


cellStats(wpop_lbr14, 'sum')

cellStats(wpop_bain14, 'sum')

clans@data$SUM_TOTAL[which(clans@data$CLNAME == "Bain")]

##############
## Contour Cluster Creation ##
##############

pop <- floor(cellStats(wpop_lbr14, 'sum'))
win <- as(clans,"owin")
set.seed(1)
bain_ppl <- rpoint(pop, f = as.im(wpop_lbr14), win = win)

## Use these two lines below to calculate appropriate sigma value for current ppp model##
bw <- bw.ppl(bain_ppl)
#bw2 <- bw.diggle(bain_ppl)

#replace the sigma=bw to bw2 if in use
bain_dens <- density.ppp(bain_ppl, sigma = bw)

plot(ecdf(bain_dens))

plot.new()
plot(bain_dens)
contour.im(bain_dens, col="red", levels=0.00025, add=T) #levels = 0.0002

Dsg <- as(bain_dens, "SpatialGridDataFrame")  # convert to spatial grid class
Dim <- as.image.SpatialGridDataFrame(Dsg)  # convert again to an image
Dcl <- contourLines(Dim, levels = .00045)  # create contour object
SLDF <- ContourLines2SLDF(Dcl, CRS(prj_lbr))
plot(SLDF)

#DF
PS1 <- SpatialLines2PolySet(SLDF)


#Spatial Polygons
SP1 <- PolySet2SpatialPolygons(PS1, close_polys=TRUE)


#DF
PS2 <- SpatialPolygons2PolySet(SP1)
SL1 <- PolySet2SpatialLines(PS2)

Polyclust <- gPolygonize(SL1)

gas <- gArea(Polyclust, byid = T)

Polyclust2 <- SpatialPolygonsDataFrame(Polyclust, data = data.frame(gas), match.ID = F)
plot(Polyclust2)
Polyclust2 <- spTransform(Polyclust2, "+proj=utm +zone=29 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")

pxy <- cbind(bain_ppl$x,bain_ppl$y)
pxy <- as.data.frame(pxy)

pxy$observation <- 1:nrow(pxy) 
pxy$observation <- as.data.frame(pxy$observation)

dfxy <- as.data.frame(cbind(pxy$V1,pxy$V2))

chocho <- SpatialPointsDataFrame(dfxy,pxy$observation)
proj4string(chocho) <- CRS("+proj=utm +zone=29 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")

cAg <- aggregate(chocho, by = Polyclust2, FUN = length)

polys <- polygons(Polyclust2)

##############
## Analysis ##
##############

r2 <- as.data.frame(extract(wpop_bain14,polys, df=T))
r2.1 <- aggregate(r2$LBR14adjv1, by=list(r2$ID), FUN=sum)
sum(r2.1$x)
plot(polys[211,])

bain_extract <- extract(wpop_bain14, Polyclust2, df=T)
plot(bain_extract)
colSums(bain_extract)

sum(wpop_bain14@data@values, na.rm=T)

Polyclust2@data["totals"] <- unlist(lapply(bain_extract, 'sum'))

plot(bain)
plot(Polyclust2, add =T)

sum(Polyclust2@data$totals) / pop

clans <- spTransform(clans, CRS = prj_wpop)
clans_f <- fortify(clans)

Polyclust <- spTransform(Polyclust, CRS = prj_wpop)
Polyclust_f <- fortify(Polyclust)

bain_extract_df <- extract(wpop_bain14, Polyclust, df = TRUE)
bain_area_totals <- aggregate(. ~ ID, bain_extract_df, sum)
bain_area_totals$ID <- as.character(bain_area_totals$ID)

Polyclust_f <- left_join(x = Polyclust_f, y = bain_area_totals, by = c("id" = "ID"))

bain_center_pts <- gCentroid(Polyclust, byid=TRUE)
bain_cluster_ttls <- cbind.data.frame(bain_center_pts@coords[,1], bain_center_pts@coords[,2], bain_area_totals[,2])
names(bain_cluster_ttls) <- c("x", "y", "totals")


##############
## Plots ##
##############

bain_map <- get_map(location = c(-9.12, 7.12, -8.94, 7.26), maptype = "watercolor")
bain_map <- ggmap(bain_map)

bain_map <- bain_map + geom_map(data=clans_f, map=clans_f, aes(x=long, y=lat, map_id=id), color ="white", fill ="orangered4", alpha = .1, lwd=.4)
bain_map <- bain_map + geom_map(data=Polyclust_f, map=Polyclust_f, aes(x=long, y=lat, map_id=id), color ="white", fill ="orangered4", alpha = .4, lwd=.25)

bain_map <- bain_map + geom_point(aes(x=x, y=y, size = totals, colour = totals), data = bain_cluster_ttls)
bain_map <- bain_map + scale_colour_gradient(low="lightblue", high="darkblue", space="Lab")

bain_map <- bain_map + labs(x = "longitude", y = "latitude")
bain_map <- bain_map + theme(legend.position=c(0.1, 0.2))
bain_map <- bain_map + guides(size=FALSE)

bain_map
ggsave("bain_map.pdf", bain_map)
