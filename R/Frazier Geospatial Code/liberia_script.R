remove(list = ls())

setwd("~/Dropbox/Research/Frazier Research/Liberia_new")

#library("xlsx")
#library(foreign)
install.packages("VIM", dependencies = TRUE)
install.packages("VIMGUI", dependencies = TRUE)

library(VIM)
library(VIMGUI)
library(sp)  # vector data
library(raster)  # raster data
library(rgdal)  # input/output, projections
library(rgeos)  # geometry ops
#library(spdep)  # spatial dependence
#install.packages("nnet")
library(nnet)

# Ghana Data

# enum_areas <- readOGR(dsn="ghana_data", layer="ghana_enums_ft", stringsAsFactors=FALSE, verbose=FALSE)
# towns <- readOGR(dsn="ghana_data", layer="ghana_twns_ft", stringsAsFactors=FALSE, verbose=FALSE)
# glss5 <- read.csv("ghana_data/sec1.csv", stringsAsFactors = FALSE, colClasses='character')
# codes <- read.csv("ghana_data/glss5_census_codes.csv", stringsAsFactors = FALSE, colClasses='character')

#west_eas <- subset(enum_areas, REG_CODE == "01")

# Liberia Data
load("cwiq10.RData")

clans <- readOGR(dsn="liberia_data", layer="clans", stringsAsFactors=FALSE, verbose=FALSE)
clans@data <- clans@data[ ,c(3,4,1,9:11,2,5:8)]
names(clans@data) <- c("county", "district", "clan", "cnty_code", "dist_code", "clan_code", "ea_count", "ttl_pop", "ttl_men", "ttl_women", "ttl_hhs")

clans@data[which(clans@data$county == "Bomi"), ]$cnty_code <- "03"
clans@data[which(clans@data$county == "Bomi" & clans@data$district == "Klay"), ]$dist_code <- "02"
clans@data[which(clans@data$county == "Bomi" & clans@data$district == "Suehn Mecca"), ]$dist_code <- "04"
clans@data[which(clans@data$county == "Bomi" & clans@data$district == "Senjeh"), ]$dist_code <- "06"
clans@data[which(clans@data$county == "Bomi" & clans@data$district == "Dowein"), ]$dist_code <- "08"
clans@data[which(clans@data$county == "Bomi" & clans@data$district == "Suehn Mecca" & clans@data$clan == "Gboor"), ]$clan_code <- NA

clans@data$code <- with(clans@data, paste(cnty_code, dist_code, clan_code, sep = ""))

View(clans@data)

sum(as.numeric(clans@data$ea_count))
sum(as.numeric(clans@data[which(clans@data$county == "Bomi"), ]$ea_count))

length(unique(clans@data$code))
table(clans@data$code)

#cwiq10_raw <- read.spss("liberia_data/cwiq10.sav", to.data.frame = TRUE, use.value.labels = FALSE)
#which( colnames(cwiq10)=="wta_hh" )

load("liberia_data/cwiq10.RData")

cwiq10$county <- substr(cwiq10$hid_mungai, 1, 2)
cwiq10$district <- substr(cwiq10$hid_mungai, 3, 4)
cwiq10$clan_town <- substr(cwiq10$hid_mungai, 5, 7)
cwiq10$ea_no <- substr(cwiq10$hid_mungai, 8, 10)
cwiq10$structure_no <- substr(cwiq10$hid_mungai, 16, 18)
cwiq10$hh_no <- substr(cwiq10$hid_mungai, 11, 12)


cwiq10$codi <- substr(cwiq10$hid_mungai, 1, 4)
cwiq10$code <- substr(cwiq10$hid_mungai, 1, 7)
#cwiq10$enum_area <- substr(cwiq10$hid_mungai, 1, 10)
cwiq10$hhid <- substr(cwiq10$hid_mungai, 1, 12)          ### <--- comment on
#cwiq10$structure <- substr(cwiq10$hid_mungai, 1, 18)

hsize <- table(cwiq10$hhid)          ### <--- add this line
cwiq10$hhid <- as.integer(hsize[as.character(cwiq10$hhid)])          ### <--- add this line

length(unique(cwiq10$codi))
table(cwiq10$codi)

length(unique(cwiq10$codi))
table(cwiq10$codi)

#length(unique(cwiq10$enum_area))
#length(unique(cwiq10$structure))
#length(unique(cwiq10$hhid))
#length(unique(cwiq10$hid_mungai))

x <- cwiq10[which(cwiq10$codi == "0104"), ]

#model

mod1 <- lm(n11 ~ n1 + n4, data = cwiq10)
null <- lm(n11 ~ 1, data = cwiq10)
summary(mod1)
AIC(mod1)
AIC(null)
# add factors, discretize age

x$n1 <- as.factor(x$n1)
#age <- x$n4
x$n11 <- as.factor(x$n11)
x$n12 <- as.factor(x$n12)

mnlr <- multinom(n11 ~ n1, x)
summary(mnlr)
head(fitted(mnlr))

mnlr <- multinom(n11 ~ 1, x)
summary(mnlr)

predictors    <- expand.grid(n1=c("1","2"))
p.fit         <- predict(mnlr, predictors, type="probs")
probabilities <- data.frame(predictors,p.fit)
probabilities