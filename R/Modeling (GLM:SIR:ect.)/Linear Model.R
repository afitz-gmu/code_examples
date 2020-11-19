#### Load packages ####
library(FD)
library(spdep)
library(MuMIn)
library(vegan)
library(nlme)
library(dplyr)

#### Read in data ####

## PAR and SHDI values for every site:
landscape = read.csv("landscape.csv")
head(landscape)

## Trait values for each species:
traits = read.csv("traits.csv", row.names = "Species")
head(traits)

## Species abundance at each site:
abundance = read.csv("abund.csv", row.names = "Sites")
head(abundance)

## LUI value for each site:
lui = read.csv("lui.csv")
head(lui)

## Coordinates of each site:
coord = read.csv("coordinates.csv")
head(coord)

#### Calculate Shannon's diversity at each site ####
H_vector = diversity(abundance)
landscape = mutate(landscape, shannon_H = H_vector)
head(landscape)

#### Calculate CWM ####
cwm = dbFD(traits, abundance, calc.FRic = F, calc.FDiv = F, calc.CWM = T)$CWM
cwm$Site = rownames(cwm)
head(cwm)

#### Join together data into final data frame ####

## Join the separate data frames by common site names
data_final = left_join(landscape, cwm, by = "Site")
data_final = left_join(data_final, lui, by = "Site")
data_final = left_join(data_final, coord, by = "Site")

## Add a column that denotes the region of a plot
data_final$Region = c(rep("A", 32), rep("H", 28), rep("S", 29))

## Log-transform PAR at all spatial scales
data_final$PAR250 = log10(data_final$PAR250)
data_final$PAR500 = log10(data_final$PAR500)
data_final$PAR750 = log10(data_final$PAR750)
data_final$PAR1000 = log10(data_final$PAR1000)
data_final$PAR1250 = log10(data_final$PAR1250)
data_final$PAR1500 = log10(data_final$PAR1500)
data_final$PAR2000 = log10(data_final$PAR2000)

## Turn HW and RW into a coordinate object
coordinates(data_final) = c("RW", "HW")

## Take a look at the final dataset
head(data_final)


###################################################
#### Mixed models with Shannon's H as response ####
###################################################

## Create the null, LUI-only model to use at all spatial scales
h_lui = lme(data= data_final, method = "ML",
            control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
            fixed = shannon_H ~ LUI, random = ~1|Region,
            correlation = corSpher(form = ~RW + HW|Region, nugget = T))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 250 Spatial scale: ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## SHDI
h_shdi_250 = lme(data= data_final, method = "ML",
                 control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                 fixed = shannon_H ~ SHDI250, random = ~1|Region,
                 correlation = corSpher(form = ~RW + HW|Region, nugget = T))

## PAR
h_par_250 = lme(data= data_final, method = "ML",
                control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                fixed = shannon_H ~ PAR250, random = ~1|Region,
                correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + SHDI
h_lui_shdi_250 = lme(data= data_final, method = "ML",
                     control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                     fixed = shannon_H ~ LUI + SHDI250, random = ~1|Region,
                     correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + PAR
h_lui_par_250 = lme(data= data_final, method = "ML",
                    control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                    fixed = shannon_H ~ LUI + PAR250, random = ~1|Region,
                    correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + SHDI + PAR
h_full_250 = lme(data= data_final, method = "ML",
                 control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                     fixed = shannon_H ~ LUI + SHDI250 + PAR250, random = ~1|Region,
                     correlation = corSpher(form = ~RW + HW, nugget = T))

## Make a table of AICs for each model at 250m 
h_aic_250 = data.frame(model = c("LUI", "LUI + SHDI", "LUI + PAR", 
                                 "SHDI", "PAR", "LUI + SHDI + PAR"),
                       AIC = c(AICc(h_lui), AICc(h_lui_shdi_250), AICc(h_lui_par_250), 
                               AICc(h_shdi_250), AICc(h_par_250), AICc(h_full_250)),
                       spatial_scale = rep(250, 6))

## Arrange the table by AIC and display it
h_aic_250 = arrange(h_aic_250, AIC)
h_aic_250

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 500 Spatial scale: ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## SHDI
h_shdi_500 = lme(data= data_final, method = "ML",
                 control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                 fixed = shannon_H ~ SHDI500, random = ~1|Region,
                 correlation = corSpher(form = ~RW + HW, nugget = T))

## PAR
h_par_500 = lme(data= data_final, method = "ML",
                control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                fixed = shannon_H ~ PAR500, random = ~1|Region,
                correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + SHDI
h_lui_shdi_500 = lme(data= data_final, method = "ML",
                     control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                     fixed = shannon_H ~ LUI + SHDI500, random = ~1|Region,
                     correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + PAR
h_lui_par_500 = lme(data= data_final, method = "ML",
                    control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                    fixed = shannon_H ~ LUI + PAR500, random = ~1|Region,
                    correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + SHDI + PAR
h_full_500 = lme(data= data_final, method = "ML",
                 control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                 fixed = shannon_H ~ LUI + SHDI500 + PAR500, random = ~1|Region,
                 correlation = corSpher(form = ~RW + HW, nugget = T))

## Make a table of AICs for each model at 500m 
h_aic_500 = data.frame(model = c("LUI", "LUI + PAR", "LUI + SHDI + PAR",
                                 "PAR", "SHDI", "LUI + SHDI"),
                       AIC = c(AICc(h_lui), AICc(h_lui_par_500), AICc(h_full_500), 
                               AICc(h_par_500), AICc(h_shdi_500), AICc(h_lui_shdi_500)),
                       spatial_scale = rep(500, 6))

## Arrange the table by AIC and display it
h_aic_500 = arrange(h_aic_500, AIC)
h_aic_500

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 750 Spatial scale: ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## SHDI
h_shdi_750 = lme(data = data_final, method = "ML",
                 control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                 fixed = shannon_H ~ SHDI750, random = ~1|Region,
                 correlation = corSpher(form = ~RW + HW, nugget = T))

## PAR
h_par_750 = lme(data = data_final, method = "ML",
                control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                fixed = shannon_H ~ PAR750, random = ~1|Region,
                correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + SHDI
h_lui_shdi_750 = lme(data= data_final, method = "ML",
                     control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                     fixed = shannon_H ~ LUI + SHDI750, random = ~1|Region,
                     correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + PAR
h_lui_par_750 = lme(data= data_final, method = "ML",
                    control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                    fixed = shannon_H ~ LUI + PAR750, random = ~1|Region,
                    correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + SHDI + PAR
h_full_750 = lme(data= data_final, method = "ML",
                 control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                 fixed = shannon_H ~ LUI + SHDI750 + PAR750, random = ~1|Region,
                 correlation = corSpher(form = ~RW + HW, nugget = T))

## Make a table of AICs for each model at 750m 
h_aic_750 = data.frame(model = c("LUI", "LUI + SHDI", "LUI + PAR", 
                                 "LUI + SHDI + PAR", "SHDI", "PAR"),
                       AIC = c(AICc(h_lui), AICc(h_lui_shdi_750), AICc(h_lui_par_750), 
                               AICc(h_full_750), AICc(h_shdi_750), AICc(h_par_750)),
                       spatial_scale = rep(750, 6))

## Arrange the table by AIC and display it
h_aic_750 = arrange(h_aic_750, AIC)
h_aic_750

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 1000 Spatial scale: ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## SHDI
h_shdi_1000 = lme(data = data_final, method = "ML",
                  control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                  fixed = shannon_H ~ SHDI1000, random = ~1|Region,
                  correlation = corSpher(form = ~RW + HW, nugget = T))

## PAR
h_par_1000 = lme(data = data_final, method = "ML",
                 control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                 fixed = shannon_H ~ PAR1000, random = ~1|Region,
                 correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + SHDI
h_lui_shdi_1000 = lme(data= data_final, method = "ML",
                      control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                  fixed = shannon_H ~ LUI + SHDI1000, random = ~1|Region,
                  correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + PAR
h_lui_par_1000 = lme(data= data_final, method = "ML",
                     control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                     fixed = shannon_H ~ LUI + PAR1000, random = ~1|Region,
                     correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + SHDI + PAR
h_full_1000 = lme(data= data_final, method = "ML",
                  control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                  fixed = shannon_H ~ LUI + SHDI1000 + PAR1000, random = ~1|Region,
                  correlation = corSpher(form = ~RW + HW, nugget = T))

## Make a table of AICs for each model at 1,000m 
h_aic_1000 = data.frame(model = c("LUI", "LUI + SHDI", "LUI + PAR", 
                                  "LUI + SHDI + PAR", "SHDI", "PAR"),
                       AIC = c(AICc(h_lui), AICc(h_lui_shdi_1000), AICc(h_lui_par_1000), 
                               AICc(h_full_1000), AICc(h_shdi_1000), AICc(h_par_1000)),
                       spatial_scale = rep(1000, 6))

## Arrange the table by AIC and display it
h_aic_1000 = arrange(h_aic_1000, AIC)
h_aic_1000

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 1250 Spatial scale: ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## SHDI
h_shdi_1250 = lme(data = data_final, method = "ML",
                  control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                  fixed = shannon_H ~ SHDI1250, random = ~1|Region,
                  correlation = corSpher(form = ~RW + HW, nugget = T))

## PAR
h_par_1250 = lme(data = data_final, method = "ML",
                 control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                 fixed = shannon_H ~ PAR1250, random = ~1|Region,
                 correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + SHDI
h_lui_shdi_1250 = lme(data= data_final, method = "ML",
                      control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                  fixed = shannon_H ~ LUI + SHDI1250, random = ~1|Region,
                  correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + PAR
h_lui_par_1250 = lme(data= data_final, method = "ML",
                     control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                 fixed = shannon_H ~ LUI + PAR1250, random = ~1|Region,
                 correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + SHDI + PAR
h_full_1250 = lme(data= data_final, method = "ML",
                  control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                  fixed = shannon_H ~ LUI + SHDI1250 + PAR1250, random = ~1|Region,
                  correlation = corSpher(form = ~RW + HW, nugget = T))

## Make a table of AICs for each model at 1,250m 
h_aic_1250 = data.frame(model = c("LUI", "LUI + SHDI", "LUI + PAR", 
                                  "LUI + SHDI + PAR", "SHDI", "PAR"),
                        AIC = c(AICc(h_lui), AICc(h_lui_shdi_1250), AICc(h_lui_par_1250), 
                                AICc(h_full_1250), AICc(h_shdi_1250), AICc(h_par_1250)),
                        spatial_scale = rep(1250, 6))

## Arrange the table by AIC and display it
h_aic_1250 = arrange(h_aic_1250, AIC)
h_aic_1250

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 1500 Spatial scale: ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## SHDI
h_shdi_1500 = lme(data = data_final, method = "ML",
                  control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                  fixed = shannon_H ~ SHDI1500, random = ~1|Region,
                  correlation = corSpher(form = ~RW + HW, nugget = T))

## PAR
h_par_1500 = lme(data = data_final, method = "ML",
                 control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                 fixed = shannon_H ~ PAR1500, random = ~1|Region,
                 correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + SHDI
h_lui_shdi_1500 = lme(data= data_final, method = "ML",
                      control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                      fixed = shannon_H ~ LUI + SHDI1500, random = ~1|Region,
                      correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + PAR
h_lui_par_1500 = lme(data= data_final, method = "ML",
                     control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                     fixed = shannon_H ~ LUI + PAR1500, random = ~1|Region,
                     correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + SHDI + PAR
h_full_1500 = lme(data= data_final, method = "ML",
                  control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                  fixed = shannon_H ~ LUI + SHDI1500 + PAR1500, random = ~1|Region,
                  correlation = corSpher(form = ~RW + HW, nugget = T))

## Make a table of AICs for each model at 1,500m 
h_aic_1500 = data.frame(model = c("LUI", "LUI + SHDI", "LUI + PAR", 
                                  "LUI + SHDI + PAR", "SHDI", "PAR"),
                        AIC = c(AICc(h_lui), AICc(h_lui_shdi_1500), AICc(h_lui_par_1500), 
                                AICc(h_full_1500), AICc(h_shdi_1500), AICc(h_par_1500)),
                        spatial_scale = rep(1500, 6))

## Arrange the table by AIC and display it
h_aic_1500 = arrange(h_aic_1500, AIC)
h_aic_1500

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 2000 Spatial scale: ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## SHDI
h_shdi_2000 = lme(data = data_final, method = "ML",
                  control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                  fixed = shannon_H ~ SHDI2000, random = ~1|Region,
                  correlation = corSpher(form = ~RW + HW, nugget = T))

## PAR
h_par_2000 = lme(data = data_final, method = "ML",
                 control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                 fixed = shannon_H ~ PAR2000, random = ~1|Region,
                 correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + SHDI
h_lui_shdi_2000 = lme(data= data_final, method = "ML",
                      control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                      fixed = shannon_H ~ LUI + SHDI2000, random = ~1|Region,
                      correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + PAR
h_lui_par_2000 = lme(data= data_final, method = "ML",
                     control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                     fixed = shannon_H ~ LUI + PAR2000, random = ~1|Region,
                     correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + SHDI + PAR
h_full_2000 = lme(data= data_final, method = "ML",
                  control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                  fixed = shannon_H ~ LUI + SHDI2000 + PAR2000, random = ~1|Region,
                  correlation = corSpher(form = ~RW + HW, nugget = T))

## Make a table of AICs for each model at 2,000m 
h_aic_2000 = data.frame(model = c("LUI", "LUI + SHDI", "LUI + PAR", 
                                  "LUI + SHDI + PAR", "SHDI", "PAR"),
                        AIC = c(AICc(h_lui), AICc(h_lui_shdi_2000), AICc(h_lui_par_2000), 
                                AICc(h_full_2000), AICc(h_shdi_2000), AICc(h_par_2000)),
                        spatial_scale = rep(2000, 6))

## Arrange the table by AIC and display it
h_aic_2000 = arrange(h_aic_2000, AIC)
h_aic_2000

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Compare models at all spatial scales: ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Combine all AIC tables into one
h_aic_all = rbind(h_aic_250, h_aic_500, h_aic_750, h_aic_1000,
                  h_aic_1250, h_aic_1500, h_aic_2000)

## Add a new column: the difference in AICc from the null AICc
h_aic_all = mutate(h_aic_all, diff_AIC = AIC - AICc(h_lui))

## Arrange the table by AIC and display it
h_aic_all = arrange(h_aic_all, diff_AIC)
h_aic_all


######################################################
#### Mixed models with larval feeding as response ####
######################################################

## Create the null, LUI-only model to use at all spatial scales
l_lui = lme(data= data_final, method = "ML",
            control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
            fixed = Phagie ~ LUI, random = ~1|Region,
            correlation = corSpher(form = ~RW + HW|Region, nugget = T))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 250 Spatial scale: ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## SHDI
l_shdi_250 = lme(data= data_final, method = "ML",
                 control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                 fixed = Phagie ~ SHDI250, random = ~1|Region,
                 correlation = corSpher(form = ~RW + HW|Region, nugget = T))

## PAR
l_par_250 = lme(data= data_final, method = "ML",
                control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                fixed = Phagie ~ PAR250, random = ~1|Region,
                correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + SHDI
l_lui_shdi_250 = lme(data= data_final, method = "ML",
                     control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                     fixed = Phagie ~ LUI + SHDI250, random = ~1|Region,
                     correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + PAR
l_lui_par_250 = lme(data= data_final, method = "ML",
                    control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                    fixed = Phagie ~ LUI + PAR250, random = ~1|Region,
                    correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + SHDI + PAR
l_full_250 = lme(data= data_final, method = "ML",
                 control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                 fixed = Phagie ~ LUI + SHDI250 + PAR250, random = ~1|Region,
                 correlation = corSpher(form = ~RW + HW, nugget = T))

## Make a table of AICs for each model at 250m 
l_aic_250 = data.frame(model = c("LUI", "LUI + SHDI", "LUI + PAR", 
                                 "SHDI", "PAR", "LUI + SHDI + PAR"),
                       AIC = c(AICc(l_lui), AICc(l_lui_shdi_250), AICc(l_lui_par_250), 
                               AICc(l_shdi_250), AICc(l_par_250), AICc(l_full_250)),
                       spatial_scale = rep(250, 6))

## Arrange the table by AIC and display it
l_aic_250 = arrange(l_aic_250, AIC)
l_aic_250

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 500 Spatial scale: ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## SHDI
l_shdi_500 = lme(data= data_final, method = "ML",
                 control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                 fixed = Phagie ~ SHDI500, random = ~1|Region,
                 correlation = corSpher(form = ~RW + HW, nugget = T))

## PAR
l_par_500 = lme(data= data_final, method = "ML",
                control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                fixed = Phagie ~ PAR500, random = ~1|Region,
                correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + SHDI
l_lui_shdi_500 = lme(data= data_final, method = "ML",
                     control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                     fixed = Phagie ~ LUI + SHDI500, random = ~1|Region,
                     correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + PAR
l_lui_par_500 = lme(data= data_final, method = "ML",
                    control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                    fixed = Phagie ~ LUI + PAR500, random = ~1|Region,
                    correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + SHDI + PAR
l_full_500 = lme(data= data_final, method = "ML",
                 control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                 fixed = Phagie ~ LUI + SHDI500 + PAR500, random = ~1|Region,
                 correlation = corSpher(form = ~RW + HW, nugget = T))

## Make a table of AICs for each model at 500m 
l_aic_500 = data.frame(model = c("LUI", "LUI + PAR", "LUI + SHDI + PAR", 
                                 "PAR", "SHDI", "LUI + SHDI"),
                       AIC = c(AICc(l_lui), AICc(l_lui_par_500), AICc(l_full_500), 
                               AICc(l_par_500), AICc(l_shdi_500), AICc(l_lui_shdi_500)),
                       spatial_scale = rep(500, 6))

## Arrange the table by AIC and display it
l_aic_500 = arrange(l_aic_500, AIC)
l_aic_500

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 750 Spatial scale: ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## SHDI
l_shdi_750 = lme(data = data_final, method = "ML",
                 control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                 fixed = Phagie ~ SHDI750, random = ~1|Region,
                 correlation = corSpher(form = ~RW + HW, nugget = T))

## PAR
l_par_750 = lme(data = data_final, method = "ML",
                control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                fixed = Phagie ~ PAR750, random = ~1|Region,
                correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + SHDI
l_lui_shdi_750 = lme(data= data_final, method = "ML",
                     control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                     fixed = Phagie ~ LUI + SHDI750, random = ~1|Region,
                     correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + PAR
l_lui_par_750 = lme(data= data_final, method = "ML",
                    control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                    fixed = Phagie ~ LUI + PAR750, random = ~1|Region,
                    correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + SHDI + PAR
l_full_750 = lme(data= data_final, method = "ML",
                 control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                 fixed = Phagie ~ LUI + SHDI750 + PAR750, random = ~1|Region,
                 correlation = corSpher(form = ~RW + HW, nugget = T))

## Make a table of AICs for each model at 750m 
l_aic_750 = data.frame(model = c("LUI", "LUI + SHDI", "LUI + PAR", 
                                 "LUI + SHDI + PAR", "SHDI", "PAR"),
                       AIC = c(AICc(l_lui), AICc(l_lui_shdi_750), AICc(l_lui_par_750), 
                               AICc(l_full_750), AICc(l_shdi_750), AICc(l_par_750)),
                       spatial_scale = rep(750, 6))

## Arrange the table by AIC and display it
l_aic_750 = arrange(l_aic_750, AIC)
l_aic_750

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 1000 Spatial scale: ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## SHDI
l_shdi_1000 = lme(data = data_final, method = "ML",
                  control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                  fixed = Phagie ~ SHDI1000, random = ~1|Region,
                  correlation = corSpher(form = ~RW + HW, nugget = T))

## PAR
l_par_1000 = lme(data = data_final, method = "ML",
                 control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                 fixed = Phagie ~ PAR1000, random = ~1|Region,
                 correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + SHDI
l_lui_shdi_1000 = lme(data= data_final, method = "ML",
                      control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                      fixed = Phagie ~ LUI + SHDI1000, random = ~1|Region,
                      correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + PAR
l_lui_par_1000 = lme(data= data_final, method = "ML",
                     control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                     fixed = Phagie ~ LUI + PAR1000, random = ~1|Region,
                     correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + SHDI + PAR
l_full_1000 = lme(data= data_final, method = "ML",
                  control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                  fixed = Phagie ~ LUI + SHDI1000 + PAR1000, random = ~1|Region,
                  correlation = corSpher(form = ~RW + HW, nugget = T))

## Make a table of AICs for each model at 1,000m 
l_aic_1000 = data.frame(model = c("LUI", "LUI + SHDI", "LUI + PAR", 
                                  "LUI + SHDI + PAR", "SHDI", "PAR"),
                        AIC = c(AICc(l_lui), AICc(l_lui_shdi_1000), AICc(l_lui_par_1000), 
                                AICc(l_full_1000), AICc(l_shdi_1000), AICc(l_par_1000)),
                        spatial_scale = rep(1000, 6))

## Arrange the table by AIC and display it
l_aic_1000 = arrange(l_aic_1000, AIC)
l_aic_1000

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 1250 Spatial scale: ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## SHDI
l_shdi_1250 = lme(data = data_final, method = "ML",
                  control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                  fixed = Phagie ~ SHDI1250, random = ~1|Region,
                  correlation = corSpher(form = ~RW + HW, nugget = T))

## PAR
l_par_1250 = lme(data = data_final, method = "ML",
                 control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                 fixed = Phagie ~ PAR1250, random = ~1|Region,
                 correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + SHDI
l_lui_shdi_1250 = lme(data= data_final, method = "ML",
                      control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                      fixed = Phagie ~ LUI + SHDI1250, random = ~1|Region,
                      correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + PAR
l_lui_par_1250 = lme(data= data_final, method = "ML",
                     control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                     fixed = Phagie ~ LUI + PAR1250, random = ~1|Region,
                     correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + SHDI + PAR
l_full_1250 = lme(data= data_final, method = "ML",
                  control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                  fixed = Phagie ~ LUI + SHDI1250 + PAR1250, random = ~1|Region,
                  correlation = corSpher(form = ~RW + HW, nugget = T))

## Make a table of AICs for each model at 1,250m 
l_aic_1250 = data.frame(model = c("LUI", "LUI + SHDI", "LUI + PAR", 
                                  "LUI + SHDI + PAR", "SHDI", "PAR"),
                        AIC = c(AICc(l_lui), AICc(l_lui_shdi_1250), AICc(l_lui_par_1250), 
                                AICc(l_full_1250), AICc(l_shdi_1250), AICc(l_par_1250)),
                        spatial_scale = rep(1250, 6))

## Arrange the table by AIC and display it
l_aic_1250 = arrange(l_aic_1250, AIC)
l_aic_1250

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 1500 Spatial scale: ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## SHDI
l_shdi_1500 = lme(data = data_final, method = "ML",
                  control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                  fixed = Phagie ~ SHDI1500, random = ~1|Region,
                  correlation = corSpher(form = ~RW + HW, nugget = T))

## PAR
l_par_1500 = lme(data = data_final, method = "ML",
                 control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                 fixed = Phagie ~ PAR1500, random = ~1|Region,
                 correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + SHDI
l_lui_shdi_1500 = lme(data= data_final, method = "ML",
                      control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                      fixed = Phagie ~ LUI + SHDI1500, random = ~1|Region,
                      correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + PAR
l_lui_par_1500 = lme(data= data_final, method = "ML",
                     control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                     fixed = Phagie ~ LUI + PAR1500, random = ~1|Region,
                     correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + SHDI + PAR
l_full_1500 = lme(data= data_final, method = "ML",
                  control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                  fixed = Phagie ~ LUI + SHDI1500 + PAR1500, random = ~1|Region,
                  correlation = corSpher(form = ~RW + HW, nugget = T))

## Make a table of AICs for each model at 1,500m 
l_aic_1500 = data.frame(model = c("LUI", "LUI + SHDI", "LUI + PAR", 
                                  "LUI + SHDI + PAR", "SHDI", "PAR"),
                        AIC = c(AICc(l_lui), AICc(l_lui_shdi_1500), AICc(l_lui_par_1500), 
                                AICc(l_full_1500), AICc(l_shdi_1500), AICc(l_par_1500)),
                        spatial_scale = rep(1500, 6))

## Arrange the table by AIC and display it
l_aic_1500 = arrange(l_aic_1500, AIC)
l_aic_1500

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 2000 Spatial scale: ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## SHDI
l_shdi_2000 = lme(data = data_final, method = "ML",
                  control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                  fixed = Phagie ~ SHDI2000, random = ~1|Region,
                  correlation = corSpher(form = ~RW + HW, nugget = T))

## PAR
l_par_2000 = lme(data = data_final, method = "ML",
                 control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                 fixed = Phagie ~ PAR2000, random = ~1|Region,
                 correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + SHDI
l_lui_shdi_2000 = lme(data= data_final, method = "ML",
                      control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                      fixed = Phagie ~ LUI + SHDI2000, random = ~1|Region,
                      correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + PAR
l_lui_par_2000 = lme(data= data_final, method = "ML",
                     control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                     fixed = Phagie ~ LUI + PAR2000, random = ~1|Region,
                     correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + SHDI + PAR
l_full_2000 = lme(data= data_final, method = "ML",
                  control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                  fixed = Phagie ~ LUI + SHDI2000 + PAR2000, random = ~1|Region,
                  correlation = corSpher(form = ~RW + HW, nugget = T))

## Make a table of AICs for each model at 2,000m 
l_aic_2000 = data.frame(model = c("LUI", "LUI + SHDI", "LUI + PAR", 
                                  "LUI + SHDI + PAR", "SHDI", "PAR"),
                        AIC = c(AICc(l_lui), AICc(l_lui_shdi_2000), AICc(l_lui_par_2000), 
                                AICc(l_full_2000), AICc(l_shdi_2000), AICc(l_par_2000)),
                        spatial_scale = rep(2000, 6))

## Arrange the table by AIC and display it
l_aic_2000 = arrange(l_aic_2000, AIC)
l_aic_2000

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Compare models at all spatial scales: ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Combine all AIC tables into one
l_aic_all = rbind(l_aic_250, l_aic_500, l_aic_750, l_aic_1000,
                  l_aic_1250, l_aic_1500, l_aic_2000)

## Add a new column: the difference in AICc from the null AICc
l_aic_all = mutate(l_aic_all, diff_AIC = AIC - AICc(l_lui))

## Arrange the table by AIC and display it
l_aic_all = arrange(l_aic_all, diff_AIC)
l_aic_all


#######################################################
#### Mixed models with forewing length as response ####
#######################################################

## Create the null, LUI-only model to use at all spatial scales
w_lui = lme(data= data_final, method = "ML",
            control=list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
            fixed = Size ~ LUI, random = ~1|Region,
            correlation = corSpher(form = ~RW + HW|Region, nugget = T))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 250 Spatial scale: ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## SHDI
w_shdi_250 = lme(data= data_final, method = "ML",
                 control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                 fixed = Size ~ SHDI250, random = ~1|Region,
                 correlation = corSpher(form = ~RW + HW|Region, nugget = T))

## PAR
w_par_250 = lme(data= data_final, method = "ML",
                control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                fixed = Size ~ PAR250, random = ~1|Region,
                correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + SHDI
w_lui_shdi_250 = lme(data= data_final, method = "ML",
                     control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                     fixed = Size ~ LUI + SHDI250, random = ~1|Region,
                     correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + PAR
w_lui_par_250 = lme(data= data_final, method = "ML",
                    control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                    fixed = Size ~ LUI + PAR250, random = ~1|Region,
                    correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + SHDI + PAR
w_full_250 = lme(data= data_final, method = "ML",
                 control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                 fixed = Size ~ LUI + SHDI250 + PAR250, random = ~1|Region,
                 correlation = corSpher(form = ~RW + HW, nugget = T))

## Make a table of AICs for each model at 250m 
w_aic_250 = data.frame(model = c("LUI", "LUI + SHDI", "LUI + PAR", 
                                 "SHDI", "PAR", "LUI + SHDI + PAR"),
                       AIC = c(AICc(w_lui), AICc(w_lui_shdi_250), AICc(w_lui_par_250), 
                               AICc(w_shdi_250), AICc(w_par_250), AICc(w_full_250)),
                       spatial_scale = rep(250, 6))

## Arrange the table by AIC and display it
w_aic_250 = arrange(w_aic_250, AIC)
w_aic_250

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 500 Spatial scale: ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## SHDI
w_shdi_500 = lme(data= data_final, method = "ML",
                 control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                 fixed = Size ~ SHDI500, random = ~1|Region,
                 correlation = corSpher(form = ~RW + HW, nugget = T))

## PAR
w_par_500 = lme(data= data_final, method = "ML",
                control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                fixed = Size ~ PAR500, random = ~1|Region,
                correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + SHDI
w_lui_shdi_500 = lme(data= data_final, method = "ML",
                     control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                     fixed = Size ~ LUI + SHDI500, random = ~1|Region,
                     correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + PAR
w_lui_par_500 = lme(data= data_final, method = "ML",
                    control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                    fixed = Size ~ LUI + PAR500, random = ~1|Region,
                    correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + SHDI + PAR
w_full_500 = lme(data= data_final, method = "ML",
                 control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                 fixed = Size ~ LUI + SHDI500 + PAR500, random = ~1|Region,
                 correlation = corSpher(form = ~RW + HW, nugget = T))

## Make a table of AICs for each model at 500m 
w_aic_500 = data.frame(model = c("LUI", "LUI + PAR", "LUI + SHDI + PAR", 
                                 "PAR", "SHDI", "LUI + SHDI"),
                       AIC = c(AICc(w_lui), AICc(w_lui_par_500), AICc(w_full_500), 
                               AICc(w_par_500), AICc(w_shdi_500), AICc(w_lui_shdi_500)),
                       spatial_scale = rep(500, 6))

## Arrange the table by AIC and display it
w_aic_500 = arrange(w_aic_500, AIC)
w_aic_500

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 750 Spatial scale: ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## SHDI
w_shdi_750 = lme(data= data_final, method = "ML",
                 control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                 fixed = Size ~ SHDI750, random = ~1|Region,
                 correlation = corSpher(form = ~RW + HW, nugget = T))

## PAR
w_par_750 = lme(data= data_final, method = "ML",
                control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                fixed = Size ~ PAR750, random = ~1|Region,
                correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + SHDI
w_lui_shdi_750 = lme(data= data_final, method = "ML",
                     control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                     fixed = Size ~ LUI + SHDI750, random = ~1|Region,
                     correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + PAR
w_lui_par_750 = lme(data= data_final, method = "ML",
                    control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                    fixed = Size ~ LUI + PAR750, random = ~1|Region,
                    correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + SHDI + PAR
w_full_750 = lme(data= data_final, method = "ML",
                 control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                 fixed = Size ~ LUI + SHDI750 + PAR750, random = ~1|Region,
                 correlation = corSpher(form = ~RW + HW, nugget = T))

## Make a table of AICs for each model at 500m 
w_aic_750 = data.frame(model = c("LUI", "LUI + PAR", "LUI + SHDI + PAR", 
                                 "PAR", "SHDI", "LUI + SHDI"),
                       AIC = c(AICc(w_lui), AICc(w_lui_par_750), AICc(w_full_750), 
                               AICc(w_par_750), AICc(w_shdi_750), AICc(w_lui_shdi_750)),
                       spatial_scale = rep(750, 6))

## Arrange the table by AIC and display it
w_aic_750 = arrange(w_aic_750, AIC)
w_aic_750

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 1000 Spatial scale: ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## SHDI
w_shdi_1000 = lme(data= data_final, method = "ML",
                 control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                 fixed = Size ~ SHDI1000, random = ~1|Region,
                 correlation = corSpher(form = ~RW + HW, nugget = T))

## PAR
w_par_1000 = lme(data= data_final, method = "ML",
                control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                fixed = Size ~ PAR1000, random = ~1|Region,
                correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + SHDI
w_lui_shdi_1000 = lme(data= data_final, method = "ML",
                     control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                     fixed = Size ~ LUI + SHDI1000, random = ~1|Region,
                     correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + PAR
w_lui_par_1000 = lme(data= data_final, method = "ML",
                    control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                    fixed = Size ~ LUI + PAR1000, random = ~1|Region,
                    correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + SHDI + PAR
w_full_1000 = lme(data= data_final, method = "ML",
                 control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                 fixed = Size ~ LUI + SHDI1000 + PAR1000, random = ~1|Region,
                 correlation = corSpher(form = ~RW + HW, nugget = T))

## Make a table of AICs for each model at 500m 
w_aic_1000 = data.frame(model = c("LUI", "LUI + PAR", "LUI + SHDI + PAR", 
                                 "PAR", "SHDI", "LUI + SHDI"),
                       AIC = c(AICc(w_lui), AICc(w_lui_par_1000), AICc(w_full_1000), 
                               AICc(w_par_1000), AICc(w_shdi_1000), AICc(w_lui_shdi_1000)),
                       spatial_scale = rep(1000, 6))

## Arrange the table by AIC and display it
w_aic_1000 = arrange(w_aic_1000, AIC)
w_aic_1000

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 1250 Spatial scale: ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## SHDI
w_shdi_1250 = lme(data= data_final, method = "ML",
                  control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                  fixed = Size ~ SHDI1250, random = ~1|Region,
                  correlation = corSpher(form = ~RW + HW, nugget = T))

## PAR
w_par_1250 = lme(data= data_final, method = "ML",
                 control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                 fixed = Size ~ PAR1250, random = ~1|Region,
                 correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + SHDI
w_lui_shdi_1250 = lme(data= data_final, method = "ML",
                      control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                      fixed = Size ~ LUI + SHDI1250, random = ~1|Region,
                      correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + PAR
w_lui_par_1250 = lme(data= data_final, method = "ML",
                     control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                     fixed = Size ~ LUI + PAR1250, random = ~1|Region,
                     correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + SHDI + PAR
w_full_1250 = lme(data= data_final, method = "ML",
                  control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                  fixed = Size ~ LUI + SHDI1250 + PAR1250, random = ~1|Region,
                  correlation = corSpher(form = ~RW + HW, nugget = T))

## Make a table of AICs for each model at 500m 
w_aic_1250 = data.frame(model = c("LUI", "LUI + PAR", "LUI + SHDI + PAR", 
                                  "PAR", "SHDI", "LUI + SHDI"),
                        AIC = c(AICc(w_lui), AICc(w_lui_par_1250), AICc(w_full_1250), 
                                AICc(w_par_1250), AICc(w_shdi_1250), AICc(w_lui_shdi_1250)),
                        spatial_scale = rep(1250, 6))

## Arrange the table by AIC and display it
w_aic_1250 = arrange(w_aic_1250, AIC)
w_aic_1250

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 1500 Spatial scale: ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## SHDI
w_shdi_1500 = lme(data= data_final, method = "ML",
                  control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                  fixed = Size ~ SHDI1500, random = ~1|Region,
                  correlation = corSpher(form = ~RW + HW, nugget = T))

## PAR
w_par_1500 = lme(data= data_final, method = "ML",
                 control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                 fixed = Size ~ PAR1500, random = ~1|Region,
                 correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + SHDI
w_lui_shdi_1500 = lme(data= data_final, method = "ML",
                      control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                      fixed = Size ~ LUI + SHDI1500, random = ~1|Region,
                      correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + PAR
w_lui_par_1500 = lme(data= data_final, method = "ML",
                     control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                     fixed = Size ~ LUI + PAR1500, random = ~1|Region,
                     correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + SHDI + PAR
w_full_1500 = lme(data= data_final, method = "ML",
                  control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                  fixed = Size ~ LUI + SHDI1500 + PAR1500, random = ~1|Region,
                  correlation = corSpher(form = ~RW + HW, nugget = T))

## Make a table of AICs for each model at 500m 
w_aic_1500 = data.frame(model = c("LUI", "LUI + PAR", "LUI + SHDI + PAR", 
                                  "PAR", "SHDI", "LUI + SHDI"),
                        AIC = c(AICc(w_lui), AICc(w_lui_par_1500), AICc(w_full_1500), 
                                AICc(w_par_1500), AICc(w_shdi_1500), AICc(w_lui_shdi_1500)),
                        spatial_scale = rep(1500, 6))

## Arrange the table by AIC and display it
w_aic_1500 = arrange(w_aic_1500, AIC)
w_aic_1500

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 2000 Spatial scale: ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## SHDI
w_shdi_2000 = lme(data= data_final, method = "ML",
                  control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                  fixed = Size ~ SHDI2000, random = ~1|Region,
                  correlation = corSpher(form = ~RW + HW, nugget = T))

## PAR
w_par_2000 = lme(data= data_final, method = "ML",
                 control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                 fixed = Size ~ PAR2000, random = ~1|Region,
                 correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + SHDI
w_lui_shdi_2000 = lme(data= data_final, method = "ML",
                      control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                      fixed = Size ~ LUI + SHDI2000, random = ~1|Region,
                      correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + PAR
w_lui_par_2000 = lme(data= data_final, method = "ML",
                     control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                     fixed = Size ~ LUI + PAR2000, random = ~1|Region,
                     correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + SHDI + PAR
w_full_2000 = lme(data= data_final, method = "ML",
                  control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                  fixed = Size ~ LUI + SHDI2000 + PAR2000, random = ~1|Region,
                  correlation = corSpher(form = ~RW + HW, nugget = T))

## Make a table of AICs for each model at 500m 
w_aic_2000 = data.frame(model = c("LUI", "LUI + PAR", "LUI + SHDI + PAR", 
                                  "PAR", "SHDI", "LUI + SHDI"),
                        AIC = c(AICc(w_lui), AICc(w_lui_par_2000), AICc(w_full_2000), 
                                AICc(w_par_2000), AICc(w_shdi_2000), AICc(w_lui_shdi_2000)),
                        spatial_scale = rep(2000, 6))

## Arrange the table by AIC and display it
w_aic_2000 = arrange(w_aic_2000, AIC)
w_aic_2000

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Compare models at all spatial scales: ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Combine all AIC tables into one
w_aic_all = rbind(w_aic_250, w_aic_500, w_aic_750, w_aic_1000,
                  w_aic_1250, w_aic_1500, w_aic_2000)

## Add a new column: the difference in AICc from the null AICc
w_aic_all = mutate(w_aic_all, diff_AIC = AIC - AICc(w_lui))

## Arrange the table by AIC and display it
w_aic_all = arrange(w_aic_all, diff_AIC)
w_aic_all

##########################################################
#### Mixed models with migratory tendency as response ####
##########################################################

## Create the null, LUI-only model to use at all spatial scales
m_lui = lme(data= data_final, method = "ML",
            control=list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
            fixed = Migration ~ LUI, random = ~1|Region,
            correlation = corSpher(form = ~RW + HW|Region, nugget = T))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 250 Spatial scale: ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## SHDI
m_shdi_250 = lme(data= data_final, method = "ML",
                 control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                 fixed = Migration ~ SHDI250, random = ~1|Region,
                 correlation = corSpher(form = ~RW + HW|Region, nugget = T))

## PAR
m_par_250 = lme(data= data_final, method = "ML",
                control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                fixed = Migration ~ PAR250, random = ~1|Region,
                correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + SHDI
m_lui_shdi_250 = lme(data= data_final, method = "ML",
                     control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                     fixed = Migration ~ LUI + SHDI250, random = ~1|Region,
                     correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + PAR
m_lui_par_250 = lme(data= data_final, method = "ML",
                    control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                    fixed = Migration ~ LUI + PAR250, random = ~1|Region,
                    correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + SHDI + PAR
m_full_250 = lme(data= data_final, method = "ML",
                 control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                 fixed = Migration ~ LUI + SHDI250 + PAR250, random = ~1|Region,
                 correlation = corSpher(form = ~RW + HW, nugget = T))

## Make a table of AICs for each model at 250m 
m_aic_250 = data.frame(model = c("LUI", "LUI + SHDI", "LUI + PAR", 
                                 "SHDI", "PAR", "LUI + SHDI + PAR"),
                       AIC = c(AICc(m_lui), AICc(m_lui_shdi_250), AICc(m_lui_par_250), 
                               AICc(m_shdi_250), AICc(m_par_250), AICc(m_full_250)),
                       spatial_scale = rep(250, 6))

## Arrange the table by AIC and display it
m_aic_250 = arrange(m_aic_250, AIC)
m_aic_250

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 500 Spatial scale: ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## SHDI
m_shdi_500 = lme(data= data_final, method = "ML",
                 control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                 fixed = Migration ~ SHDI500, random = ~1|Region,
                 correlation = corSpher(form = ~RW + HW|Region, nugget = T))

## PAR
m_par_500 = lme(data= data_final, method = "ML",
                control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                fixed = Migration ~ PAR500, random = ~1|Region,
                correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + SHDI
m_lui_shdi_500 = lme(data= data_final, method = "ML",
                     control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                     fixed = Migration ~ LUI + SHDI500, random = ~1|Region,
                     correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + PAR
m_lui_par_500 = lme(data= data_final, method = "ML",
                    control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                    fixed = Migration ~ LUI + PAR500, random = ~1|Region,
                    correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + SHDI + PAR
m_full_500 = lme(data= data_final, method = "ML",
                 control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                 fixed = Migration ~ LUI + SHDI500 + PAR500, random = ~1|Region,
                 correlation = corSpher(form = ~RW + HW, nugget = T))

## Make a table of AICs for each model at 250m 
m_aic_500 = data.frame(model = c("LUI", "LUI + SHDI", "LUI + PAR", 
                                 "SHDI", "PAR", "LUI + SHDI + PAR"),
                       AIC = c(AICc(m_lui), AICc(m_lui_shdi_500), AICc(m_lui_par_500), 
                               AICc(m_shdi_500), AICc(m_par_500), AICc(m_full_500)),
                       spatial_scale = rep(500, 6))

## Arrange the table by AIC and display it
m_aic_500 = arrange(m_aic_500, AIC)
m_aic_500

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 750 Spatial scale: ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## SHDI
m_shdi_750 = lme(data= data_final, method = "ML",
                 control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                 fixed = Migration ~ SHDI750, random = ~1|Region,
                 correlation = corSpher(form = ~RW + HW|Region, nugget = T))

## PAR
m_par_750 = lme(data= data_final, method = "ML",
                control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                fixed = Migration ~ PAR750, random = ~1|Region,
                correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + SHDI
m_lui_shdi_750 = lme(data= data_final, method = "ML",
                     control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                     fixed = Migration ~ LUI + SHDI750, random = ~1|Region,
                     correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + PAR
m_lui_par_750 = lme(data= data_final, method = "ML",
                    control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                    fixed = Migration ~ LUI + PAR750, random = ~1|Region,
                    correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + SHDI + PAR
m_full_750 = lme(data= data_final, method = "ML",
                 control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                 fixed = Migration ~ LUI + SHDI750 + PAR750, random = ~1|Region,
                 correlation = corSpher(form = ~RW + HW, nugget = T))

## Make a table of AICs for each model at 250m 
m_aic_750 = data.frame(model = c("LUI", "LUI + SHDI", "LUI + PAR", 
                                 "SHDI", "PAR", "LUI + SHDI + PAR"),
                       AIC = c(AICc(m_lui), AICc(m_lui_shdi_750), AICc(m_lui_par_750), 
                               AICc(m_shdi_750), AICc(m_par_750), AICc(m_full_750)),
                       spatial_scale = rep(750, 6))

## Arrange the table by AIC and display it
m_aic_750 = arrange(m_aic_750, AIC)
m_aic_750

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 1000 Spatial scale: ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## SHDI
m_shdi_1000 = lme(data= data_final, method = "ML",
                 control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                 fixed = Migration ~ SHDI1000, random = ~1|Region,
                 correlation = corSpher(form = ~RW + HW|Region, nugget = T))

## PAR
m_par_1000 = lme(data= data_final, method = "ML",
                control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                fixed = Migration ~ PAR1000, random = ~1|Region,
                correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + SHDI
m_lui_shdi_1000 = lme(data= data_final, method = "ML",
                     control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                     fixed = Migration ~ LUI + SHDI1000, random = ~1|Region,
                     correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + PAR
m_lui_par_1000 = lme(data= data_final, method = "ML",
                    control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                    fixed = Migration ~ LUI + PAR1000, random = ~1|Region,
                    correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + SHDI + PAR
m_full_1000 = lme(data= data_final, method = "ML",
                 control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                 fixed = Migration ~ LUI + SHDI1000 + PAR1000, random = ~1|Region,
                 correlation = corSpher(form = ~RW + HW, nugget = T))

## Make a table of AICs for each model at 250m 
m_aic_1000 = data.frame(model = c("LUI", "LUI + SHDI", "LUI + PAR", 
                                 "SHDI", "PAR", "LUI + SHDI + PAR"),
                       AIC = c(AICc(m_lui), AICc(m_lui_shdi_1000), AICc(m_lui_par_1000), 
                               AICc(m_shdi_1000), AICc(m_par_1000), AICc(m_full_1000)),
                       spatial_scale = rep(1000, 6))

## Arrange the table by AIC and display it
m_aic_1000 = arrange(m_aic_1000, AIC)
m_aic_1000

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 1250 Spatial scale: ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## SHDI
m_shdi_1250 = lme(data= data_final, method = "ML",
                  control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                  fixed = Migration ~ SHDI1250, random = ~1|Region,
                  correlation = corSpher(form = ~RW + HW|Region, nugget = T))

## PAR
m_par_1250 = lme(data= data_final, method = "ML",
                 control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                 fixed = Migration ~ PAR1250, random = ~1|Region,
                 correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + SHDI
m_lui_shdi_1250 = lme(data= data_final, method = "ML",
                      control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                      fixed = Migration ~ LUI + SHDI1250, random = ~1|Region,
                      correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + PAR
m_lui_par_1250 = lme(data= data_final, method = "ML",
                     control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                     fixed = Migration ~ LUI + PAR1250, random = ~1|Region,
                     correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + SHDI + PAR
m_full_1250 = lme(data= data_final, method = "ML",
                  control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                  fixed = Migration ~ LUI + SHDI1250 + PAR1250, random = ~1|Region,
                  correlation = corSpher(form = ~RW + HW, nugget = T))

## Make a table of AICs for each model at 250m 
m_aic_1250 = data.frame(model = c("LUI", "LUI + SHDI", "LUI + PAR", 
                                  "SHDI", "PAR", "LUI + SHDI + PAR"),
                        AIC = c(AICc(m_lui), AICc(m_lui_shdi_1250), AICc(m_lui_par_1250), 
                                AICc(m_shdi_1250), AICc(m_par_1250), AICc(m_full_1250)),
                        spatial_scale = rep(1250, 6))

## Arrange the table by AIC and display it
m_aic_1250 = arrange(m_aic_1250, AIC)
m_aic_1250

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 1500 Spatial scale: ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## SHDI
m_shdi_1500 = lme(data= data_final, method = "ML",
                  control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                  fixed = Migration ~ SHDI1500, random = ~1|Region,
                  correlation = corSpher(form = ~RW + HW|Region, nugget = T))

## PAR
m_par_1500 = lme(data= data_final, method = "ML",
                 control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                 fixed = Migration ~ PAR1500, random = ~1|Region,
                 correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + SHDI
m_lui_shdi_1500 = lme(data= data_final, method = "ML",
                      control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                      fixed = Migration ~ LUI + SHDI1500, random = ~1|Region,
                      correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + PAR
m_lui_par_1500 = lme(data= data_final, method = "ML",
                     control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                     fixed = Migration ~ LUI + PAR1500, random = ~1|Region,
                     correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + SHDI + PAR
m_full_1500 = lme(data= data_final, method = "ML",
                  control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                  fixed = Migration ~ LUI + SHDI1500 + PAR1500, random = ~1|Region,
                  correlation = corSpher(form = ~RW + HW, nugget = T))

## Make a table of AICs for each model at 250m 
m_aic_1500 = data.frame(model = c("LUI", "LUI + SHDI", "LUI + PAR", 
                                  "SHDI", "PAR", "LUI + SHDI + PAR"),
                        AIC = c(AICc(m_lui), AICc(m_lui_shdi_1500), AICc(m_lui_par_1500), 
                                AICc(m_shdi_1500), AICc(m_par_1500), AICc(m_full_1500)),
                        spatial_scale = rep(1500, 6))

## Arrange the table by AIC and display it
m_aic_1500 = arrange(m_aic_1500, AIC)
m_aic_1500

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 2000 Spatial scale: ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## SHDI
m_shdi_2000 = lme(data= data_final, method = "ML",
                  control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                  fixed = Migration ~ SHDI2000, random = ~1|Region,
                  correlation = corSpher(form = ~RW + HW|Region, nugget = T))

## PAR
m_par_2000 = lme(data= data_final, method = "ML",
                 control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                 fixed = Migration ~ PAR2000, random = ~1|Region,
                 correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + SHDI
m_lui_shdi_2000 = lme(data= data_final, method = "ML",
                      control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                      fixed = Migration ~ LUI + SHDI2000, random = ~1|Region,
                      correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + PAR
m_lui_par_2000 = lme(data= data_final, method = "ML",
                     control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                     fixed = Migration ~ LUI + PAR2000, random = ~1|Region,
                     correlation = corSpher(form = ~RW + HW, nugget = T))

## LUI + SHDI + PAR
m_full_2000 = lme(data= data_final, method = "ML",
                  control = list(opt = "optim", msMaxIter = 1000, returnObject = TRUE),
                  fixed = Migration ~ LUI + SHDI2000 + PAR1000, random = ~1|Region,
                  correlation = corSpher(form = ~RW + HW, nugget = T))

## Make a table of AICs for each model at 250m 
m_aic_2000 = data.frame(model = c("LUI", "LUI + SHDI", "LUI + PAR", 
                                  "SHDI", "PAR", "LUI + SHDI + PAR"),
                        AIC = c(AICc(m_lui), AICc(m_lui_shdi_2000), AICc(m_lui_par_2000), 
                                AICc(m_shdi_2000), AICc(m_par_2000), AICc(m_full_2000)),
                        spatial_scale = rep(2000, 6))

## Arrange the table by AIC and display it
m_aic_2000 = arrange(m_aic_2000, AIC)
m_aic_2000

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Compare models at all spatial scales: ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Combine all AIC tables into one
m_aic_all = rbind(m_aic_250, m_aic_500, m_aic_750, m_aic_1000,
                     m_aic_1250, m_aic_1500, m_aic_2000)

## Add a new column: the difference in AICc from the null AICc
m_aic_all = mutate(m_aic_all, diff_AIC = AIC - AICc(m_lui))


## Arrange the table by AIC and display it
m_aic_all = arrange(m_aic_all, AIC)
m_aic_all

#################################
#### Recreating Figure 2 (a) ####
#################################

## First, change data_final into a data frame and save it as a new variable
data_plots = as.data.frame(data_final)

## Create a new data frame of values with which to predict model values
pred_frame_a = data.frame(SHDI250 = seq(min(data_plots$SHDI250), max(data_plots$SHDI250), 
                                      length.out = nrow(data_plots)),
                        LUI = rep(mean(data_plots$LUI), nrow(data_plots)),
                        Region = rep("A", nrow(data_plots)))

## Pass in data
fig2_a = ggplot(data_plots,
                aes(x = SHDI250,
                    y = shannon_H)) +
  
  ## Add a line of predicted points in red
  geom_line(data = pred_frame, aes(y = predict(h_lui_shdi_250, pred_frame_a)),
            color = "red", size = 1.5) +
  
  ## Plot the individual values of Shannon's H
  geom_point(alpha = 0.5, size = 1.5) +
  
  labs(x = "SHDI", y = "Taxonomic Diversity (H`)") +
  theme_minimal()

#################################
#### Recreating Figure 2 (b) ####
#################################

## Create a new data frame of values with which to predict model values
pred_frame_b = data.frame(PAR750 = seq(min(data_plots$PAR750), max(data_plots$PAR750), 
                                        length.out = nrow(data_plots)),
                          Region = rep("A", nrow(data_plots)))

## Pass in data
fig2_b = ggplot(data_plots,
                aes(x = PAR750,
                    y = Phagie)) +
  
  ## Add a line of predicted points in red
  geom_line(data = pred_frame_b, aes(y = predict(l_par_750, pred_frame_b)),
            color = "red", size = 1.5) +
  
  ## Plot the individual values of Larval Feeding
  geom_point(alpha = 0.5) +
  
  labs(x = "PAR (log-transform)", y = "Larval Feeding") +
  theme_minimal()

#################################
#### Recreating Figure 2 (c) ####
#################################

## Create a new data frame of values with which to predict model values
pred_frame_c = data.frame(PAR250 = seq(min(data_plots$PAR250), max(data_plots$PAR250), 
                                       length.out = nrow(data_plots)),
                          Region = rep("A", nrow(data_plots)))

## Pass in data
fig2_c = ggplot(data_plots,
                aes(x = PAR250,
                    y = Size)) +
  
  ## Add a line of predicted points in red
  geom_line(data = pred_frame_c, aes(y = predict(w_par_250, pred_frame_c)),
            color = "red", size = 1.5) +
  
  ## Plot the individual values of Forewing length
  geom_point(alpha = 0.5) +
  
  ## Add labels and theme
  labs(x = "PAR (log-transform)", y = "Forewing Length") +
  theme_minimal()

#################################
#### Recreating Figure 2 (d) ####
#################################

## Create a new data frame of values with which to predict model values
pred_frame_d = data.frame(PAR250 = seq(min(data_plots$PAR250), max(data_plots$PAR250), 
                                       length.out = nrow(data_plots)),
                          LUI = rep(mean(data_plots$LUI), nrow(data_plots)),
                          Region = rep("A", nrow(data_plots)))

## Pass in data
fig2_d = ggplot(data_plots,
                aes(x = PAR250,
                    y = Migration)) +
  
  ## Add a line of predicted points in red
  geom_line(data = pred_frame_d, aes(y = predict(m_lui_par_250, pred_frame_d)),
            color = "red", size = 1.5) +
  
  ## Plot the individual values of Shannon's H
  geom_point(alpha = 0.5) +
  
  labs(x = "PAR (log-transform)", y = "Migratory Tendency") +
  theme_minimal()
