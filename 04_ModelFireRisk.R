#----------------------------------------------------------------------------------------------------#
# Name:     03_ModelFireRisk                                                                         #
# Purpose:                                                                                           #
# Created:  2019-06-04                                                                               #
# Updated:  2019-06-06                                                                               #
# Author:   Sandro Groth                                                                             #
# Based on: [March 2016](www.animove.org, www.ecosens.org), B. Reineking, M. Wegmann:                #
#           Species distribution models (SDM) using buffalo data from GBIF                           #
#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#
# SETTINGS

study_area_path <- "C:\\Users\\sandr\\Documents\\EAGLE_Data\\EuropeanStates_Clip.shp"
burn_occ_path <- "C:\\Users\\sandr\\Documents\\EAGLE_Data\\Burn_Acc_Vec.shp"

windspeed_path <- "C:\\Users\\sandr\\Documents\\EAGLE_Data\\SS2019_2nd_Term\\MET1_Spatial_modeling_and_prediction\\Final_Project\\R\\FireriskSDM\\wc2.0_2.5m_wind"
glcc_landcov_path <- "C:\\Users\\sandr\\Documents\\EAGLE_Data\\MODIS_GLCC\\GLCC_Majority_Type_1.tif"

#----------------------------------------------------------------------------------------------------#
# IMPORTS

library(rms) 
library(raster)
library(mgcv)
library(randomForest)
library(dismo)
library(rgdal)
library(ellipse)
library(randomForest)
library(rJava)
library(XML)
source("varImpBiomod.R")

#----------------------------------------------------------------------------------------------------#

study_area <- readOGR(study_area_path)
burn_occ <- readOGR(burn_occ_path)
glcc_landcov <- raster(glcc_landcov_path)
windspeed_layer_paths <- list.files(windspeed_path, pattern = ".tif", full.names = T, no.. = T)

# create mean windspeed layer
wind_jan <- raster(windspeed_layer_paths[1])
wind_feb <- raster(windspeed_layer_paths[2])
wind_mar <- raster(windspeed_layer_paths[3])
wind_apr <- raster(windspeed_layer_paths[4])
wind_may <- raster(windspeed_layer_paths[5])
wind_jun <- raster(windspeed_layer_paths[6])
wind_jul <- raster(windspeed_layer_paths[7])
wind_aug <- raster(windspeed_layer_paths[8])
wind_sep <- raster(windspeed_layer_paths[9])
wind_oct <- raster(windspeed_layer_paths[10])
wind_nov <- raster(windspeed_layer_paths[11])
wind_dec <- raster(windspeed_layer_paths[12])
wind_avg <- (wind_jan + wind_feb + wind_mar + wind_apr + wind_may + wind_jun +
             wind_jul + wind_aug + wind_sep + wind_oct + wind_nov + wind_dec) / 12

bio <- raster::getData('worldclim', var = "bio", res = 2.5)
bio_2050 <- raster::getData('CMIP5', var = "bio", model = 'AC', rcp = 45, year = 50,  res = 2.5)
bio_2070 <- raster::getData('CMIP5', var = "bio", model = 'AC', rcp = 45, year = 70,  res = 2.5)

biocrop <- crop(bio, extent(study_area) + 10)
biocrop_2050 <- crop(bio_2050, extent(study_area) + 10)
biocrop_2070 <- crop(bio_2070, extent(study_area) + 10)
wind_crop <- crop(wind_avg, extent(study_area) + 10)
glcc_crop <- crop(glcc_landcov, extent(study_area) + 10)

biocrop <- stack(biocrop, glcc_crop, wind_crop)
biocrop_2050 <- stack(biocrop_2050, glcc_crop, wind_crop)
biocrop_2070 <- stack(biocrop_2070, glcc_crop, wind_crop)

names(biocrop_2050) <- names(biocrop)
names(biocrop_2070) <- names(biocrop)

########################################
#' Collinearity
#' -----------------------------
#' ### Visual inspection of collinearity ###
cm <- cor(getValues(biocrop), use = "complete.obs")
plotcorr(cm, col=ifelse(abs(cm) > 0.7, "red", "grey"))
########################################

#' Select species records for which environmental information is available
#' -------------------------------
burn_occ <- burn_occ[complete.cases(extract(biocrop, burn_occ)), ]

#' ### Select an uncorrelated subset of environmental variables ###
env <- subset(biocrop, c( "bio4", "bio9", "bio14", 
                          "bio15", "bio17", "bio18", "GLCC_Majority_Type_1", "layer"))
env_2050 <- subset(biocrop_2050, c("bio4", "bio9", "bio14", 
                   "bio15", "bio17", "bio18", "GLCC_Majority_Type_1", "layer"))
env_2070 <- subset(biocrop_2070, c( "bio4", "bio9", "bio14", 
                                    "bio15", "bio17", "bio18", "GLCC_Majority_Type_1", "layer"))

#' Selecting 2000 random background points, excluding cells where
#' the species is present
set.seed(2)
background <- randomPoints(env, 2000, burn_occ)
#' Select only one presence record in each cell of the environmental layer
presence <- gridSample(burn_occ, env, n = 1)

#' 
#' Now we combine the presence and background points, adding a 
#' column "species" that contains the information about presence (1)
#' and background (0)
fulldata <- SpatialPointsDataFrame(rbind(presence, background),
                                   data = data.frame("burn_occ" = rep(c(1,0), 
                                                                     c(nrow(presence), nrow(background)))),
                                   match.ID = FALSE,
                                   proj4string = CRS(projection(env)))


#' 
# Split data set into a training and test data set
set.seed(2)
fold <- kfold(fulldata, k = 5)
traindata <- fulldata[fold != 1, ]
testdata <- fulldata[fold == 1, ]

#' We can now use a range of statistical methods to estimate the
#' probability of species occurrence.
#' Unfortunately, there are often subtle differences in how the models
#' are specified and in which data formats are useable

varnames <- c("bio4", "bio9", "bio14", 
              "bio15", "bio17", "bio18", "GLCC_Majority_Type_1", "layer")

gammodel <- gam(burn_occ ~ s(bio4) + s(bio9) + s(bio14) + 
                s(bio15) + s(bio17) + s(bio18) + s(GLCC_Majority_Type_1) + s(layer),
                family="binomial", data=traindata)
summary(gammodel)

plot(gammodel[6])


gamtest <- predict(gammodel, newdata = testdata, type = "response")
# b) Calculate performance indices
val.prob(gamtest, testdata[["burn_occ"]])

# Variable importance
gamimp <- varImpBiomod(gammodel, varnames,
                       traindata)
barplot(100 * gamimp/sum(gamimp), ylab = "Variable importance (%)")

# Response functions
plot(gammodel, pages = 1)

# png("gammodel_resp.png", 800, 800)
# plot(gammodel, pages = 1)
# dev.off()

# Prediction map
gammap <- predict(env, gammodel, type = "response")
gammap_2050 <- predict(env_2050, gammodel, type = "response")
gammap_2070 <- predict(env_2070, gammodel, type = "response")

par(mfrow=c(2,2))
plot(gammap)
plot(gammap_2050)
plot(gammap_2070)
