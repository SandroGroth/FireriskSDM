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

study_area_path <- "C:\\Users\\sandr\\Documents\\EAGLE_Data\\EuropeanStates.shp"
burn_occ_path <- "C:\\Users\\sandr\\Documents\\EAGLE_Data\\Burn_Acc_Vec.shp"

corine_landcov <- "C:\\Users\\sandr\\Documents\\EAGLE_Data\\CORINE\\CLC2018_CLC2018_V2018_20b2.tif"

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

bio <- raster::getData('worldclim', var = "bio", res = 2.5)

biocrop <- crop(bio, extent(study_area) + 10)

plot(burn_occ, add = TRUE)

########################################
#' Collinearity
#' -----------------------------
#' ### Visual inspection of collinearity ###
cm <- cor(getValues(bio), use = "complete.obs")
plotcorr(cm, col=ifelse(abs(cm) > 0.7, "red", "grey"))
########################################

#' Select species records for which environmental information is available
#' -------------------------------
burn_occ <- burn_occ[complete.cases(extract(biocrop, burn_occ)), ]

#' ### Select an uncorrelated subset of environmental variables ###
env <- subset(biocrop, c( "bio3", "bio5", "bio14", "bio15"))

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
#' Add information of environmental conditions at point locations
fulldata@data <- cbind(fulldata@data, extract(env, fulldata))

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

varnames <- c("bio3", "bio5", "bio14", "bio15")

gammodel <- gam(burn_occ ~ s(bio3) + s(bio5) + s(bio14)+ s(bio15),
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

plot(gammap)

