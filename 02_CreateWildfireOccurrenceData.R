# 02_CreateWildfireOccurenceData

library(getSpatialData)
library(raster)
library(rgdal)
library(sp)

burn_acc_path <- "C:\\Users\\sandr\\Documents\\EAGLE_Data\\Burn_Acc_Reproj.tif"

burn_acc <- raster(burn_acc_path)
study_area <- readOGR("C:\\Users\\sandr\\Documents\\EAGLE_Data\\SS2019_2nd_Term\\MET1_Spatial_modeling_and_prediction\\Final_Project\\R\\FireriskSDM\\data\\europe_boundaries.shp")

clipped <- mask(burn_acc, study_area)
raster::writeRaster(writeRaster(clipped, paste0(burn_acc_out_path, "\\", "Burn_Acc_StudyArea.tif"), format = 'GTiff', 
                        overwrite = T, progress = 'text'))

values(clipped)[values(clipped) < 200] <- NA
plot(clipped)
burn_acc_vec <- rasterToPoints(clipped, function(x){x>=200}, spatial = TRUE)
plot(burn_acc_vec)

writeOGR(burn_acc_vec, "C:\\Users\\sandr\\Documents\\EAGLE_Data\\Burn_Acc_Vec.shp", layer = "burn_acc_vec", driver = "ESRI Shapefile")
