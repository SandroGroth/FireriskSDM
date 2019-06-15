library(getSpatialData)
library(raster)
library(rgdal)
library(sp)
library(gdalUtils)

#----------------------------------------------------------------------------------------------------#
# SETTINGS

study_area_path <- "C:\\Users\\sandr\\Documents\\EAGLE_Data\\EuropeanStates.shp"
GLCC_out_dir <- "C:\\Users\\sandr\\Documents\\EAGLE_Data\\MODIS_GLCC"
time_range <- c('2017-01-01', '2017-12-31')
GLCC_tif_out <- "C:\\Users\\sandr\\Documents\\EAGLE_Data\\MODIS_GLCC\\GLCC_Majority_Type_1.tif"

#----------------------------------------------------------------------------------------------------#

study_area <- readOGR(study_area_path)

set_aoi()

login_USGS("Sandro.Groth")

# Query the LAADS DAAC
query <- getMODIS_query(time_range, "MODIS_MCD12C1_V6")

# download landcover product
query_data <- getMODIS_data(query, dir_out = GLCC_out_dir)


sds <- get_subdatasets(query_data)
name <- sds[1]
gdal_translate(name, dst_dataset = GLCC_tif_out)

# reproject
glcc_tif <- raster(GLCC_tif_out)
glcc_reproj <- projectRaster(glcc_tif, crs = study_area@proj4string)
glcc_resample <- resample(glcc_reproj, biocrop)

raster::writeRaster(glcc_resample, paste0(GLCC_out_dir, "\\", "GLCC_Majority_Type_1.tif"), fromat = "GTiff", overwrite = T, progress = "text")
