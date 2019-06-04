#----------------------------------------------------------------------------------------------------#
# Name:     01_GetMODISBurnIntensity                                                                 #
# Purpose:  Downloads a monthly MODIS time series for a given period of time and creates a burn      #
#           accumulation raster representing the number of wildfires occured within the specified    #
#           time. The result of this script is used to create wildfire occurence data.               #
# Created:  2019-06-01                                                                               #
# Updated:  2019-06-04                                                                               #
# Author:   Sandro Groth                                                                             #
#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#
# SETTINGS

MODIS_out_dir <- "C:\\Users\\sandr\\Documents\\EAGLE_Data\\MODIS_Temp"
mosaic_out_path <- "C:\\Users\\sandr\\Documents\\EAGLE_Data\\Burn_Mosaic"
burn_acc_out_path <- "C:\\Users\\sandr\\Documents\\EAGLE_Data"
time_range <- c('2015-01-01', '2018-12-31')

#----------------------------------------------------------------------------------------------------#
# IMPORTS

library(getSpatialData)
library(raster)
library(rgdal)
library(sp)
library(gdalUtils)

#----------------------------------------------------------------------------------------------------#

# read study area
study_area <- readOGR("data\\europe_boundaries.shp")
# Sources: Esri, Global Mapping International, U.S. Central Intelligence Agency (The World Factbook)

# set getSpatialData aoi
set_aoi()

# Login to USGS
login_USGS("Sandro.Groth")

# Query the LAADS DAAC
query <- getMODIS_query(time_range, "MODIS_MCD64A1_V6")

# Build up a sequence of product dates
date_start <- gsub("-", "/", time_range[1])
date_end <- gsub("-", "/", time_range[2])
product_dates <- seq.Date(from = as.Date(date_start), to = as.Date(date_end), by = "month")
product_dates <- as.character(product_dates)

# create a burn acc dummy
burn_acc <- NULL

# loop thourgh all monthly records in query and add burned pixels to burn acc
i <- 1
while (i <= length(product_dates)) {
  message(paste0("Starting processing for date: ", as.character(product_dates[i])))
  
  # check if the mosaiced burned area product already exists
  if (length(list.files(mosaic_out_path, pattern = as.character(product_dates[i])) > 0)) {
    # import existing burned area mosaic
    raster_mosaic <- raster(paste0(mosaic_out_path, "\\", "Burn_Mosaic_", 
                                   as.character(product_dates[i]), ".tif"))
  } else {
    # select only query records for current month
    query_sel <- query[query$acquisitionDate == product_dates[i],]
    # download all records for current month
    query_data <- getMODIS_data(query_sel, dir_out = MODIS_out_dir)
    # initialize an empty raster mosaic
    raster_mosaic <- NULL
    # loop through all downloaded records and add them to the mosaic
    j <- 1
    while (j <= length(query_data)) {
      out_file <- gsub("hdf", "tif", query_data[j])
      # check if a .tif version of the tile already exists
      if (!file.exists(out_file)) {
        # convert .hdf to .tif
        sds <- get_subdatasets(query_data[j])
        name <- sds[1]
        gdal_translate(name, dst_dataset = out_file)
        message(paste0("TIF created: ", out_file))
      } else {
        message(paste0("TIF already exists: ", out_file))
      }
      # add tile to mosaic
      if (is.null(raster_mosaic)) {
        raster_mosaic <- raster(out_file)
        message(paste0("New Mosaic created for date: ", product_dates[i]))
      } else {
        message(paste0("Dataset added to mosaic for date: ", product_dates[i]))
        curr_ras <- raster(out_file)
        raster_mosaic <- mosaic(raster_mosaic, curr_ras, fun = mean)
      }
      j <- j + 1
    }
    # export raster mosaic
    message(paste0("Writing raster mosaic for date: ", product_dates[i]))
    mosaic_name <- paste0("Burn_Mosaic_", as.character(product_dates[i]))
    writeRaster(raster_mosaic, paste0(mosaic_out_path, "\\", mosaic_name, ".tif"), 
                format = 'GTiff', overwrite = T, progress = 'text')
    
  }
  message("Reclassifying raster mosaic...")
  # add values to burn accumulation
  raster_mosaic[is.na(raster_mosaic[])] <- 0 
  values(raster_mosaic)[values(raster_mosaic < 0)] <- 0
  values(raster_mosaic)[values(raster_mosaic > 0)] <- 1
  message("Plotting Raster Mosaic...")
  plot(raster_mosaic)
  if (is.null(burn_acc)) {
    message("Creating Burn Accumulation...")
    burn_acc <- raster_mosaic
    message("Burn accumulation created.")
  } else {
    message("Adding burned pixels to burn accumulation...")
    values(burn_acc) <- values(burn_acc) + values(raster_mosaic)
    message(paste0("Burned areas added to Burn Accumulation for date: ", product_dates[i]))
  }
  i <- i + 1
  # Savie burn accumulation
  writeRaster(burn_acc, paste0(burn_acc_out_path, "\\", "Burn_Acc.tif"), format = 'GTiff', 
              overwrite = T, progress = 'text')
  message("Plotting Burn Accumulation...")
  plot(burn_acc)
}
