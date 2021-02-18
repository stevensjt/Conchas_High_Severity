##Author: Jens Stevens; stevensjt@gmail.com
##Purpose: Process Sean Parks CBI_bc layers, into thresholded high-severity pixels CBI >= 225, then into polygons.
##Complete as of 12/10/20, don't need to run again.

####0. Read Libraries and Import Data####
library(sf)
library(raster)
library(stars)
library(dplyr)


fire_name_table <- #Get fire names as part of Fire_Atlas shapefile
  read_sf("./SpatialInput/1.Fire_Atlas", layer = "fires_that_intersect_LC")
lc_file <- dir("./SpatialInput/2.CBI_Raw_LasConchas/")
other_hs <- dir("./SpatialInput/2.CBI_Raw_OtherHS/") #Identify other raw CBI rasters w/ substantial HS
hs_list <- list() #Set up list object to hold high severity polygons
perim <- read_sf("../../../GIS/Fires/Las Conchas/Perimeter/", layer="Las_Conchas_Perim")

####1. Process Las Conchas####
#Las Conchas goes first because we want to save it separately in the ./SpatialInput/Scenarios folder
r <- #Read in Las Conchas CBI raster
  raster(paste0("./SpatialInput/2.CBI_Raw_LasConchas/",lc_file))
#plot(r, main = fire_name) #Option to plot fire and make sure it looks good.
cbi.binary <- r #New object to store the binary high severity pixels as "1"
cbi.binary[r >= 2.25] <- 1 #Set high severity pixels to 1
cbi.binary[r < 2.25] <- NA #Set other pixels to NA
#plot(cbi.binary, main = "Las Conchas") #Option to plot high severity fire (raster) and make sure it looks good
#writeRaster(cbi.binary,"./SpatialInput/4.Spatial_Scenarios/lc_raw.tif") Option to write Las Conchas raster to file
hs_vector <- cbi.binary %>% #Steps to efficiently convert raster to vector
  st_as_stars() %>% #convert raster to a "stars" 
  st_as_sf(as_points = FALSE, merge = TRUE) %>% #"sf" object is a multipart polygon w/ data frame
  st_buffer(dist = 0) %>% #This corrects self-intersecting rings (per Stack Overflow)
  st_transform(crs = 5070) %>% #Set CRS as NAD83 / Conus Albers (EPSG 5070)
  st_intersection(st_geometry(perim)) #Only take what intersects the official Las Conchas perim.
  
names(hs_vector)[1] <- c("Name") #Treat first column as "id" for merging purposes (it's constant)
hs_vector$Name = "lc" #All polygons here belong to Las Conchas "lc"
hs_vector <- hs_vector %>% #make singlepart and check validity
  group_by(Name) %>% #dissolve step 1: Common variable for all polygons (merge_id = 1)
  summarize() #dissolve step 2
#st_is_valid(hs_vector) #Option to check validity

#plot(hs_vector$geometry, add = T) #Option to overlay high severity fire (vector) and make sure it looks good

#write_sf( #Write sf object to shapefile
#  hs_vector, 
#  ("./SpatialInput/4.Spatial_Scenarios/lc_raw.shp"),
#  overwrite = TRUE)


####2. Process other fires####
all_hs <- #Read the Las-Conchas-only CBI high severity layer, 
  #which becomes the template for merging high severity patches from other fires
  read_sf("./SpatialInput/4.Spatial_Scenarios/lc_raw.shp") #if running from file
  #hs_vector #or store from section 1 above, if just completed
all_hs$Name <- "all" #This becomes the template for adding on the other high-severity fires.
st_is_valid(all_hs) #Check validity.
#all_hs <- st_make_valid(all_hs) #Fix topology errors, can happen if reading from file.

for(f in 1:length(other_hs)){ #takes about 1 minute
  #For every fire within Las Conchas, here we read the raw CBI file, select pixels >= 225,
    #Mask out everything else, convert to shapefile, and save shapefile to file.
  r <- #Read in the given fire
    raster(paste0("./SpatialInput/2.CBI_Raw_OtherHS/",other_hs[f]))
  fire_name <- #Temporary object to hold that fire's name
    fire_name_table$Fire_Name[which(!is.na(pmatch(fire_name_table$Fire_ID,other_hs[f])))]
  plot(r, main = fire_name) #Option to plot fire and make sure it looks good.
  
  cbi.binary <- r #New object to store the binary high severity pixels as "1"
  cbi.binary[r >= 2.25] <- 1 #Set high severity pixels to 1
  cbi.binary[r < 2.25] <- NA #Set other pixels to NA
  plot(cbi.binary, main = fire_name) #Option to plot fire and make sure it looks good
  #CHECKME eventually can add script here to merge these rasters onto the Las Conchas and save the "all" raster to file.
  
 hs_vector <- cbi.binary %>% #Steps to efficiently convert raster to vector
   st_as_stars() %>% #"stars" object is a multipart polygon (?)
   st_as_sf(as_points = FALSE, merge = TRUE) %>% #"sf" object is a singlepart polygon
   st_buffer(dist = 0) %>% #This corrects self-intersecting rings (per Stack Overflow)
   st_transform(crs = 5070) %>% #Set CRS as NAD83 / Conus Albers (EPSG 5070)
   st_intersection(st_geometry(perim))  #Only take what intersects the official Las Conchas perimeter.
 #Warning is ok
 
 names(hs_vector)[1] <- c("Name") #Treat first column as "id" for merging purposes (it's constant)
 hs_vector$Name = "all" #Regardless of which fire this is, we label it "all" for merging.
 plot(hs_vector$geometry, add = T) #Option to plot fire and make sure it looks good
 #write_sf( #Option to write sf object to shapefile to inspect individual fires.
 #  hs_vector, 
 #  paste0("./SpatialInput/3.CBI_Vector_OtherHS/", fire_name, ".shp"),
 #  overwrite = TRUE)
 
 all_hs <- #Merge to a multipart polygon (only one "object") + dissolve boundaries ("summarize" does this)
   rbind(all_hs,hs_vector) %>% #merge
   group_by(Name) %>% #dissolve step 1: Common variable for all polygons (Name = "all")
   summarize() #dissolve step 2
 #Move on to the next fire
}

plot(all_hs$geometry, col = "darksalmon") #Option to plot and make sure everything looks ok
#write_sf( #Option to write sf object to shapefile to inspect individual fires.
#  all_hs, 
#  ("./SpatialInput/4.Spatial_Scenarios/all_raw.shp"),
#  overwrite = TRUE)

##END