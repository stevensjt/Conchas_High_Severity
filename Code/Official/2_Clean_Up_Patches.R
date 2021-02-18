##Author: Jens Stevens; stevensjt@gmail.com
##Purpose: Implement a series of simplifying operations on the vectorized high-severity patch layers 
#for three different scenarios within the Las Conchas footprint. These exercises produce the "Clean_Cut" layers
#which then in turn serve as the null models for tuning and comparing the simulations.


####0. Read Libraries and Import Data###
library(sf) #version 0.9-6 for read_sf()
library(smoothr) #version 0.1.2 for fill_holes (used in 0_Functions)
library(lwgeom) #version 0.2-5 for st_split (used in 0_Functions)
library(ggplot2) #version 3.2.1 for ggplot()
source("./Code/Official/0_Functions.R") #read in customized functions necessary for processing.

hs_fire <- #read raw high severity patch layer in sf format for one of three scenarios
  #read_sf("./SpatialInput/4.Spatial_Scenarios", layer="lc_raw") #Scenario 1: Las Conchas
  #read_sf("./SpatialInput/4.Spatial_Scenarios", layer="all_raw")
  read_sf("./SpatialInput/4.Spatial_Scenarios", layer="treeless_raw") #Skip to intermediate step in #2

####Step 1: Process high severity patch layer and fill holes, to prepare for cutting####
hs_fire <- #Set as singlepart (one row per polygon, each polygon has a "single part"); suppress warning about applying attributes
  st_cast(hs_fire, to ="POLYGON", warn = FALSE)
hs_fire <- #Simplify and standardize columns for patch layer, just keeping geometry column
  hs_fire[,which(names(hs_fire)=="geometry")]
hs_fire$Area_m2 <- #Calculate area of each patch; assumes units are meters.
  st_area(hs_fire)
hs_fire$Area_ha <- #Calculate area in ha of each patch
  as.vector(hs_fire$Area_m2)*0.0001
par(mar = c(0.5, 0.1, 0.9, 0.1)) #Shrink margins for efficient use of plot space

hs_fill <- #Fill holes <= 1 ha; delete crumbs <= 1 ha. Takes 30 seconds.
  donut_holes_and_crumbs(hs_fire)
hs_fill$PatchName <-
  seq(1:nrow(hs_fill))
hs_fill$RemoveLater <- NA #Placeholder for patches that are eventually split and need to be removed.

#Histogram of original patch distribution
ggplot(hs_fill, aes(x = Area_ha)) +
  geom_histogram(binwidth = 0.1) +
  scale_x_log10(breaks = c(1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192)) + 
  theme(axis.text.x = element_text(angle = 45, size = 14))

####Step 1 modified: Process high severity patch layer as raster - remove crumbs and hoes - and prepare for cutting####
#CHECKME option to do this and save time but not sure what priority is.

####2. split polygon chunks separated by narrow bridges####
hs_pinch <- #temporarily drop crumbs <= 1,000,000 m2 = 100ha, minimum crumb size retained is 10,001 m2
  #here we are only looking to chop up large (>100 ha) polygons
  hs_fill[hs_fill$Area_ha > 100,]
hs_fill_small <- #Store the polygons <= 100 ha, for subsequent re-integration
  hs_fill[hs_fill$Area_ha <= 100,]


#If working with Las Conchas RdNBR layer: START HERE
Sys.time()
lc_clean_cut <- #Conduct splits via "breath increments" and "lungs" functions. SLOW.
  breath_increments(hs_pinch,increments = c(-15, -30, -60, -120))
lc_clean_cut_full <- #Combine small unsplit polygons (hs_fill_small) with larger split polygons
  #into a single shapefile
  rbind(lc_clean_cut, hs_fill_small)
Sys.time()

#If working with Las Conchas - Dome - Cerro[grande] RdNBR layer (LCDC):
Sys.time()
all_clean_cut <- #Conduct splits via "breath increments" and "lungs" functions. SLOW.
  breath_increments(hs_pinch,increments = c(-15, -30, -60, -120))
all_clean_cut_full <- #Combine small unsplit polygons (hs_fill_small) with larger split polygons
  rbind(all_clean_cut, hs_fill_small)
Sys.time()

#If working with Las Conchas Treeless layer:
#Intermediate steps to split along Hwy 4 and save mega time:
#st_write(st_make_valid(hs_pinch), "./SpatialOutput/LasConchas_Treeless_Intermediate.shp",delete_layer = T)
#Then process in QGIS to split polygon manually along Hwy 4 and near crater rim
#hs_pinch <- read_sf("./SpatialOutput", layer="LasConchas_Treeless_Intermediate") #Increase from 9 to 13 patches
#names(hs_pinch)[c(3,4)] <- c("PatchName", "RemoveLater") #Fix formatting error in column names

Sys.time() #Processing time here is SLOW and could probably be sped up by splitting the giant chunk along Hwy 4 before anything. Done.
LasConchasTreeless_clean <-  #Conduct splits via "breath increments" and "lungs" functions. SLOW.
  breath_increments(hs_pinch2,increments = c(-15, -30, -60, -120))
LasConchasTreeless_clean_full <- #Combine small unsplit polygons (hs_fill_small) with larger split polygons
  rbind(LasConchasTreeless_clean, hs_fill_small)
Sys.time()

#Option to write products to file
#st_write(lc_clean_cut_full, "./SpatialOutput/lc_clean_cut.shp",delete_layer=T)
#st_write(all_clean_cut_full, "./SpatialOutput/all_clean_cut.shp",delete_layer=T)
#st_write(LasConchasTreeless_clean_full, "./SpatialOutput/treeless_clean_cut.shp",delete_layer=T)
