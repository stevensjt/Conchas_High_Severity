library(sf) #version 0.9-3 for read_sf()
library(lwgeom)
library(ggplot2) #version 3.2.1 for ggplot()
source("./Code/Official/Lungs2.R") #read in customized functions necessary for processing.

####Step 1: Read and clean up high severity patch layer####
hs_fire <- #read raw high severity patch layer in sf format for one of three scenarios
  #read_sf("./SpatialInput/4.Spatial_Scenarios", layer="lc_raw") #Scenario 1: Las Conchas
  #read_sf("./SpatialInput/4.Spatial_Scenarios", layer="all_raw")
  read_sf("./SpatialInput/4.Spatial_Scenarios", layer="treeless_raw") #Skip to intermediate step in #2

hs_fire <- #Set as singlepart; suppress warning about applying attributes
  st_cast(hs_fire, to ="POLYGON", warn = FALSE)
hs_fire <- #Simplify and standardize columns for patch layer
  hs_fire[,which(names(hs_fire)=="geometry")]
hs_fire$Area_m2 <- #Calculate area of each patch; assumes units are meters.
  st_area(hs_fire)
hs_fire$Area_ha <- #Calculate area in ha of each patch
  as.vector(hs_fire$Area_m2)*0.0001
par(mar = c(0.5, 0.1, 0.9, 0.1)) #Shrink margins for efficient use of plot space


hs_fill <- #Fill holes <= 1 ha; delete crumbs <= 1 ha. Takes 30 seconds.
  #Number of discrete patches drops from 3916 to 674 for LC; 4478 to 720 for All, 5928 to 359 for Treeless
  donut_holes_and_crumbs(hs_fire)
hs_fill$PatchName <-
  seq(1:nrow(hs_fill))
hs_fill$RemoveLater <- NA #Placeholder for patches that are eventually split and need to be removed.

#Histogram of original patch distribution
ggplot(hs_fill, aes(x = Area_ha)) +
  geom_histogram(binwidth = 0.1) +
  scale_x_log10(breaks = c(1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192)) + 
  theme(axis.text.x = element_text(angle = 45, size = 14))

####2. split polygon chunks separated by narrow bridges####
hs_pinch <- #temporarily drop crumbs <= 1,000,000 m2 = 100ha, minimum crumb size retained is 10,001 m2
  #here we are looking at which large (>100 ha) polygons have pinch points that we can chop up
  #n = 23 for Las Conchas
  #n = n for All
  #n = 9 for Treeless (but subsequent manual splitting brings to 13, see below); 
  hs_fill[hs_fill$Area_ha > 100,]
hs_fill_small <- #Store the polygons <= 100 ha, for subsequent re-integration
  hs_fill[hs_fill$Area_ha <= 100,]


#If working with Las Conchas RdNBR layer: START HERE
Sys.time()
LasConchasRDNBR_clean <- #Conduct splits via "breath increments" and "lungs" functions. SLOW.
  breath_increments(hs_pinch,increments = c(-15, -30, -60, -120))
LasConchasRDNBR_clean_full <- #Combine small unsplit polygons (hs_fill_small) with larger split polygons
  #into a single shapefile
  rbind(LasConchasRDNBR_clean, hs_fill_small)
Sys.time()

#If working with Las Conchas - Dome - Cerro[grande] RdNBR layer (LCDC):
Sys.time()
LCDC_large_clean <- #Conduct splits via "breath increments" and "lungs" functions. SLOW.
  breath_increments(hs_pinch,increments = c(-15, -30, -60, -120))
LCDC_large_clean_full <- #Combine small unsplit polygons (hs_fill_small) with larger split polygons
  rbind(LCDC_large_clean, hs_fill_small)
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
#st_write(LasConchasRDNBR_clean_full, "./SpatialOutput/LasConchas_643_RdNBR_Clean.shp",delete_layer=T)
#st_write(LCDC_large_clean_full, "./SpatialOutput/LasConchas_Dome_Cerro_643_RdNBR_Clean.shp",delete_layer=T)
st_write(LasConchasTreeless_clean_full, "./SpatialOutput/LasConchas_Treeless_Clean.shp",delete_layer=T)
