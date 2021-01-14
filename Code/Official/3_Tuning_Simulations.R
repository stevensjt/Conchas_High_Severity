#Author: Jens Stevens; stevensjt@gmail.com
#This code tweaks the parameters of the gaussian field simulations, trying to iteratively match to the values of the observed "null model" layer.

####0. Libraries####
source("./Code/Official/0_Functions.R") #read in customized functions necessary for processing.
library(lwgeom) #version 0.2-5 for st_split (used in 0_Functions)
library(ggplot2) #version 3.2.1 for ggplot()
library(tidyverse)

####1. Read Data####
lc <- read_sf("./SpatialOutput/", layer = "lc_clean_cut")
perim <- read_sf("../../../GIS/Fires/Las Conchas/Perimeter/", layer="Las_Conchas_Perim")
#perim_sp <- readOGR("../../../GIS/Fires/Las Conchas/Perimeter/", layer="Las_Conchas_Perim") #Need to access projection. deprecated?
scars <- read_sf("../../../GIS/Fires/Las Conchas/FireScarLocations/", layer="fs_in_las_conchas")


####2. Define original parameters, to compare tuning against####
#Initialize comparison data frame with observed original patch size distribution
df_original <- data.frame(sim = 0, Area_ha = sort(lc$Area_ha, decreasing = T))

####3a. Find optimal values####
sims_list_spatial <- list()

nugget <- 1
magvar <- 100
pct <- 5
prop_hs <- as.vector(round(sum(st_area(lc))/st_area(perim),3)[[1]])
parms_text <- paste0("pct_", pct,"_nug_",nugget,"_mv_",magvar)

#Run one to check it out
s = 1 #simulation number
gc()
out = random_cheese_patchdist(pct = pct, nugget = nugget, magvar = magvar, prop_hs = prop_hs)
plot(out$geometry,col = out$Area_ha, border = "black", pal = rainbow(50), main = parms_text)

#Run ten and store the results
df_compare <- df_original
for(s in c(13:20)){
  gc()
  out <- random_cheese_patchdist(pct = pct, nugget = nugget, magvar = magvar, prop_hs = prop_hs)
  
  #Test if the random cookie cutter of the Las Conchas perimeter picked up an accurate proportion of high severity
  test_prop_hs <- abs( as.vector(round(sum(st_area(out))/st_area(perim),3)[[1]]) - prop_hs ) < 0.03
  if(test_prop_hs){#store results if proportion high severity within 3% of true value
    #double check but this test can probably be deprecated
    df_out <- data.frame(sim = s, Area_ha = sort(out$Area_ha, decreasing = T)) #store patch distribution
    df_compare <- rbind(df_compare, df_out) #add patch distribution to comparison data frame
    sims_list_spatial[[s]] <- out #store spatial layer in list
    print(paste("sim", s, "stored"))
    }
  else {print(paste("sim", s, "not stored"))}
  
}
write_csv(df_compare,"test.csv")

#df_plot_c <- dv_compare %>%
  

a <- ggplot() + 
  geom_sf(data = perim, fill = NA) + 
  geom_sf(data = lc, aes(fill = Area_ha), col = NA) +
  labs(title = "Las Conchas high severity")

b_plot <- sims_list_spatial[[2]]
b <- ggplot() + 
  geom_sf(data = perim, fill = NA) + 
  geom_sf(data = b_plot, aes(fill = Area_ha), col = NA) +
  labs(title = "Las Conchas high severity")

c <- ggplot(df