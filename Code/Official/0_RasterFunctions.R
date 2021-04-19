##Author: Jens Stevens; stevensjt@gmail.com
##Purpose: Define a series of spatial processing functions that are the workhorses for simplifying
#and simultaing high severity patch layers, working more extensively in raster than before.

raster_holes_and_crumbs <- function(hs_fire){
  #This function takes a vector layer of high-severity patches and 
  #1) 2) drops crumbs (small patches) <= 1ha, and fills holes <= 1ha
  
  #1) drop crumbs <= 10,000 m2 = 1 ha = 11.1 pixels
  rc <- clump(hs_fire, directions = 4) #pretty fast
  #plot(rc) #deprecate
  clump <- data.frame(freq(rc)) #number of cells in each unique clump
  clump <- clump[ clump$count > 11, ] #keep clumps with > 11 cells (> 1 ha)
  clump_ids = as.vector(na.exclude(clump$value)) #store ID's of clumps to keep

  rc[which(! getValues(rc)%in%clump_ids)] <- NA #set pixel values NOT in a "large" clump to NA
  rc[which(!is.na(getValues(rc)))] <- 1 #set pixel values in "large" clumps to 1
  #plot(rc)
  
  #2) fill holes <= 10,000 m2 = 1 ha = 11.1 pixels
  r_inverse <- rc
  #re-cast raster as inverse
  r_inverse[which(is.na(getValues(r_inverse)))] <- 2
  r_inverse[which(getValues(r_inverse)==1)] <- NA
  r_inverse[which(getValues(r_inverse)==2)] <- 1
  
  #re-run the "drop crumbs" operation
  rc <- clump(r_inverse, directions = 4) #pretty fast
  #plot(rc) #deprecate
  clump <- data.frame(freq(rc)) #number of cells in each unique hole
  clump <- clump[ clump$count > 11, ] #keep holes with > 11 cells (> 1 ha)
  clump_ids = as.vector(na.exclude(clump$value)) #store ID's of holes to keep
  
  rc[which(! getValues(rc)%in%clump_ids)] <- NA #set pixel values NOT in a "hole to keep" to NA
  rc[which(is.na(getValues(rc)))] <- -1 #set NA's (high-severity fire) to -1 temporarily
  rc[which(getValues(rc) != -1)] <- NA #set anything not high-severity to NA
  rc[which(getValues(rc) == -1)] <- 1 #set high-severity (including small holes) to 1 to finish.
  #plot(rc)
  
  #START HERE check area inside or outside of function?
  
  rc #return the raster without small holes and crumbs.
  
}

raster_cheese_patchdist <- function(pct,nugget,magvar,prop_hs, seed){
  #This function generates one random patch layer for a given set of variogram parameters
  test_prop_hs <- FALSE
  prop_hs_adj <- 0.01
  print("check1")
  while(!test_prop_hs){#Simulate surfaces, clip by perim, and delete crumbs/holes
    #until a processed simulation falls within 2% ('prop_hs_adj') of the target area burned at HS
    rng <- #set the variogram range as the radius of a circle with an area equal to 'pct' percent of the grid area
      round(sqrt(((pct/100) * xmax * ymax)/pi), 0) #define rng based on pct
    # Generate a continuous Gaussian random field based on the variogram input parameters
    gf <- #Slow, ~20 seconds w 30 m resolution (less now)?
      nlm_gaussianfield(ncol = xmax, nrow = ymax, autocorr_range = rng, mag_var = magvar, nug = nugget, user_seed = seed)
    #Reproject the random grid to be on the same CRS as Las Conchas perimeter.
    extent(gf) <- round(extent(perim),0)
    gf@crs <- crs(perim)
    gf.mask <- raster::mask(gf, perim) #mask by Las Conchas
    #plot(gf.mask)
    print("check2")
    # Convert the continuous Gaussian field into a binary grid
    convert.to.binary <- function(x) { ifelse(x <  cutoff, NA, 1) }
    # convert grid to a vector and sort in decreasing order to find a cutoff point where the proportion of cells with values 
    # above the cutoff point is equal to 'prop_hs'
    gfv <- sort(as.vector(gf.mask), decreasing = TRUE) #gfv = gaussian field vector; length = n cells in grid
    cutoff <- #prop_hs percent cells of total grid, everything larger is high severity (when sorted by descending value)
      #NOTE here we adjust the target proportion high severity up by 1% the first time, because the "holes and crumbs" operation below 
      #generally reduces the total HS area following the convert.to.binary operation. So increase the target by 1%, and
      #if it still fails the check below, we iteratively increase until it is within 1% of target.
      gfv[length(gfv) * (prop_hs + prop_hs_adj)] 
    gf.binary <- calc(gf.mask, fun = convert.to.binary)
    #plot(gf.binary)

    #"Simplify" raster by removing holes and crumbs < 1 ha
    gf.simple <- raster_holes_and_crumbs(gf.binary) #Now takes < 10 seconds
    
    area_sim <- # number of cells = 1, times area of a cell 900.5704 m2 
      freq(gf.simple, value = 1) * 900.5704 #in m2
    prop_hs_sim <- as.vector(round(area_sim/st_area(perim),3)[[1]])
    
    test_prop_hs <- #Ensure simulation is still within 1% of target high-severity
      abs( prop_hs_sim - prop_hs ) < 0.01 
    print(paste("target proportion hs =",prop_hs, "; simulated proportion =", prop_hs_sim))
    #here we adjust target proportion if we fail the check
    if(prop_hs_sim < (prop_hs-prop_hs_adj)) {
      #if our simulation was outside the target proportion high-severity, adjust the target
      prop_hs_adj <- prop_hs - prop_hs_sim
      print(paste("adjusting target proportion hs upward by ", prop_hs_adj))
    }
    if(prop_hs_sim > (prop_hs)) {
      #if our simulation was outside the target proportion high-severity, adjust the target
      prop_hs_adj <- prop_hs - prop_hs_sim
      print(paste("adjusting target proportion hs downward by ", prop_hs_adj))
    }
  }

  #Now we polygonize the patches
  writeRaster(gf.simple,"./SpatialOutput/sim.tif", overwrite = TRUE)
  #Polygonize using "stars" from https://cran.r-project.org/web/packages/stars/vignettes/stars5.html
  hs_fill <- read_stars("./SpatialOutput/sim.tif")
  hs_fill <- st_as_sf(hs_fill, as_points = FALSE, merge = TRUE)
  
  #Now proceed with the polygon cleaning and cutting operations
  hs_fill$Area_m2 <- #Calculate area of each patch after deleting holes and crumbs
    st_area(hs_fill)
  hs_fill$Area_ha <- #Calculate area in ha of each patch
    as.vector(hs_fill$Area_m2)*0.0001
  hs_fill$PatchName <-
    seq(1:nrow(hs_fill))
  hs_fill$RemoveLater <- NA #Placeholder for patches that are eventually split and need to be removed.
    
  hs_pinch <- #drop crumbs <= 1,000,000 m2 = 100ha, minimum crumb size retained is 10,001 m2
    #here we are looking at which large (>100 ha) polygons have pinch points that we can chop up
    hs_fill[hs_fill$Area_ha > 100,]
  hs_fill_small <- #Store the polygons <= 100 ha, for subsequent re-integration
    hs_fill[hs_fill$Area_ha <= 100,]
  
  print("got a good sim; starting breath increments")
  sim_clean <- #Conduct splits via "breath increments" and "lungs" functions. SLOW.
    breath_increments(hs_pinch,increments = c(-15, -30, -60, -120))
  sim_clean_full <- #Combine small unsplit polygons (hs_fill_small) with larger split polygons
    rbind(sim_clean, hs_fill_small)
  print(paste(Sys.time(), "finished sim ", s))
  
  #sim_clean_full$Area_ha #return the vector of patch size values
  sim_clean_full #return the polygons
  
}

raster_cheese_patchdist_multipart <- function(pct,nugget,magvar,prop_hs, seed, cutoffs){
  #This function generates one random patch layer for a given set of variogram parameters
  #But the variogram runs separately for large and small patches
  test_prop_hs <- FALSE
  prop_hs_adj <- 0 #For later adjustments
  print("check1")
  while(!test_prop_hs){#Simulate surfaces, clip by perim, and delete crumbs/holes
    #until a processed simulation falls within 2% ('prop_hs_adj') of the target area burned at HS
    rng <- #set the variogram range as the radius of a circle with an area equal to 'pct' percent of the grid area
      round(sqrt(((pct/100) * xmax * ymax)/pi), 0) #define rng based on pct
    # Generate a continuous Gaussian random field based on the variogram input parameters
    gf <- #Slow, ~20 seconds w 30 m resolution (less now)?
      nlm_gaussianfield(ncol = xmax, nrow = ymax, 
                        autocorr_range = rng, mag_var = magvar, nug = nugget, user_seed = seed)
    #Reproject the random grid to be on the same CRS as Las Conchas perimeter.
    extent(gf) <- round(extent(perim),0)
    gf@crs <- crs(perim)
    gf.mask <- raster::mask(gf, perim) #mask by Las Conchas
    #plot(gf.mask)
    print("check2")

    # Convert the continuous Gaussian field into a binary grid
    convert.range.to.binary  <- function(x) { 
      v = cutoff_vals
       if(length(v)==1){ #Standard case, only take large values
        ifelse(x > v[1], 1, NA)
      } else if(length(v) == 3){ #Take the large patches from large rnorm values
        #And take smaller patches from middling rnorm values
        ifelse(between(x,v[1],v[2]) | x > v[3], 1, NA)
        #ifelse(between(x,v[1],v[2]) , 1, ifelse( x > v[3] ,2, NA)) #option for color-coding
        
      }
    }
    gfv <- sort(as.vector(gf.mask), decreasing = FALSE) #gfv = gaussian field vector
    cutoff_vals <- 
      c(gfv[length(gfv) * cutoffs[1]], 
        gfv[length(gfv) * cutoffs[2]],
        gfv[length(gfv) * cutoffs[3]])
    cutoff_vals <- cutoff_vals[!is.na(cutoff_vals)] # mod for single threshold
    
    gf.binary <- calc(gf.mask, fun = convert.range.to.binary)
    #plot(gf.binary, col = c("darkred", "darkblue"))
    
    #"Simplify" raster by removing holes and crumbs < 1 ha
    gf.simple <- raster_holes_and_crumbs(gf.binary) #Now takes < 10 seconds
    #plot(gf.simple)
    
    area_sim <- # number of cells = 1, times area of a cell 900.5704 m2 
      freq(gf.simple, value = 1) * 900.5704 #in m2
    prop_hs_sim <- as.vector(round(area_sim/st_area(perim),3)[[1]])
    
    test_prop_hs <- #Ensure simulation is still within 1% of target high-severity
      abs( prop_hs_sim - prop_hs ) < 0.01 
    print(paste("target proportion hs =",prop_hs, "; simulated proportion =", prop_hs_sim))
    #here we adjust target proportion if we fail the check
    if(prop_hs_sim < (prop_hs)) {
      #if our simulation was below the target proportion high-severity, bump up the lower range by increment
      inc <- (prop_hs - prop_hs_sim)/2 #halve the distance; ensures you don't overshoot
      cutoffs[2] = cutoffs[2] +inc #add more area to the simulation
      print(paste("adjusting target proportion hs upward by ", inc))
    }
    if(prop_hs_sim > (prop_hs)) { #this statement probably deprecated because redundant except for "adjusting downward" text
      #if our simulation was above the target proportion high-severity, reduce the lower range by increment
      inc <- (prop_hs - prop_hs_sim)/2 #will be negative value in this case
      cutoffs[2] = cutoffs[2] + inc #trim more area from the simulation
      print(paste("adjusting target proportion hs downward by ", inc))
    }
  }
  
  #Now we polygonize the patches
  writeRaster(gf.simple,"./SpatialOutput/sim.tif", overwrite = TRUE)
  #Polygonize using "stars" from https://cran.r-project.org/web/packages/stars/vignettes/stars5.html
  hs_fill <- read_stars("./SpatialOutput/sim.tif")
  hs_fill <- st_as_sf(hs_fill, as_points = FALSE, merge = TRUE)
  
  #Now proceed with the polygon cleaning and cutting operations
  hs_fill$Area_m2 <- #Calculate area of each patch after deleting holes and crumbs
    st_area(hs_fill)
  hs_fill$Area_ha <- #Calculate area in ha of each patch
    as.vector(hs_fill$Area_m2)*0.0001
  hs_fill$PatchName <-
    seq(1:nrow(hs_fill))
  hs_fill$RemoveLater <- NA #Placeholder for patches that are eventually split and need to be removed.
  
  hs_pinch <- #drop crumbs <= 1,000,000 m2 = 100ha, minimum crumb size retained is 10,001 m2
    #here we are looking at which large (>100 ha) polygons have pinch points that we can chop up
    hs_fill[hs_fill$Area_ha > 100,]
  hs_fill_small <- #Store the polygons <= 100 ha, for subsequent re-integration
    hs_fill[hs_fill$Area_ha <= 100,]
  
  print("got a good sim; starting breath increments")
  sim_clean <- #Conduct splits via "breath increments" and "lungs" functions. SLOW.
    breath_increments(hs_pinch,increments = c(-15, -30, -60, -120))
  sim_clean_full <- #Combine small unsplit polygons (hs_fill_small) with larger split polygons
    rbind(sim_clean, hs_fill_small)
  print(paste(Sys.time(), "finished sim ", s))
  
  #sim_clean_full$Area_ha #return the vector of patch size values
  sim_clean_full #return the polygons
  
}
