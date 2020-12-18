##Author: Jens Stevens; stevensjt@gmail.com
##Purpose: Define a series of spatial processing functions that are the workhorses for simplifying
#and simultaing high severity patch layers.


donut_holes_and_crumbs <- function(hs_fire){
  #This function takes a vector layer of high-severity patches and 
  #1) fills holes <= 1ha, and 2) drops crumbs (small patches) <= 1ha
  
  #1) fill holes
  #plot(hs_fire$geometry, col = "darksalmon") #option to plot
  hs_fire$geometry = #fill holes <= 10,000 m2 = 1ha, minimum holesize retained is 10,001 m2
    smoothr::fill_holes(hs_fire$geometry, threshold = 10000) #Takes ~25 seconds
  hs_fire$Area_m2 <- #Re-Calculate area of each patch after filling holes; assumes metric
    st_area(hs_fire)
  hs_fire$Area_ha <- #Calculate area in ha of each patch
    as.vector(hs_fire$Area_m2)*0.0001
  #plot(hs_fire$geometry, col = "darksalmon") #option to plot
  
  #2) drop crumbs
  hs2 = #drop crumbs <= 10,000 m2 = 1ha, minimum crumb size retained is 10,001 m2
    #hs2 is second temporary object holding high severity polygons
    #need to create new object here because "drop_crumbs" changes the object length
    hs_fire[hs_fire$Area_ha > 1,]
  #plot(hs2$geometry, col = "darkred") #option to plot
  

  #return the polygons with holes filled and crumbs deleted, in decreasing size order
  hs2[order(hs2$Area_ha, decreasing = TRUE),]
}

test_chunks <- function(pol, sdist){
  #
  bparts = #if the value of sdist severs pinch points, bparts is a sfc_MULTIPOLYGON, 
    #bparts[[1]][[1]] is a list of coords for the first pinched polygon, which may or may not have holes
    #but bparts[[1]][[1]] is NOT an sfc_POLYGON, it is *just* a list.
    #-0.1 ensures polygons aren't connected by a "thread"
    st_buffer(pol, sdist - 0.1, joinStyle="MITRE", mitreLimit=2) 
  features = #if bparts is a sfc_MULTIPOLYGON, st_cast turns it into a POLYGON
    #features[[1]][[1]] is equal to bparts[[1]][[1]][[1]]; 
    #contains matrix of coords for outer wall or (more likley) a hole
    #st_cast(st_sfc(bparts), "POLYGON") #original line from web, mod below
    st_cast(bparts, "POLYGON")
  inhale1 = #here we only take the new "pinched" polygons that are > 10 ha
    features[which(as.vector(st_area(features))*0.0001 > 10)]
  
  output = length(inhale1) > 1
  output
}

lungs <- function(pol, sdist){
  #This function operates on a single polygon. It looks for narrow pinch points less than sdist (in meters)
  #If does this by "inhaling" (STEP 1) and seeing if that creates multiple polygons
  #If multiple polygons are created, STEP 2 "exhales", 
  #and the difference between the inhale and the exhale is used to identify where to split the polygon.
  
  output = list() #object to store output
  
  ##STEP 1: Inward Buffer ("Inhale")
  #pol = st_make_valid(pol) #Corrects rings that self-overlap. Maybe not needed; issue for l = 10
  bparts = #if the value of sdist severs pinch points, bparts is a sfc_MULTIPOLYGON, 
    #bparts[[1]][[1]] is a list of coords for the first pinched polygon, which may or may not have holes
    #but bparts[[1]][[1]] is NOT an sfc_POLYGON, it is *just* a list.
    #-0.1 ensures polygons aren't connected by a tiny "thread"
    st_buffer(pol, sdist - 0.1, joinStyle="MITRE", mitreLimit=2)
  features = #if bparts is a sfc_MULTIPOLYGON, st_cast turns it into a POLYGON
    #features[[1]][[1]] is equal to bparts[[1]][[1]][[1]]; 
    #contains matrix of coords for outer wall or (more likley) a hole
    #st_cast(st_sfc(bparts), "POLYGON") #original line from web, mod below
    st_cast(bparts, "POLYGON")
  inhale1 = #here we only take the new "pinched" polygons that are > 10 ha
    features[which(as.vector(st_area(features))*0.0001 > 10)]
  #plot(pol) #original
  ##plot(bparts, col = rainbow(8), add = T) #default plot for the sfc_MULTIPOLYGON doesn't separate colors
  ##plot(features, col = rainbow(8), add = T) #default plot for the sfc_POLYGON separates colors
  #plot(inhale1, col = rainbow(8), add = F) #The larger pinched polygons only (>10 ha)
  
  if(length(inhale1) > 1){ #Proceed if there is a split to make (if we have inhaled >=2 pinched polygons > 10 ha)
    ##STEP 2: Outward Buffer a bit beyond the original("Exhale")
    small_chunk = #If only working with 2 chunks after the inhale split, pick the small one
      inhale1[which.min(st_area(inhale1))]
    large_chunk = #Identify the "other" (larger) chunk for use during splitting process.
      inhale1[which.max(st_area(inhale1))]
    exhale1 = #This approximates the original small chunk prior to the inhale
      st_make_valid(
        st_cast(st_sfc(
          st_buffer(small_chunk,
                    -sdist * 1.5, #+ 1 is to push a bit past the original polygon
                    #had exhaled by twice as much (*-1.95) but don't need to? TBD
                    joinStyle="MITRE",
                    mitreLimit=2)
          ), "POLYGON")
        )
    
    #plot(exhale1,border="goldenrod",add=F)
    #plot(small_chunk, border = "blue", add = T)
    #plot(pol,add=T)
    
    ##STEP 3: Split and process
    exhale1_mls = #Convert the exhale polygon into a multilinestring
      st_cast(exhale1,"MULTILINESTRING")
    first_split = #Split original polygon by the exhale.
      #This is an intermediate step because the exhale likely missed some "danglers" off the small chunk,
      #In addition to severing the main "bridge"
      st_cast(st_split(pol,exhale1_mls)) 
    #plot(first_split, col = rainbow(8), add = F) #Tiny multicolored polygons ("danglers") are to be discarded
    second_split = #Modify the polygon for the second split by discarding any new danglers < 10000m2 (1 ha)
      st_cast(first_split[which(as.vector(st_area(first_split)) > 10000)])
    bridge = #The only polygon we want to split off is the largest one from the first split; 
      #this is the "mainland" whose bridge we want to sever. Buffering ensures a clean split.
      st_buffer(first_split[which.max(st_area(first_split))],1)
    #Temporary CHECKME; if there is a "within" step to properly identify the bridge, it happens here
    #plot(bridge,col="darkred",add=T)
    splitter = #Here we find the critical line segment(s) to do the splitting from the mainland
      st_intersection(exhale1_mls,bridge)
    #plot(splitter,col = "red", lwd = 3, add = T)
    splitpol = #This is the key split: Sever the original polygon at the splitter. Should be length 2.
      st_cast(st_split(pol,splitter))
    splitpol = #any remaining tiny polygons created by this are the result of "kitty-corner" splits and are discarded.
      splitpol[which(as.vector(st_area(splitpol)) > 10000)]
    #plot(splitpol, col = rainbow(8), add = T)
    
    output[[1]] = "Split polygon" 
    output[[2]] = splitpol
  } else if (length(inhale1) == 1) { 
    output[[1]] = "Not split" 
  }
   
  output
    
}

breath_increments <- function(hs, increments, to_plot = FALSE){
  #hs is the full set of patches to be (potentially) pinched
  #increments is a vector of negative values arranged from small to large buffer
  hs = hs[order(hs$Area_ha, decreasing = FALSE),] #Order polygons, smallest to largest. Only done once.
  for(sdist in increments){
    keep_splitting = TRUE ; l = 1
    print(paste("starting to buffer by ", sdist))
    #l is polygon index; resets after split polygons have been removed and sdist increases
    while(keep_splitting){ #For the Las Conchas MTBS layer, takes 
      #pol = hs$geometry[l] #Pick the polygon to test and split #Deprecate because passed to Test_chunks and lungs below
      
      if(test_chunks(hs$geometry[l],sdist)){
        #Proceed with modifying original polygons if the polygon in question has chunks to split:
        oldname = hs$PatchName[l]
        splitpol2 = lungs(hs$geometry[l], sdist)[[2]] #Store the split polygons
        #plot(pol, main = oldname) #deprecate
        if(to_plot){
          plot(splitpol2, col = rainbow(8), main = oldname)
        }
        print(paste("cut geometries for original polygon", oldname, "row", l, "by", sdist)) #optional counter
        hs$RemoveLater[l] = 1 #Flag the split polygon (row l) for removal
        if(length(splitpol2)>2){
          #This is a failsafe if a split happened to produce 3 polygons instead of 2. Not the most elegant.
          min.holder = which.min(st_area(splitpol2))
          splitpol3 = splitpol2[min.holder,] #Pick the smallest one to split
          splitpol4 = st_union(splitpol2[-min.holder,])
          splitpol2 = splitpol3
          splitpol2[[2]] = splitpol4[[1]]
        }
        hs[nrow(hs)+c(1,2),"geometry"]$geometry = st_geometry(splitpol2) #add split geometries
        hs$PatchName[c(nrow(hs)-1,nrow(hs))] = rep(oldname, 2) #split geometries get OG PatchName
        #hs$Area_m2 = st_area(hs) #Calculate area (m2) for new split polygons, deprecate?
        #hs$Area_ha = as.vector(hs$Area_m2)*0.0001 #Calculate area (ha) for new split polygons, deprecate?
        #hs = hs[order(hs$Area_ha, decreasing = FALSE),] #Big polygons at the end now, to keep pinching em. deprecate?
        if(is.na(hs[l,]$Area_ha) & to_plot == TRUE) {mtext("repeat split\n", 1)} #optional plot label to track splits of splits
      }
      if(l == nrow(hs)){
        #the check to end the while loop: if we've reached the end of the data frame
        print(paste("reached the end at row", l))
        keep_splitting = FALSE
      } else {l = l+1}
      
    } #Works well for Las Conchas fire, though l=24 is slow to process.
    print(paste("finished buffering by ", sdist))
    gc() #free memory
    hs = hs[-which(hs$RemoveLater == 1),]
    hs$Area_m2 = st_area(hs) 
    hs$Area_ha = as.vector(hs$Area_m2)*0.0001
  }
  hs
}

random_cheese_patchdist <- function(pct,nugget,magvar,pct_hs){
  rng <- round(sqrt(((pct/100) * xmax * ymax)/pi), 0)
  
  # Alternatively, set the variogram range in number of cells
  # rng <- 20
  
  # Generate a continuous Gaussian random field based on the variogram input parameters
  gf <- #Slow, ~20 seconds w 30 m resolution
    nlm_gaussianfield(ncol = xmax, nrow = ymax, autocorr_range = rng, mag_var = magvar, nug = nugget)
  
  #Reproject the random grid to be on the same CRS as Las Conchas perimeter.
  extent(gf) <- round(extent(perim),0)
  gf@crs <- perim_sp@proj4string
  gf <- projectRaster(gf, crs = crs(perim))
  
  #### Convert the continuous Gaussian field into a binary grid ###
  # Function to convert continuous Gaussian random field to binary data
  # The binary grid will have Values of 1 = burned and 0 = unburned
  convert.to.binary <- function(x) 
  {
    ifelse(x <  cutoff, 0, 1)
  }
  
  prop.burned = pct_hs
  
  # convert grid to a vector and sort in decreasing order to find a cutoff point where the proportion of cells with values 
  # above the cutoff point is equal to 'prop.burned'
  gfv <- sort(as.vector(gf), decreasing = TRUE) #gfv = gaussian field vector; length = n cells in grid
  cutoff <- #prop.burned percent cells of total grid, everything larger is burned (when sorted by descending value)
    gfv[length(gfv) * prop.burned] 
  gf.binary <- calc(gf, fun = convert.to.binary)
  gf.binary1 <- mask(gf.binary,perim)
  gf.binary1[gf.binary != 1] <- NA
  writeRaster(gf.binary1,"./SpatialOutput/sim.tif", overwrite = TRUE)
  
  #Polygonize using "stars" from https://cran.r-project.org/web/packages/stars/vignettes/stars5.html
  hs_fire = read_stars("./SpatialOutput/sim.tif")
  hs_fire = st_as_sf(hs_fire, as_points = FALSE, merge = TRUE)
  
  #Now proceed with the polygon cleaning and cutting operations
  hs_fire$Area_m2 <- #Calculate area of each patch; assumes metric
    st_area(hs_fire)
  hs_fire$Area_ha <- #Calculate area in ha of each patch
    as.vector(hs_fire$Area_m2)*0.0001
  hs_fill <- #Fill holes <= 1 ha; delete crumbs <= 1 ha.
    donut_holes_and_crumbs(hs_fire)
  hs_fill$PatchName <-
    seq(1:nrow(hs_fill))
  hs_fill$RemoveLater <- NA #Placeholder for patches that are eventually split and need to be removed.
  hs_pinch <- #drop crumbs <= 1,000,000 m2 = 100ha, minimum crumb size retained is 10,001 m2
    #here we are looking at which large (>100 ha) polygons have pinch points that we can chop up
    hs_fill[hs_fill$Area_ha > 100,]
  hs_fill_small <- #Store the polygons <= 100 ha, for subsequent re-integration
    hs_fill[hs_fill$Area_ha <= 100,]
  
  print(Sys.time())
  sim_clean <- #Conduct splits via "breath increments" and "lungs" functions. SLOW.
    breath_increments(hs_pinch,increments = c(-15, -30, -60, -120))
  sim_clean_full <- #Combine small unsplit polygons (hs_fill_small) with larger split polygons
    rbind(sim_clean, hs_fill_small)
  print(Sys.time()); print(paste("finished sim ", s))
  
  sim_clean_full$Area_ha #return the vector of patch size values
  
}

random_cheese_spatial <- function(pct,nugget,magvar,pct_hs){
  rng <- round(sqrt(((pct/100) * xmax * ymax)/pi), 0)
  
  # Alternatively, set the variogram range in number of cells
  # rng <- 20
  
  # Generate a continuous Gaussian random field based on the variogram input parameters
  gf <- #Slow, ~20 seconds w 30 m resolution
    nlm_gaussianfield(ncol = xmax, nrow = ymax, autocorr_range = rng, mag_var = magvar, nug = nugget)
  
  #Reproject the random grid to be on the same CRS as Las Conchas perimeter.
  extent(gf) <- round(extent(perim),0)
  gf@crs <- perim_sp@proj4string
  gf <- projectRaster(gf, crs = crs(perim))
  
  #### Convert the continuous Gaussian field into a binary grid ###
  # Function to convert continuous Gaussian random field to binary data
  # The binary grid will have Values of 1 = burned and 0 = unburned
  convert.to.binary <- function(x) 
  {
    ifelse(x <  cutoff, 0, 1)
  }
  
  prop.burned = pct_hs
  
  # convert grid to a vector and sort in decreasing order to find a cutoff point where the proportion of cells with values 
  # above the cutoff point is equal to 'prop.burned'
  gfv <- sort(as.vector(gf), decreasing = TRUE) #gfv = gaussian field vector; length = n cells in grid
  cutoff <- #prop.burned percent cells of total grid, everything larger is burned (when sorted by descending value)
    gfv[length(gfv) * prop.burned] 
  gf.binary <- calc(gf, fun = convert.to.binary)
  gf.binary1 <- mask(gf.binary,perim)
  gf.binary1[gf.binary != 1] <- NA
  writeRaster(gf.binary1,"./SpatialOutput/sim.tif", overwrite = TRUE)
  
  #Polygonize using "stars" from https://cran.r-project.org/web/packages/stars/vignettes/stars5.html
  hs_fire = read_stars("./SpatialOutput/sim.tif")
  hs_fire = st_as_sf(hs_fire, as_points = FALSE, merge = TRUE)
  
  #Now proceed with the polygon cleaning and cutting operations
  hs_fire$Area_m2 <- #Calculate area of each patch; assumes metric
    st_area(hs_fire)
  hs_fire$Area_ha <- #Calculate area in ha of each patch
    as.vector(hs_fire$Area_m2)*0.0001
  hs_fill <- #Fill holes <= 1 ha; delete crumbs <= 1 ha.
    donut_holes_and_crumbs(hs_fire)
  hs_fill$PatchName <-
    seq(1:nrow(hs_fill))
  hs_fill$RemoveLater <- NA #Placeholder for patches that are eventually split and need to be removed.
  hs_pinch <- #drop crumbs <= 1,000,000 m2 = 100ha, minimum crumb size retained is 10,001 m2
    #here we are looking at which large (>100 ha) polygons have pinch points that we can chop up
    hs_fill[hs_fill$Area_ha > 100,]
  hs_fill_small <- #Store the polygons <= 100 ha, for subsequent re-integration
    hs_fill[hs_fill$Area_ha <= 100,]
  
  Sys.time()
  sim_clean <- #Conduct splits via "breath increments" and "lungs" functions. SLOW.
    breath_increments(hs_pinch,increments = c(-15, -30, -60, -120))
  sim_clean_full <- #Combine small unsplit polygons (hs_fill_small) with larger split polygons
    rbind(sim_clean, hs_fill_small)
  Sys.time()
  
  sim_clean_full #return the vector of patch size values ONLY CHANGE from similar function above; should adjust that one.
  
}
