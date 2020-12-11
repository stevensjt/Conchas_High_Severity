#Generate patches function (from Alan Tepley)

generate_patches <- function(xmax, ymax, nugget, magvar, pct=NA, rng=NA, perim, prop.burned){
  #if using the percentage approach:
  if(!is.na(pct)){
    rng <- round(sqrt(((pct/100) * xmax * ymax)/pi), 0)
  }

  # Generate a continuous Gaussian random field based on the variogram input parameters
  gf <- #Slow, ~20 seconds w 30 m resolution
    nlm_gaussianfield(ncol = xmax, nrow = ymax, autocorr_range = rng, mag_var = magvar, nug = nugget)
  
  #Reproject the random grid to be on the same CRS as Las Conchas perimeter.
  extent(gf) <- round(extent(perim),0)
  gf@crs <- perim@proj4string
  gf <- projectRaster(gf, crs = crs(perim))
  
  #### Convert the continuous Gaussian field into a binary grid
  
  # Function to convert continuous Gaussian random field to binary data
  # The binary grid will have Values of 1 = burned and 0 = unburned
  
  convert.to.binary <- function(x) 
  {
    ifelse(x <  cutoff, 0, 1)
  }

  # convert grid to a vector and sort in decreasing order to find a cutoff point 
  # where the proportion of cells with values 
  # above the cutoff point is equal to 'prop.burned'
  gfv <- #length = n cells in grid
    sort(as.vector(gf), decreasing = TRUE) 
  cutoff <- #prop.burned percent cells of total grid;
    # everything larger is burned (when sorted by descending value)
    gfv[length(gfv) * prop.burned] 
  gf.binary <- #Everything above cutoff value is converted to 1
    calc(gf, fun = convert.to.binary)
  
  gf.binary.mask <- #Mask the binary grid by the fire perimeter
    mask(gf.binary,perim)

  # copy the binary grid and convert the zeros to NAs  
  gf.binary.mask[gf.binary.mask != 1] <- 
    NA
  return(gf.binary.mask)
}