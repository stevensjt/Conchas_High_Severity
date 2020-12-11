#Jens Stevens; stevensjt@gmail.com
#Code modified from Alan Tepley

####1. Load libraries####
library(NLMR)
library(raster)
library(igraph)
library(rgdal) #for readOGR(); version 1.3-6
library(sf) #for st_as_sf(); version 0.7-3
library(smoothr)
library(lwgeom)

source("./Code/SpatialFunctions.R")
source("./Code/generate_patches.R")
par(mar = c(0.1, 0.1, 0.9, 0.1)) #Shrink margins for efficient use of plot space


####2. Process Las Conchas layers for comparison####
##2a Load Las Conchas specific files
# All spatial projections unless otherwise noted: EPSG:5070 - NAD83 / Conus Albers - Projected
perim <- readOGR("../../../GIS/Fires/Las Conchas/Perimeter/", layer="Las_Conchas_Perim")
scars <- readOGR("../../../GIS/Fires/Las Conchas/FireScarLocations/", layer="fs_in_las_conchas")

hs_lc <- #Las Conchas burn severity with 643 RdNBR threshold, from Sean Parks
  readOGR("../../../GIS/Fires/Las Conchas/Severity", layer="LasConchas_643_RdNBR")
  #Better for plotting holes in ggplot framework:
  #st_read("../../../GIS/Fires/Las Conchas/Severity", layer="LasConchas_643_RdNBR")
#Process hs_lc:
hs_lc@data <- data.frame(Name = "NA", Area_m2 = NA, Area_ha = NA)

hs_treeless <- #Cumulative treeless area over multiple fires, from Jon Coop
  readOGR("../../../GIS/Fires/Las Conchas/Severity", layer="LasConchas_Treeless")
  #Better for plotting holes in ggplot framework:
  #st_read("../../../GIS/Fires/Las Conchas/Severity", layer="LasConchas_Treeless")

##2b. Fill holes of low-severity that are less than 9 pixels (8100m2, 0.81 ha) large
lc_fill <-#Fill holes of low-severity that are less than 9 pixels (8100m2, 0.81 ha) large
  fill_holes(hs_lc) 
treeless_fill <-#Fill holes of low-severity that are less than 9 pixels (8100m2, 0.81 ha) large
  fill_holes(hs_treeless) 
gf.binary1[gf.binary != 1] <- NA


writeOGR(treeless_fill, "../../../GIS/Fires/Las Conchas/Severity", layer="LasConchas_Treeless_fill",driver="ESRI Shapefile")

####3c. Rasterize####
r_fill <- rasterize(treeless_fill, gf.binary.mask)
#Stats:
#Area of Las Conchas footprint:
lc.area <- area(perim)*0.0001 #61056.63 ha
(length(which(!is.na(getValues(r_fill)))) *0.09) / #Area in HS, in ha (1 pixel = 0.09 ha)
  lc.area
#0.7624347 percent high severity

####3d. Get stats for "treeless" layer####
clumps <- clump(r_fill, directions = 8)
# Exclude NAs and sort patch sizes from largest to smallest
patch.size <- na.omit(data.frame(freq(clumps)))
patch.size <- patch.size[order(-patch.size$count),]
# calculate patch sizes as a proportion of the grid area
patch.size$prop <- (patch.size$count / (xmax[i] * ymax[i]))
#head(patch.size)
#hist(patch.size$prop[c(1:20)])

# check if the patch sizes sum to the desired proportion of the landscape burned (value specified earlier)
treeless_output <- list()
treeless_output$obs_pct_hs[1] <- 
  sum(patch.size$count)*900 / #puts in m2 units to compare to area(perim)
  #If only considering patches within the LC perimeter, this value represents the
  #proportion of Las Conchas perimeter in high-severity (should be ~= prop.burned)
  area(perim)
treeless_output$obs_pct_hs[1] <- round(treeless_output$obs_pct_hs[1],2)
# prop.burned #original prop.burned parameter

#Percent of scars intersected
scar_hits <- extract(r_fill,scars)
treeless_output$pct_scars_burned_hs[1] <- 
  round(sum(scar_hits, na.rm = TRUE)/length(scar_hits), 2)

#output$plot[i] <-  #eventually use ggplot object
plot(r_fill, col = "#df5757", legend = FALSE)
plot(perim, add=TRUE, border = "darkgray")
plot(scars, add=TRUE, pch = ".")
title(paste(
  "Actual Las Conchas footprint- treeless", 
  "\n prop treeless, obs =", treeless_output$obs_pct_hs[1], 
  "\n proportion of scars burned @ HS = ", treeless_output$pct_scars_burned_hs[1]),
  cex.main = 0.9
)
dev.copy2pdf(file = "Las Conchas Treeless.pdf")


#### 2.1 Create the simulated landscape###
#Set the grid size
# xmax and ymax set the number of cells in the x and y direction
# increasing the grid size slows down the 'nlm_gaussianfied' function
xmax <- #30 converts from UTM 1m resolution to LANDSAT 30 m resolution
  round((extent(perim)@xmax - extent(perim)@xmin)/30,0) #Grid is 964 pixels wide
ymax <- 
  round((extent(perim)@ymax - extent(perim)@ymin)/30,0) #Grid is 1447 pixels tall
#Grid is 964*1447 = 1394908 pixels large = (*900 m2) 1255417200 m2 large = 125542 ha 1255 km2


#### 2.2 Generate spatially continuous Gaussian random field####
#set variogram parameters
nugget <- 0
magvar <- 100

# set the variogram range as the radius of a circle with an area equal to 'pct' percent of the grid area
# increasing the variogram range will make the 'nlm_gaussianfield' function run slower
pct <- 10 #Key parameter to vary autocorrelation (JTS)
rng <- round(sqrt(((pct/100) * xmax * ymax)/pi), 0)

# Alternatively, set the variogram range in number of cells
# rng <- 20

# Generate a continuous Gaussian random field based on the variogram input parameters
gf <- #Slow, ~20 seconds w 30 m resolution
  nlm_gaussianfield(ncol = xmax, nrow = ymax, autocorr_range = rng, mag_var = magvar, nug = nugget)

#Reproject the random grid to be on the same CRS as Las Conchas perimeter.
extent(gf) <- round(extent(perim),0)
gf@crs <- perim@proj4string
gf <- projectRaster(gf, crs = crs(perim))

plot(gf, col = terrain.colors(25))
plot(perim, add=TRUE)


#### Convert the continuous Gaussian field into a binary grid ###

# Function to convert continuous Gaussian random field to binary data
# The binary grid will have Values of 1 = burned and 0 = unburned

convert.to.binary <- function(x) 
{
  ifelse(x <  cutoff, 0, 1)
}

# set the desired proportion burned
prop.burned <- 0.20

# convert grid to a vector and sort in decreasing order to find a cutoff point where the proportion of cells with values 
# above the cutoff point is equal to 'prop.burned'
gfv <- sort(as.vector(gf), decreasing = TRUE) #length = n cells in grid
cutoff <- gfv[length(gfv) * prop.burned] #prop.burned percent cells of total grid, everything larger is burned (when sorted by descending value)
gf.binary <- calc(gf, fun = convert.to.binary)

#Visualize binary output
bin_out <- na.exclude(sample(getValues(gf),10000)) #Gaussian output without NA's
#plot(ecdf(bin_out),main = "cumulative density")
#abline(h = 1-prop.burned, v = quantile(bin_out,probs = (1-prop.burned)))
color_cutoff <- #Color value associated with the prop.burned color
  terrain.colors(25)[round(25*quantile(bin_out,probs = (1-prop.burned)))]

plot(gf.binary, col = c("white",color_cutoff)) 
plot(perim, add = TRUE)
gf.binary.mask <- mask(gf.binary,perim)
plot(gf.binary.mask, col = c("white",color_cutoff))
plot(perim, add = TRUE)


### Calculate patch size distribution ###

# These steps will create a data frame called 'patch.size' with columns "value" and "count"
# "value" is the patch number (arbitrary values)
# "count" is the number of cells in the patch, which will be sorted from largest to smallest
# unburned patches (cell values of 0 in the binary grid) are excluded from the patch.size data frame

# copy the binary grid and convert the zeros to NAs  
gf.binary1 <- gf.binary.mask #Specify whether to calculate clipped to perim ("gf.binary.mask") or not("gf.binary").
gf.binary1[gf.binary != 1] <- NA

# group all contiguous pixels with a value of 1 into clumps based on neighboring cells in 8 directions
# set "directions = 4" to generate clumps only using neighboring cells in the 4 cardinal directions
clumps <- clump(gf.binary1, directions = 8)

# Exclude NAs and sort patch sizes from largest to smallest
patch.size <- na.omit(data.frame(freq(clumps)))
patch.size <- patch.size[order(-patch.size$count),]

# calculate patch sizes as a proportion of the grid area
patch.size$prop <- (patch.size$count / (xmax * ymax))

# check if the patch sizes sum to the desired proportion of the landscape burned (value specified earlier)
sum(patch.size$prop) #Proportion of entire simulated landscape in high-severity
sum(patch.size$count)*900 / #puts in m2 units to compare to area(perim)
  #If only considering patches within the LC perimeter, this value represents the
  #proportion of Las Conchas perimeter in high-severity (should be ~= prop.burned)
  area(perim) 
prop.burned

#head(patch.size)
#hist(patch.size$prop[c(1:20)])
hist(patch.size$count[c(1:20)])

#START HERE; raster to vector and then do the overlapping scars piece.
#Skipping this now
#hs_lc <- rasterToPolygons(gf.binary1) #large and unwieldy... find an alternative.
#lc_fill <-#Fill holes of low-severity that are less than 9 pixels (8100m2, 0.81 ha) large
#  fill_holes(hs_lc) 

####3. Compare random to Conchas####


####3e. Compare: Set up parameters####
(length(which(!is.na(getValues(gf.binary1)))) *0.09) / #Area in HS, in HA
  lc.area #only 17% vs 76%

xmax <- #30 converts from UTM 1m resolution to LANDSAT 30 m resolution
  round((extent(perim)@xmax - extent(perim)@xmin)/30,0) #Grid is 964 pixels wide
ymax <- 
  round((extent(perim)@ymax - extent(perim)@ymin)/30,0) #Grid is 1447 pixels tall
nugget <- 0
magvar <- 100
pct <- 1 #Key parameter to vary autocorrelation (JTS)
rng <- round(sqrt(((pct/100) * xmax * ymax)/pi), 0)
perim <- perim
prop.burned <- (length(which(!is.na(getValues(r_fill)))) *0.09) / #Area in HS, in HA
  lc.area

pars <- list(
  xmax = rep(xmax,16),
  ymax = rep(ymax,16),
  nugget = c(0,0,0,0,0,0,0,0,10,10,10,10,10,10,10,10),
  magvar = rep(c(10,10,10,10,100,100,100,100), 2),
  rng = rep(c(50,50,100,100),4),
  prop.burned = rep(c(0.75,0.25),8)
)
sims <- list()
i <- 1
output <- list(
  obs_pct_hs <- NA,
  pct_scars_burned_hs <- NA,
  plot <- NA
)
#dev.off() #in case previous graphics device had been set up
par(mar=c(1.5, 1.5, 4, 0))

####3f: Run loop of simulations with different configurations####
for(i in 1:1){
  sims[[i]] <- 
    generate_patches(xmax=pars$xmax[i], ymax=pars$ymax[i], nugget=pars$nugget[i], magvar=pars$magvar[i],
                     pct=NA, rng=pars$rng[i], perim=perim, prop.burned=pars$prop.burned[i])
  
  clumps <- clump(sims[[i]], directions = 8)
  # Exclude NAs and sort patch sizes from largest to smallest
  patch.size <- na.omit(data.frame(freq(clumps)))
  patch.size <- patch.size[order(-patch.size$count),]
  # calculate patch sizes as a proportion of the grid area
  patch.size$prop <- (patch.size$count / (xmax[i] * ymax[i]))
  #head(patch.size)
  #hist(patch.size$prop[c(1:20)])
  
  # check if the patch sizes sum to the desired proportion of the landscape burned (value specified earlier)
  output$obs_pct_hs[i] <- 
    sum(patch.size$count)*900 / #puts in m2 units to compare to area(perim)
    #If only considering patches within the LC perimeter, this value represents the
    #proportion of Las Conchas perimeter in high-severity (should be ~= prop.burned)
    area(perim)
  output$obs_pct_hs[i] <- round(output$obs_pct_hs[i],2)
  # prop.burned #original prop.burned parameter
  
  #Percent of scars intersected
  scar_hits <- extract(sims[[i]],scars)
  output$pct_scars_burned_hs[i] <- 
    round(sum(scar_hits, na.rm = TRUE)/length(scar_hits), 2)
  
  #output$plot[i] <-  #eventually use ggplot object
  plot(sims[[i]], col = "darksalmon", legend = FALSE)
  plot(perim, add=TRUE, border = "darkgray")
  plot(scars, add=TRUE, pch = ".")
  title(paste(
    "magnitude of autocorrelation variation =", pars$nugget[i],
    "\n magnitude of landscape variation = ", pars$magvar[i],
    "\n max range of spatial autocorrelation (px) = ", pars$rng[i],
    "\n prop high sev, obs =", output$obs_pct_hs[i], 
    "; prop high sev, set =", pars$prop.burned[i],
    "\n proportion of scars burned @ HS = ", output$pct_scars_burned_hs[i]),
    cex.main = 0.6
  )
  dev.copy2pdf(file = paste0("plot",i,".pdf"))
}

####3g: Run scenario 1 iteratively####
#Using scenario 1, run 10 times for each percent high sev from 75% down to 5%, look at scar hits distribution

pars <- list(
  xmax = xmax,
  ymax = xmax,
  nugget = rep(0,8),
  magvar = rep(0,8),
  rng = rep(50,8),
  prop.burned = seq(from = 75, to = 5, by = -10)
)
sims <- list()
i <- 1
output <- list(
  obs_pct_hs <- vector(),
  pct_scars_burned_hs <- NA,
  plot <- NA
)
for(i in 1:2){
  for(n in 1:3){
    sim <- generate_patches(xmax=pars$xmax, ymax=pars$ymax, nugget=pars$nugget[i], magvar=pars$magvar[i],
                            pct=NA, rng=pars$rng[i], perim=perim, prop.burned=pars$prop.burned[i])
    clumps <- clump(sims, directions = 8)
    patch.size <- na.omit(data.frame(freq(clumps)))
    patch.size <- patch.size[order(-patch.size$count),]
    output$obs_pct_hs[i][n] <- 
      sum(patch.size$count)*900 / area(perim)
    scar_hits <- extract(sims[[i]],scars)
    output$pct_scars_burned_hs[i] <- 
      round(sum(scar_hits, na.rm = TRUE)/length(scar_hits), 2)
  }
} #START HERE and figure out how to incorporate multiple sims into a vector.
  
    


  #Percent of scars intersected
  scar_hits <- extract(sims[[i]],scars)
  output$pct_scars_burned_hs[i] <- 
    round(sum(scar_hits, na.rm = TRUE)/length(scar_hits), 2)
  
  #output$plot[i] <-  #eventually use ggplot object
  plot(sims[[i]], col = "darksalmon", legend = FALSE)
  plot(perim, add=TRUE, border = "darkgray")
  plot(scars, add=TRUE, pch = ".")
  title(paste(
    "magnitude of autocorrelation variation =", pars$nugget[i],
    "\n magnitude of landscape variation = ", pars$magvar[i],
    "\n max range of spatial autocorrelation (px) = ", pars$rng[i],
    "\n prop high sev, obs =", output$obs_pct_hs[i], 
    "; prop high sev, set =", pars$prop.burned[i],
    "\n proportion of scars burned @ HS = ", output$pct_scars_burned_hs[i]),
    cex.main = 0.6
  )
  dev.copy2pdf(file = paste0("plot",i,".pdf"))
#}
