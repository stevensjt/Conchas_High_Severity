library(sf) #for st_as_sf(); version 0.7-3
library(rgdal) #for readOGR(); version 1.3-6
library(raster)
library(NLMR)
library(stars)
#This will replace Random Cheese

####1. Read in the three observed layers####
lc <- read_sf("./SpatialOutput/", layer = "lc_clean_cut")
#lcdc <- read_sf("./SpatialOutput/", layer = "LasConchas_Dome_Cerro_643_RdNBR_Clean")
#treeless <- read_sf("./SpatialOutput/", layer = "LasConchas_Treeless_Clean")
perim <- read_sf("../../../GIS/Fires/Las Conchas/Perimeter/", layer="Las_Conchas_Perim")
#perim_sp <- readOGR("../../../GIS/Fires/Las Conchas/Perimeter/", layer="Las_Conchas_Perim")
scars <- read_sf("../../../GIS/Fires/Las Conchas/FireScarLocations/", layer="fs_in_las_conchas")
#scars_sp <- readOGR("../../../GIS/Fires/Las Conchas/FireScarLocations/", layer="fs_in_las_conchas")

#CHECKME need to run stats, can pull from Random Cheese
#sum(st_area(lc))/st_area(perim)

#### 2.2 Generate spatially continuous Gaussian random field: EXAMPLE####
#Set the grid size
# xmax and ymax set the number of cells in the x and y direction
# increasing the grid size slows down the 'nlm_gaussianfied' function
xmax <- #30 converts from UTM 1m resolution to LANDSAT 30 m resolution
  round((extent(perim)@xmax - extent(perim)@xmin)/30,0)[[1]] #Grid is 964 pixels wide
ymax <- 
  round((extent(perim)@ymax - extent(perim)@ymin)/30,0)[[1]] #Grid is 1447 pixels tall
#Grid is 964*1447 = 1394908 pixels large = (*900 m2) 1255417200 m2 large = 125542 ha 1255 km2

#set variogram parameters
nugget <- 0
magvar <- 100

# set the variogram range as the radius of a circle with an area equal to 'pct' percent of the grid area
# increasing the variogram range will make the 'nlm_gaussianfield' function run slower
pct <- 1 #Key parameter to vary autocorrelation (JTS)
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

#par(mar = c(0.1, 0.1, 0.9, 0.1)) #Shrink margins for efficient use of plot space
plot(gf, col = terrain.colors(25))
plot(perim_sp, add=TRUE)


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
plot(perim_sp, add = TRUE)
gf.binary.mask <- mask(gf.binary,perim)
plot(gf.binary.mask, col = c("white",color_cutoff))
plot(perim_sp, add = TRUE)


####2.3 Calculate patch size distribution using sf/stars; EXAMPLE####
# These steps will create a data frame called 'patch.size' with columns "value" and "count"
# "value" is the patch number (arbitrary values)
# "count" is the number of cells in the patch, which will be sorted from largest to smallest
# unburned patches (cell values of 0 in the binary grid) are excluded from the patch.size data frame

# copy the binary grid and convert the zeros to NAs  
gf.binary1 <- gf.binary.mask #Specify whether to calculate clipped to perim ("gf.binary.mask") or not("gf.binary").
gf.binary1[gf.binary != 1] <- NA
writeRaster(gf.binary1,"./SpatialOutput/sim.tif")

#Polygonize using "stars" from https://cran.r-project.org/web/packages/stars/vignettes/stars5.html
x = read_stars("./SpatialOutput/sim.tif")
hs_fire = st_as_sf(x, as_points = FALSE, merge = TRUE)
plot(hs_fire,add=T,col=rainbow(8))

hs_fire$Area_m2 <- #Calculate area of each patch; assumes metric
  st_area(hs_fire)
hs_fire$Area_ha <- #Calculate area in ha of each patch
  as.vector(hs_fire$Area_m2)*0.0001

hs_fill <- #Fill holes <= 1 ha; delete crumbs <= 1 ha.
  donut_holes_and_crumbs(hs_fire)
hs_fill$PatchName <-
  seq(1:nrow(hs_fill))
hs_fill$RemoveLater <- NA #Placeholder for patches that are eventually split and need to be removed.


####2. split polygon chunks separated by narrow bridges####
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

####3. Simulate random patches to dial optimal parameter configuration####
##3.1: Initialize data frame with observed data
#e2 <- function(x){2^x}
#b <- c(1,e2(c(1:13)))
#h <- hist(lc$Area_ha, breaks = b) #deprecate
#Histogram of original patch distribution
og_hist = 
  ggplot(lc, aes(x = log2(Area_ha))) +
  geom_histogram(binwidth = 0.5) +
  #scale_x_log10(breaks = c(1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192)) + 
  theme(axis.text.x = element_text(angle = 45, size = 14))
og_data <- ggplot_build(og_hist)
og_data$data[[1]]$y
og_data$data[[1]]$x

df_sim <- #Set up data frame with original patch stats to compare sims against
  data.frame(sim_num = 0, bins = og_data$data[[1]]$x, 
             original_patches = og_data$data[[1]]$y)


#Set the grid size
# xmax and ymax set the number of cells in the x and y direction
# increasing the grid size slows down the 'nlm_gaussianfied' function
xmax <- #30 converts from UTM 1m resolution to LANDSAT 30 m resolution
  round((extent(perim)@xmax - extent(perim)@xmin)/30,0)[[1]] #Grid is 964 pixels wide
ymax <- 
  round((extent(perim)@ymax - extent(perim)@ymin)/30,0)[[1]] #Grid is 1447 pixels tall
#Grid is 964*1447 = 1394908 pixels large = (*900 m2) 1255417200 m2 large = 125542 ha 1255 km2

####3.1 Simulation 1: This is the original standard against which other parameter values will be compared####
nugget <- 0
magvar <- 100
pct <- 1
pct_hs <- 0.424
#sum(st_area(lc))/st_area(perim) #42.4% high severity for Las Conchas RdNBR


for(s in c(1:10)){
  out = random_cheese_patchdist(pct = 1, nugget = 0, magvar = 100, pct_hs = 0.424)
  sim_hist = 
    ggplot(data = data.frame(out), aes(x = log2(out))) +
    geom_histogram(binwidth = 0.5) +
    #scale_x_log10(breaks = c(1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192)) + 
    theme(axis.text.x = element_text(angle = 45, size = 14))
  sim_data <- ggplot_build(sim_hist)
  y = sim_data$data[[1]]$y
  x = sim_data$data[[1]]$x
  
  rows = c((nrow(df_sim)+1):(nrow(df_sim) + length(y)))
  df_sim[rows,"sim_num"] = rep(s, length(rows))
  df_sim[rows,"bins"] = x
  df_sim[rows,"original_patches"] = y
}

#Compare patch distribution in Sim 1 vs observed
p_resid =
  ggplot(data = df_sim[df_sim$sim_num == 0,], aes(x = bins, y = original_patches)) + 
  geom_point() + 
  geom_smooth()
p_resid_data <- ggplot_build(p_resid)
#View(p_resid_data$data[[1]])
#View(p_resid_data$data[[2]])
obs_values <- p_resid_data$data[[1]]
pred_values <- p_resid_data$data[[2]]

p_full <-
  p_resid +
  geom_point(data = df_sim[df_sim$sim_num!=0,], aes(x = bins, y = original_patches, col = sim_num)) #+
  #geom_smooth(data = df_sim[df_sim$sim_num!=0,], aes(x = bins, y = original_patches)) #Option to smooth the simulated data
p_full_data <- ggplot_build(p_full)
View(p_full_data$data[[3]])

resid_data_long <- p_full_data$data[[3]][,c("x","y")]
resid_data_long$y_obs <- #get the actual patch count for this bin from the observed fire in question
  p_resid_data$data[[1]][pmatch(resid_data_long$x,obs_values$x, duplicates.ok = T),"y"]
resid_data_long$y_obs[which(is.na(resid_data_long$y_obs))] = 0 #Patches larger than largest observed patch get predicted abundance value of 0
resid_data_long$y_pred = 0

#get the predicted patch count for this bin from the loess smooth.
for(i in 1:nrow(resid_data_long)){
  index <- which.min(abs(pred_values$x-resid_data_long$x[i]) ) #identify the nearest x value.
  if(resid_data_long$x[i]<= max(pred_values$x)){
    #We only replace 0 values for patch sizes that are less than the maximum patch in the observed fire. 
    #Otherwise the predicted value stays 0.
    print(paste("replacing row", i))
    resid_data_long$y_pred[i] <- pred_values$y[index]
  }
}

resid_data_long$resid_obs = resid_data_long$y - resid_data_long$y_obs 
resid_data_long$resid_pred = resid_data_long$y - resid_data_long$y_pred 
ggplot(resid_data_long)+
  geom_point(aes(x = x, y = resid_obs))
#START HERE test the slope of this line, should be intercept 0 and slope 0.

####3.2 Tuning Simulations: nugget = 0, magvar = 100, what is optimal pct?####
#TAKES A LONG TIME FOR Nug = 10, magvar = 100, mct = 1
#And for Nug = 1, magvar = 1, pct = 1
nugget <- 10 #0 5 10 
magvar <- 100 #0/1 50 100
pct <- 1 #0.3, 0.6, 1, 1.3
pct_hs <- 0.424
#sum(st_area(lc))/st_area(perim) #42.4% high severity for Las Conchas RdNBR

df_sim <- #Set up data frame with original patch stats to compare sims against
  data.frame(sim_num = 0, bins = og_data$data[[1]]$x, 
             original_patches = og_data$data[[1]]$y)

for(s in c(3:10)){ #about 20 mins to run 10 tuning sims
  out = random_cheese_patchdist(pct = pct, nugget = nugget, magvar = magvar, pct_hs = pct_hs)
  sim_hist = 
    ggplot(data = data.frame(out), aes(x = log2(out))) +
    geom_histogram(binwidth = 0.5) +
    #scale_x_log10(breaks = c(1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192)) + 
    theme(axis.text.x = element_text(angle = 45, size = 14))
  sim_data <- ggplot_build(sim_hist)
  y = sim_data$data[[1]]$y
  x = sim_data$data[[1]]$x
  
  rows = c((nrow(df_sim)+1):(nrow(df_sim) + length(y)))
  df_sim[rows,"sim_num"] = rep(s, length(rows))
  df_sim[rows,"bins"] = x
  df_sim[rows,"original_patches"] = y
}

#Write df_sim to csv for the given set of parameters. TUNING SIMULATION OUTPUT
tuning_sim_pathname <- paste0("./Output/Tuning_Sims/sim_pct_",pct,"_nug_",nugget,"_mv_",magvar,".csv")
write.csv(df_sim,tuning_sim_pathname)

####3.3 Take tuning simulations and figure out the optimal pct####
#Read the files
df_sim_list <- list(); d <- "./Output/Tuning_Sims/"
tuning_test <- data.frame(pct = NA, nug = NA, mv = NA)

for(file in dir(d)){

df_sim <- df_sim_list[[file]] <- read.csv(paste0(d,file))

#Compare patch distribution in Sim 1 vs observed
p_resid =
  ggplot(data = df_sim[df_sim$sim_num == 0,], aes(x = bins, y = original_patches)) + 
  geom_point() + 
  geom_smooth()
p_resid_data <- ggplot_build(p_resid)
#View(p_resid_data$data[[1]])
#View(p_resid_data$data[[2]])
obs_values <- p_resid_data$data[[1]]
pred_values <- p_resid_data$data[[2]]

p_full <-
  p_resid +
  geom_point(data = df_sim[df_sim$sim_num!=0,], aes(x = bins, y = original_patches, col = sim_num)) #+
#geom_smooth(data = df_sim[df_sim$sim_num!=0,], aes(x = bins, y = original_patches)) #Option to smooth the simulated data
p_full_data <- ggplot_build(p_full)
#View(p_full_data$data[[3]])

resid_data_long <- p_full_data$data[[3]][,c("x","y")]
resid_data_long$y_obs <- #get the actual patch count for this bin from the observed fire in question
  p_resid_data$data[[1]][pmatch(resid_data_long$x,obs_values$x, duplicates.ok = T),"y"]
resid_data_long$y_obs[which(is.na(resid_data_long$y_obs))] = 0 #Patches larger than largest observed patch get predicted abundance value of 0
resid_data_long$y_pred = 0

#get the predicted patch count for this bin from the loess smooth.
for(i in 1:nrow(resid_data_long)){
  index <- which.min(abs(pred_values$x-resid_data_long$x[i]) ) #identify the nearest x value.
  if(resid_data_long$x[i]<= max(pred_values$x)){
    #We only replace 0 values for patch sizes that are less than the maximum patch in the observed fire. 
    #Otherwise the predicted value stays 0.
    print(paste("replacing row", i))
    resid_data_long$y_pred[i] <- pred_values$y[index]
  }
}

resid_data_long$resid_obs <- resid_data_long$y - resid_data_long$y_obs 
resid_data_long$resid_pred <- resid_data_long$y - resid_data_long$y_pred 
ggplot(resid_data_long)+
  geom_point(aes(x = x, y = resid_obs))
m <- glm(resid_obs ~ x, data= resid_data_long)

#Populate tuning test table
tuning_test[file,"pct"] <- strsplit(file,c("_"))[[1]][3]
tuning_test[file,"nug"] <- strsplit(file,c("_"))[[1]][5]
tuning_test[file,"mv"] <- strsplit(file,c("_"))[[1]][7]
tuning_test[file,"resid_slope"] <- summary(m)$coefficients[2,1] #coef
tuning_test[file,"resid_pval"] <- round(summary(m)$coefficients[2,4],3) #P value

}
write.csv(tuning_test[-1,],paste0(d,"tuning_test.csv")) #minus 1 because first row is NA
#See README in the Tuning_Sims/Beta folder for lessons learned from the dry run about initial parameterization.

####4. Simulate random patches to look for intersection with scars####
hits_df <- #set up data frame to track scar hits among sims
  data.frame(sim = NA, n_hits = NA, pct_hits = NA) 
#Set the grid size
xmax <- #30 converts from UTM 1m resolution to LANDSAT 30 m resolution
  round((extent(perim)@xmax - extent(perim)@xmin)/30,0)[[1]] #Grid is 964 pixels wide
ymax <- 
  round((extent(perim)@ymax - extent(perim)@ymin)/30,0)[[1]] #Grid is 1447 pixels tall


for(s in c(1:10)){
  out = random_cheese_spatial(pct = 1, nugget = 0, magvar = 100, pct_hs = 0.424)
  print("cutting of random cheese complete for sim", s)
  out = st_transform(out, crs = st_crs(scars)) #at some point the crs didn't transfer to the simulation so we re-assign
  hits = st_intersects(scars, out)
  hits = which(apply(hits, 1, any)) #vector of which scars intersected HS patches
  
  #option to plot
  #scar_hits_sim = scars #option to plot
  #scar_hits_sim$hit = "no"
  #scar_hits_sim$hit[hits] = "yes"
  #plot(out$geometry, col = "lightblue", border = "lightblue", main = paste0(round(length(hits)/nrow(scars),2)*100, "% scars hit a patch"))
  #plot(perim$geometry,add=T)
  #plot(scar_hits_sim[,"hit"], add=T, pch = 20)
  
  hits_df[s,"sim"] <- s
  hits_df[s,"n_hits"] <- length(hits)
  hits_df[s,"pct_hits"] = round(100* length(hits)/nrow(scars),2)
  print(paste("finished sim", s))
  
  
}
Sys.time()
#hist(hits_df)
hits_df_hist = 
  ggplot(data = hits_df, aes(x = log2(out))) +
  geom_histogram(binwidth = 0.5) +
  #scale_x_log10(breaks = c(1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192)) + 
  theme(axis.text.x = element_text(angle = 45, size = 14))