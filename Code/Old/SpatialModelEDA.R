###Swiss cheese exploratory analysis
###Jens Stevens; jtstevens@usgs.gov

####0. Load libraries and other code####
library(rgdal) #for readOGR(); version 1.3-6
library(raster) #for extent(); version 2.8-4
library(rgeos) #for gBuffer(); version 0.4-2
library(pgirmess) #for MergeTrackObs(); version 1.6.9
library(tidyverse) #for ggplot() etc; version 1.3.0
library(sf) #for st_as_sf(); version 0.7-3
library(NLMR) #for random patches
library(igraph) #for random patches
source("./Code/SpatialFunctions.R")

####1. Read and Process Data####
##1a. Read Spatial Data
scars <- readOGR("../../../GIS/Fires/Las Conchas/FireScarLocations/", layer="fs_in_las_conchas")
perim <- readOGR("../../../GIS/Fires/Las Conchas/Perimeter/", layer="Las_Conchas_Perim")
hs_lc <- #Las Conchas burn severity
  readOGR("../../../GIS/Fires/Las Conchas/Severity", layer="LasConchas_643_RdNBR")
  #Better for plotting holes in ggplot framework:
  #st_read("../../../GIS/Fires/Las Conchas/Severity", layer="LasConchas_643_RdNBR")

hs_total <- #Cumulative treeless area over multiple fires, from Jon Coop
  readOGR("../../../GIS/Fires/Las Conchas/Severity", layer="LasConchas_Treeless")
  #Better for plotting holes in ggplot framework:
  #st_read("../../../GIS/Fires/Las Conchas/Severity", layer="LasConchas_Treeless")

##1b. Fill holes
lc_fill <-#Fill holes of low-severity that are less than 9 pixels (8100m2, 0.81 ha) large
  fill_holes(hs_lc) 
treeless_fill <-#Fill holes of low-severity that are less than 9 pixels (8100m2, 0.81 ha) large
  fill_holes(hs_total) 


##1c. Visualize
#sf method:
scars_sf <- as(scars,"sf")
perim_sf <- as(perim,"sf")
treeless_sf <- as(treeless_fill,"sf")
p<-
  ggplot() +
  geom_sf(data = treeless_sf, fill="red", color=NA) +
  geom_sf(data = perim_sf, fill=NA, color = "black") +
  geom_sf(data = scars_sf, shape = 1) +
  theme_bw()

####2. Calculate data statistics####
overlap_lc <- #for each scar, does it overlap a high-severity patch? NA if no.
  over(scars,lc_fill)
length(na.exclude(overlap_lc))/length(overlap_lc) 
#28% of scars are in Las Conchas high-severity areas
overlap_treeless <- #for each scar, does it overlap a cumulative high-severity patch? NA if no.
  over(scars,treeless_fill)
length(na.exclude(overlap_treeless))/length(overlap_treeless) 
#61% of scars are in cumulative high-severity areas

lc_fill_disag <- #Make each distinct polygon its own feature ("multipart")
  disaggregate(lc_fill) 
lc_fill_disag <- #Remove hs patches < 0.81 ha ("small blowouts")
  lc_fill_disag[-which(area(lc_fill_disag)<8101)] 
treeless_fill_disag <- #Make each distinct polygon its own feature ("multipart")
  disaggregate(treeless_fill) 
treeless_fill_disag <- #Remove hs patches < 0.81 ha
  treeless_fill_disag[-which(area(treeless_fill_disag)<8101)] 

overlap_lc <- #for each scar, does it overlap a high-severity patch >0.81 ha? NA if no.
  over(scars,lc_fill_disag)
length(na.exclude(overlap_lc))/length(overlap_lc) 
#28% of scars are in Las Conchas high-severity areas > 0.81 ha (no change from above)
overlap_treeless <- #for each scar, does it overlap a cumulative high-severity patch? NA if no.
  over(scars,treeless_fill_disag)
length(na.exclude(overlap_treeless))/length(overlap_treeless) 
#59% of scars are in cumulative high-severity areas > 0.81 ha

#Visualize initial conditions
plot(perim)
plot(scars,add=T)
plot(perim)
plot(lc_fill_disag,col="coral3",add=T)
plot(perim)
plot(treeless_fill_disag,col="coral3",add=T)


####2b. SDC for both Las Conchas files####
hs_patches <- list(lc_fill,treeless_fill); fires.to.sample <- c("Las Conchas", "Cumulative Fires")

#for(f in c(1:length(fires.to.sample))){
#  hs_fire=hs_patches[[f]]
  #This step below is necessary to buffer from internal "holes" within the patch.
#  hs_fire2 <- createSPComment(hs_fire)
  #plot(hs_fire,col="darkred",border="transparent")
#  Sys.time()
#  decay.table=decay(hs_fire=hs_fire2,buf_max=1000,buf_inc=100,name=fires.to.sample[f])
#  Sys.time()
#  ifelse(f==1, fires_long <- decay.table, fires_long <- rbind(fires_long,decay.table) )
#  print(paste(fires.to.sample[f],Sys.time()))
#  gc()
#}

#fires_long <- calculate.sdc(fires_long)
#summary_fires <- 
#  fires_long %>%
#  group_by(name) %>%
#  summarise(sdc = mean(sdc),
#            log_sdc = log(mean(sdc)))

####3a. Play with Swiss Cheese simulation model growing patches####
#3.1: 
centroids <- SpatialPoints(getSpPPolygonsLabptSlots(lc_fill_disag), 
                           proj4string = CRS(proj4string(lc_fill_disag)))
proj4string(centroids) <- proj4string(lc_fill_disag)

#Convert perimeter to points:
lc_line  <- as(perim,"SpatialLines")
lc_pts <- as(lc_line, "SpatialPointsDataFrame")
lc_pts <- lc_pts[,-c(2,3)]; 
names(lc_pts) <- #Add marker to the perimeter points to identify the boundary type.
  "perim" 
scars$scar <-1 #Add marker to the scar points to identify the boundary type.
scars <- #Delete extra variables from scars data frame
  scars[,-c(1:5)] 
boundaries <- #Single layer combining plot boundaries and fire scars
  bind(lc_pts,scars, makeUniqueIDs = TRUE) 
proj4string(boundaries) <- CRS(proj4string(lc_pts))
centroids_random <- #Randomly position the initial centroid for 698 hs patches
  spsample(perim,n=698, type = "random")

#Visualize:
plot(perim)
plot(lc_fill_disag,col="coral3",add=T)
points(scars, col = "cyan")
plot(perim); points(centroids, col = "coral3")
plot(perim)
points(centroids_random, col = "cyan4")

#Simulation
keep_going <- #to start, we want to buffer every point
  rep(TRUE,length(centroids_random)) 
Sys.time()
sim_patches <- #Buffer every patch by 1m to start
  gBuffer(centroids_random,width = 1, byid = TRUE) 
for(w in seq(10,1000,by=10)){ #~10 seconds
  intersections <- #has a patch hit a scar or the perimeter?
    over(sim_patches,boundaries) 
  done <- #TRUE if the patch has hit a scar or the perimeter
    apply(intersections,2,is.na) 
  done <- #TRUE if the patch has hit a scar or the perimeter
    apply(!done,1,any) 
  keep_going <- #TRUE if the patch has not hit a scar or the perimeter
    !done 
  done <- as.vector(which(done)) #formatting
  keep_going <- as.vector(which(keep_going)) #formatting
  hold_patches <- #save patches that are done and expand the rest
    sim_patches[done]
  new_patches <- #Expand patches that haven't hit anything yet
    gBuffer(centroids_random[keep_going],width = w, byid = TRUE)
  sim_patches <- 
    bind(hold_patches,new_patches,keepnames=FALSE)
  centroids_random <- #reorder the centroids for consistency
    bind(centroids_random[done],centroids_random[keep_going],keepnames=FALSE)
  #print(w); #print(keep_going)
  #print(sim_patches[1])
  #plot(sim_patches[1:10])
}
Sys.time()

#Visualize simulation results
hs_merge<- st_as_sf(sim_patches)
hs_merge<- st_union(hs_merge)
plot(perim)
plot(hs_merge, col = "coral3", add = T)
plot(scars,add=T)

####3b. Simulation of largest patches####
lc_largest <- #20 largest patches in Las Conchas, in ha
  sort(area(lc_fill_disag),decreasing = TRUE)[1:20]
#20336ha/26429ha total = 77%
lc_largest_ha <- #20 largest patches in Las Conchas, in ha
  sort(area(lc_fill_disag),decreasing = TRUE)[1:20]*0.0001
treeless_largest_ha <- #20 largest patches in Las Conchas, in ha
  sort(area(treeless_fill_disag),decreasing = TRUE)[1:20]*0.0001
lc_radii <- sqrt(lc_largest/pi)  #convert from ha to m2
ggplot()+
  geom_histogram(aes(lc_largest_ha))+
  theme_bw() + 
  labs(x="patch size (ha)")
ggplot()+
  geom_histogram(aes(lc_largest_ha[3:20]))+
  theme_bw() + 
  lims(x = c(0,800)) +
  labs(x="patch size (ha)")

c.df <- data.frame(c.rad = lc_radii, prop_w_scars=NA, x=NA, y=NA)
c.rad <- c.polys <- list()
for(r in 1:nrow(c.df)){
  c.rad[[r]] <- c.df[r,"c.rad"]
  c.df[r,"area_ha"] <- pi*(c.df[r,"c.rad"]^2) * 0.0001
  coords_sp <- #1000 random centroid points within the perimeter
    spsample(perim,n=100, type = "random")
  c.poly <- gBuffer(coords_sp,width=c.rad[[r]],quadsegs = 40,byid=TRUE)
  #plot(c.poly)
  c.df[r,"prop_w_scars"] <- 
    nrow(na.exclude(over(c.poly,scars)))/nrow(over(c.poly,scars))
  #c.poly.decay <- decay(c.poly,buf_inc=10,buf_max=1000,name=as.character(c.rad[[r]]))
  #c.df[r,"sdc"] <- calculate.sdc(c.poly.decay)[,"sdc.name"] %>% unique()
  #c.polys[[r]] <- c.poly
  
}
library(splines)
ggplot(c.df,aes(x = area_ha, y = 1-prop_w_scars)) +
  stat_smooth(method="glm", method.args = list(family="quasipoisson"), 
              formula = y ~ ns(x, 3)) +
  geom_point() +
  theme_bw() + 
  lims(y = c(0,1)) +
  labs(x = "patch size (ha)", y = "likelihood of not \nencountering a scar")


#Try to place all circles simultaneously
sims <- c()
for(s in 1:400){
c.rad <- c.polys <- list();
c.df <- data.frame(c.rad = lc_radii, prop_w_scars=NA, x=NA, y=NA)
for(r in 1:nrow(c.df)) {
  c.rad[[r]] <- c.df[r,"c.rad"]
  c.df[r,"area_ha"] <- pi*(c.df[r,"c.rad"]^2) * 0.0001
  c.df[r,c("x","y")] <- #1 random centroid point within the perimeter
    coordinates(spsample(perim,n=1, type = "random", iter = 30))
  c.polys[[r]] <- gBuffer(SpatialPoints(c.df[r,c("x","y")]),
                          width=c.rad[[r]],quadsegs = 40,byid=TRUE)
  c.polys.comb <- bind(c.polys)
  if(r>1){
    
    while(any(gOverlaps(unlist(c.polys.comb),byid = TRUE)) | 
          any(gContainsProperly(unlist(c.polys.comb),byid = TRUE))){
      c.df[r,c("x","y")] <- #1 random centroid point within the perimeter
        coordinates(spsample(perim,n=1, type = "random"))
      c.polys[[r]] <- gBuffer(SpatialPoints(c.df[r,c("x","y")]),
                              width=c.rad[[r]],quadsegs = 40,byid=TRUE)
      c.polys.comb <- bind(c.polys)
      print(paste("overlap ", r))
    }
  }
  print(r)
}
#plot(c.polys.comb)
proj4string(c.polys.comb) <- proj4string(scars)
sims[s] <- sum(na.exclude(over(c.polys.comb,scars)))/20
}
tmp<-c(tmp,sims)
ggplot()+
  geom_histogram(aes(tmp),binwidth = 0.05)+
  theme_bw() + 
  lims(x = c(0,1)) +
  labs(x="fraction of 20 largest patches encountering a scar")

plot(perim)
plot(c.polys.comb,add=T,col="red")
points(scars, pch = "+")


####4. Scar distribution within treeless areas####
sim_index_lc<-c()
Sys.time()
for(s in c(1:1000)){ #3 minutes
  scars_random <- spsample(perim,n=length(scars), type = "random")
  overlap <- #for each scar, does it overlap a high-severity patch? NA if no.
    over(scars_random,lc_fill_disag)
  sim_index_lc[s] <- length(na.exclude(overlap))/length(overlap) 
  #28% of scars are in high-severity areas <0.81 ha
}
Sys.time()

#hist(sim_index_lc)
ggplot()+
  geom_histogram(aes(x = sim_index_lc), binwidth = 0.01, fill = "white", col = "black") +
  geom_vline(aes(xintercept = length(na.exclude(overlap_lc))/length(overlap_lc)),
             col = "red") +
  labs(title = "proportion random fire scars falling \nwithin Las Conchas high-severity") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

sim_index_treeless<-c()
Sys.time()
for(s in c(1:1000)){ #4 minutes
  scars_random <- spsample(perim,n=length(scars), type = "random")
  overlap <- #for each scar, does it overlap a high-severity patch? NA if no.
    over(scars_random,treeless_fill_disag)
  sim_index_treeless[s] <- length(na.exclude(overlap))/length(overlap) 
  #28% of scars are in high-severity areas <0.81 ha
}
Sys.time()

ggplot()+
  geom_histogram(aes(x = sim_index_treeless), binwidth = 0.01, fill = "white", col = "black") +
  geom_vline(aes(xintercept = length(na.exclude(overlap_treeless))/length(overlap_treeless)),
             col = "red") +
  labs(title = "proportion random fire scars falling \nwithin cumulative high-severity") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

