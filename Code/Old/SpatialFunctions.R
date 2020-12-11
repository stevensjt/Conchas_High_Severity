##Fill holes that are less than or equal to 11 pixels large (11*900m2=9900m2, or ~1 ha)
fill_holes <- 
  function(hs_fire){
    #hs_fire is a multipart SpatialPolygonsDataFrame; 
    #a single polygon (e.g. #3) is pol= hs_fire@polygons[[1]]@Polygons[[3]]
    #to view: plot(SpatialPolygons(list(Polygons(list(pol),ID="3"))))
    hs_fire_p = #store each polygon as a list element
      slot(hs_fire, "polygons") 
    holes = #list with one element, a vector TRUE if the polygon is a hole
      lapply(hs_fire_p, function(x) sapply(slot(x, "Polygons"), slot, "hole"))
    areas = 
      lapply(hs_fire_p, function(x) sapply(slot(x, "Polygons"), slot, "area"))
  res <- lapply(1:length(hs_fire_p), function(i) slot(hs_fire_p[[i]], "Polygons")[!(holes[[1]]&areas[[1]]<8100)])#Select the polygons that are not holes. The "i" here is an artifact of the example code; the fires here only have one polygon ID so i=1 (it's a multipart polygon).
  IDs <- row.names(hs_fire)
  hs_fire_fill <- SpatialPolygons(lapply(1:length(res), function(i)
    Polygons(res[[i]], ID=IDs[i])), proj4string=CRS(proj4string(hs_fire)))
  return(hs_fire_fill)
  ##One consequence of this is that it's no longer a sp data frame, and there are some warnings. 
}

##
##First converting to a singlepart SpatialPolygonsDataFrame using 'disaggregate'
##Then Fill holes that are less than or equal to 11 pixels large (11*900m2=9900m2, or ~1 ha)
##Then remove crumbs that are less than or equal to 11 pixels large (11*900m2=9900m2, or ~1 ha)
##Both using the awesome smoothr https://cran.r-project.org/web/packages/smoothr/vignettes/smoothr.html
##Then split narrow points of polygons using splitnarrow from
##https://gis.stackexchange.com/questions/333817/splitting-polygons-at-narrowest-part-using-r?rq=1
clean_up_1 <- 
  function(hs_fire){
    #if using "sp" approach:
    #hs_fire is a multipart SpatialPolygonsDataFrame; 
    #a single polygon (e.g. #3) is pol= hs_fire@polygons[[1]]@Polygons[[3]]
    #to view: plot(perim); plot(SpatialPolygons(list(pol)),add=T) 
    hs1 = #hs1 is first temporary object holding high severity polygons
      sp::disaggregate(hs_fire) #Convert to singlepart; takes ~20-60 seconds depending on file size
    #hs_fire_fill is now a singlepart SpatialPolygonsDataFrame
    #a single polygon (e.g. #1) is pol= hs_fire_fill@polygons[[1]], including holes
    #holes = #list with one element, a vector TRUE if the polygon is a hole
    #  lapply(hs_fire_fill[[1]], function(x) sapply(slot(x, "Polygons"), slot, "hole"))
    
    #convert to "sf" object to use special functions from smoothr:
    hs1 = st_as_sf(hs1) 
    #plot(hs1$geometry, col = "darkred") #Plot whole fire
    #plot(hs1$geometry[51], col="darksalmon") #Polygon 51 is the second-largest
    hs1$geometry = #fill holes <= 10,000 m2 = 1ha, minimum holesize retained is 10,001 m2
      smoothr::fill_holes(hs1$geometry, threshold = 10000) #Takes ~25 seconds
    
    hs2 = #drop crumbs <= 10,000 m2 = 1ha, minimum crumb size retained is 10,001 m2
      #hs2 is second temporary object holding high severity polygons
      #need to create new object here because "drop_crumbs" changes the object length
      drop_crumbs(hs1,threshold = 10000) #Takes ~40 seconds
    #Number of discrete patches drops from 2948 to 607 for Las Conchas
    #plot(hs2$geometry, col = "darkred") #Plot whole fire
    #rm(hs1) #Clean up working environment
    
    hs2$Area_m2 = st_area(hs2)
    hs2$Area_ha = as.vector(hs2$Area_m2)*0.0001
    #st_write(hs2, "./SpatialOutput/LasConchas_643_RdNBR_CleanUp1.shp",delete_layer=T)
    
    hs3 = #drop crumbs <= 1,000,000 m2 = 100ha, minimum crumb size retained is 10,001 m2
      #here we are looking at which large (>100 ha) polygons have pinch points that we can chop up
      drop_crumbs(hs2,threshold = 1000000) #Takes ~15 seconds
    
    
  }

splitnarrow <- function(pol, sdist, eps){
  ###https://gis.stackexchange.com/questions/333817/splitting-polygons-at-narrowest-part-using-r?rq=1
  ### split a polygon at its narrowest point.
  ### modified by Jens Stevens
  
  sdist = -15
  eps = 1
  ### sdist is the smallest value for internal buffering that splits the
  ### polygon into a MULTIPOLYGON and needs computing before running this.
  
  ### eps is another tolerance that is needed to get the points at which the
  ### narrowest point is to be cut.
  
  ## split the polygon into multiple separate polygons based on sdist
  bparts = st_buffer(pol, sdist)
  features = st_cast(st_sfc(bparts), "POLYGON")
  
  ##Find which separate polygons are still "large" (> 100 ha)
  sub_areas = features[which(as.vector(st_area(features))*0.0001 > 100)]
  #plot(sub_areas,col=rainbow(8)) #option to plot
  
  
  ## find where the two separate polygons are closest, this is where
  ## the internal buffering pinched off into two polygons.
  
  pinch = st_nearest_points(sub_areas[1],sub_areas[2]) #Can take ~10 or more seconds
  
  ## buffering the pinch point by a slightly larger buffer length should intersect with
  ## the polygon at the narrow point. 
  cast = st_cast(pol,"MULTILINESTRING")
  bparts2 = st_buffer(pinch,-(sdist-(eps)))
  inter = st_intersection(cast, bparts2) #single feature
  join = st_cast(st_sfc(inter), "LINESTRING") #same as inter but multiple features
  
  #CHECKME: Will eventually want to modify this to cut polygons multiple times. And calibrate the buffer size to the polygon size somehow (bigger buffer for bigger polygon but maybe need to simplify first, though the large polygon filter above might solve it.)
  ## join is now two small line segments of the polygon across the "waist".
  ## find the line of closest approach of them:
  splitline = st_nearest_points(join[1], join[2])
  #plot(splitline, add=T, col = "darkred")
  
  ## that's our cut line. Now put that with the polygon and make new polygons:
  splitpol = st_split(pol,splitline)
  splitpol2 = st_cast(splitpol) #Cast as sfc_POLYGON
  #plot(splitpol2, col = rainbow(8))
  #Old code not used:
  #mm = st_union(splitline, st_cast(pol, "LINESTRING"))
  #parts = st_collection_extract(st_polygonize(mm))
  #parts
  
  splitpol2
  
}

