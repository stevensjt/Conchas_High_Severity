library(sf)

pol <- read_sf("/Users/jtstevens/Documents/New Mexico/USGS/GIS/Fires/Las Conchas/Severity/Treeless_mods_postfire/LasConchas_Treeless_1ha.shp")
pol <- pol[,"Area_ha"]
plot(pol)
sdist = -29
ppol = splitnarrow(pol, sdist, 1e-3)
plot(ppol, col=1:2)






splitnarrow <- function(pol, sdist, eps){
  ###
  ### split a polygon at its narrowest point.
  ###
  
  ### sdist is the smallest value for internal buffering that splits the
  ### polygon into a MULTIPOLYGON and needs computing before running this.
  
  ### eps is another tolerance that is needed to get the points at which the
  ### narrowest point is to be cut.
  
  ## split the polygon into two separate polygons
  bparts = st_buffer(pol, sdist)
  features = st_cast(st_sfc(bparts), "POLYGON")
  
  ## find where the two separate polygons are closest, this is where
  ## the internal buffering pinched off into two polygons.
  
  pinch = st_nearest_points(features[1],features[2])
  
  ## buffering the pinch point by a slightly larger buffer length should intersect with
  ## the polygon at the narrow point. 
  inter = st_intersection(
    st_cast(pol,"MULTILINESTRING"),
    st_buffer(pinch,-(sdist-(eps))
    )
  )
  join = st_cast(st_sfc(inter), "LINESTRING")
  
  ## join is now two small line segments of the polygon across the "waist".
  ## find the line of closest approach of them:
  splitline = st_nearest_points(join[1], join[2])
  
  ## that's our cut line. Now put that with the polygon and make new polygons:
  mm = st_union(splitline, st_cast(pol, "LINESTRING"))
  parts = st_collection_extract(st_polygonize(mm))
  parts
}
