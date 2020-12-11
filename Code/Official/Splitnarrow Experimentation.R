#Splitnarrow experimentatino

splitnarrow <- function(pol){
  ###https://gis.stackexchange.com/questions/333817/splitting-polygons-at-narrowest-part-using-r?rq=1
  ### split a polygon at its narrowest point.
  ### modified by Jens Stevens
  #pol = hs2$geometry[51] #pol is the test polygon, coming out of SpatialFunctions. 51, 117
  
  sdist = -15
  eps = 1
  possible_connections = list()
  bparts = st_buffer(pol, sdist)
  features = st_cast(st_sfc(bparts), "POLYGON")
  
  ##Find which separate polygons are still "large" (> 100 ha)
  sub_areas = features[which(as.vector(st_area(features))*0.0001 > 100)]
  #plot(sub_areas,col=rainbow(8)) #option to plot
  
  while(sdist >= -180){
    if(length(sub_areas)>1){
      #Do something
    }
    
    #Beta testing polygon 51.
    #Split into three units and correctly identify where to insert the cut lines
    sdist = -50  
    bparts = st_buffer(pol, sdist)
    features = st_cast(st_sfc(bparts), "POLYGON")
    sub_areas = features[which(as.vector(st_area(features))*0.0001 > 100)]
    lost_areas = st_cast(st_difference(pol,bparts),"POLYGON")
    #plot(lost_areas, col = "darksalmon")
    #plot(sub_areas,col=rainbow(8), add = TRUE) #option to plot
    possible_connections[[1]] = st_buffer(sub_areas[[1]],50)
    
    ## find where the two separate polygons are closest, this is where
    ## the internal buffering pinched off into two polygons.
    pinch = st_nearest_points(sub_areas[1],sub_areas[2]) #Can take ~10 or more seconds
    plot(pinch,add=TRUE,col="pink",lwd = 4)
    
    ## buffering the pinch point by a slightly larger buffer length should intersect with
    ## the polygon at the narrow point. 
    cast = st_cast(pol,"MULTILINESTRING")
    bparts2 = st_buffer(pinch,-(sdist-(eps)))
    inter = st_intersection(cast, bparts2) #single feature
    join = st_cast(st_sfc(inter), "LINESTRING") #same as inter but multiple features
  }
  
}

#Sequential buffering

tmp = st_buffer(pol,dist = -30, joinStyle = 'MITRE', mitreLimit = 2)
plot(pol,col="darksalmon")
plot(tmp,col="darkgreen",add=T)
tmp2 = st_voronoi(pol)
plot(tmp2,col="darkgreen",add=T)

buffer_seq <- function(pol){
  #START HERE: sequentially buffer inward at 15 m increments and sequentially pull out resulting polygons > 100 ha for export
  sdist = -15
  eps = 1
  bparts = st_buffer(pol, sdist)
  features = st_cast(st_sfc(bparts), "POLYGON")
  
  ##Find which separate polygons are still "large" (> 100 ha)
  sub_areas = features[which(as.vector(st_area(features))*0.0001 > 100)]
  #plot(sub_areas,col=rainbow(8),main = sdist) #option to plot
  previous_iteration = sub_areas #initial previous iteration is the 15 m buffer
  index = 1
  
  while(sdist >= -180){
        sdist = sdist - 15
        #bparts = st_buffer(pol, sdist) #original, if operating on original polygon
        bparts = st_buffer(previous_iteration,-15) #if buffereng previous version by 15 more. Might not matter.
        features = st_cast(st_sfc(bparts), "POLYGON")
        
      
        ##Find which separate polygons are still "large" (> 100 ha)
        sub_areas = features[which(as.vector(st_area(features))*0.0001 > 100)]
        print(paste0("sdist = ", sdist))
        print(paste0("number of sub areas = ",length(sub_areas)))
        #plot(sub_areas,col = rainbow(10), main = sdist)
        
        if(any(st_overlaps(previous_iteration,sub_areas,sparse = FALSE))){
          #If an individual polygon in the previous iteration no longer overlaps the current iteration
          #That is a storage polygon to be "cut". It goes into out_storage
          out_storage[[index]] = which(st_overlaps(previous_iteration,sub_areas,sparse = FALSE))
          index = index + 1
          
        }
        previous_iteration = sub_areas
        #plot(previous_iteration, col = rainbow(10))
        #START HERE: Can we do a "difference" of the buffered and unbuffered versions and then make that the "slicer" which continously expands, splits the original polygon, and every time a new polygon is split off, it is stored elsewhere?
    }
  
}

buffer_split <- function(pol){
  #Generate a new (buffered) and old (unbuffered) singlepart polygon
  large_buf = st_buffer(pol, 20)
  small_buf = st_buffer(pol, 10)
  large_buf = st_make_valid(st_cast(st_sfc(st_buffer(pol, 20)), "POLYGON"))
  small_buf = st_make_valid(st_cast(st_sfc(st_buffer(pol, 10)), "POLYGON"))
  #plot(large_buf,col="darksalmon")
  #plot(small_buf,col="blue",add=TRUE)
  slicer = st_cast(st_difference(large_buf,small_buf),"POLYGON")
  #plot(slicer,col="darksalmon")
  #plot(pol,add=TRUE,col="darkred")

  slicer_buf = st_make_valid(st_cast(st_sfc(st_buffer(slicer, 30)), "POLYGON"))
  plot(slicer_buf,col="purple")
  plot(slicer,col="darksalmon",add=TRUE)
  
}

nearest_point <- function(pol){
  st_nearest_points()
}

lungs <- function(pol,sdist){
  ###https://gis.stackexchange.com/questions/333817/splitting-polygons-at-narrowest-part-using-r?rq=1
  ### split a polygon at its narrowest point.
  ### modified by Jens Stevens
  #pol = hs2$geometry[25] #pol is the test polygon, coming out of all polys > 1 ha. 51, 117. 25 to test.
  #pol = hs3$geometry[5] #pol is the test polygon, coming out of all polys >10 ha. 5 to test.
  #pol is a sfc_POLYGON - single 2D geometry with holes. First list element is the outer polygon wall.
  #an sfc_POLYGON is a list of length n; if n>1 (as with hs3$geometry), each list element is a different outer wall
  #so pol is an sfc_POLYGON; pol[[1]] is essentially the same as pol but it is a generic list
  #pol[[1]][[1]] is a simple matrix of coordinates (x = first column, y = second) of outer wall
  #pol[[1]][[2:n]] are simple matrices of coordinates for each of n holes.
  
  output = list() #object to store output
  
  #Step 1: inward buffer on pol by sdist units (generall meters).
    #This may split pol into multiple polygons, creating a MULTIPOLYGON, if pinch points exist at 2 * sdist
  #sdist = -15 #user-defined; sdist should be half the pixel size (e.g. LANDSAT 30m, sdist = -15)
  pol = st_make_valid(pol) #Corrects rings that self-overlap. Maybe not needed; issue for l = 10
  bparts = #if the value of sdist severs pinch points, bparts is a sfc_MULTIPOLYGON, 
    #bparts[[1]][[1]] is a list of coords for the first pinched polygon, which may or may not have holes
    #but bparts[[1]][[1]] is NOT an sfc_POLYGON, it is *just* a list.
    st_buffer(pol, sdist - 0.1, joinStyle="MITRE", mitreLimit=2) #had cut by sdist - 0.1 but deprecated?
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
  #plot(inhale1, col = rainbow(8), add = T) #The larger pinched polygons only (>10 ha)
  
  #Identify split points
  n_chunks = length(inhale1)
  if(n_chunks >= 2){#split off the smallest qualifying polygon
    diff1 = #diff1 is the difference between the original polygon and the inhale
      st_cast(st_difference(pol,bparts),"POLYGON")
    diff1 = #Need this for hs3[9], maybe more, to get rid of the unconnected "holes"
      #CHECKME not sure if I'm keeping this, might mess up polygons != 9 but doesn't look like it.
      #an alternative might be to st_cast diff1 as a single sfc_POLYGON of length 1.
      diff1[which.max(st_area(diff1))]
    #plot(pol) #original
    #plot(inhale1, col = rainbow(8), add = T) #The larger pinched polygons only
    #plot(diff1,border="darksalmon", col = "darksalmon", add = T)
    small_chunk = inhale1[which.min(st_area(inhale1))]
    exhale1 = #This approximates the original small chunk
      st_make_valid(st_cast(st_sfc(st_buffer(small_chunk,
                                             -sdist + 0.1, 
                                             #had exhaled by twice as much (*-1.95) but don't need to? TBD
                                             joinStyle="MITRE",
                                             mitreLimit=2)
                                   ), "POLYGON")
                    )
    #plot(small_chunk, border = "blue"); plot(pol,add=T)
    #plot(exhale1,border="goldenrod",add=T)
    #plot(diff1,col="darksalmon",border=NA,add=T)
    buffer1 = #buffer1 is an sfc_POLYGON of the "ring" buffer around the small chunk from the inhale pinch
      st_cast(st_intersection(diff1,exhale1),"POLYGON")
    buffer1 = #here, any "holes" in buffer2 are removed, 
    #since it's only the exterior polygon which might touch diff1
      buffer1[which.max(st_area(buffer1))] #testing; taking the holes out from diff might make this redundant.
    #update: Yes, fixing "diff1" by removing unneccessary holes makes this step redundant.
    #update2: Actually need this step, l = 10 is an example of why (kitty corner chunks).
    
    buffer2 = #buffer2 is an sfc_POLYGON of the "ring" around the larger chunk 
      #(remove buffer1 from the full diff1)
      #st_difference(diff1, buffer1) #temp, deprecate
      #warning here if diff1 or buffer1 have multiple parts; should be removed by which.max() operations above
      st_cast(st_difference(diff1,buffer1),"POLYGON") #warning ok (warns here for l = 8 but not 9)
    buffer2 = #here, any "small" parts of buffer2 are likely small protrusions off the small chunk captured in diff1
      #these can be removed so the pinch point remains a single contact point to the "main buffer 2"
      buffer2[which(as.vector(st_area(buffer2))*0.0001 > 1)] #retaining only the large "main buffer2"
    #problem here because for larger complex polygons there are many that meet this criteria...
    #plot(pol, main = "buffer2 large only")
    #plot(buffer1,col="blue",border="blue",add = T) 
    #plot(buffer2,col="darkred",border="darkred",add = T)
    
    #Pinch it off
    chunks_line_intersect = st_intersection(buffer1,buffer2)
    #Sonofabitch I need to merge the multiple line segments into a single line I think.
    #Example l = 13
    #CHECKME Start Here, it's close but st_
    #May need "if" statement here to caveat this if the original CLI string is only one line...
    mls_pts = st_cast(chunks_line_intersect,"POINT")
    chunks_line_intersect = #generate MULTILINESTRING from segments, first breaking to points above
      st_cast(st_linestring(st_coordinates(mls_pts)),"MULTILINESTRING")
    #Ah it's close but not quite there. A problem with the format of chunks_line_intersect now.
    #Above is based on https://github.com/r-spatial/sf/issues/321
    
    
    
    if(any(lapply(chunks_line_intersect[[1]],length)==2)){
      #This tests for "single point" intersections at kitty corner, and removes them
      #This is because chunks_line_intersect must be a MULTILINESTRING for st_split, below, to work
      remove_points = 
        chunks_line_intersect[[1]][-which(lapply(chunks_line_intersect[[1]],length)==2)]
      chunks_line_intersect = st_cast(st_sfc(remove_points),"MULTILINESTRING")
      #Not the most elegant but it works (e.g. polygon 10 as an example).
      #chunks_line_intersect now lacks a CRS
    }

    #plot(pol)
    #plot(buffer1,col="blue",border="blue",add = T) 
    #plot(buffer2,col="darkred",border="darkred",add = T)
    #plot(chunks_line_intersect, col = "pink", lwd = 3, add = T)
    splitpol <- st_split(pol,chunks_line_intersect) 
    splitpol2 <- st_cast(splitpol[1]) #Cast as sfc_POLYGON to finish the split; should be length 2...
    #length(splitpol2)
    #plot(splitpol2, col = rainbow(8), main = l,add=T)
    output[[1]] = "Split polygon" 
    output[[2]] = splitpol2
    
  } else if (n_chunks == 1) { 
    output[[1]] = "Not split" 
    } else { #START HERE
      #Need to modularize chunks = 2 in to a single function and call it multiple times here
      #ultimately storing the output in a 3+ feature polygon in output[[2]]
      output[[1]] = "Check me"
    }
  output
}

for (l in 1:13 ){ #eventually 1:length(hs3$geometry)
  pol = hs3$geometry[l] #breaking on l = 9
  lungs_output = lungs(pol, sdist = -15)
  if(lungs_output[[1]] == "Split polygon"){
    splitpol2 = lungs_output[[2]]
    plot(pol, main = l)
    plot(splitpol2, col = rainbow(8), add = T)
    #START HERE store splitpol2 and remove pol[l]
    #maybe need to find closest feature here:
    #https://stackoverflow.com/questions/49200458/find-nearest-features-using-sf-in-r
  } else if(lungs_output[[1]] == "Check me") {
    plot(pol, main = l)
    mtext("Check me", side = 3)
  } else {print(c(l, " not split"))}
}
