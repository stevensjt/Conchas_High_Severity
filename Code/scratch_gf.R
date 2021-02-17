#scratch
pct <- 5 #must be such that rng >= 1
nugget <- 0
magvar <- 10
rng <- #set the variogram range as the radius of a circle with an area equal to 'pct' percent of the grid area
  round(sqrt(((pct/100) * 100 * 100)/pi), 0) #define rng based on pct
parms_text <- paste0("pct_", pct,"rng_",rng,"_nug_",nugget,"_mv_",magvar)

####Hypothetical landscape####
#Lessons learned
#1. If nugget = 0, magvar doesn't matter
#2. Nugget of even 1 starts to create scatter and it creates more crumbs than holes.
#3. When nugget is high, reducing magvar creates even more scatter, to the point of basically removing AC
#4. Nugget = 1 and magvar = 10 is the EXACT SAME as Nugget = 10 and magvar = 100.
#4b Likewise Nugget = 10 and magvar = 10 is the SAME as nugget = 100 and magvar = 100
#4c. When nugget is greater than magvar, you get almost perfect scatter; when smaller, clusters are tighter.

gaussian_field <- nlm_gaussianfield(ncol = 100, nrow = 100,
                                    user_seed = 2,
                                    resolution = 1, #Can set to 30 but maybe doesn't matter?
                                    autocorr_range = rng,
                                    mag_var = magvar,
                                    nug = nugget)


convert.to.binary <- function(x) { ifelse(x <  cutoff, NA, 1) }
convert.range.to.binary  <- function(x) { 
  v = cutoffs
  l = length(v) # v = vector of values to test between-ness 
  iters = l/2 + 0.5
  if(iters == 2){
    ifelse(between(x,v[1],v[2]) | x > v[3], 1, NA)
  }
}
gfv <- sort(as.vector(gaussian_field), decreasing = TRUE) #gfv = gaussian field vector
#cutoff <- gfv[length(gfv) * (0.35)] 
cutoffs <- c(gfv[length(gfv) * 0.99],
             gfv[length(gfv) * 0.90],
             gfv[length(gfv) * 0.25] 
)
#gaussian_field <- calc(gaussian_field, fun = convert.to.binary)
gaussian_field <- calc(gaussian_field, fun = convert.range.to.binary)

landscapetools::show_landscape(gaussian_field) + labs(title = parms_text)

####Apply to actual landscape
pct <- 0.25 #must be such that rng >= 1
nugget <- 5
magvar <- 50
rng <- #set the variogram range as the radius of a circle with an area equal to 'pct' percent of the grid area
  round(sqrt(((pct/100) * 100 * 100)/pi), 0) #define rng based on pct
parms_text <- paste0("pct_", pct,"rng_",rng,"_nug_",nugget,"_mv_",magvar)

gf <- #Slow, ~20 seconds w 30 m resolution (less now)?
  nlm_gaussianfield(ncol = xmax, nrow = ymax, 
                    user_seed = 1,
                    resolution = 1,
                    autocorr_range = rng, 
                    mag_var = magvar, 
                    nug = nugget)

extent(gf) <- round(extent(perim),0) #This plus next step put it on a 30 m resolution
gf@crs <- crs(perim)
gf.mask <- raster::mask(gf, perim) #mask by Las Conchas
#convert.to.binary <- function(x) { ifelse(x <  cutoff, NA, 1) }
gfv <- sort(as.vector(gf.mask), decreasing = TRUE) #gfv = gaussian field vector
cutoff <- gfv[length(gfv) * (prop_hs+0.03)] 
gf.binary <- calc(gf.mask, fun = convert.to.binary)
gf.simple <- raster_holes_and_crumbs(gf.binary) #Now takes < 10 seconds
area_sim <- # number of cells = 1, times area of a cell 900.5704 m2 
  freq(gf.simple, value = 1) * 900.5704 #in m2
prop_hs_sim <- as.vector(round(area_sim/st_area(perim),3)[[1]])

png(paste0("./Output/Tuning_Sims/old_bad/scratch/",parms_text,".png"), height = 12, width = 8, units = "in", res = 300)
landscapetools::show_landscape(gf.simple) + labs(title = parms_text)
dev.off()

#Start here compare this to the actual landscape