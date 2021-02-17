##Other_nlmr##
####1. Percolate randome raster with base tools####
r = raster(nrows = c(0:xmax), ncols = ymax, vals = 1)
plot(r)
ymax

####2. random_cluster####
#SLOW at scale
#minimal control?
Sys.time() #16s for 100x100
random_cluster <- nlm_randomcluster(ncol = 100, nrow = 100,
                                    p = 0.4,
                                    ai = c(0.65, 0.35))
Sys.time()
landscapetools::show_landscape(random_cluster)

####3. random neighborhood model (neutral landscape) ####
#might be something here
neigh_raster <- nlm_neigh(ncol = 100, nrow = 100, p_neigh = 0.9, p_empty = 0.01,
                          categories = 2, neighbourhood = 4,
                          proportions = c(0.65, 0.35))
## Not run:
# visualize the NLM
landscapetools::show_landscape(neigh_raster)
## End(Not run)
#p_neigh = 0.9 and p_empty = 0.01 is about as clumpy as I can get it, but clumps still not random enough.



####5. Another idea: just play with the cutoff in the gaussian? e.g. select a high cutoff with a highly correlated model and make that the largest patch? And/or try to match the distribution for the cutoffs >= 100 ha, do the clean-cutting, and then set another cutoff range for all patches <100 ha > 1 ha ? Like a mid-range value?

####6. Another idea: Explore landscape metrics to figure out patch area?