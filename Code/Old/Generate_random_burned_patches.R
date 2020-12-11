########################
##### Load packages ####
########################

library(NLMR)
library(raster)
library(igraph)


###########################
#### Set the grid size ####
###########################

# xmax and ymax set the number of cells in the x and y direction
# increasing the grid size slows down the 'nlm_gaussianfied' function
xmax <- 200
ymax <- 200


##################################
#### set variogram parameters ####
##################################

nugget <- 0
magvar <- 100

# set the variogram range as the radius of a circle with an area equal to 'pct' percent of the grid area
# increasing the variogram range will make the 'nlm_gaussianfield' function run slower
pct <- 5
rng <- round(sqrt(((pct/100) * xmax * ymax)/pi), 0)

# Alternatively, set the variogram range in number of cells
# rng <- 20

# Generate a continuous Gaussian random field based on the variogram input parameters
gf <- nlm_gaussianfield(ncol = xmax, nrow = ymax, autocorr_range = rng, mag_var = magvar, nug = nugget)

par(bty="n") 
plot(gf, axes = FALSE) 


##################################################################
#### Convert the continuous Gaussian field into a binary grid ####
##################################################################

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
gfv <- sort(as.vector(gf), decreasing = TRUE)
cutoff <- gfv[length(gfv) * prop.burned]
gf.binary <- calc(gf, fun = convert.to.binary)

plot(gf.binary, axes = FALSE) 


###########################################
#### Calculate patch size distribution ####
###########################################

# These steps will create a data frame called 'patch.size' with columns "value" and "count"
# "value" is the patch number (arbitrary values)
# "count" is the number of cells in the patch, which will be sorted from largest to smallest
# unburned patches (cell values of 0 in the binary grid) are excluded from the patch.size data frame

# copy the binary grid and convert the zeros to NAs  
gf.binary1 <- gf.binary
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
sum(patch.size$prop)
prop.burned

head(patch.size)





