#Author: Jens Stevens; stevensjt@gmail.com
#This code tweaks the parameters of the gaussian field simulations, trying to iteratively match to the values of the observed "null model" layer.

####0. Libraries####
source("./Code/Official/0_Functions.R") #read in customized functions necessary for processing.
source("./Code/Official/0_RasterFunctions.R") #read in customized raster functions necessary for processing.
library(lwgeom) #version 0.2-5 for st_split (used in 0_Functions)
library(ggplot2) #version 3.2.1 for ggplot()
library(sf) #version 0.9-6 for read_sf() and many other things
library(tidyverse)
library(raster) #version 3.4-5 for extent()
library(NLMR) #version 1.0 for nlm_gaussianfield() to simulate
library(stars) #version 4.0.3 for stars() to vectorize
library(gridExtra)
library(ggpubr)
library(viridisLite) #version 0.3.0 for magma()

####1. Read Data####
#lc <- read_sf("./SpatialOutput/", layer = "lc_clean_cut") #"Las Conchas" only sim
#lc <- read_sf("./SpatialOutput/", layer = "all_clean_cut") #"all fires" sim
lc <- read_sf("./SpatialOutput/", layer = "lc_clean_cut") #"treeless" sim
#lc_raw <- read_sf("./SpatialInput/4.Spatial_Scenarios/", layer = "lc_raw")
perim <- read_sf("../../../GIS/Fires/Las Conchas/Perimeter/", layer="Las_Conchas_Perim")
scars <- read_sf("../../../GIS/Fires/Las Conchas/FireScarLocations/", layer="fs_in_las_conchas")


####2. Define original parameters, to compare tuning against####
#Set the simulated grid size
# xmax and ymax set the number of cells in the x and y direction
# increasing the grid size slows down the 'nlm_gaussianfied' function
xmax <- #30 converts from UTM 1m resolution to LANDSAT 30 m resolution
  round((extent(perim)@xmax - extent(perim)@xmin)/30,0)[[1]] #Grid is 964 pixels wide
ymax <- 
  round((extent(perim)@ymax - extent(perim)@ymin)/30,0)[[1]] #Grid is 1447 pixels tall
#Grid is 964*1447 = 1394908 pixels large = (*900 m2) 1255417200 m2 large = 125542 ha 1255 km2

#Initialize comparison data frame with observed original patch size distribution (after having filled holes)
prop_hs <- as.vector(round(sum(st_area(lc))/st_area(perim),3)[[1]])
df_original <- data.frame(sim = 0, Area_ha = sort(lc$Area_ha, decreasing = T), prop_hs = prop_hs)

####3. Run simulation for given set of parameters####
sims_list_spatial <- list()
df_list <- list()

#So far best combo for Las Conchas is 0.25/5/50
#Best "lc" combo with the cutoffs method is 0.6/0.6/10, USING THIS ONE FOR SUBMISSION
#Best "all" combo with the cutoffs method is 0.6/0.6/10, USING THIS ONE FOR SUBMISSION
#Best "treeless" combo with the cutoffs method is 
pct <- 0.25
nugget <- 5
magvar <- 50
#cutoffs: Original option is to set cutoff at ~prop_hs
#cutoffs = (1-(prop_hs + 0.02)) #everything above the threshold is "in"
#cutoffs: If nugget >0, set cutoff #2 to ~10% over the minimum needed
#0.354/0.6/0.6/10: c(0.01, 0.14, 0.75) *using for "lc"
#0.428/0.6/0.6/10: c(0.01, 0.18, 0.71) *using for "all"
#0.655/0.6/0.6/10: c(0.01, 0.18,0.50) 
#0.655/1,5/0/10: c(0.02, 0.18,0.50) 
#0.655/7/1/10: c(0.01, 0.19,0.50) 
#0.655/7/0.2,0.3/10: c(0.02, 0.18,0.50) *using 0.3 for "treeless"
cutoffs <- c(0.01,0.14,0.75) #needs to be defined outside of function below to work properly.
parms_text <- paste0("pct_", pct,"_nug_",nugget,"_mv_",magvar)

#Run one to check it out
#s = 1 #simulation number
#gc()
#out = random_cheese_patchdist(pct = pct, nugget = nugget, magvar = magvar, prop_hs = prop_hs)
#plot(out$geometry,col = out$Area_ha, border = "black", pal = rainbow(50), main = parms_text)

#Run ten and store the results
df_compare <- df_original
for(s in c(1:10)){
  gc()
  print(Sys.time())
  out <- raster_cheese_patchdist_multipart(pct = pct, nugget = nugget, magvar = magvar, prop_hs = prop_hs, 
                                           seed = s, cutoffs = cutoffs)
  
  #Test if the random cookie cutter of the Las Conchas perimeter picked up an accurate proportion of high severity
  area_sim <- sum(st_area(out))
  prop_hs_sim <- as.vector(round(area_sim/st_area(perim),3)[[1]])
  test_prop_hs <- abs( prop_hs_sim - prop_hs ) < 0.02
  if(test_prop_hs){#store results if proportion high severity within 2% of true value
    #double check but this test can probably be deprecated CHECKME if this doesn't produce any errors delete the test
    df_out <- data.frame(sim = s, Area_ha = sort(out$Area_ha, decreasing = T), prop_hs = prop_hs_sim) #store patch distribution
    df_compare <- rbind(df_compare, df_out) #add patch distribution to comparison data frame
    sims_list_spatial[[s]] <- out #store spatial layer in list
    print(paste(prop_hs_sim - prop_hs, ";  sim", s, "stored"))
    }
  else {print(paste("sim", s, "not stored"))}
  
}
df_list[[parms_text]] <- df_compare
write_csv(df_compare,paste0("./Output/Tuning_Sims/treeless/Theta1/",parms_text,".csv"))
write_rds(sims_list_spatial,paste0("./Output/Tuning_Sims/treeless/Theta1/",parms_text,".rds"))

####4. Dataframe mods and plotting####

csv_list <-  
  dir("./Output/Tuning_Sims/treeless/Theta1/") [grep("csv", dir("./Output/Tuning_Sims/treeless/Theta1/"))]
parms_list <- sub("\\.c.*", "", csv_list)

for(parms_text in parms_list){
#if doing manually: 
#parms_text <- parms_list[2]
df_compare <- read.csv(paste0("./Output/Tuning_Sims/treeless/Theta1/", parms_text, ".csv"))
sims_list_spatial <- read_rds(paste0("./Output/Tuning_Sims/treeless/Theta1/", parms_text, ".rds"))

df_compare$bin <- ceiling(log2(df_compare$Area_ha))

a <- ggplot() + 
  geom_sf(data = perim, fill = NA) + 
  geom_sf(data = lc, aes(fill = log2(Area_ha)), col = NA) +
  scale_fill_gradientn(colours=magma(7), breaks=seq(1:13)[c(TRUE, FALSE)], labels=(2^seq(1:13))[c(TRUE, FALSE)], limits = c(1, 13) ) +
  labs(title = "Las Conchas obs high severity", fill = "patch area (ha)") + 
  theme(axis.text.x = element_text(angle = 45), legend.position = "none")

b_plot <- sims_list_spatial[[1]] #pick a random simulation to plot
b <- tryCatch( ggplot() + 
  geom_sf(data = perim, fill = NA) + 
  geom_sf(data = b_plot, aes(fill = log2(Area_ha)), col = NA) +
    scale_fill_gradientn(colours=magma(7), breaks=seq(1:13)[c(TRUE, FALSE)], labels=(2^seq(1:13))[c(TRUE, FALSE)], limits = c(1, 13) ) +
    labs(title = "Las Conchas sim high severity", fill = "patch area (ha)") + 
    theme(axis.text.x = element_text(angle = 45), legend.position = "none")
) #There is an occasional error in plot generation; tryCatch just means "keep trying till it works"

c_plot1 <- df_compare[df_compare$sim == 0,] %>%
  count(bin)
c_plot2 <- df_compare[df_compare$sim > 0 & df_compare$sim <11,] %>%
  count(sim,bin)

c <- ggplot() + 
  #geom_bar(data = c_plot2, aes(x = bin, y = n), stat = "summary", fun = mean ) +
  stat_summary(geom = "crossbar", data = c_plot2, 
               aes(x = bin, y = n, group = bin, col = bin)) + 
  geom_point(data = c_plot1, aes(x = bin, y = n) ) + 
  scale_colour_gradientn(colours=magma(7), breaks=seq(1:13)[c(TRUE, FALSE)], labels=(2^seq(1:13))[c(TRUE, FALSE)], limits = c(1, 13) ) +
  #geom_crossbar(data = c_plot2, aes(x = bin, y = n, group = sim), stat = "summary", fun = mean_se ) +
  theme_bw() +
  scale_x_continuous(breaks = seq(1:13), minor_breaks = NULL, labels = paste(2^seq(1:13))) +
  labs(x = "patch area (ha)", y = "n patches", col = "patch area (ha)") +
  theme(legend.position = c(0.8,0.6))

d_plot1 <- c_plot2 %>%
  group_by(bin) %>%
  summarise(mean_n_sim = round(mean(n),0))

#manipulate data frame for residuals plot
d_plot2 <- merge.data.frame(d_plot1,c_plot1, all = TRUE)
d_plot2[is.na(d_plot2)] <- 0
d_plot2$resid <- d_plot2$mean_n_sim - d_plot2$n
d_plot2$prop_resid <- d_plot2$resid / d_plot2$n
d_plot2$prop_resid[!is.finite(d_plot2$prop_resid)] = 
  d_plot2$resid[!is.finite(d_plot2$prop_resid)]
d_plot2$residual <- ifelse(d_plot2$resid > 0, "overestimate", "underestimate")

m_resid <- lm(resid~bin,data = d_plot2)
m_summary <- data.frame(coef(summary(m_resid)))
names(m_summary) = c("est", "se", "t", "P")
m_summary$P = round(m_summary$P,3)
m_summary$P <- ifelse(m_summary$P == 0, "<0.001",m_summary$P)


d <- ggplot(d_plot2,aes(x = bin, y = prop_resid)) +
  geom_point(aes(col = residual))+
  scale_color_manual(values=c("blue","red"))+
  geom_smooth(method = lm, se = TRUE, col = "black") +
  annotate("text", x = 0, y = 1, hjust = 0,
           label = paste("int. P = ", m_summary$P[1], 
                           "\nslp. P = ", m_summary$P[2])
           ) +
  theme_bw() +
  scale_x_continuous(breaks = seq(1:13), minor_breaks = NULL, labels = paste(2^seq(1:13))) +
  labs(x = "patch area (ha)", y = "n patches") +
  theme(legend.position = c(0.8,0.2))
  
final <-
  tryCatch(ggarrange(
    ggarrange(a,b,ncol = 2, labels = c("a","b")),
    c, d, nrow = 3, labels = c("","c", "d")
  ) 
  ) #There is an occasional error in plot generation; tryCatch just means "keep trying till it works"

png(paste0("./Output/Tuning_Sims/treeless/Theta1/Figures/",parms_text,".png"), height = 12, width = 8, units = "in", res = 300)
print({
  tryCatch(

  annotate_figure(final,
                  top = text_grob(parms_text))

) #There is an occasional error in plot generation; tryCatch just means "keep trying till it works"
})
dev.off()
}