#Author: Jens Stevens; stevensjt@gmail.com
#This code looks at the intersection of fire scars and clusters with real and simulated data


####0. Libraries####
source("./Code/Official/0_Functions.R") #read in customized functions necessary for processing.
source("./Code/Official/0_RasterFunctions.R") #read in customized raster functions necessary for processing.
library(sf)
library(tidyverse)

####1. Read Data####
#scars <- read_sf("../../../GIS/Fires/Las Conchas/FireScarLocations/", layer="fs_in_las_conchas")
scars <- read_sf("./SpatialInput/5.Fire_Scars/", layer="fs_in_las_conchas")
clusters <- read_sf("../../../GIS/Fires/Las Conchas/FireScarLocations/", layer="fire_scar_clusters")
#lc <- read_sf("./SpatialOutput/", layer = "treeless_clean_cut")
#scenario <- "treeless"


####2. Study observed patterns first####
##2a. Spatial Layers
hsm1_single_raw <- #high severity model 1, Las Conchas severity only, unprocessed
  read_sf("./SpatialInput/4.Spatial_Scenarios", layer = "lc_raw")
hsm1_single_cleancut <- #high severity model 1, Las Conchas severity only, processed (removed holes+crumbs)
  read_sf("./SpatialOutput/", layer = "lc_clean_cut")
hsm2_single_raw <- #high severity model 2, all fires since '84, unprocessed
  read_sf("./SpatialInput/4.Spatial_Scenarios", layer = "all_raw")
hsm2_single_cleancut <- #high severity model 2, all fires since '84, processed
  read_sf("./SpatialOutput/", layer = "all_clean_cut")
hsm3_single_raw <- #high severity model 3, treeless area, unprocessed
  read_sf("./SpatialInput/4.Spatial_Scenarios", layer = "treeless_raw")
hsm3_single_cleancut <- #high severity model 3, treeless area, processed
  read_sf("./SpatialOutput/", layer = "treeless_clean_cut")
hs_models <- rbind(hsm1_single_raw,hsm2_single_raw,hsm3_single_raw)

fire_name_table <- #Get fire names as part of Fire_Atlas shapefile
  read_sf("./SpatialInput/1.Fire_Atlas", layer = "fires_that_intersect_LC") %>%
  st_transform(crs = st_crs(hsm1_single_cleancut)) %>%
  filter(Fire_Name %in% c("LAS CONCHAS", "CERRO GRANDE", "OSO", "DOME")) %>%
  select(Fire_Name, Year, geometry) %>%
  arrange(desc(Year)) %>%
  add_row(Fire_Name = "Cerro overlap", .before = 3) %>%
  add_row(Fire_Name = "Oso overlap", .before = 5) %>%
  add_row(Fire_Name = "Dome overlap", .before = 7) 

fire_name_table$geometry[3] <- #Intersect Cerro Grande with Las Conchas
  st_geometry(st_intersection(fire_name_table[1,],fire_name_table[2,]))
fire_name_table$geometry[5] <- #Intersect Oso with Las Conchas
  st_geometry(st_intersection(fire_name_table[1,],fire_name_table[4,]))
fire_name_table$geometry[7] <- #Intersect Dome with Las Conchas
  st_geometry(st_intersection(fire_name_table[1,],fire_name_table[6,]))

##2b. Statistical calculations for individual fires and overlap (Table 1)
fire_stats <- #initialize area table
  tibble(Fire = c("Las Conchas", "Cerro Grande","Cerro overlap", 
                  "Oso", "Oso overlap", "Dome", "Dome overlap"),
         Year = c("2011", "2000","", "1998","", "1996",""),
         Area = NA, Area_HS = NA, Prop_HS = NA)

fire_stats$Area <- #add fire perimeter area
  fire_name_table %>%
  st_area %>%
  units::set_units(ha) %>%
  as.vector %>%
  round(0)

fire_stats$Area_HS[1] <- #Las Conchas high severity area
  hsm1_single_raw %>%
  st_area %>% units::set_units(ha) %>% as.vector %>% round(0)
fire_stats$Area_HS[2] <- #Cerro Grande high severity area
  read_sf("../../../GIS/Fires/Parks Layers For Swiss Cheese/Vector", layer = "CERRO GRANDE") %>%
  st_area %>% units::set_units(ha) %>% as.vector %>% sum %>% round(0) 
fire_stats$Area_HS[3] <- #Cerro overlap high severity area
  read_sf("./SpatialInput/3.CBI_Vector_OtherHS", layer = "CERRO GRANDE") %>%
  st_area %>% units::set_units(ha) %>% as.vector %>% sum %>% round(0) 
fire_stats$Area_HS[4] <- #Cerro Grande high severity area
  read_sf("../../../GIS/Fires/Parks Layers For Swiss Cheese/Vector", layer = "OSO") %>%
  st_area %>% units::set_units(ha) %>% as.vector %>% sum %>% round(0) 
fire_stats$Area_HS[5] <- #Cerro overlap high severity area
  read_sf("./SpatialInput/3.CBI_Vector_OtherHS", layer = "OSO") %>%
  st_area %>% units::set_units(ha) %>% as.vector %>% sum %>% round(0) 
fire_stats$Area_HS[6] <- #Cerro Grande high severity area
  read_sf("../../../GIS/Fires/Parks Layers For Swiss Cheese/Vector", layer = "DOME") %>%
  st_area %>% units::set_units(ha) %>% as.vector %>% sum %>% round(0) 
fire_stats$Area_HS[7] <- #Cerro overlap high severity area
  read_sf("./SpatialInput/3.CBI_Vector_OtherHS", layer = "DOME") %>%
  st_area %>% units::set_units(ha) %>% as.vector %>% sum %>% round(0) 

fire_stats$Prop_HS <- round(fire_stats$Area_HS / fire_stats$Area, 2)
#write_csv(fire_stats,"./Tables/table1_firestats.csv")

##2c. Model stats (independent of simulations; precursor to Table 2)
scars_core <- scars[scars$Dewar == "n",] #option to refine the analysis and exclude Dewar plots in Valles Caldera

model_stats = data.frame(single = rep(NA,9), multi = NA, treeless = NA)
rownames(model_stats) <- c("hs_area","prop_hs","prop_dead_pred","n_dead_pred", "false_mortality", "false_survival", 
                   "classification_error", "Parameter_1", "Parameter_2")

model_stats["hs_area",] <- #Get area and percentage of each high severity model (hsm)
  hs_models %>% st_area %>% units::set_units(ha) %>% as.vector %>% round(0)
model_stats["prop_hs",] <- round(model_stats["hs_area",] / 61057, 3) #61057 (ha) is the area of Las Conchas

hits_1 <- grep(1,st_intersects(scars_core, hs_models[1,])) #which scars intersect model 1 (single burn)
hits_2 <- grep(1,st_intersects(scars_core, hs_models[2,])) #which scars intersect model 2 (multi burn)
hits_3 <- grep(1,st_intersects(scars_core, hs_models[3,])) #which scars intersect model 3 (treeless)

validate_hits <- scars_core %>% select(SampleID, cluster_id, hs_valid) 
validate_hits$hsm1 <- 0
validate_hits$hsm1[hits_1] <-1
validate_hits$hsm2 <- 0
validate_hits$hsm2[hits_2] <-1
validate_hits$hsm3 <- 0
validate_hits$hsm3[hits_3] <-1

hits <- as.data.frame(validate_hits)[,c("hsm1","hsm2","hsm3")]
hits$hs_valid <- scars_core$hs_valid
hits[,c("errors_1","errors_2","errors_3")] <-
  hits[,c("hsm1","hsm2","hsm3")] - hits$hs_valid

model_stats["n_dead_pred",] <- apply(hits[,c(1:3)], 2, sum)
model_stats["prop_dead_pred",] <- round(model_stats["n_dead_pred",] / nrow(scars_core), 3)
model_stats["false_mortality","single"] <- sum(hits$errors_1 == 1)
model_stats["false_survival","single"] <- sum(hits$errors_1 == -1)
model_stats["false_mortality","multi"] <- sum(hits$errors_2 == 1)
model_stats["false_survival","multi"] <- sum(hits$errors_2 == -1)
model_stats["false_mortality","treeless"] <- sum(hits$errors_3 == 1)
model_stats["false_survival","treeless"] <- sum(hits$errors_3 == -1)
model_stats["classification_error",] <- 
  round(colSums(model_stats[c("false_mortality","false_survival"),])/ nrow(scars_core),3)
  
#write.csv(model_stats,"./Tables/table2_modelstats_core.csv", row.names = TRUE)


####3. Iterate simulation of choice####
#So far best combo for Las Conchas is 0.25/5/50
#Update: 0.6/1.6/10 for lc (cutoffs 0.01, 0.14, 0.75) and all (cutoffs 0.01, 0.18, 0.71)
#7/0.3/10 for treeless (0.02, 0.18, 0.50)
pct <- 0.6
nugget <- 0.6
magvar <- 10
#cutoffs: 
cutoffs <- c(0.01,0.14,0.75) #needs to be defined outside of function below to work properly.

parms_text <- paste0(scenario,"_pct_", pct,"_nug_",nugget,"_mv_",magvar)

#Run one to check it out
#s = 1 #simulation number
#gc()
#out = random_cheese_patchdist(pct = pct, nugget = nugget, magvar = magvar, prop_hs = prop_hs)
#plot(out$geometry,col = out$Area_ha, border = "black", pal = rainbow(50), main = parms_text)

#Initialize comparison data frame with observed original patch size distribution (after having filled holes)
prop_hs <- as.vector(round(sum(st_area(lc))/st_area(perim),3)[[1]])
df_original <- data.frame(sim = 0, Area_ha = sort(lc$Area_ha, decreasing = T), prop_hs = prop_hs)

#Run one hundred and store the results
sims_list_spatial <- list()
df_list <- list()
df_compare <- df_original
for(s in c(1:100)){
  gc()
  print(Sys.time())
  out <- raster_cheese_patchdist_multipart(pct = pct, nugget = nugget, magvar = magvar, prop_hs = prop_hs, 
                                           seed = s, cutoffs = cutoffs)
  #out <- random_cheese_patchdist(pct = pct, nugget = nugget, magvar = magvar, prop_hs = prop_hs)
  
  
  #Test if the random cookie cutter of the Las Conchas perimeter picked up an accurate proportion of high severity
  area_sim <- sum(st_area(out))
  prop_hs_sim <- as.vector(round(area_sim/st_area(perim),3)[[1]])
  test_prop_hs <- abs( prop_hs_sim - prop_hs ) < 0.01
  if(test_prop_hs){#store results if proportion high severity within 1% of true value
    #CHECKME this test can probably be deprecated if runs are producing 100/100 within 1%
    df_out <- data.frame(sim = s, Area_ha = sort(out$Area_ha, decreasing = T), prop_hs = prop_hs_sim) #store patch distribution
    df_compare <- rbind(df_compare, df_out) #add patch distribution to comparison data frame
    sims_list_spatial[[s]] <- out #store spatial layer in list
    print(paste("sim", s, "stored"))
  }
  else {print(paste("sim", s, "not stored"))}
  
}
#df_list[[parms_text]] <- df_compare #Deprecated when running for analysis
write_csv(df_compare,paste0("./Output/Analysis_Ready/",parms_text,".csv"))
write_rds(sims_list_spatial,paste0("./large_files/",parms_text,".rds"))

####4. Assemble plot data####
df_compare <- read_csv("./Output/Analysis_Ready/treeless_pct_7_nug_0.3_mv_10.csv")
sims_list_spatial <- read_rds("./large_files/treeless_pct_7_nug_0.3_mv_10.rds")
#all_pct_0.6_nug_0.6_mv_10 [same for lc]
#treeless_pct_7_nug_0.3_mv_10

df_compare$bin <- ceiling(log2(df_compare$Area_ha))
baseline <- lc

sim_comparison <- df_compare %>%
  group_by(sim) %>%
  summarise(prop_scars = NA, prop_patch_scars = NA, prop_patch_scars_20 = NA,
            prop_clusters = NA, prop_cluster_scars = NA,prop_cluster_scars_20 = NA)

##CHECKME the premise is wrong here for intersecting the clusters, because a cluster can be intersected by multiple patches which artificially inflates the numbers. To solve this:
#lc_simple_intersect <- st_make_valid(st_combine(lc)) 
#below where lc etc are used as second term of st_intersection. Not implemented yet 2/18.

for(s in c(1:100)){ #about 2 minutes (more for treeless)
  sim <- st_make_valid(sims_list_spatial[[s]])
  sim <- sim[order(sim$Area_ha, decreasing = TRUE),]
  sim$sim = s
  sim$unique_patch_id = c(1:nrow(sim))
  sim_20 <- sim[c(1:20),]
  scar_hits <- st_intersection(scars,sim)
  scar_hits_20 <- st_intersection(scars,sim_20)
  patch_scar_hits <- length(unique(scar_hits$unique_patch_id))
  patch_scar_hits_20 <- length(unique(scar_hits_20$unique_patch_id))
  cluster_hits <- st_intersection(clusters,sim) #CHECKME for entirely within; this just does partial intersection
  cluster_hits_20 <- st_intersection(clusters, sim_20)
  patch_cluster_hits <- length(unique(cluster_hits$unique_patch_id))
  patch_cluster_hits_20 <-  length(unique(cluster_hits_20$unique_patch_id))
  
  sim_comparison$prop_scars[s] <- nrow(scar_hits)/nrow(scars)
  sim_comparison$prop_patch_scars[s] <- patch_scar_hits/nrow(sim)
  sim_comparison$prop_patch_scars_20[s] <- patch_scar_hits_20/nrow(sim_20)
  sim_comparison$prop_clusters[s] <- nrow(cluster_hits)/nrow(clusters)
  sim_comparison$prop_patch_clusters[s] <- patch_cluster_hits/nrow(sim)
  sim_comparison$prop_patch_clusters_20[s] <- patch_cluster_hits_20/nrow(sim_20)  
}

#Actual values
lc <- lc[order(lc$Area_ha, decreasing = TRUE),]
lc <- st_make_valid(lc) #treeless layer has a self-intersecting problem
lc$unique_patch_id = c(1:nrow(lc))
lc_20 <- lc[c(1:20),]
scar_hits_lc <- st_intersection(scars,lc)
scar_hits_20_lc <- st_intersection(scars,lc_20)
patch_scar_hits_lc <- length(unique(scar_hits_lc$unique_patch_id))
patch_scar_hits_20_lc <- length(unique(scar_hits_20_lc$unique_patch_id))
cluster_hits_lc <- st_intersection(clusters,lc) #CHECKME for entirely within; this just does partial intersection and furthermore it does it incorrectly if multiple patches intersect
cluster_hits_lc2 <- st_intersection(clusters,lc_simple_intersect) 
#Per Sean suggestion, just pick a tree at random from a cluster first?
cluster_hits_20_lc <- st_intersection(clusters, lc_20)
patch_cluster_hits_lc <- length(unique(cluster_hits_lc$unique_patch_id))
patch_cluster_hits_20_lc <-  length(unique(cluster_hits_20_lc$unique_patch_id))

observed <- data.frame(scenario = scenario)
observed$prop_scars_lc <- nrow(scar_hits_lc)/nrow(scars)
observed$prop_patch_scars_lc <- patch_scar_hits_lc/nrow(lc)
observed$prop_patch_scars_20_lc <- patch_scar_hits_20_lc/nrow(lc_20)
observed$prop_clusters_lc <- nrow(cluster_hits_lc)/nrow(clusters)
observed$prop_patch_clusters_lc <- patch_cluster_hits_lc/nrow(lc)
observed$prop_patch_clusters_20_lc <- patch_cluster_hits_20_lc/nrow(lc_20)  
observed$prop_hs <- as.vector(round(sum(st_area(lc))/st_area(perim),3)[[1]])
observed$observed <- "treeless"

####4. Build plots####

p1 <-
  ggplot(sim_comparison) +
  geom_histogram(aes(prop_scars), fill = NA,col = "black") +
  geom_vline(data = observed, aes(xintercept = prop_scars_lc, lty = observed)) +
  scale_linetype_manual(values = "dashed") +
  lims(x = c(0,1)) +
  labs(x = "proportion scars \nburned at HS", y = "n simulations")+
  theme_bw() +
  theme(legend.position = c(0.8,0.8), axis.text=element_text(size=14),
        axis.title=element_text(size=16))

p2 <- 
  ggplot(sim_comparison) +
  geom_histogram(aes(prop_patch_scars_20), fill = NA,col = "black", binwidth = 0.05) +
  geom_vline(data = observed, aes(xintercept = prop_patch_scars_20_lc, lty = observed)) +
  scale_linetype_manual(values = "dashed") +
  lims(x = c(0,1)) +
  labs(x = "proportion of 20 largest HS patches \nburning a scar", y = "n simulations")+
  theme_bw() +
  theme(legend.position = c(0.8,0.8), axis.text=element_text(size=14),
        axis.title=element_text(size=16))

p3 <- 
  ggplot(sim_comparison) +
  geom_histogram(aes(prop_clusters), fill = NA,col = "black") +
  geom_vline(data = observed, aes(xintercept = prop_clusters_lc, lty = observed)) +
  scale_linetype_manual(values = "dashed") +
  lims(x = c(0,1)) +
  labs(x = "proportion clusters \nburned at HS", y = "n simulations")+
  theme_bw() +
  theme(legend.position = c(0.2,0.8), axis.text=element_text(size=14),
        axis.title=element_text(size=16))

p4 <- 
  ggplot(sim_comparison) +
  geom_histogram(aes(prop_patch_clusters_20), fill = NA,col = "black", binwidth = 0.05) +
  geom_vline(data = observed, aes(xintercept = prop_patch_clusters_20_lc, lty = observed)) +
  scale_linetype_manual(values = "dashed") +
  lims(x = c(0,1)) +
  labs(x = "proportion of 20 largest HS patches \nburning a cluster of scars", y = "n simulations")+
  theme_bw() +
  theme(legend.position = c(0.8,0.8), axis.text=element_text(size=14),
        axis.title=element_text(size=16))

final <-
    ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2, labels = c("a","c", "b", "d"))

png(paste0("./Output/Analysis_Ready/final_treeless_v2.png"), height = 8, width = 10, units = "in", res = 300)
final
dev.off()


####5. Simulate variation in prop_hs####
#So far best combo for Las Conchas is 0.25/5/50
pct <- 0.25
nugget <- 5
magvar <- 50
prop_hs <- 0.5
prop_range <- c(0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1, 0.08, 0.06, 0.04)
prop_range_full <- rep(prop_range, each = 10)
#rev(prop_range)

#Run ten and store the results
df_compare <- df_original
df_compare$prop_hs_sim <- 0.354

for(s in 1:length(prop_range_full)){
  gc()
  prop_hs = prop_range_full[s]
  out <- 
    random_cheese_patchdist(pct = pct, nugget = nugget, magvar = magvar, prop_hs = prop_hs)
  
  #Test if the random cookie cutter of the Las Conchas perimeter picked up an accurate proportion of high severity
  area_sim <- sum(st_area(out))
  prop_hs_sim <- as.vector(round(area_sim/st_area(perim),3)[[1]])
  test_prop_hs <- abs( prop_hs_sim - prop_hs ) < 0.035
  print(paste("simulated prop_hs", prop_hs_sim, "; target prop_hs", prop_hs))
  if(test_prop_hs){#store results if proportion high severity within 3% of true value
    #double check but this test can probably be deprecated
    df_out <- data.frame(sim = s, Area_ha = sort(out$Area_ha, decreasing = T), prop_hs = prop_hs, prop_hs_sim = prop_hs_sim) #store patch distribution
    df_compare <- rbind(df_compare, df_out) #add patch distribution to comparison data frame
    sims_list_spatial_prop_hs[[s]] <- out #store spatial layer in list
    print(paste("sim", s, "stored"))
  }
  else {print(paste("sim", s, "not stored"))}
  
}
#df_list[[parms_text]] <- df_compare #Deprecated when running for analysis
write_csv(df_compare,paste0("./Output/Analysis_Ready/prop_hs_comparison.csv"))
write_rds(sims_list_spatial_prop_hs,paste0("./Output/Analysis_Ready/prop_hs_comparison.rds"))

#df_compare <- read_csv("./Output/Analysis_Ready/prop_hs_comparison.csv")
#df_compare$bin <- ceiling(log2(df_compare$Area_ha)) #deprecated
#baseline <- lc #deprecatrd

sim_comparison_prop_hs <- df_compare[df_compare$sim>0,] %>%
  group_by(sim) %>%
  summarise(prop_hs_bin = mean(prop_hs), prop_hs = mean(prop_hs_sim), prop_scars = NA)

for(s in c(1:nrow(sim_comparison_prop_hs))){
  sim <- st_make_valid(sims_list_spatial_prop_hs[[s]])
  sim <- sim[order(sim$Area_ha, decreasing = TRUE),]
  sim$sim = s
  sim$prop_hs = sim_comparison_prop_hs$prop_hs[s]
  #sim_20 <- sim[c(1:20),]
  scar_hits <- st_intersection(scars,sim)
  #scar_hits_20 <- st_intersection(scars,sim_20)
  #patch_scar_hits <- length(unique(scar_hits$unique_patch_id))
  #patch_scar_hits_20 <- length(unique(scar_hits_20$unique_patch_id))
  #cluster_hits <- st_intersection(clusters,sim) #CHECKME for entirely within; this just does partial intersection
  #cluster_hits_20 <- st_intersection(clusters, sim_20)
  #patch_cluster_hits <- length(unique(cluster_hits$unique_patch_id))
  #patch_cluster_hits_20 <-  length(unique(cluster_hits_20$unique_patch_id))
  
  sim_comparison_prop_hs$prop_scars[s] <- nrow(scar_hits)/nrow(scars)
  #sim_comparison$prop_patch_scars[s] <- patch_scar_hits/nrow(sim)
  #sim_comparison$prop_patch_scars_20[s] <- patch_scar_hits_20/nrow(sim_20)
  #sim_comparison$prop_clusters[s] <- nrow(cluster_hits)/nrow(clusters)
  #sim_comparison$prop_patch_clusters[s] <- patch_cluster_hits/nrow(sim)
  #sim_comparison$prop_patch_clusters_20[s] <- patch_cluster_hits_20/nrow(sim_20)  
}

prop_hs_plot <-
  ggplot(sim_comparison_prop_hs, aes(x = prop_hs, y = prop_scars)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  geom_smooth( method = "lm") +
  geom_point(data = observed, aes(x = prop_hs, y = prop_scars_lc), col = "red") +
  geom_text(data = observed, aes(x = prop_hs, y = prop_scars_lc, label=observed),hjust=0, vjust=0) +
  labs(y = "proportion scars burned at high severity", x = "proportion high severity") +
  theme_bw() +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16))

png(paste0("./Output/Analysis_Ready/prop_hs.png"), height = 6, width = 6, units = "in", res = 300)
prop_hs_plot
dev.off()

####6. How many fires to cover all the scars####
target_rows <- which (sim_comparison_prop_hs$prop_hs_bin==0.35)
cumulative_df <- sim_comparison_prop_hs[target_rows,]
for(s in target_rows){
  sim <- st_make_valid(sims_list_spatial_prop_hs[[s]])
  if(s == target_rows[1]){
    cumulative = sim
    r = 1
  } else { cumulative = rbind(cumulative, sim); r = r+1 }
  
  scar_hits_all <- st_intersection(scars,cumulative)
  cumulative_df$total_prop_scars[r] <- length(unique(scar_hits_all$SampleID))/nrow(scars)
 
}

cumulative_plot <-
  ggplot(cumulative_df, aes(x = sim-20, y = total_prop_scars)) +
  geom_col(fill = NA, col = "black") +
  scale_x_continuous(breaks = c(1:10)) +
  labs(y = "proportion scars burned at high severity", x = "n successive fires", title = "35% high severity") + 
  theme_bw() +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16))

png(paste0("./Output/Analysis_Ready/cumulative.png"), height = 6, width = 6, units = "in", res = 300)
cumulative_plot
dev.off()
