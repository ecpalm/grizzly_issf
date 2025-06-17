
################################################################################
#
# Journal: Conservation Science and Practice
#
# Article title: "Changing grizzly bear space use and functional connectivity in
#                response to human disturbance in the southern Canadian Rocky 
#                Mountains"
#
# Article authors: Eric Palm, Clayton Apps, Tal Avgar, Melanie Dickie, 
#                  Bruce McLellan, Joseph Northrup, Michael Sawaya, 
#                  Julie Turner, Jesse Whittington, Erin Landguth, 
#                  Katherine Zeller, Clayton Lamb
#
# Script description: Simulate utilization distribution from a fitted 
#                     mixed-effects iSSF model (spring season model)
#
# Script author: Eric Palm
#
# Portions of this script were adapted from Jesse Whittington's simulation code
# in Whittington et al., "Towns and trails drive carnivore movement behaviour,
# resource selection, and connectivity". 2022. Movement Ecology. 
# https://doi.org/10.1186/s40462-022-00318-5
#
################################################################################

# Load packages
sapply(
  c("amt", "caret", "data.table", "doParallel", "dplyr", "foreach", "glmmTMB",
    "parallel", "raster", "Rfast", "sf", "sfheaders", 
    "stringr", "terra", "tibble", "tidyr", "tictoc"), 
  require, character.only = T)

# Import function for simulating coefficients from full variance-covariance
# matrix of the model, which only simulates from fixed-effects for the 
# movement covariates because incorporating random effects variance for those
# variables was resulting in nonsensical values (e.g., negative concentration
# parameters for von Mises distribution) in some instances. This is a known
# issue with iSSFs.
source("scripts/functions/simulate_coefs_random.R")

# Define season for model and simulation
m_season <- "spring"

# Import model fixed effects coefficients
m_fixef <- readRDS(stringr::str_c("outputs/models/fixef_", m_season, ".rds")) 

# Import input issf data frame
df_issf <- readRDS("data/df_issf.rds")

# Calculate means of all variables to use for centering/scaling the raster values so the 
# raster values are scaled the same way as the variables were in the models.
# Make sure not to select the binary variables or the movement variables.
# Save the values as a list for later.
# Manually include the "_start" variables because they only appear in interactions
# and not on their own. Probably a more automated way to do this step.
means <- df_issf %>%
  dplyr::filter(season == m_season) %>% 
  dplyr::select(dplyr::any_of(rownames(data.frame(m_fixef))),
                greenness_broad_end,
                canopy_cover_start, greenness_start, tri_start, elev_start,
                -mine_end, -crossing_lake_reservoir, -crossing_hwy, -crossing_town, -crossing_rock_ice,
                -log_sl_km, -cos_ta) %>% 
  dplyr::summarize(dplyr::across(dplyr::everything(), ~mean(.))) %>% 
  as.list()

# Repeat the same step as above but calculate the standard deviation for scaling the raster values.
sds <- df_issf %>%
  dplyr::filter(season == m_season) %>% 
  dplyr::select(dplyr::any_of(rownames(data.frame(m_fixef))),
                greenness_broad_end,
                canopy_cover_start, greenness_start, tri_start, elev_start,
                -mine_end, -crossing_lake_reservoir, -crossing_hwy, -crossing_town, -crossing_rock_ice,
                -log_sl_km, -cos_ta) %>% 
  dplyr::summarize(dplyr::across(dplyr::everything(), ~sd(.))) %>% 
  as.list()

# This saves the scaling parameters of the distance to spring start variable 
# from the input dataset, so we can apply these to center and scale this variable
# while simulating.
d_spring_start_scaled <- df_issf %>%
  dplyr::filter(season == m_season) %>% 
  dplyr::select(d_spring_start_end_sim = d_spring_start_end) %>% 
  caret::preProcess(., method = c("center", "scale"))

###################### CENTER AND SCALE RASTERS ###################
# Import and scale all rasters in model and save to temporary folder.
# This way we can call these rasters using terra package directly within the 
# foreach loop below, because we can't import terra spatRaster objects into the 
# foreach loop (terra uses C++ which doesn't work with parallelization)
# Probably could automate this step to avoid hard coding.
terra::rast("spatial/raster/greenness.tif") %>% 
  terra::scale(., center = means$greenness_end, scale = sds$greenness_end) %>% 
  `names<-`("greenness_end") %>% 
  terra::writeRaster(., "spatial/raster/temp_sim/greenness_end.tif", datatype = "FLT4S", overwrite = T)

terra::rast("spatial/raster/greenness.tif") %>% 
  terra::scale(., center = means$greenness_start, scale = sds$greenness_start) %>% 
  `names<-`("greenness_start") %>% 
  terra::writeRaster(., "spatial/raster/temp_sim/greenness_start.tif", datatype = "FLT4S", overwrite = T)

terra::rast("spatial/raster/elev.tif") %>% 
  terra::scale(., center = means$elev_end, scale = sds$elev_end) %>% 
  `names<-`("elev_end") %>% 
  terra::writeRaster(., "spatial/raster/temp_sim/elev_end.tif", datatype = "FLT4S", overwrite = T)

terra::rast("spatial/raster/elev.tif") %>% 
  terra::scale(., center = means$elev_start, scale = sds$elev_start) %>% 
  `names<-`("elev_start") %>% 
  terra::writeRaster(., "spatial/raster/temp_sim/elev_start.tif", datatype = "FLT4S", overwrite = T)

terra::rast("spatial/raster/tri.tif") %>% 
  terra::scale(., center = means$tri_end, scale = sds$tri_end) %>% 
  `names<-`("tri_end") %>% 
  terra::writeRaster(., "spatial/raster/temp_sim/tri_end.tif", datatype = "FLT4S", overwrite = T)

terra::rast("spatial/raster/tri.tif") %>% 
  terra::scale(., center = means$tri_start, scale = sds$tri_start) %>% 
  `names<-`("tri_start") %>% 
  terra::writeRaster(., "spatial/raster/temp_sim/tri_start.tif", datatype = "FLT4S", overwrite = T)

terra::rast("spatial/raster/canopy_cover.tif") %>% 
  terra::scale(., center = means$canopy_cover_end, scale = sds$canopy_cover_end) %>% 
  `names<-`("canopy_cover_end") %>% 
  terra::writeRaster(., "spatial/raster/temp_sim/canopy_cover_end.tif", datatype = "FLT4S", overwrite = T)

terra::rast("spatial/raster/canopy_cover.tif") %>% 
  terra::scale(., center = means$canopy_cover_start, scale = sds$canopy_cover_start) %>% 
  `names<-`("canopy_cover_start") %>% 
  terra::writeRaster(., "spatial/raster/temp_sim/canopy_cover_start.tif", datatype = "FLT4S", overwrite = T)

terra::rast("spatial/raster/northness.tif") %>% 
  terra::scale(., center = means$northness_end, scale = sds$northness_end) %>% 
  `names<-`("northness_end") %>% 
  terra::writeRaster(., "spatial/raster/temp_sim/northness_end.tif", datatype = "FLT4S", overwrite = T)

terra::rast("spatial/raster/greenness_broad.tif") %>% 
  terra::scale(., center = means$greenness_broad_end, scale = sds$greenness_broad_end) %>% 
  `names<-`("greenness_broad_end") %>% 
  writeRaster(., "spatial/raster/temp_sim/greenness_broad_end.tif", datatype = "FLT4S", overwrite = T)

# This applies the negative exponential decay for d_highways before scaling.
terra::rast("spatial/raster/d_hwy.tif") %>% 
  terra::app(., function(x) 1 - exp(x*-.01)) %>% 
  terra::scale(., center = means$d_hwy_decay_end, scale = sds$d_hwy_decay_end) %>% 
  `names<-`("d_hwy_decay_end") %>% 
  terra::writeRaster(., "spatial/raster/temp_sim/d_hwy_decay_end.tif", datatype = "FLT4S", overwrite = T)

# This applies the negative exponential decay for d_roads before scaling.
terra::rast("spatial/raster/d_road.tif") %>% 
  terra::app(., function(x) 1 - exp(x*-.01)) %>% 
  terra::scale(., center = means$d_road_decay_end, scale = sds$d_road_decay_end) %>% 
  `names<-`("d_road_decay_end") %>% 
  terra::writeRaster(., "spatial/raster/temp_sim/d_road_decay_end.tif", datatype = "FLT4S", overwrite = T)

####################################################################################################
# Import a study area raster (1 = in study, 0 = out)
r_study <- raster::raster("spatial/raster/r_study.tif")

# Creates a mask for TRI so we don't start animals in the steepest terrain
# which we define as > 2 standard deviations from the mean raster value
terra::rast("spatial/raster/temp_sim/tri_end.tif") %>% 
  terra::classify(., cbind(2, Inf, NA)) %>% 
  terra::writeRaster(., "spatial/raster/temp_sim/tri_mask.tif", datatype = "FLT4S", overwrite = T)

###############################
# Import semi permeable crossing shapefiles

sf_crossing_hwy <- sf::read_sf("spatial/shp/hwy.shp") 
sf_crossing_lake_reservoir <- sf::read_sf("spatial/shp/barrier_lake.shp") 
sf_crossing_town <- sf::read_sf("spatial/shp/barrier_town.shp") 
sf_crossing_rock_ice <- sf::read_sf("spatial/shp/barrier_rock_ice.shp") 
###############################

###############################
# Extract parameters for turning angle (von Mises concentration parameter, kappa) and
# step length (gamma scale and shape) from the input dataset and correct season.
# Use these to create random (proposed) steps in simulation.

params_move <-
  df_issf %>%
  dplyr::filter(season == m_season,
                case_ == T)  %>% 
  dplyr::summarize(sl_km_scale = amt::fit_distr(.$sl_km, "gamma") %>% .$params %>% .$scale,
         sl_km_shape = amt::fit_distr(.$sl_km, "gamma") %>% .$params %>% .$shape,
         ta_kappa = amt::fit_distr(.$ta_, "vonmises") %>% .$params %>% .$kappa) %>% 
  dplyr::mutate(i_strata = 1)


###########################################################
# Define the simulation start extent and start zones.
# This layer is the top 40% (bins 7, 8, 9, 10) of predicted values
# from a second-order RSF fit to spring grizzly bear data, and then open water, 
# mines, snow and ice, and towns are masked out so
# we don't start simulated bears in these areas. We also make sure not to start
# them on extremely steep slopes, using the TRI mask created above. 

hab_good <- terra::rast("spatial/raster/start_zones_2nd_order_spring_masked_top_4.tif") %>%
  terra::mask(., rast("spatial/raster/temp_sim/tri_mask.tif")) %>%
  raster::raster()


# This is the simulation start extent. It is the full study area with a -10km
# buffer, so simulated bears can move outwards and still be in the study area.
# This step might not really be necessary.
start_extent <- sf::read_sf("spatial/shp/ssf_simulation_start_extent.shp")

# We don't allow simulated bears to enter pixels with terrain ruggedness above 26
# (i.e., we set TRI for those pixels to NA) because otherwise, with the negative 
# coefficient for log_sl:tri_end, they can get "stuck" in cliff areas, 
# as the model says that step lengths are nearly zero if terrain is that rugged.
# So here we save the scaled equivalent of TRI = 26 for later.
tri_scaled <- 
  df_issf %>%
  dplyr::filter(season == m_season) %>% 
  dplyr::select(tri_end) %>% 
  caret::preProcess(., method = c("center", "scale")) 

tri_max <- predict(tri_scaled,
                   tibble::tibble(tri_end = 26)) %>%
  dplyr::pull(tri_end)

###########################################
####      Simulate movement paths      ####
###########################################

# Set simulation parameters. Values in this example are small so the simulation
# will run quickly and (hopefully) successfully.

n_mega <- 5    # Number of "mega" loops, equal to the number of cores. Used 100 for publication.      
n_start_locs <- 500   # Number of random paths (animals) within one mega loop. 
                      # Used 15000 for publication. Must be even.
n_rand <- 5             # Number of random locations to generate for each step. Used 50 for publication.
n_days <- 10            # Number of days to simulate for each path. Used 68, 61, and 47 days for spring,
                        # summer and fall seasons, respectively, for publication.
i_hours <- rep(seq(0, 18, by = 6), n_days) # One location every six hours for 'n_days' days.
nt <- length(i_hours)   # number of steps to simulate for each path.

# Sequence variable to be used for the mega loops
iseq <- 1:n_mega

# Define the number of cores for the foreach loop and initiate the cluster
cl <- parallel::makeCluster(n_mega)
doParallel::registerDoParallel(cl)

# Create a log text file so you can track progress 
writeLines(c(" "), paste0("outputs/log_sim/log_sim_", m_season, ".txt"))

# Use tictoc to record the overall time for the simulation.
tictoc::tic()

# Start the foreach loop; include all packages needed within loop. Need to call 
# terra SpatRaster files inside the loop from 'temp_sim' folder created above.
path_mega <- foreach::foreach(yyy=iseq, .combine = rbind,
                                 .packages = c("sf", "Rfast", "dplyr", "tidyr", "data.table", 
                                               "sfheaders", "terra", "raster", "caret", "stringr")) %dopar% 
  {
    # Create an empty data frame with the necessary columns, including path id, stratum,
    # previous and current locations.
    df_path <- dplyr::tibble(path_id = 1:n_start_locs,
                             i_strata = 1,
                             x_prev = as.numeric(NA),
                             y_prev = as.numeric(NA),
                             x = as.numeric(NA),
                             y = as.numeric(NA)) 
    
    # Create an excess of random start locations and subsample desired amount 
    # from inside the start extent.
    xy <- raster::sampleRandom(hab_good, n_start_locs*3, method = "random", na.rm = T, xy = T)[,-3] %>% 
      sfheaders::sf_point() %>% 
      sf::st_set_crs(., 26911) %>% 
      dplyr::filter(lengths(sf::st_intersects(., start_extent)) > 0) %>% 
      dplyr::slice_sample(., n = n_start_locs)
    
    # Create an initial location for each simulated path using the start locations above
    # Create a random initial bearing
    # Then simulate a set of coefficients using the 'sim_coef_ran' function, which 
    # draws from multivariate normal distributions of each coefficient in the model,
    # and includes random slope variation (except for movement variables in this case).
    # Each path (or simulated animal) has its own set of coefficients throughout the simulation.
    tmp_0 <- df_path %>% 
      dplyr::mutate(i_strata = 1,
                    x_prev = st_coordinates(xy)[,1],
                    y_prev = st_coordinates(xy)[,2],
                    direction_prev = round(runif(n_start_locs, -pi, pi), 3),
                    sex_01 = c(rep(0, n_start_locs*.5), rep(1, n_start_locs*.5)),
                    good = 1) %>% 
      dplyr::bind_cols(., sim_coef_ran(m_season, n = nrow(.))) %>% 
      dplyr::inner_join(., params_move)
    
    # start simulations by each individual step across all animals
    path_list <- vector(mode = 'list', length = nt + 1)
    
    for (i in 1:nt){
      
      # Add each step to the log file created earlier to monitor progress 
      cat(paste("Starting iteration: mega loop", yyy, "smaller loop",  i, lubridate::now(), "\n"),
          file = paste0("outputs/log_sim/log_sim_", m_season, ".txt"), append=TRUE)
      
      # Only keeps animals that are still in the study area
      tmp <- tmp_0 %>% dplyr::filter(good == 1)
      
      # Start the n_keep loop (won't run if animal has left the study area)
      n_keep <- nrow(tmp)
      if(n_keep > 0){
        
        # Add the simulated strata number (starting at 1 and going to nt)
        tmp <- tmp %>% mutate(i_strata = i)
        
        # At all the start locations, extract the covariates that are used in interactions with the 
        # movement parameters (log_sl and cos_ta), so the distribution of step lengths and turning angles
        # depends on the habitat conditions at the start of steps.
        # Gamma distribution for the step lengths has shape and scale parameters.
        # Only using log_sl_km (modifies the gamma distribution shape parameter), but one could also use sl_km.
        # The code for modifying sl_km (scale) is a bit more complicated but can be found in amt's source code on github.
        # This code chunk generates all the proposed (candidate) steps.
        tmp_2 <- tmp %>%
          dplyr::bind_cols(., terra::extract(list.files("spatial/raster/temp_sim", full.names = T, pattern = "start") %>% terra::rast(), 
                                             cbind(.$x_prev, .$y_prev))) %>% 
          dplyr::rowwise() %>% 
          # proposed step lengths
          dplyr::mutate(step = list(rgamma(n = n_rand,
                                           scale = sl_km_scale,
                                           shape = sl_km_shape + `log_sl_km` +
                                             tri_start*`log_sl_km:tri_start` +
                                             canopy_cover_start*`log_sl_km:canopy_cover_start` +
                                             greenness_start*`log_sl_km:greenness_start` +
                                             elev_start*`log_sl_km:elev_start` +
                                             `log_sl_km:sex_01`*sex_01 +
                                             tri_start*`log_sl_km:tri_start:sex_01`*sex_01 +
                                             canopy_cover_start*`log_sl_km:canopy_cover_start:sex_01`*sex_01 +
                                             greenness_start*`log_sl_km:greenness_start:sex_01`*sex_01 +
                                             elev_start*`log_sl_km:elev_start:sex_01`*sex_01))) %>% 
          tidyr::unnest(step) %>% 
          dplyr::filter(!is.na(step), step > 0) %>%
          # proposed concentration parameters for von Mises distribution
          dplyr::mutate(conc = ta_kappa + `cos_ta` + log(step)*`cos_ta:log_sl_km` + 
                          `cos_ta`*sex_01 + log(step)*`cos_ta:log_sl_km:sex_01`*sex_01,
                        mu = dplyr::if_else(conc < 0, pi, 0)) %>% 
          dplyr::rowwise() %>% 
          dplyr::mutate(angle = as.numeric(Rfast::rvonmises(1, m=mu, k = abs(conc)))) %>% 
          dplyr::ungroup() %>% 
          dplyr::mutate(angle = if_else(angle > pi, angle - (2 * pi), angle),
                        direction = direction_prev + angle,
                        direction = dplyr::if_else(direction < -2*pi, direction + 2*pi,
                                                   dplyr::if_else(direction > 2*pi, direction - 2*pi, direction)),
                        dx = step*1000 * cos(direction), # *1000 converts meters to km
                        dy = step*1000 * sin(direction), # *1000 converts meters to km
                        x = round(x_prev + dx), 
                        y = round(y_prev + dy),
                        # determines whether candidate locations are within study area
                        in_study = terra::extract(terra::rast("spatial/raster/r_study.tif"), cbind(x, y))[,1]) %>%  
          dplyr::group_by(path_id) %>% 
          # Calculates proportion of locs are within study area
          dplyr::mutate(prop_in_study = sum(in_study) / n_rand, 
                        # Calculates the distance to first location in spring for that path
                        d_spring_start_end_sim = sqrt((x - dplyr::first(x))^2 + (y - dplyr::first(y))^2)) %>% 
          dplyr::ungroup() %>% 
          dplyr::filter(in_study == 1) %>% # only keeps locations within the study area
          # centers/scales the distance to den variable based on saved parameters from before
          predict(d_spring_start_scaled, .) 
        
        ############# CREATE LINE SEGMENTS FOR CROSSING VARIABLES #########
        # Use the 'data.table' and 'sfheaders' packages to quickly create line segments for crossing variables.
        dt <- data.table::as.data.table(tmp_2 %>% dplyr::select(x_prev:y) %>% dplyr::mutate(id = dplyr::row_number()))
        
        ## To use `sfheaders` the data needs to be in long form.
        dt1 <- dt[, .(id, lon = x_prev, lat = y_prev)]
        dt2 <- dt[, .(id, lon = x, lat = y)]
        
        ## Add on a 'sequence' variable so we know which one comes first.
        dt1[, seq := 1L ]
        dt2[, seq := 2L ]
        
        ## put back together
        dt <- data.table::rbindlist(list(dt1, dt2), use.names = TRUE)
        data.table::setorder(dt, id, seq)
        
        # Create an sf line object, where every individual step is a line    
        line_sf <- sfheaders::sf_linestring(
          obj = dt
          , x = "lon"
          , y = "lat"
          , linestring_id = "id"
        ) %>% 
          sf::st_set_crs(26911) # UTM Zone 11 North
        ########################################################## 
        
        # Get probability of each step, using the "end" covariates and the line segment intersections.
        # This step should be more automated to avoid manual hard coding.
        # Then calculate the exponentiated linear predictor by multiplying extracted covariate values
        # by the simulated coefficients for that covariate.
        # Then calculate the relative probability of each candidate step being chosen 
        tmp_3 <- 
          tmp_2 %>%
          dplyr::bind_cols(., terra::extract(list.files("spatial/raster/temp_sim", full.names = T, pattern = "end") %>% terra::rast(), 
                                             cbind(.$x, .$y)) %>%
                             dplyr::rename_with(~stringr::str_c(., "_sim"), dplyr::everything())) %>%
          # Determine whether candidate steps cross (intersect with) crossings
          dplyr::mutate(crossing_hwy_sim = as.numeric(lengths(sf::st_intersects(line_sf, sf_crossing_hwy)) > 0),
                        crossing_lake_reservoir_sim = as.numeric(lengths(sf::st_intersects(line_sf, sf_crossing_lake_reservoir)) > 0),
                        crossing_town_sim = as.numeric(lengths(sf::st_intersects(line_sf, sf_crossing_town)) > 0),
                        crossing_rock_ice_sim = as.numeric(lengths(sf::st_intersects(line_sf, sf_crossing_rock_ice)) > 0),
                        mine_end_sim = terra::extract(terra::rast("spatial/raster/mine.tif"), cbind(tmp_2$x, tmp_2$y))[,1],
                        lake_reservoir_end_sim = terra::extract(terra::rast("spatial/raster/barrier_lake_reservoir.tif"), cbind(tmp_2$x, tmp_2$y))[,1],
                        tri_end_sim = dplyr::if_else(tri_end_sim > tri_max, NA_real_, tri_end_sim), # set extracted value to NA if ends on a cliff
                        crossing_lake_reservoir_sim = dplyr::na_if(lake_reservoir_end_sim, 1), # set extracted value to NA if ends in a lake
                        
                        # Calculate the relative probability of selection by
                        # exponentiating the linear predictor
                        exp_lp = exp(canopy_cover_end_sim*canopy_cover_end + canopy_cover_end_sim^2*`I(canopy_cover_end^2)` +
                                       d_spring_start_end_sim*d_spring_start_end +
                                       canopy_cover_end_sim*elev_end_sim*`canopy_cover_end:elev_end` +
                                       greenness_end_sim*canopy_cover_end_sim*`greenness_end:canopy_cover_end` +
                                       greenness_end_sim*elev_end_sim*`greenness_end:elev_end` +
                                       d_hwy_decay_end_sim*d_hwy_decay_end +
                                       d_road_decay_end_sim*d_road_decay_end +
                                       d_hwy_decay_end_sim*greenness_broad_end_sim*`d_hwy_decay_end:greenness_broad_end` +
                                       d_road_decay_end_sim*greenness_broad_end_sim*`d_road_decay_end:greenness_broad_end` +
                                       greenness_end_sim*greenness_end + 
                                       northness_end_sim*northness_end +
                                       elev_end_sim*elev_end + elev_end_sim^2*`I(elev_end^2)` +
                                       mine_end_sim*mine_end +
                                       tri_end_sim*tri_end + tri_end_sim^2*`I(tri_end^2)` +
                                       crossing_hwy_sim*crossing_hwy +
                                       crossing_hwy_sim*greenness_broad_end_sim*`crossing_hwy:greenness_broad_end` +
                                       crossing_lake_reservoir_sim*crossing_lake_reservoir +
                                       crossing_town_sim*crossing_town +
                                       crossing_rock_ice_sim*crossing_rock_ice +
                                       crossing_rock_ice_sim*greenness_broad_end_sim*`crossing_rock_ice:greenness_broad_end` +
                                       
                                       # same variables as above but with sex (sex_01) interactions (male = 1, female = 0)
                                       canopy_cover_end_sim*sex_01*`canopy_cover_end:sex_01` + 
                                       canopy_cover_end_sim^2*sex_01*`I(canopy_cover_end^2):sex_01` +
                                       d_spring_start_end_sim*sex_01*`d_spring_start_end:sex_01` +
                                       canopy_cover_end_sim*elev_end_sim*sex_01*`canopy_cover_end:elev_end:sex_01` +
                                       greenness_end_sim*canopy_cover_end_sim*sex_01*`greenness_end:canopy_cover_end:sex_01` +
                                       greenness_end_sim*elev_end_sim*sex_01*`greenness_end:elev_end:sex_01` +
                                       d_hwy_decay_end_sim*sex_01*`d_hwy_decay_end:sex_01` +
                                       d_road_decay_end_sim*sex_01*`d_road_decay_end:sex_01` +
                                       d_hwy_decay_end_sim*greenness_broad_end_sim*sex_01*`d_hwy_decay_end:greenness_broad_end:sex_01` +
                                       d_road_decay_end_sim*greenness_broad_end_sim*sex_01*`d_road_decay_end:greenness_broad_end:sex_01` +
                                       greenness_end_sim*sex_01*`greenness_end:sex_01` + 
                                       northness_end_sim*sex_01*`northness_end:sex_01` +
                                       elev_end_sim*sex_01*`elev_end:sex_01` + 
                                       elev_end_sim^2*sex_01*`I(elev_end^2):sex_01` +
                                       mine_end_sim*sex_01*`mine_end:sex_01` +
                                       tri_end_sim*sex_01*`tri_end:sex_01` + 
                                       tri_end_sim^2*sex_01*`I(tri_end^2):sex_01` +
                                       crossing_hwy_sim*sex_01*`crossing_hwy:sex_01` +
                                       crossing_hwy_sim*greenness_broad_end_sim*sex_01*`crossing_hwy:greenness_broad_end:sex_01` +
                                       crossing_lake_reservoir_sim*sex_01*`crossing_lake_reservoir:sex_01` +
                                       crossing_town_sim*sex_01*`crossing_town:sex_01` +
                                       crossing_rock_ice_sim*sex_01*`crossing_rock_ice:sex_01` +
                                       crossing_rock_ice*sex_01*greenness_broad_end_sim*`crossing_rock_ice:greenness_broad_end:sex_01`))  %>%
          dplyr::filter(exp_lp > 0 & is.finite(exp_lp))
        
        # For each path, use the slice_sample function (weighted by relative probability of selection)
        # to probabilistically select one location from the candidate set of locations. 
        # If more than 60% of the candidate locations are outside of the study area, terminate the path. 
        if (nrow(tmp_3) > 0){
          tmp_4 <- tmp_3 %>% 
            dplyr::group_by(path_id) %>% 
            dplyr::slice_sample(n = 1, weight_by = exp_lp) %>% 
            dplyr::ungroup() %>% 
            dplyr::arrange(path_id) %>% 
            dplyr::mutate(good = dplyr::if_else(prop_in_study >= 0.6, 1, 0))   
          
          # Quick rounding to make cleaner (and smaller file sizes) output tables
          tmp_4 <- tmp_4 %>% 
            dplyr::mutate(step = round(step),
                          angle = round(angle, 3),
                          direction = round(direction, 3))
          
          # Only retain the necessary information for saving location data frame to file
          tmp_save <- tmp_4 %>% 
            dplyr::select(path_id, i_strata, x, y, sex_01) 
          
          # Store the selected steps for this iteration in a list
          path_list[[i]] <- tmp_save
          
          # Then set these x,y locations to x_prev and y_prev for use in the next iteration
          # Do the same for the previous direction
          tmp_5 <- tmp_4 %>% 
            dplyr::mutate(x_prev = x,
                          y_prev = y,
                          x = NA,
                          y = NA,
                          direction_prev = direction) %>% 
            dplyr::select(dplyr::all_of(names(tmp_0)))
          
          # use tmp_0 in next iteration at top of loop
          tmp_0 <- tmp_5  
          
        } 
        
      }  # END of if n_keep > 0
      
    }  # End of i (individual steps)
    
    # Take the list with all saved steps and add the i_mega sequence
    dplyr::bind_rows(path_list) %>% 
      dplyr::mutate(i_mega = yyy, .before = path_id)
  }

# END OF MEGA LOOP
tictoc::toc()

# Stop the cluster
parallel::stopCluster(cl)

# Delete temporary rasters used in simulation
unlink("spatial/raster/temp_sim/*")


###### SAVE OUTPUTS #########

# Save the path_mega data frame, which has all the actual simulated locations
saveRDS(path_mega, stringr::str_c("outputs/paths_sim/path_all_", m_season, ".rds"))

# Rasterize the paths so pixel values are the number of simulated locations that
# fell within the pixel across all simulated animals.

# All animals
terra::rasterize(cbind(path_mega$x, path_mega$y), terra::rast("spatial/raster/r_study.tif"), fun = length) %>% 
  terra::classify(., cbind(NA, NA, 0)) %>% 
  terra::crop(., sf::read_sf("spatial/shp/extent_canada.shp"), mask = T) %>% 
  terra::writeRaster(., stringr::str_c("outputs/maps/raw/ud", m_season, "all_current.tif", sep = "_"), 
                     overwrite = T, datatype = "INT2U")

# Females only
terra::rasterize(cbind(path_mega$x[path_mega$sex_01 == 0], path_mega$y[path_mega$sex_01 == 0]), 
                 terra::rast("spatial/raster/r_study.tif"), fun = length) %>% 
  terra::classify(., cbind(NA, NA, 0)) %>% 
  terra::crop(., sf::read_sf("spatial/shp/extent_canada.shp"), mask = T) %>% 
  terra::writeRaster(., stringr::str_c("outputs/maps/raw/ud", m_season, "female_current.tif", sep = "_"), 
                     overwrite = T, datatype = "INT2U")

# Males only
terra::rasterize(cbind(path_mega$x[path_mega$sex_01 == 1], path_mega$y[path_mega$sex_01 == 1]), 
                 terra::rast("spatial/raster/r_study.tif"), fun = length) %>% 
  terra::classify(., cbind(NA, NA, 0)) %>% 
  terra::crop(., sf::read_sf("spatial/shp/extent_canada.shp"), mask = T) %>% 
  terra::writeRaster(., stringr::str_c("outputs/maps/raw/ud", m_season, "male_current.tif", sep = "_"), 
                     overwrite = T, datatype = "INT2U")
