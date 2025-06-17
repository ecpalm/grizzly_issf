
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
# Script description: Create final seasonal and annual utilization distributions
#
# Script author: Eric Palm
#
################################################################################

# Load packages
sapply(
  c("terra", "dplyr", "lubridate", "purrr", "stringr"), 
  require, character.only = T)

# Functions to find 99% quantile of raster pixel values and to normalize raster values between 0 and 1
q99 <- function(x) as.numeric(global(x, fun = quantile, probs = .99, na.rm = T))
norm <- function(x) (x - min(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm = T))

# Get season lengths in days from issf data frame by calculating median start and
# end dates by animal, then median start and end dates across animals
lengths_season <-
  readRDS("data/df_issf.rds") %>%
  dplyr::mutate(year = lubridate::year(t2_)) %>% 
  dplyr::group_by(season, year, id) %>%
  dplyr::summarize(start = min(lubridate::yday(t2_)),
                   end = max(lubridate::yday(t2_))) %>%
  dplyr::group_by(season, year) %>%
  dplyr::summarize(med_start_year = median(start),
                   med_end_year = median(end)) %>%
  dplyr::summarize(med_start_season = median(med_start_year),
                   med_end_season = median(med_end_year)) %>%
  dplyr::mutate(length_season = med_end_season - med_start_season + 1) %>% 
  dplyr::dplyr::select(season, length_season)

# Import all raw simulated UD rasters and set extreme pixel values 
# (those above 99th quantile) to the 99th quantile value
uds_seasonal_q99 <-
  dplyr::tibble(sex = rep(c("all", "female", "male"), each = 3, times = 3),
                season = rep(c("fall", "spring", "summer"), each = 9),
                period = rep(c("current", "future", "past"), times = 9),
                name = stringr::str_c("ud", season, sex, period, "q99", sep = "_")) %>% 
  dplyr::mutate(ud_q99 = list.files("outputs/maps/raw", full.names = T, pattern = ".tif") %>%
                  purrr::map(., ~terra::rast(.) %>%
                               terra::classify(., cbind(q99(.), Inf, q99(.))))) 

# Write seasonal rasters to file for males, females, and all bears
uds_seasonal_q99 %>% 
  dplyr::pull(ud_q99) %>% 
  terra::rast() %>% 
  `names<-`(uds_seasonal_q99$name) %>% 
  terra::writeRaster(., stringr::str_c("outputs/maps/final/seasonal/", names(.), ".tif"), datatype = "INT2U")

# Normalize each seasonal raster and then create annual rasters for each of the
# three disturbance scenarios and by sex (male, female, all):
# 1) past (disturbance-free), 2) current, 3) future (additional disturbance) 
# Each season weighted by its corresponding duration in days
uds_annual_q99_norm <- 
  uds_seasonal_q99 %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(ud_q99_norm = list(terra::app(ud_q99, norm))) %>% 
  dplyr::inner_join(., lengths_season) %>% 
  dplyr::group_by(sex, period) %>%
  dplyr::summarize(ud_annual = list(terra::weighted.mean(terra::rast(ud_q99_norm), length_season) %>% 
                                   terra::app(., norm))) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(name = stringr::str_c("ud_annual", sex, period, "q99", sep = "_"))

# Write annual rasters to file for males, females, and all bears
uds_annual_q99_norm %>% 
  dplyr::pull(ud_annual) %>% 
  terra::rast() %>% 
  `names<-`(uds_annual_q99_norm$name) %>% 
  terra::writeRaster(., stringr::str_c("outputs/maps/final/annual/", names(.), ".tif"), datatype = "FLT4S")

