
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
# Script description: Run mixed-effect iSSF model in glmmTMB for spring season
#
# Script author: Eric Palm
#
################################################################################

# Install 'glmmTMB' using the following steps:
# 1) install 'TMB' and 'remotes' using install.packages()
# 2) install 'glmmTMB' with remotes::install_github("glmmTMB/glmmTMB/glmmTMB")

# Load packages
sapply(
  c("glmmTMB", "dplyr"), 
  require, character.only = T)

# Load input data frame (location coordinates omitted)
df_issf <- readRDS("data/df_issf.rds") 

# define season
m_season <- "spring"

# Center and scale all continuous covariates separately by season, 
# but don't center/scale movement parameters (log_sl_km, sl_km, cos_ta).
df_issf_scaled <- 
  df_issf %>% 
  dplyr::group_by(season) %>% 
  dplyr::mutate(dplyr::across(c(canopy_cover_end:d_road_end, 
                                d_hwy_decay_end:greenness_broad_end, d_spring_start_end),
                              ~as.numeric(scale(.)[,1]))) %>% 
  dplyr::ungroup()

# Number of "cores" for running models in glmmTMB
nt <- 10

# Run the mixed effects model in glmmTMB with random slopes on (almost) every variable.
# We didn't include random slopes on GFR interactions because there's very little 
# variation within an animal's home range (GFRs use home range-level averages).
system.time(
  
  m <- glmmTMB::glmmTMB(
    
    # drop the fixed-effect intercept
    case_ ~ -1 +
      
      # habitat covariates
      greenness_end + canopy_cover_end + I(canopy_cover_end^2) +
      tri_end  + I(tri_end^2) + northness_end + elev_end + I(elev_end^2) + 
      d_hwy_decay_end + d_road_decay_end +
      mine_end +
      
      # habitat interactions 
      canopy_cover_end:elev_end + 
      greenness_end:elev_end +
      greenness_end:canopy_cover_end +
      
      # distance to spring starting location
      d_spring_start_end +
      
      # crossing variables
      crossing_lake_reservoir + crossing_hwy + crossing_rock_ice + crossing_town +
      
      # GFR interactions
      crossing_rock_ice:greenness_broad_end +
      crossing_hwy:greenness_broad_end + 
      d_hwy_decay_end:greenness_broad_end + d_road_decay_end:greenness_broad_end + 
      
      # movement parameters
      cos_ta + log_sl_km + 
      
      # movement interactions
      # use the 'start' habitat variables for these
      # The cos_ta:log_sl_km tests for straighter movements when there are longer steps
      (canopy_cover_start + tri_start + greenness_start + elev_start):log_sl_km +
      cos_ta:log_sl_km + 
      
      # Repeat all terms above but include interactions with binary sex variable
      # habitat covariates and sex
      greenness_end:sex_01 + canopy_cover_end:sex_01 + I(canopy_cover_end^2):sex_01 +
      tri_end:sex_01  + I(tri_end^2):sex_01 + northness_end:sex_01 + elev_end:sex_01 + I(elev_end^2):sex_01 + 
      d_hwy_decay_end:sex_01 + d_road_decay_end:sex_01 +
      mine_end:sex_01 +
      
      # habitat interactions and sex
      canopy_cover_end:elev_end:sex_01 + 
      greenness_end:elev_end:sex_01 +
      greenness_end:canopy_cover_end:sex_01 +
      
      # distance to spring start and sex
      d_spring_start_end:sex_01:sex_01 +
      
      # crossing variables and sex
      crossing_lake_reservoir:sex_01 + crossing_hwy:sex_01 + crossing_rock_ice:sex_01 + crossing_town:sex_01 +
      
      # GFR interactions
      crossing_rock_ice:greenness_broad_end:sex_01 +
      crossing_hwy:greenness_broad_end:sex_01 + 
      d_hwy_decay_end:greenness_broad_end:sex_01 + d_road_decay_end:greenness_broad_end:sex_01 +  
      
      # movement parameters and sex
      cos_ta:sex_01 + log_sl_km:sex_01 + 
      
      # movement interactions and sex
      (canopy_cover_start:sex_01 + tri_start:sex_01 + greenness_start:sex_01 + elev_start:sex_01):log_sl_km:sex_01 +
      cos_ta:log_sl_km:sex_01 + 
      
      # random effects
      # random intercept (always put this first of random effects)
      (1 | stratum) + 
      
      # habitat random slopes at the animal level
      (0 + canopy_cover_end | id) + (0 + I(canopy_cover_end^2) | id) + 
      (0 + greenness_end | id) +
      (0 + d_hwy_decay_end | id) + (0 + northness_end | id) +
      (0 + d_road_decay_end | id) +
      (0 + elev_end | id) + (0 + I(elev_end^2) | id) +
      (0 + tri_end | id) + (0 + I(tri_end^2)| id) +
      (0 + mine_end | id) + 
      
      # habitat interaction random slopes
      (0 + greenness_end:canopy_cover_end |id) +
      (0 + greenness_end:elev_end | id) +
      (0 + canopy_cover_end:elev_end | id)  +
      
      # distance to spring start random slope
      (0 + d_spring_start_end | id)  +
      
      # movement parameter random slopes
      (0 + cos_ta | id) + 
      (0 + log_sl_km | id) +
      (0 + cos_ta:log_sl_km | id) + 
      
      # crossing random slopes (couldn't fit these on towns or lakes
      # because there were too few used steps that crossed these rare features)
      (0 + crossing_hwy | id) +
      (0 + crossing_rock_ice | id) +
      
      # movement interaction random slopes
      (0 + log_sl_km:tri_start | id) + 
      (0 + log_sl_km:canopy_cover_start | id) +
      (0 + log_sl_km:greenness_start | id) + 
      (0 + log_sl_km:elev_start | id)
    ,
    
    # The NA in the 'map' argument tells model not to estimate the variance of the random intercept
    # (the first random effect), but to do so for the remaining 24 random slopes.
    # The NA in the 'start' argument tells it to fix the variance of that random intercept 
    # to a high value (log(1e3)), while all 24 of the random slopes get a 0. This is from Muff et al. 2020. 
    # Adjust the 24s in both arguments to whatever number of random slopes you have in the model.
    data = filter(df_issf_scaled, season == m_season), family = poisson,
    map = list(theta = factor(c(NA, 1:24))), start = list(theta = c(log(1e3), rep(0, 24))),
    control = glmmTMBControl(parallel = nt))
)

# Save the model - very large file size
# saveRDS(m, stringr::str_c("outputs/models/m_", m_season, ".rds"))

# Save fixed effects coefficients
glmmTMB::fixef(m)$cond %>% 
  saveRDS(stringr::str_c("outputs/models/fixef_", m_season, ".rds"))

# Save full variance-covariance matrix
vcov(m, full = T) %>% 
  saveRDS(stringr::str_c("outputs/models/vcov_", m_season, ".rds"))

