
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
# Script description: Function for simulating iSSF coefficients from a glmmTMB 
#                     mixed-effects model and including random slope variance
#
# Script author: Eric Palm
#
# This version of the function does not incorporate random slope 
# variance into the movement covariates (e.g., log_sl_, sl_, cos_ta_)
#
################################################################################

## Load packages
sapply(
  c("stringr", "tidyr", "dplyr"), 
  require, character.only = T)

# Main function
sim_coef_ran <-
  function(m_season, n) {
    
    # Import the saved fixed effects coefficients from the fitted iSSF
    mu_fixed <- readRDS(stringr::str_c("outputs/models/fixef_", m_season, ".rds"))
    
    # Import the saved full variance-covariance matrix from the fitted iSSF
    sig_temp <- readRDS(stringr::str_c("outputs/models/vcov_", m_season, ".rds"))
    
    # Extract the names of all the random effects
    names_rand <- attr(sig_temp, "dimnames")[[1]][-c(1:length(mu_fixed))] %>% 
      as.character()
    
    # Create a vector of mu values where mu is set to 0 for all random effects
    # Random slopes should be offsets from fixed effect coefficients
    mu_all <- c(mu_fixed, rep(0, length(names_rand)) %>% `names<-`(names_rand))
    
    # Find the index of the random intercept, which we will drop
    index_int <- which(names(mu_all) == "theta_1|stratum.1")
    
    # Drop the random intercept from the vector of mus
    mu_final <- mu_all[-index_int]
    
    # Now drop the random intercept from the vcov matrix
    sig <- sig_temp[-index_int, -index_int]
    
    # Simulate n sets of coefficients from a multivariate normal distribution
    coefs_raw <-
      MASS::mvrnorm(n = n, mu = mu_final, Sigma = sig) %>% 
      dplyr::as_tibble()
    
    # Rename the random slopes to match the fixed effects
    names_clean <- 
      coefs_raw %>% 
      dplyr::select(dplyr::contains("theta")) %>% 
      colnames() %>% 
      stringr::str_sub(., 9, -6)  
    
    # Sum the coefficients for fixed effects and random slopes for each variable
    # within each set (simulated animal) of coefficients 
  coefs_raw %>% 
    dplyr::select(!dplyr::contains("theta")) %>% 
    dplyr::mutate(index = dplyr::row_number()) %>% 
    tidyr::pivot_longer(cols = -index, names_to = "covariate") %>%
    dplyr::bind_rows(., coefs_raw %>%
                dplyr::select(dplyr::contains("theta")) %>%
                dplyr::rename_with(~ names_clean) %>% 
                dplyr::mutate(index = dplyr::row_number()) %>% 
                tidyr::pivot_longer(cols = -index, names_to = "covariate") %>% 
                  # Sets the value of random slope coefficients from movement covariates to 0
                  dplyr::mutate(value = dplyr::if_else(stringr::str_detect(covariate, "log_sl|cos_ta"), 0, value))) %>% 
    dplyr::mutate(dplyr::across(dplyr::everything(), ~tidyr::replace_na(., 0))) %>% 
    dplyr::group_by(index, covariate) %>% 
    dplyr::summarize(coef = sum(value)) %>% 
    tidyr::pivot_wider(names_from = "covariate", values_from = "coef") %>% 
    dplyr::ungroup() %>% 
    dplyr::select(-index)
}
