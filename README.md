# Habitat selection and functional connectivity analysis for grizzly bears in the Rocky Mountains in southeast British Columbia and southwest Alberta, Canada

This repository includes data and code for reproducing analyses in the following article, accepted for publication in *Conservation Science and Practice*:

**"Changing grizzly bear space use and functional connectivity in response to human disturbance in the southern Canadian Rocky Mountains"**
by Eric Palm, Clayton Apps, Tal Avgar, Melanie Dickie, Bruce McLellan, Joseph Northrup, Michael Sawaya, Julie Turner, Jesse Whittington, Erin Landguth, Katherine Zeller, Clayton Lamb.


Please open the `grizzly_issf.Rproj` file to start RStudio before opening individual code files to ensure that relative file paths in the code work correctly.


Here is a list of script files with descriptions: 

`1a. fit_issf_spring.R` – Fit the final integrated step-selection function (iSSF) for the spring season. 

`1b. fit_issf_summer.R` – Fit the final iSSF for the summer season. 

`1c. fit_issf_fall.R` – Fit the final iSSF for the fall season. 

`2. simulate_from_issf.R` – Simulate a utilization distribution from the fitted spring season iSSF model. 

`3. save_final_UDs.R` – Save final utilization distributions by season, sex and disturbance scenario.


**FINAL UTILIZATION DISTRIBUTION RASTERS ARE IN** `outputs/maps/final`