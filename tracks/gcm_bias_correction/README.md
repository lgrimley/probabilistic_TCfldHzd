**1. Bias Correcting Intensity Distributions**
GCM TC outputs need to be bias-corrected because they over/under predict the TC intensity (vamx) distribution when compared to historical storms. This is done using delta-quantile mapping which essentially provides "TC weights" that are used to shift the CDF when calculating probabilities.

Script: `calc_storm_weights.m`
- Overview: 
  - This code was originally produced in Gori et al. (2022) and was modified for use in this study. 
- Dependencies:
  - `fit_gpd.m` (function for selecting threshold and fitting gpd distribution)
  - `QDM.m` (functino for quantile  delta matching)
- Inputs:
  - Landfall wind speeds (vmax) for NCEP Reanalysis  (observed historical)
  - Landfall wind speeds (vmax) for GCM historical simulations (modeled historical)
  - Landfall wind speeds (vmax) for GCM projected simulations (modeled projected future)
- Outputs:
  - A table of normalized weights for each GCM TC that can be used when calculating probabilities.

**2. Weighted sampling from Intensity Distributions**
This code can be used to sample a specificed number of storms from the bias-corrected intensity distribution. The goal is to sample them in such a way that the difference between the cdf of the full suite (>1k) is similar to the cdf of the sample set (100-500). 
Script: `weighted_selection_and_cdf.py`
