

1. **calc_storm_weights.m**
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
    
2. **weighted_selection_and_cdf.py**
   - Overview:
       - This code is used to compare the CDF of a subset of storms using weighted selection from the full suite.
