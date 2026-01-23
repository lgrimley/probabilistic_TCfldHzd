# Synthetic TC Track Datasets

Below is a list of the scripts in this folder and there general usages.  

1. **Getting landfall location and info for synthetic TC tracks**
   - These scripts clip the storm tracks (polyline) to a polygon (buffer of the SFINCS model domain) and then calculate the distance from the center of the TC tracks (2 hours) to the ADCIRC storm surge point locations. The TC center that is closest to one of these offshore gage locations (<1km off the coastline) is considered the "landfall" timestep. The TC track info (minimum pressure, max wind speed, radius to max wind) are assigned as the landfall characteristics. This information is interesting to keep track of since it is often used to categorize storm intensity. We also use this info to bias correct the GCM model.
   - Scripts:
     - calculate_landfall_vmax.py --> getting the max wind speeds at landfall
     - tracks_select_by_point.py and tracks_select_by_polygon.py --> spatial selection of TC tracks
     - tc_track_info_to_table.py --> get all TC characteristics at landfall


2. **Estimate the TC durations using the storm tracks**
   - Script: storm_durations_table.py
   - Overview: we wanted to know how long the TC tracks were over the model domain to estimate approximately how long the SFINCS model runs would be to estimate simulation time. 

3. **Global climate model Bias Correction**
   - Overview: GCM TC outputs need to be bias-corrected because they over/under predict the TC intensity (vamx) distribution when compared to historical storms. This is done using delta-quantile mapping which essentially provides "TC weights" that are used to shift the CDF when calculating probabilities. The scripts to do this are in this folder.