# Synthetic TC Tracks 

## Data Requirements
1. TC Dataset from Gori et al. 2022 (https://www.nature.com/articles/s41558-021-01272-7)
   - **Track files** (MAT-files) containing the output from the MIT TC Model: TC tracks with variables `lat100`, `lon100`, `vstore100`, `pstore100`, `rmw100`, `ro100`, `uinc100`, `vinc100`
   - **ADCIRC Stormtide Locations***
2. SFINCS model domain or study area (polygon)
   
---
## TC Track Information
These scripts return the TC storm characteristics near landfall. First, the code clips all the intersecting tracks to a specified area (buffered of the SFINCS model domain). Of the portions of the track remaining, the code calculates the distance between the center of the TC tracks (updated every 2-hours) to the modeled ADCIRC storm surge locations which are offshore (<1km off the coastline). The gage that has the minimum distance to the TC center is selected as the nearest landfall gage and the TC characteristics at this timestep are returned (max wind speed, radius to max wind, minimum pressure, etc.). This information is used in downstream analysis.
### Scripts
- **tracks_select_by_point.py** and **tracks_select_by_polygon.py**
  - Usage: selecting the TC tracks that come within a distance of a point or intersect with a polygon
  - Input:
     - point or polygon shapefile
     - track files (.mat)
  - Output:
     - list of TC IDS corresponding to the track files that intersect input geometry (buffered point or polygon)
- **calculate_landfall_vmax.py**
  - Usage: this script returns the max wind speeds at landfall for the select TCs
  - Input:
     - track files (.mat)
     - list of TC IDs (i.e., generated from scripts above)
  - Output:
     - table of landfall vmax for each TC ID
- **tc_track_info_to_table.py**
  - Usage: this script returns all of the variable stored in the TC track files for the select TC IDs at landfall
  - Input:
     - track files (.mat)
     - list of TC IDs (i.e., generated from scripts above)
  - Output:
     - table of all TC track variables for each TC ID
- **storm_durations_table.py**
  - Usage: this script estimates the duration of the TC over the model domain to estimate approximately how long the SFINCS model runs would be to estimate simulation time.
  - Input:
     - track files (.mat)
     - list of TC IDs (i.e., generated from scripts above)
---
## Bias Correcting Intensity Distributions
The `gcm_bias_correction` folder has scripts to bias correct the TC intensity at landfall. GCM TC outputs need to be bias-corrected because they over/under predict the TC intensity (vamx) distribution when compared to historical storms. This is done using delta-quantile mapping which essentially provides "TC weights" that are used to shift the CDF when calculating probabilities. The scripts to do this are in this folder.


