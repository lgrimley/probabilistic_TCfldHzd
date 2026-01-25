# Synthetic TC Tracks 
Overview XYZ
MIT model (cite)

## Data Requirements
1. **Track files** (MAT-files) containing TC information:
   - Variables: `lat100`, `lon100`, `vstore100`, `pstore100`, `rmw100`, `ro100`, `uinc100`, `vinc100`
   - Location: `.\tracks\UScoast6_AL_<GCM>cal_roEst1rmEst1_trk100`
2. Points (i.e., coastline_lores.csv)
3. Polygon (i.e., model domain)
---
## TC Track Information
These scripts return the TC storm characteristics near landfall. First, the code clips all the intersecting tracks to a specified area (buffered of the SFINCS model domain). Of the portions of the track remaining, the code calculates the distance between the center of the TC tracks (updated every 2-hours) to the modeled ADCIRC storm surge locations which are offshore (<1km off the coastline). The gage that has the minimum distance to the TC center is selected as the nearest landfall gage and the TC characteristics at this timestep are returned (max wind speed, radius to max wind, minimum pressure, etc.). This information is used in downstream analysis.
### Scripts
- **tracks_select_by_point.py** and **tracks_select_by_polygon.py**
  - Usage: this scripts read in the TC track files (.mat) and selects the TC IDs that intersect a point/polygon given a buffer
  - Output: a list of TC IDs 
- **calculate_landfall_vmax.py**
  - Usage: this script returns the max wind speeds at landfall for the selected TC IDs
  - Output: a table of the landfall Vmax by TC ID
- **tc_track_info_to_table.py**
  - Usage: this script returns all the track variables (vmax, pstore, etc.) for a select TC ID
  - Output:
- **storm_durations_table.py**
  - Usage: This script estimates the duration of select TC tracks (by TC ID) over the model domain. The output is used to estimate simulation time.
  - Output:
---
## Bias Correcting Intensity Distributions
This folder has scripts to bias correct the TC intensity at landfall. GCM TC outputs need to be bias-corrected because they over/under predict the TC intensity (vamx) distribution when compared to historical storms. This is done using delta-quantile mapping which essentially provides "TC weights" that are used to shift the CDF when calculating probabilities. The scripts to do this are in this folder.


