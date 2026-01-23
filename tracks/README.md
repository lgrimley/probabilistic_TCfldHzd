# Synthetic TC Tracks 
Overview XYZ
MIT model (cite)

## Data Requirements
1. **Track files** (MAT-files) containing TC information:
   - Variables: `lat100`, `lon100`, `vstore100`, `pstore100`, `rmw100`, `ro100`, `uinc100`, `vinc100`
   - Location: `.\tracks\UScoast6_AL_<GCM>cal_roEst1rmEst1_trk100`
2. coastline_lores.csv
3. Model domain??
---
## TC Track Information
These scripts return the TC storm characteristics near landfall. First, the code clips all the intersecting tracks to a specified area (buffered of the SFINCS model domain). Of the portions of the track remaining, the code calculates the distance between the center of the TC tracks (updated every 2-hours) to the modeled ADCIRC storm surge locations which are offshore (<1km off the coastline). The gage that has the minimum distance to the TC center is selected as the nearest landfall gage and the TC characteristics at this timestep are returned (max wind speed, radius to max wind, minimum pressure, etc.). This information is used in downstream analysis.
### Scripts
- **tracks_select_by_point.py**
  - Usage:
  - Output:
  - Configuration:
- **tracks_select_by_polygon.py**
  - Usage:
  - Output:
  - Configuration:
- **calculate_landfall_vmax.py**
  - Usage: this script returns the max wind speeds at landfall for the select TCs
  - Output:
  - Configuration:
- **tc_track_info_to_table.py**
  - Usage: this script returns all TC info from the tracks dataset at landfall for the select TCs
  - Output:
  - Configuration:
- **storm_durations_table.py**
  - Usage: this script estimates the duration of the TC over the model domain to estimate approximately how long the SFINCS model runs would be to estimate simulation time.
  - Output:
  - Configuration:
---
## Bias Correcting Intensity Distributions
This folder has scripts to bias correct the TC intensity at landfall. GCM TC outputs need to be bias-corrected because they over/under predict the TC intensity (vamx) distribution when compared to historical storms. This is done using delta-quantile mapping which essentially provides "TC weights" that are used to shift the CDF when calculating probabilities. The scripts to do this are in this folder.


