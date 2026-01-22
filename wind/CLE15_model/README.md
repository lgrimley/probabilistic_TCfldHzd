# Tropical Cyclone Wind Processing Script
## Script 
**getmaxwind_county_SynStorms_grid_LG.m**


## Overview
The script uses the **CLE15 wind model** (Chavas et al., 2015) to generate surface wind profiles for each storm along its track.
This MATLAB script calculates the **maximum wind speed** and **U/V wind components** for synthetic tropical cyclones (TCs) tracks.
The outputs are stored as **NetCDF files**, ready for input into SFINCS.

https://journals.ametsoc.org/view/journals/atsc/72/9/jas-d-15-0014.1.xml
---

## Input Data

1. **Track files** (MAT-files) containing TC information:
   - Variables: `lat100`, `lon100`, `vstore100`, `pstore100`, `rmw100`, `ro100`, `uinc100`, `vinc100`
   - Location: `.\tracks\UScoast6_AL_<GCM>cal_roEst1rmEst1_trk100`

2. **Selected TC IDs** (CSV):
   - Columns: `tc_id`, gage counts
   - Location: `.\stormTide\stormTide_TCIDs_and_gageCounts_<GCM>cal.csv`

---

## Outputs

- **NetCDF files** for each tropical cyclone:
  - Variables:
    - `x` : longitude (degrees)
    - `y` : latitude (degrees)
    - `time` : time steps (arbitrary units)
    - `wind_speed` : total wind speed (m/s)
    - `wind10_u` : U-component (m/s)
    - `wind10_v` : V-component (m/s)
  - Location: `.\wind\CLE15_Gridded_<GCM>cal\`

- **Master array of maximum winds** (optional for analysis):
  - `maxwindgrid_all`

---

## Configuration
- **Study region**: Latitude: 31.95 – 36.81° N; Longitude: -83.67 – -75.11° E
- **Grid spacing**: `inc = 0.017°` (~1-2 km)
- **GCMs processed**: `canesm_ssp585, cnrm6_ssp585, ecearth6_ssp585, ipsl6_ssp585, miroc6_ssp585`

---

## Usage

1. Set working directory to the location of CMIP6 data
2. Run the script in MATLAB.
3. Outputs will be saved in NetCDF format for each TC in the corresponding GCM folder.

---

## Notes
- Ensure BoundingGrid.m and ER11E04_nondim_rmaxinput.m functions are in the MATLAB path.
- Background winds (uinc, vinc) are added to the cyclonic wind.
- Optional plots for intermediate steps can be enabled by uncommenting sections.