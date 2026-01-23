## Generate Synthetic TC Wind Fields
### Script 
**getmaxwind_county_SynStorms_grid_LG.m**


### Overview
The script uses the **CLE15 wind model** (Chavas et al., 2015) to generate surface wind profiles for each storm along its track.
This MATLAB script calculates the **maximum wind speed** and **U/V wind components** for synthetic tropical cyclones (TCs) tracks.
The outputs are stored as **NetCDF files**, ready for input into SFINCS.

https://journals.ametsoc.org/view/journals/atsc/72/9/jas-d-15-0014.1.xml

**External Functions** (must be on MATLAB path)
- `ER11E04_nondim_rmaxinput` (CLE15 wind model)
- `BoundingGrid`
- `m_lldist_L` (Mapping Toolbox)

---

### Input Data

1. **Track files** (MAT-files) containing TC information:
   - Variables: `lat100`, `lon100`, `vstore100`, `pstore100`, `rmw100`, `ro100`, `uinc100`, `vinc100`
   - Location: `.\tracks\UScoast6_AL_<GCM>cal_roEst1rmEst1_trk100`

2. **Selected TC IDs** (CSV):
   - Columns: `tc_id`, gage counts
   - Location: `.\stormTide\stormTide_TCIDs_and_gageCounts_<GCM>cal.csv`

---

### Outputs

- Individual NetCDF of the wind per synthetic TC
- Master array of maximum winds (optional for analysis): `maxwindgrid_all`

Each NetCDF file contains:

| Variable       | Description                  | Units |
|----------------|------------------------------|-------|
| `x`            | Longitude                    | deg   |
| `y`            | Latitude                     | deg   |
| `time`         | Timestep index               | —     |
| `wind_speed`  | 10-m wind speed              | m/s   |
| `wind10_u`    | Zonal wind component         | m/s   |
| `wind10_v`    | Meridional wind component    | m/s   |

---

### Configuration
- **Study region**: Latitude: 31.95 – 36.81° N; Longitude: -83.67 – -75.11° E
- **Grid spacing**: `inc = 0.017°` (~1-2 km)
- **GCMs processed**: `canesm_ssp585, cnrm6_ssp585, ecearth6_ssp585, ipsl6_ssp585, miroc6_ssp585`

---

### Usage

1. Set working directory to the location of CMIP6 data
2. Run the script in MATLAB.
3. Outputs will be saved in NetCDF format for each TC in the corresponding GCM folder.

---

### Notes
- Background winds (uinc, vinc) are added to the cyclonic wind.
- Optional plots for intermediate steps can be enabled by uncommenting sections.
- Wind speeds are converted from knots to m/s and adjusted from surface to gradient winds. 
- Time is stored as an index (not physical units).