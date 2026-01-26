"""
Script Name: build_sfincs_observation_locations.py

Description:
    This script generates a combined SFINCS observation location file (.obs)
    by merging coastal gauge locations with river gauge locations. The resulting
    file is formatted according to SFINCS observation input requirements.

    The script:
    - Loads an existing SFINCS observation gauge file
    - Loads river gauge locations from a CSV file
    - Converts river gauges to a GeoDataFrame and extracts coordinates
    - Standardizes station naming conventions
    - Combines coastal and river gauges into a single table
    - Writes a formatted .obs file for use in SFINCS simulations

Notes:
    - Coordinate reference systems are assumed to be compatible
    - Output formatting (spacing, quoting) is critical for SFINCS
    - Input column names and file paths must remain consistent
"""

import os
import datetime
import hydromt
import pandas as pd
import rasterio.merge
from hydromt import DataCatalog
from hydromt_sfincs import SfincsModel, utils
import geopandas as gpd

# Path to HydroMT data catalog
yml = os.path.join(r'Z:\Data-Expansion\users\lelise\data', 'data_catalog_SFINCS_Carolinas.yml')
cat = hydromt.DataCatalog(yml)

# Change working directory to SFINCS observation file location
os.chdir(
    r'Z:\Data-Expansion\users\lelise\projects\Carolinas_SFINCS\Chapter3_SyntheticTCs'
    r'\SFINCS_mod_setup\base_model\obsfile'
)

# ---------------------------------------------------------------------
# Load existing SFINCS observation gauge locations
# ---------------------------------------------------------------------
obs_gage = pd.read_csv('sfincs_gages.obs', delim_whitespace=True, header=None)
obs_gage.columns = ['x', 'y', 'name']

# ---------------------------------------------------------------------
# Load river gauge locations and convert to GeoDataFrame
# ---------------------------------------------------------------------
riv_gage = pd.read_csv('river_obs_locs_v3.csv')
riv_gage = gpd.GeoDataFrame(
    riv_gage,
    geometry=gpd.points_from_xy(
        x=riv_gage['xcoord'],
        y=riv_gage['ycoord'],
        crs=32617
    )
)

# Extract projected x/y coordinates
riv_gage['x'] = riv_gage.geometry.x
riv_gage['y'] = riv_gage.geometry.y

# Construct standardized station numbers with zero padding
staNum = [str(s).zfill(4) for s in riv_gage['ORIG_SEQ']]

# Build station name from HYDRAID and padded station number
riv_gage['name'] = riv_gage['HYDRAID'] + '_' + staNum

# Retain only columns matching the observation gauge format
riv_gage = riv_gage[obs_gage.columns.tolist()]

# ---------------------------------------------------------------------
# Combine coastal and river gauges
# ---------------------------------------------------------------------
combined = pd.concat(objs=[obs_gage, riv_gage], axis=0, ignore_index=True)

# ---------------------------------------------------------------------
# Write combined observation locations to SFINCS .obs file
# ---------------------------------------------------------------------
with open(r'sfincs_obsLocs_gauge_riv_v3.obs', mode='w') as obs_file:
    for i in range(len(combined)):
        x = combined['x'].loc[i].round(2)
        y = combined['y'].loc[i].round(2)
        id = combined['name'].loc[i]

        # Fixed-width formatting required by SFINCS
        line_entry = f'{x:<011}    {y:<011}    "{id}"\n'
        obs_file.write(line_entry)

obs_file.close()
