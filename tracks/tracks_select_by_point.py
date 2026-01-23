import os

import scipy.io as sio
import geopandas as gpd
import pandas as pd
import numpy as np
from src.utils import track_points_to_linestring


"""
This script identifies tropical cyclone (TC) tracks that intersect a buffered
ADCIRC gage location and compares those storms against TC IDs available in
ADCIRC storm tide simulations.

The workflow:
1. Load ADCIRC gage locations and apply a spatial buffer.
2. Convert TC track points into LineString geometries.
3. Identify TC tracks that intersect the buffered gage location.
4. Compare selected storms with those included in ADCIRC simulations.
5. Compare spatially selected storms against storms within a 200 km buffer.
"""

# -------------------------------------------------------------------------
# Load ADCIRC gage locations (points)
# -------------------------------------------------------------------------
# CSV is assumed to contain at least:
#   - 'x': longitude
#   - 'y': latitude
#   - 'gage_id': unique identifier for each ADCIRC gage
points = pd.read_csv(
    r'Z:\Data-Expansion\users\lelise\Chapter3\NCEP_Reanalysis\stormTide\coastline_lores_NCSC.csv'
)

# Convert point data to a GeoDataFrame using WGS84 coordinates (EPSG:4326)
points_gdf = gpd.GeoDataFrame(
    points,
    geometry=gpd.points_from_xy(
        x=points['x'].values,
        y=points['y'].values,
        crs=4326
    )
)

# Apply a spatial buffer around each gage location
# NOTE: Distance units depend on CRS (degrees for EPSG:4326)
points_gdf_buffer = gpd.GeoDataFrame(
    points_gdf,
    geometry=points_gdf.geometry.buffer(distance=2)
)

# -------------------------------------------------------------------------
# Load tropical cyclone track data
# -------------------------------------------------------------------------
# TC track data are stored in MATLAB format and contain latitude/longitude
# points for each storm
fname = (
    r'Z:\Data-Expansion\users\lelise\Chapter3\NCEP_Reanalysis\tracks'
    r'\UScoast6_AL_ncep_reanal_roEst1rmEst1_trk100'
)
tc_tracks = sio.loadmat(f'{fname}.mat')

# Convert TC track points to LineString geometries
# Each LineString represents the full track of a single TC
tc_tracks_polyline = track_points_to_linestring(
    tc_tracks,
    output_shpfile=None
)

# -------------------------------------------------------------------------
# Select TC tracks that intersect a specific ADCIRC gage
# -------------------------------------------------------------------------
gage = 194

# Select buffered geometry for the target gage
pt = points_gdf_buffer[points_gdf_buffer['gage_id'] == gage]

# Identify TC tracks that intersect the buffered gage region
select_trks = tc_tracks_polyline.overlay(pt, how='intersection')

# Extract unique TC IDs for intersecting tracks
trk_ids = select_trks['tc_id'].unique()

# -------------------------------------------------------------------------
# Compare with TC IDs used in ADCIRC storm tide simulations
# -------------------------------------------------------------------------
# Load table containing TC IDs used in ADCIRC simulations by gage
adc = pd.read_csv(
    r'Z:\Data-Expansion\users\lelise\Chapter3\NCEP_Reanalysis\stormTide\WL_inds_subset.csv',
    index_col=0
)

# Extract TC IDs associated with the selected gage
adc_pt = adc[f'{gage}'].unique()

# Remove NaN values (no storm assigned)
adc_pt = adc_pt[~np.isnan(adc_pt)]

# Identify TCs that both intersect the gage buffer and were simulated in ADCIRC
common_tcs = select_trks[select_trks['tc_id'].isin(adc_pt)]

# -------------------------------------------------------------------------
# Compare intersecting storms with storms inside a 200 km buffer
# -------------------------------------------------------------------------
os.chdir(r'Z:\Data-Expansion\users\lelise\Chapter3\NCEP_Reanalysis\tracks')

# TC IDs that were modeled in ADCIRC
adcirc_tcs = pd.read_csv('adcirc_modeled_TCs_all.csv')

# TC IDs that pass within a 200 km buffer of the study region
buffer_200km_tcs = pd.read_csv(
    'UScoast6_AL_ncep_reanal_roEst1rmEst1_trk100_200km.csv'
)

# TCs that are both modeled in ADCIRC and within the 200 km buffer
common_tcs = adcirc_tcs[
    adcirc_tcs['tc_id'].isin(buffer_200km_tcs['tc_id'])
]

# TCs within the 200 km buffer that were not modeled in ADCIRC
no_surge_tcs = buffer_200km_tcs[
    ~buffer_200km_tcs['tc_id'].isin(adcirc_tcs['tc_id'])
]

# Summary output
print(f'{len(common_tcs)} of the ADCIRC TCs are also within the 200km buffer.')
