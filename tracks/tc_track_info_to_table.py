import os
import pandas as pd
import scipy.io as sio
from src.utils import get_track_info_in_df, get_landfall_info, get_track_datetime
from hydromt_sfincs import SfincsModel

"""
This script estimates approximate landfall times and wind speeds for a set
of tropical cyclones (TCs) from CMIP6 synthetic tracks using a buffered SFINCS
model domain and ADCIRC gage locations.

Landfall is operationally defined as the intersection of the TC track with a
buffered domain around the model area in proximity to ADCIRC gages. 
"""

# -------------------------------------------------------------------------
# Set working directory
# -------------------------------------------------------------------------
os.chdir(r'Z:\Data-Expansion\users\lelise\projects\Carolinas_SFINCS')

# -------------------------------------------------------------------------
# Load SFINCS model domain for landfall clipping
# -------------------------------------------------------------------------
mod = SfincsModel(root=r'.\03_MODEL_RUNS\sfincs_base_mod', mode='r')

# Convert the model domain to geographic coordinates and buffer it by 5 degrees
# FLAG: Buffer is large (~500 km) and may overestimate the area for landfall detection
region = mod.region.to_crs(4326).buffer(5)

# -------------------------------------------------------------------------
# Load ADCIRC station locations used as reference points for landfall
# -------------------------------------------------------------------------
gage_locs = pd.read_csv(
    r'.\02_DATA\NCEP_Reanalysis\stormTide\coastline_lores_NCSC.csv'
)

# -------------------------------------------------------------------------
# Load TC track data
# -------------------------------------------------------------------------
# Uncomment the following line to use NCEP reanalysis tracks
# matfile = r'.\02_DATA\NCEP_Reanalysis\tracks\UScoast6_AL_ncep_reanal_roEst1rmEst1_trk100.mat'

# Current run uses CMIP6 synthetic tracks
matfile = r'.\02_DATA\CMIP6_585\tracks\UScoast6_AL_canesm_ssp585cal_roEst1rmEst1_trk100.mat'

tc_tracks = sio.loadmat(f'{matfile}')

# -------------------------------------------------------------------------
# Load the list of TC IDs that were modeled in ADCIRC
# -------------------------------------------------------------------------
# Uncomment for NCEP TCs
# modeled_tcs = pd.read_csv(r'.\02_DATA\NCEP_Reanalysis\stormTide\gage_peaks_ZerosRemoved_ncep.csv', index_col=0)

# Current run uses CMIP6 synthetic TCs
modeled_tcs = pd.read_csv(
    r'.\02_DATA\CMIP6_585\stormTide\gage_peaks_ZerosRemoved_canesm_ssp585.csv',
    index_col=0
)
modeled_tcs['tc_id'] = modeled_tcs.index.tolist()
tc_ids = modeled_tcs['tc_id'].tolist()

# -------------------------------------------------------------------------
# Initialize empty dataframe to store landfall info
# -------------------------------------------------------------------------
landfall_info_df = pd.DataFrame()
issue_tcs = []

# -------------------------------------------------------------------------
# Loop through each TC to calculate approximate landfall time and wind speed
# -------------------------------------------------------------------------
for tc_id in tc_ids:
    # Load track data for the TC
    tc_df = get_track_info_in_df(tc_id, tc_tracks)
    tc_df = get_track_datetime(tc_df)
    
    # Remove points with zero wind speed
    tc_df_clean = tc_df[tc_df['vt100'] != 0]

    try:
        # Estimate landfall using buffered domain and gage locations
        landfall_df = get_landfall_info(
            tc_df=tc_df_clean,
            gage_locs=gage_locs,
            clip_gdf=region
        )
        landfall_df['tc_id'] = tc_id

        # Append to master dataframe
        landfall_info_df = pd.concat(
            objs=[landfall_info_df, landfall_df],
            axis=0,
            ignore_index=True
        )

    except Exception:
        # Track failed landfall detection
        issue_tcs.append(tc_id)
        print(f'Issue with {tc_id}')
        break  # FLAG: Stops the loop on first error; may not process remaining TCs

# -------------------------------------------------------------------------
# Save landfall info to CSV
# -------------------------------------------------------------------------
landfall_info_df.to_csv(
    r'Z:\Data-Expansion\users\lelise\projects\Carolinas_SFINCS\Chapter3_SyntheticTCs'
    r'\02_DATA\CMIP6_585\tracks\canesm_landfall_track_info.csv'
)

# -------------------------------------------------------------------------
# FLAGS / NOTES
# -------------------------------------------------------------------------
# 1. The region buffer of 5 degrees is very large and may overestimate landfall.
# 2. Landfall is detected based on proximity to gages and the buffered domain,
#    not exact coastline intersection.
# 3. The 'break' in the loop stops processing on the first TC that causes an error.
#    For production runs, this should be removed to process all TCs.
# 4. Zero wind points are excluded before landfall detection.
# 5. Results are saved to CSV for downstream analysis.
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
# METHODS PARAGRAPH (for documentation/publication)
# -------------------------------------------------------------------------
"""
Landfall Definition and Methods:

Landfall is operationally defined as the first instance when a tropical cyclone 
track intersects a buffered SFINCS model domain and comes within proximity to 
ADCIRC gage locations. 

Steps:
1. TC tracks are filtered to remove zero wind-speed points.
2. Tracks are clipped to a 5-degree buffered SFINCS domain to focus on coastal impacts.
3. ADCIRC gages are used as reference points to estimate landfall location.
4. The time of landfall and maximum wind speed at landfall (vstore100 or vt100) 
   are extracted for each TC.
5. Results are saved to CSV for further analysis such as wind-speed statistics, 
   storm surge modeling, or bias correction.
"""
