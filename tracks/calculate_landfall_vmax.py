import os
import pandas as pd
import scipy.io as sio

from src.utils import (
    get_track_info_in_df,
    get_landfall_info,
    get_track_datetime
)

from hydromt_sfincs import SfincsModel


"""
This script estimates tropical cyclone (TC) characteristics near landfall using track data.

Inputs: 
- TC track files (.mat)
- ADCIRC storm surge gage locations x,y (.csv)
- SFINCS model domain (geopandas dataframe extracted from SFINCS model OR read in as a shapefile)

Outputs:
TC track landfall max wind speed (vmax) are saved in a table for each synthetic TC (TC_ID) 
and output as a CSV file for downstream analysis (e.g., bias correction). 
Note, landfall is approximated:
 - we clip the tracks to a buffered polygon of the domain (user defined)
 - we calculate the distance from the ADCIRC points to the TC track center (lat,lon)
 - we return the TC track info for the center that is closest to one of the ADCIRC points

Assumptions / Flags:
- The SFINCS domain is buffered by 1.8 degrees for landfall detection, but is not a precise coastline-based measure.
- Some TCs may fail the landfall detection step (logged in `issue_tcs`).

"""

# -------------------------------------------------------------------------
# Load SFINCS model domain and define analysis region
# -------------------------------------------------------------------------
wdir = r'Z:\Data-Expansion\users\lelise\projects\Carolinas_SFINCS\Chapter3_SyntheticTCs'
os.chdir(wdir)

# Load the SFINCS model (to get the region)
mod = SfincsModel(root=r'.\03_MODEL_RUNS\sfincs_base_mod',mode='r')

# Convert model region to geographic coordinates and buffer
buffer_range_deg = 1.8
region = mod.region.to_crs(4326).buffer(buffer_range_deg)

# -------------------------------------------------------------------------
# Load ADCIRC storm surge gage locations (x,y) used for landfall estimation
# -------------------------------------------------------------------------
gage_locs = pd.read_csv(r'.\02_DATA\NCEP_Reanalysis\stormTide\coastline_lores_NCSC.csv')

# -------------------------------------------------------------------------
# Process synthetic TC tracks
# -------------------------------------------------------------------------
# Update the directory depending on where the track files are stored (.mat)
climate_scenario = 'CMIP6_585' # or CMIP6_585
tc_dir = os.path.join(os.getcwd(), '02_DATA',f'{climate_scenario}')

for file in os.listdir(os.path.join(tc_dir,'tracks')):

    if file.endswith('.mat'):
        print(file)
        # Parse GCM and scenario from filename
        gcm = file.split('_')[2]
        ssp = file.split('_')[3].removesuffix('cal')

        # Load TC track data for the different GCM ensembles for the selected climate scenario
        filepath = os.path.join(tc_dir,'tracks', file)
        tc_tracks = sio.loadmat(filepath)

        # Load TC IDs modeled in ADCIRC (./stormTide/stormTide_process_adcircWL.py)
        modeled_tcs = pd.read_csv(os.path.join(os.getcwd(), tc_dir, 'stormTide',rf'gage_peaks_{gcm}_{ssp}.csv'),index_col=0)
        modeled_tcs['tc_id'] = modeled_tcs.index.tolist()
        tc_ids = modeled_tcs['tc_id'].tolist()

        # -----------------------------------------------------------------
        # Estimate landfall characteristics for each TC
        # -----------------------------------------------------------------
        landfall_df = pd.DataFrame()
        counter = 0
        tot = len(tc_ids)

        for tc_id in tc_ids:
            tc_df = get_track_info_in_df(tc_id, tc_tracks)
            tc_df = get_track_datetime(tc_df)
            lf_df = get_landfall_info(
                tc_df=tc_df,
                gage_locs=gage_locs,
                clip_gdf=region
            )
            lf_df['tc_id'] = tc_id
            landfall_df = pd.concat(
                objs=[landfall_df, lf_df],
                axis=0,
                ignore_index=True
            )
            print(f'{counter} out of {tot}')
            counter += 1

        out = landfall_df[['tc_id', 'vstore100']]
        out.to_csv(os.path.join(tc_dir, fr'{gcm}_{ssp}_landfall_vmax.csv'))

# -------------------------------------------------------------------------
# Process NCEP reanalysis tropical cyclones
# -------------------------------------------------------------------------
os.chdir(
    r'Z:\Data-Expansion\users\lelise\projects\Carolinas_SFINCS'
    r'\Chapter3_SyntheticTCs\02_DATA\NCEP_Reanalysis\tracks'
)
# There is just one track file for the reanalysis
matfile = 'UScoast6_AL_ncep_reanal_roEst1rmEst1_trk100.mat'
tc_tracks = sio.loadmat(matfile)

modeled_tcs = pd.read_csv(
    r'../stormTide/gage_peaks_ZerosRemoved_ncep.csv',
    index_col=0
)
modeled_tcs['tc_id'] = modeled_tcs.index.tolist()
tc_ids = modeled_tcs['tc_id'].tolist()

landfall_df = pd.DataFrame()
counter = 0
tot = len(tc_ids)
issue_tcs = []

for tc_id in tc_ids:
    tc_df = get_track_info_in_df(tc_id, tc_tracks)
    tc_df = get_track_datetime(tc_df)
    tc_df2 = tc_df[tc_df['vt100'] != 0]  # remove zero wind values

    try:
        lf_df = get_landfall_info(
            tc_df=tc_df2,
            gage_locs=gage_locs,
            clip_gdf=region
        )
        lf_df['tc_id'] = tc_id
        landfall_df = pd.concat(
            objs=[landfall_df, lf_df],
            axis=0,
            ignore_index=True
        )
    except Exception:
        issue_tcs.append(tc_id)
        print(f'Issue with {tc_id}')

    print(f'{counter} out of {tot}')
    counter += 1

out = landfall_df
out.to_csv('ncep_landfall_vmax.csv')
