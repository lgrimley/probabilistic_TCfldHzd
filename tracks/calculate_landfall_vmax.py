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
This script estimates approximate tropical cyclone (TC) landfall characteristics
(time and maximum wind speed) for both:

1. Synthetic CMIP6 tropical cyclones, and
2. Historical NCEP reanalysis tropical cyclones.

Landfall is approximated based on proximity to ADCIRC gage locations and
intersection with the SFINCS model domain. Results are saved as CSV files
for downstream analysis (e.g., wind-speed bias correction or scenario design).

---------------------------
Important Assumptions / Flags
---------------------------
- The SFINCS domain is buffered by 1.8 degrees for landfall detection.
  This expands the effective area used to detect landfall, but is not a precise
  coastline-based measure.
- In the CMIP6 processing loop, there is a `break` statement after the first TC.
  This is likely for testing/debugging purposes and currently limits processing
  to one TC per file.
- Landfall is proxy-based: we detect when a TC comes within proximity of
  ADCIRC gages inside the buffered model domain. This is not a physical landfall
  (coastline intersection) but an operational definition useful for surge/wind
  impact studies.
- Some TCs may fail the landfall detection step (logged in `issue_tcs`).
"""

# -------------------------------------------------------------------------
# Load SFINCS model domain and define analysis region
# -------------------------------------------------------------------------
mod = SfincsModel(
    root=r'Z:\Data-Expansion\users\lelise\projects\Carolinas_SFINCS\Chapter3_SyntheticTCs\03_MODEL_RUNS\sfincs_base_mod',
    mode='r'
)
# Convert model region to geographic coordinates and buffer by 1.8 degrees
region = mod.region.to_crs(4326).buffer(1.8)  # FLAG: buffered domain, not exact coastline

# -------------------------------------------------------------------------
# Load ADCIRC gage locations used for landfall estimation
# -------------------------------------------------------------------------
gage_locs = pd.read_csv(
    r'Z:\Data-Expansion\users\lelise\projects\Carolinas_SFINCS\Chapter3_SyntheticTCs'
    r'\02_DATA\NCEP_Reanalysis\stormTide\coastline_lores_NCSC.csv'
)

# -------------------------------------------------------------------------
# Process CMIP6 synthetic TC tracks
# -------------------------------------------------------------------------
os.chdir(
    r'Z:\Data-Expansion\users\lelise\projects\Carolinas_SFINCS'
    r'\Chapter3_SyntheticTCs\02_DATA\CMIP6_245\tracks'
)

for file in os.listdir(os.getcwd()):

    if file.endswith('.mat'):
        print(file)

        # Load TC track data
        tc_tracks = sio.loadmat(file)

        # Parse GCM and scenario from filename
        gcm = file.split('_')[2]
        ssp = file.split('_')[3].removesuffix('cal')

        # Load TC IDs modeled in ADCIRC
        modeled_tcs = pd.read_csv(
            fr'../stormTide/gage_peaks_{gcm}_{ssp}.csv',
            index_col=0
        )
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

            break  # FLAG: stops after one TC, for debugging/testing

        out = landfall_df[['tc_id', 'vstore100']]
        out.to_csv(fr'{gcm}_{ssp}_landfall_vmax.csv')

# -------------------------------------------------------------------------
# Process NCEP reanalysis tropical cyclones
# -------------------------------------------------------------------------
os.chdir(
    r'Z:\Data-Expansion\users\lelise\projects\Carolinas_SFINCS'
    r'\Chapter3_SyntheticTCs\02_DATA\NCEP_Reanalysis\tracks'
)

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
out.to_csv('ncep_landfall_vmax_v2.csv')

# -------------------------------------------------------------------------
# METHODS PARAGRAPH (for documentation/publication)
# -------------------------------------------------------------------------
"""
Landfall Definition:

Landfall is operationally defined as the first instance when a tropical cyclone
approaches or intersects the buffered SFINCS model domain within proximity
of ADCIRC gage locations. Specifically:

1. The SFINCS domain is projected to geographic coordinates (EPSG:4326) and
   buffered by 1.8 degrees (~200 km) to account for model resolution and
   potential storm track uncertainties.

2. ADCIRC gage locations serve as reference points to identify coastal impact.

3. TC track data are filtered for valid timestamps and wind speeds
   (vstore100 or vt100 > 0).

4. Landfall time and maximum wind speed are extracted using a spatial overlay
   between the buffered domain, gage locations, and TC tracks.

This proxy-based definition provides a practical estimate of landfall for
hydrodynamic modeling (storm surge, flooding, wind) but does not necessarily
correspond to the exact point where the storm center crosses the coastline.
"""
