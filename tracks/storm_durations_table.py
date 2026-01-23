import matplotlib.pyplot as plt
import scipy.io as sio
from src.utils import get_track_info_in_df

import os
import pandas as pd
import numpy as np


"""
This script analyzes the temporal characteristics of selected synthetic tropical cyclones (TCs) 

Specifically, it:
    Extracts start and end times for each selected TC to a table to compute storm durations.
    We then visualize the distribution of TC durations using a histogram and
    estimate avg computational cost for hydrodynamic simulations.
"""

# -------------------------------------------------------------------------
# Set working directory and define track data location
# -------------------------------------------------------------------------
climate_scenario = 'CMIP6_585'
os.chdir(
    rf'Z:\Data-Expansion\users\lelise\projects\Carolinas_SFINCS\Chapter3_SyntheticTCs\02_DATA\{climate_scenario}'
)

track_dir = os.path.join(os.getcwd(), 'tracks')

# -------------------------------------------------------------------------
# Load synthetic TC track data
# -------------------------------------------------------------------------
# Here we focus on a single emission scenario (SSP5-8.5) and GCM member (CanESM)
gcm_ensemble = 'canesm'
ssp = 'ssp585'
track_file = os.path.join(
    track_dir,
    f'UScoast6_AL_{gcm_ensemble}_{ssp}cal_roEst1rmEst1_trk100'
)
tc_tracks = sio.loadmat(f'{track_file}.mat')

# -------------------------------------------------------------------------
# Load ADCIRC-modeled peak water levels and extract TC IDs
# -------------------------------------------------------------------------
# Each row corresponds to a TC that was modeled in ADCIRC and includes the storm surge peaks for each TC
gage_peaks = pd.read_csv(
    fr'.\stormTide\gage_peaks_{gcm_ensemble}_{ssp}.csv',
    index_col=0
)

# Get the TC IDs used for analysis
selected_tcs = gage_peaks.index.tolist()

# -------------------------------------------------------------------------
# Extract start and end times for each TC from the track file (.mat)
# -------------------------------------------------------------------------
tc_time_info = pd.DataFrame(
    index=selected_tcs,
    columns=['tstart', 'tend']
)

for tc_id in selected_tcs:

    # Retrieve track information with datetime for the TC
    track_df = get_track_info_in_df(tc_id=tc_id, tc_tracks=tc_tracks)

    # Remove zero-padded datetime entries
    track_df = track_df[track_df['datetime'] != 0]

    # Store start and end times of the TC
    tc_time_info.loc[tc_id, 'tstart'] = track_df['datetime'].iloc[0]
    tc_time_info.loc[tc_id, 'tend']   = track_df['datetime'].iloc[-1]

# Print the temporal range of TC start times
print(tc_time_info['tstart'].min())
print(tc_time_info['tstart'].max())

# -------------------------------------------------------------------------
# Compute TC durations
# -------------------------------------------------------------------------
tc_time_info['duration'] = (
    tc_time_info['tend'] - tc_time_info['tstart']
)

# Summary statistics of duration
mean   = np.mean(tc_time_info['duration'].values)
maxval = np.max(tc_time_info['duration'].values)
median, pct95 = np.percentile(
    tc_time_info['duration'],
    q=[50, 95]
)

# Convert duration to days
tc_time_info['duration_days'] = [
    x.total_seconds() / 86400.0
    for x in tc_time_info['duration']
]

# -------------------------------------------------------------------------
# Plot histogram of TC durations
# -------------------------------------------------------------------------
fig, ax = plt.subplots(
    nrows=1,
    ncols=1,
    figsize=(6, 4),
    tight_layout=True
)

tc_time_info['duration_days'].plot.hist(
    ax=ax,
    legend=False
)

ax.set_xlabel('TC Duration (days)')

plt.savefig(r'.\tracks\TC_duration_histogram.png')
plt.close()

# -------------------------------------------------------------------------
# Estimate computational cost for compound simulations
# -------------------------------------------------------------------------
# Assumption:
# - A reference 8-day TC simulation takes ~1 hour using ~32 CPUs
# - Cost is scaled linearly with storm duration
# - A factor of 3 is applied to account for compound or ensemble runs

tc_time_info['cpu_time_hrs_x3'] = (
    tc_time_info['duration_days'] * (1 / 8)
) * 3

# Total estimated CPU time (rounded up and doubled for safety margin)
total_cpu_time = np.ceil(
    tc_time_info['cpu_time_hrs_x3'].sum()
) * 2
