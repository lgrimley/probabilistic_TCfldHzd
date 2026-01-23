import matplotlib.pyplot as plt
import xarray as xr
import scipy.io as sio
import os
import warnings

# Suppress future warnings (e.g., from xarray or pandas)
warnings.simplefilter(action='ignore', category=FutureWarning)

# Project-specific utility and core functions
from src.utils import *
from src.core import *


"""
This script converts gridded tropical cyclone (TC) wind or rainfall data
associated with synthetic CMIP6 storms into NetCDF format with proper
time coordinates and spatial metadata.

The workflow:
- Loops over multiple CMIP6 GCMs (SSP5-8.5 scenario)
- Loads TC track metadata from MATLAB files
- Assigns datetime information to gridded TC data
- Writes standardized NetCDF files for downstream modeling (e.g., SFINCS)
"""

%% ------------------------------------------------------------------------
# Set working directory containing CMIP6 synthetic TC data
# ------------------------------------------------------------------------
os.chdir(
    r'Z:\Data-Expansion\users\lelise\projects\Carolinas_SFINCS\Chapter3_SyntheticTCs\02_DATA\CMIP6_585'
)

# List of CMIP6 GCMs included in the analysis
gcms = [
    'canesm_ssp585cal',
    'cnrm6_ssp585cal',
    'ecearth6_ssp585cal',
    'ipsl6_ssp585cal',
    'miroc6_ssp585cal'
]

%% ------------------------------------------------------------------------
# Loop through each GCM
# ------------------------------------------------------------------------
for gcm in gcms:

    # Directory containing original gridded TC data
    input_dir = os.path.join(os.getcwd(), 'wind', f'CLE15_Gridded_{gcm}')

    # Directory to store re-gridded / reformatted NetCDF output
    output_dir = os.path.join(os.getcwd(), 'wind', f'CLE15_ReGridded_{gcm}')

    # Create output directory if it does not exist
    if os.path.exists(output_dir) is False:
        os.makedirs(output_dir)

    # Load TC track metadata (MATLAB format)
    fname = rf'.\tracks\UScoast6_AL_{gcm}_roEst1rmEst1_trk100'
    tc_tracks = sio.loadmat(f'{fname}.mat')

    %% --------------------------------------------------------------------
    # Loop through individual TC NetCDF files
    # --------------------------------------------------------------------
    for file in os.listdir(input_dir):

        # Extract TC ID from filename (assumes filename is "<tc_id>.nc")
        tc_id = int(file.strip('.nc'))

        out_name = os.path.join(output_dir, file)

        # Only process files that do not already exist
        if os.path.exists(out_name) is False:
            print(tc_id)

            # --------------------------------------------------------------
            # Retrieve TC track information including datetime
            # --------------------------------------------------------------
            df = get_track_info_in_df(tc_id=tc_id, tc_tracks=tc_tracks)

            # Remove zero-padded datetime entries
            time = df['datetime'][df['datetime'] != 0]

            # --------------------------------------------------------------
            # Load gridded TC data
            # --------------------------------------------------------------
            data = xr.open_dataset(os.path.join(input_dir, file))

            # Assign proper time coordinate from TC track metadata
            data['time'] = time.values

            # Assign geographic coordinate reference system (WGS84)
            data = data.rio.write_crs("epsg:4326", inplace=True)

            # Optional: resample to hourly resolution (not currently used)
            # data = data.resample(time='1H').ffill()

            # --------------------------------------------------------------
            # Save standardized NetCDF file
            # --------------------------------------------------------------
            data.to_netcdf(out_name)
