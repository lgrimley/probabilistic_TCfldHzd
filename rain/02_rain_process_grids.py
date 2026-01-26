# ==========================================================================
# Formatting TCR Rainfall for Input to SFINCS: Hourly Regridding
#
# Purpose:
#   Converts storm-based, gridded TCR rainfall NetCDF files from 2-hourly
#   resolution to hourly resolution and aligns rainfall time coordinates
#   with storm track datetimes. Outputs are written as hourly NetCDF files
#   suitable for use as rainfall forcing in SFINCS.
#
# Key Features:
#   - Associates rainfall fields with storm track timestamps
#   - Downscales 2-hour rainfall rates to hourly resolution
#   - Preserves total precipitation via consistency checks
#
# Inputs:
#   - TCR gridded rainfall NetCDF files (2-hourly)
#   - Storm track MAT files
#   - CSV files listing selected storm IDs
#
# Outputs:
#   - Hourly gridded rainfall NetCDF files (one per storm)
#
# Assumptions:
#   - TCR rainfall represents 2-hour accumulations or rates
#   - First two track timesteps contain zero rainfall and are discarded
#   - Rainfall grids are in EPSG:4326
#
# Notes:
#   - Total precipitation before and after resampling must match within
#     a tolerance; storms failing this check are skipped
#
# ==========================================================================

# Suppress future warnings for cleaner output
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

# Import core and utility functions from project source
from src.core import *
from src.utils import *

# Change working directory to CMIP6 SSP585 data location
os.chdir(r'Z:\Data-Expansion\users\lelise\projects\Carolinas_SFINCS\Chapter3_SyntheticTCs\02_DATA\CMIP6_585')

# List of GCM scenarios to process
gcms = ['canesm_ssp585cal', 'cnrm6_ssp585cal', 'ecearth6_ssp585cal',
        'ipsl6_ssp585cal', 'miroc6_ssp585cal']

# Track missing files or storms
missing = []

# Loop through each GCM
for gcm in gcms:
    missing.append(gcm)

    # Load storm track MAT file for this GCM
    fname = f'UScoast6_AL_{gcm}_roEst1rmEst1_trk100.mat'
    tc_tracks = sio.loadmat(fr'.\tracks\{fname}')

    # Load CSV containing selected storm IDs
    tc_select = pd.read_csv(
        rf'Z:\Data-Expansion\users\lelise\projects\Carolinas_SFINCS\Chapter3_SyntheticTCs'
        rf'\02_DATA\CMIP6_585\stormTide\stormTide_TCIDs_and_gageCounts_{gcm}.csv'
    )

    ''' RAINFALL to netcdf '''

    # Input directory containing 2-hourly gridded TCR rainfall
    input_dir = fr'.\rain\TCR_Gridded_{gcm}'

    # Output directory for hourly rainfall NetCDFs
    output_dir = fr'.\rain\TCR_Gridded_{gcm}_hourly'

    # Create output directory if it does not exist
    if os.path.exists(output_dir) is False:
        os.makedirs(output_dir)

    # Loop through each selected storm ID
    for tc_id in tc_select['tc_id'].tolist():
        tc_id = int(tc_id)

        # Construct rainfall NetCDF filename
        file = f'{str(tc_id).zfill(4)}.nc'

        # Check that the raw gridded rainfall file exists
        if os.path.exists(os.path.join(input_dir, file)) is False:
            missing.append(tc_id)
            print(f'{tc_id} raw rainfall has not been processed!')
            continue

        # Output hourly rainfall filename
        out_name = os.path.join(output_dir, file)

        # Skip storm if hourly file already exists
        if os.path.exists(out_name) is False:
            print(f'Regridding rainfall for {tc_id}')

            # Extract storm track information with datetime
            df = get_track_info_in_df(tc_id=tc_id, tc_tracks=tc_tracks)

            # Remove zero-valued datetime entries
            time = df['datetime'][df['datetime'] != 0]

            # Load gridded TCR rainfall dataset
            rain_data = xr.open_dataset(os.path.join(input_dir, file))

            # The first two track timesteps have zero rainfall
            # Align rainfall timesteps with track datetimes
            time_rain = time.iloc[2:(2+len(rain_data['time']))]
            rain_data['time'] = time_rain.values

            # Assign geographic coordinate reference system
            rain_data = rain_data.rio.write_crs("epsg:4326", inplace=True)

            # Convert 2-hour rainfall to hourly-equivalent rate
            rain_data = rain_data * 2.0

            # Downscale to hourly resolution via temporal resampling
            rain_data_1hr = rain_data.resample(time='1H').mean()

            # Compute total precipitation before and after resampling
            totP_2hr = rain_data.sum()['precip'].item()
            totP_1hr = rain_data_1hr.sum()['precip'].item()

            # Verify precipitation conservation
            if abs(totP_2hr - totP_1hr) > 100:
                print(f'Check {tc_id}: precip difference {abs(totP_2hr - totP_1hr)}')
                continue
            else:
                # Write hourly rainfall dataset to NetCDF
                rain_data_1hr.to_netcdf(out_name)

                # Optional diagnostic plotting (commented out)
                # track_gdf = track_as_gdf(df)
                # tot_rain = rain_data.sum(dim='time')['precip']
                # fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6, 4), tight_layout=True)
                # tot_rain.plot(ax=ax, cmap='jet')
                # track_gdf.geometry.plot(ax=ax, color='black')
                # usa.geometry.plot(ax=ax, color='none', edgecolor='black')
                # plt.savefig(os.path.join(figs_dir, file.replace('.nc', '_total_precip.png')))
                # plt.close()
