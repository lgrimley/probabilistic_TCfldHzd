# ==========================================================================
# Aggregation and Summary Statistics of TCR Rainfall
#
# Purpose:
#   Computes storm-level rainfall summary statistics and total precipitation
#   grids from gridded TCR rainfall NetCDF files.
#   Outputs are organized by storm selection category and include both
#   tabular statistics and spatial NetCDF products.
#
# Key Features:
#   - Processes multiple storm subsets (e.g., within/outside buffer)
#   - Computes total precipitation over storm lifetime
#   - Extracts maximum rainfall metrics for each storm
#   - Aggregates storm-total precipitation grids into a single NetCDF
#
# Inputs:
#   - Gridded rainfall NetCDF files (per storm)
#   - CSV files listing storm IDs by category
#
# Outputs:
#   - CSV file containing per-storm rainfall statistics
#   - NetCDF file containing storm-total precipitation grids
#
# Metrics Computed:
#   - max_total_precip_cell : Maximum storm-total precipitation at any grid cell
#   - cumm_precip_domain   : Total precipitation summed over space and time
#   - max_rain_rate        : Maximum instantaneous rainfall rate
#
# Assumptions:
#   - Rainfall units are consistent across all NetCDF files
#   - Time dimension represents storm duration
#   - Spatial grids are identical across storms
#
# ==========================================================================

import pandas as pd
import os
import xarray as xr

# Directory containing gridded TCR rainfall NetCDF files
rain_dir = r'Z:\Data-Expansion\users\lelise\Chapter3\NCEP_Reanalysis\rain'

# Directory containing storm track CSV files
tracks_dir = r'Z:\Data-Expansion\users\lelise\Chapter3\NCEP_Reanalysis\tracks'

# List of storm-selection categories to process
track_files = ['UScoast6_AL_ncep_reanal_roEst1rmEst1_trk100_200km',
               'adcirc_modeled_TCs_all',
               'adcirc_modeled_TCs_outside200kmBuffer',
               'noSurge_TCs_within200kmBuffer']

# Loop through each storm-selection category
for f in track_files:

    # Load storm track information for this category
    tc_d = pd.read_csv(os.path.join(tracks_dir, f'{f}.csv'))

    # List to store storm-total precipitation grids
    total_precip_grids = []

    # DataFrame to store per-storm rainfall statistics
    tc_rain_info = pd.DataFrame()

    # Number of storms in this category
    n_tc = len(tc_d['tc_id'])

    # Loop through each storm
    for i in range(n_tc):
        tc_id = tc_d['tc_id'][i]

        # Construct rainfall NetCDF filename
        grid = f'{str(tc_id).zfill(4)}.nc'

        # Load gridded rainfall dataset
        d = xr.open_dataset(os.path.join(rain_dir, '03_TCR_RainOutput_Gridded', grid))

        # Compute total precipitation over the storm duration
        total_precip = d.sum(dim='time')

        # Store total precipitation grid
        total_precip_grids.append(total_precip)

        # ------------------------------------
        # Compute storm-level rainfall metrics
        # ------------------------------------

        # Maximum storm-total precipitation at any grid cell
        max_total_precip_cell = round(d.sum(dim='time').max()['precip'].item(), 3)

        # Cumulative precipitation summed over all grid cells and time
        cumm_precip_domain = round(d.sum()['precip'].item(), 3)

        # Maximum instantaneous rainfall rate anywhere in the domain
        max_rain_rate = round(d.max(dim='time').max()['precip'].item(), 3)

        # Store metrics in a temporary DataFrame
        x = pd.DataFrame([tc_id,
                          max_total_precip_cell,
                          cumm_precip_domain,
                          max_rain_rate]).T

        # Append to summary table
        tc_rain_info = pd.concat(objs=[tc_rain_info, x], axis=0)

        # Progress update
        print(f'{i} out of {n_tc}')

    # ------------------------------------
    # Write outputs for this storm category
    # ------------------------------------

    # Create output directory if it does not exist
    outdir = os.path.join(rain_dir, f)
    if os.path.exists(outdir) is False:
        os.mkdir(outdir)

    # Assign column names to summary statistics table
    tc_rain_info.columns = ['tc_id',
                            'max_total_precip_cell',
                            'cumm_precip_domain',
                            'max_rain_rate']

    # Write storm-level rainfall statistics to CSV
    tc_rain_info.to_csv(os.path.join(outdir, 'tc_rain_info.csv'), index=False)

    # Concatenate all storm-total precipitation grids along a new dimension
    da = xr.concat(total_precip_grids, dim='run')

    # Assign storm IDs as the coordinate for the run dimension
    da['run'] = xr.IndexVariable('run', tc_rain_info['tc_id'].values)

    # Write storm-total precipitation grids to NetCDF
    da.to_netcdf(os.path.join(outdir, 'TC_TotPrecip_grids.nc'))
