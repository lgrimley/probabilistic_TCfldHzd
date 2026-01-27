"""
===============================================================================
Description:
------------
This script processes SFINCS model outputs for synthetic tropical cyclone (TC)
simulations over the Carolinas. It:
- Loads a base SFINCS model configuration
- Identifies which synthetic storms have valid ADCIRC/SFINCS outputs
- Groups storms by landfall maximum wind speed
- Extracts flood hazard and attribution results using TCFloodHazard
- Writes grouped results to NetCDF files for downstream analysis

Key Outputs:
------------
- NetCDF files of maximum water level (zsmax) per TC group
- NetCDF files of attribution results per TC group

------------
Before running, confirm that the following directories exist and are accessible:
- SFINCS base model directory
- Completed TC run directories
- Input CSV files for storm metadata
- Output directory

Author: Lauren Grimley
Created: 01/23/2025
===============================================================================
"""

import os

# Core hydrologic and geospatial modeling libraries
import hydromt
import xarray as xr
from hydromt_sfincs import SfincsModel

# Standard utilities
from pathlib import Path
import pandas as pd
import numpy as np
import sys
import time

# Add local project path so custom modules can be imported
sys.path.append(r'C:\Users\lelise\Documents\GitHub\flood_model_carolinas\syntheticTCs_cmpdfld')

# Custom class for extracting TC flood hazard results
from src.core import TCFloodHazard

# ------------------------------------------------------------------------------
# Load SFINCS base model (used for reading outputs consistently)
# ------------------------------------------------------------------------------

# Path to HydroMT data catalog
yml_base = r'Z:\Data-Expansion\users\lelise\data\data_catalog_BASE_Carolinas.yml'

# Root directory for the initialized SFINCS model
base_root = r'Z:\Data-Expansion\users\lelise\projects\Carolinas_SFINCS\Chapter3_SyntheticTCs\03_MODEL\sfincs_initcond_mod'

# Initialize SFINCS model in read-only mode
sfincs_mod = SfincsModel(root=base_root, mode='r', data_libs=yml_base)

# Read model configuration and static data
# NOTE: This is intentionally called twice due to a known quirk
sfincs_mod.read()

# ------------------------------------------------------------------------------
# Set working directory to synthetic TC project
# ------------------------------------------------------------------------------

os.chdir(r'Z:\Data-Expansion\users\lelise\projects\Carolinas_SFINCS\Chapter3_SyntheticTCs')

# Scenario and grouping metadata
run_group = 'canesm_ssp585'
foldname = 'CMIP6_585'
group_size = 20

# Root directory containing completed SFINCS TC runs
root = Path(
    fr'Z:\Data-Expansion\users\lelise\projects\Carolinas_SFINCS\Chapter3_SyntheticTCs\03_MODEL\{run_group}_runs\completed_runs_SLR112cm'
)

# ------------------------------------------------------------------------------
# Load and filter storm metadata
# ------------------------------------------------------------------------------

# CSV containing all storms and their landfall maximum wind speeds
tc_ids_vmax = pd.read_csv(
    fr'.\02_DATA\{foldname}\tracks\landfall_vmax\adcirc_tcs\{run_group}_landfall_vmax.csv',
    index_col=0
)

# CSV listing storms that actually have ADCIRC storm tide outputs
correct = pd.read_csv(
    fr'.\02_DATA\{foldname}\stormTide\gage_peaks_ZerosRemoved_{run_group}.csv',
    index_col=0
)

# Select only storms that:
# 1) Appear in the landfall vmax file
# 2) Have valid ADCIRC outputs
selected = tc_ids_vmax[tc_ids_vmax['tc_id'].isin(correct.index)]

# Sort storms by landfall wind speed (ascending)
tc_ids_vmax_ordered = selected.sort_values(by='vstore100', ascending=True)

# Split storms into fixed-size groups for batch NetCDF writing
groups = [
    tc_ids_vmax_ordered[i:i + group_size]
    for i in range(0, len(tc_ids_vmax_ordered), group_size)
]

# ------------------------------------------------------------------------------
# Prepare output directories
# ------------------------------------------------------------------------------

outputdir = r'Z:\Data-Expansion\users\lelise\projects\Carolinas_SFINCS\Chapter3_SyntheticTCs\04_RESULTS'
outdir = os.path.join(outputdir, run_group)

# Create output directory if it does not exist
if os.path.exists(outdir) is False:
    os.makedirs(outdir)

# ------------------------------------------------------------------------------
# Loop over TC groups and process results
# ------------------------------------------------------------------------------

counter = 0
for group in groups:

    # Define group label based on min/max landfall wind speed
    vmin, vmax = (int(group['vstore100'].min()), int(group['vstore100'].max()))
    group_info = f'VmaxLf_{vmin}_{vmax}'
    start_time = time.time()

    # Containers for model outputs
    scenario_results = []
    attribution_results = []
    keep_ids = []

    # Loop over individual tropical cyclones in the group
    for tc_index in group['tc_id'].values:

        # Directory containing SFINCS outputs for this TC
        tc_root = os.path.join(root, f'TC_{str(tc_index).zfill(4)}')

        # Skip storms that did not produce runoff outputs
        if os.path.exists(os.path.join(tc_root, 'runoff')) is False:
            print(f'skipping {tc_index}')
            continue
        else:
            try:
                # Load flood hazard results for this TC
                r = TCFloodHazard(
                    tc_root=Path(tc_root),
                    sfincs_mod=sfincs_mod,
                    tc_index=tc_index,
                    tracks=run_group,
                    zsmax_threshold=0.05
                )

                # Store outputs for later concatenation
                scenario_results.append(r.scenario_results)
                attribution_results.append(r.attribution_results)
                keep_ids.append(tc_index)

                print(f'Successfully processed {tc_index}')

            except Exception as e:
                # Catch and report failures without stopping the batch
                print(e)
                print(f'Issue with {tc_index}')
                continue

    # ------------------------------------------------------------------------------
    # Write grouped outputs to NetCDF
    # ------------------------------------------------------------------------------

    attrs = {
        'description': f'Modeled SFINCS ouputs for {group_info}',
        'track_model': run_group,
        'author': 'Lauren Grimley',
        'date_created': '1/23/2025'
    }

    # Concatenate scenario results across TCs
    ds1 = xr.concat(objs=scenario_results, dim='tc_id')
    ds1['tc_id'] = xr.IndexVariable(dims='tc_id', data=keep_ids)

    # Write maximum water level (zsmax) output
    var = 'zsmax'
    outfilepath = os.path.join(outputdir, run_group, f'{var}_{group_info}.nc')

    if os.path.exists(outfilepath) is False:
        print(f'Writing {var}...')
        dsout = ds1[var]
        dsout.attrs = attrs
        dsout.to_netcdf(outfilepath)

    # Write attribution results
    outfilepath = os.path.join(outputdir, run_group, f'attribution_{group_info}.nc')

    if os.path.exists(outfilepath) is False:
        ds2 = xr.concat(objs=attribution_results, dim='tc_id')
        ds2['tc_id'] = xr.IndexVariable(dims='tc_id', data=keep_ids)
        ds2.attrs = attrs
        print(f'Writing attribution...')
        ds2.to_netcdf(outfilepath)

    # ------------------------------------------------------------------------------
    # Timing and progress reporting
    # ------------------------------------------------------------------------------

    print(f'Finished writing output for {group_info}')
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time writing {group_info}: {elapsed_time} seconds")

    counter += 1
    print(f'{counter} out of {len(groups)} groups processed...')
