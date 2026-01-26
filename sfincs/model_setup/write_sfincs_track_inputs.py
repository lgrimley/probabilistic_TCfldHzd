"""
Description:
    This script generates SFINCS boundary condition and meteorological forcing
    files for a single synthetic tropical cyclone track. It is intended to be run 
    from the command line with arguments specifying the TC index (TC ID) and output directory.

    The script:
    - Loads a synthetic TC track and associated forcing datasets
    - Updates SFINCS configuration timing based on storm duration
    - Generates water level, precipitation, and wind forcing inputs
    - Fills missing storm tide locations with tidal reanalysis data
    - Writes SFINCS forcing and configuration files
    - Copies outputs to a track-specific directory
    - Generates tide-only boundary condition files for runoff-only scenarios

Usage:
    python write_sfincs_track_inputs.py <tc_index> <track_dir>

Notes:
    - Paths and data catalogs are hard-coded
    - Output formatting and file naming follow SFINCS requirements
    - Script is designed for batch execution across multiple TC indices
"""

import sys
import os
import shutil
import hydromt
from hydromt_sfincs import SfincsModel
from src.core import cansem_ssp585_DataPaths, SyntheticTrack

# Add local repository path for custom modules
sys.path.append(
    r'C:\Users\lelise\Documents\GitHub\flood_model_carolinas\syntheticTCs_cmpdfld'
)

# ---------------------------------------------------------------------
# Parse command-line arguments
# ---------------------------------------------------------------------
arg1 = sys.argv[1]
tc_index = int(arg1)           # Synthetic tropical cyclone index
track_dir = sys.argv[2]        # Output directory for this TC

# Assign data paths and data catalog
DataPaths = cansem_ssp585_DataPaths
DatCat = hydromt.DataCatalog(
    fr'Z:\Data-Expansion\users\lelise\projects\Carolinas_SFINCS'
    r'\Chapter3_SyntheticTCs\02_DATA\CMIP6_585\data_catalog_canesm_ssp585.yml'
)

'''
Code below doesn't need to be edited
'''

# ---------------------------------------------------------------------
# Setup working directories and base SFINCS model
# ---------------------------------------------------------------------
cat_dir = r'Z:\Data-Expansion\users\lelise\data'
yml_base_CONUS = os.path.join(cat_dir, 'data_catalog_BASE_CONUS.yml')
yml_base_Carolinas = os.path.join(cat_dir, 'data_catalog_BASE_Carolinas.yml')
yml_sfincs_Carolinas = os.path.join(cat_dir, 'data_catalog_SFINCS_Carolinas.yml')
yml_sfincs_Carolinas_Ch3 = os.path.join(cat_dir, 'data_catalog_SFINCS_Carolinas_Ch3.yml')

# Initialize SFINCS model in read-only mode
mod = SfincsModel(
    root=r'Z:\Data-Expansion\users\lelise\projects\Carolinas_SFINCS'
         r'\Chapter3_SyntheticTCs\03_MODEL\sfincs_base_mod',
    mode='r',
    data_libs=[
        yml_base_CONUS,
        yml_base_Carolinas,
        yml_sfincs_Carolinas,
        yml_sfincs_Carolinas_Ch3
    ]
)

try:
    # -----------------------------------------------------------------
    # Load synthetic TC track and associated datasets
    # -----------------------------------------------------------------
    track = SyntheticTrack(
        DataPaths=DataPaths,
        DatCat=DatCat,
        tc_index=tc_index
    )

    # -----------------------------------------------------------------
    # Update SFINCS configuration for storm timing
    # -----------------------------------------------------------------
    mod.setup_config(
        **{
            "tref": track.time_range[0],
            "tstart": track.time_range[0],
            "tstop": track.time_range[1]
        }
    )
    mod.write_config(config_fn='sfincs.inp')

    # Setup primary water level forcing from ADCIRC storm tide
    mod.setup_waterlevel_forcing(
        geodataset=track.merged_waterlevel,
        offset='lmsl_to_navd88',
        buffer=100000,
        merge=False
    )

    # -----------------------------------------------------------------
    # Fill missing storm tide locations with tidal reanalysis data
    # -----------------------------------------------------------------
    tide_stations = track.tide_offset_table[
        track.tide_offset_table['st_ht_dt'].isna()
    ]['Point'].tolist()

    print(
        f'Tides used where ADCIRC stormtide output is not available at gages: '
        f'{tide_stations}'
    )

    if len(tide_stations) > 0:
        tide_stations_int = [int(s.split('P')[1]) for s in tide_stations]

        tide_filler_data = track.reanalysis_data_offset.sel(
            index=tide_stations
        )
        tide_filler_data = tide_filler_data.assign_coords(
            index=tide_stations_int
        )

        mod.setup_waterlevel_forcing(
            geodataset=tide_filler_data,
            offset='lmsl_to_navd88',
            buffer=100000,
            merge=True
        )

    x = len(mod.forcing['bzs'].index.values)
    print(f'Total number of gage inputs for the track: {x}')

    # -----------------------------------------------------------------
    # Setup meteorological forcing (rainfall and wind)
    # -----------------------------------------------------------------
    mod.setup_precip_forcing_from_grid(
        precip=track.precip,
        aggregate=False
    )
    mod.setup_wind_forcing_from_grid(wind=track.wind)
    print('Wind and rainfall SFINCS inputs created.')

    # Write all forcing files
    print('Writing the boundary condition files.')
    mod.write_forcing()

    # -----------------------------------------------------------------
    # Copy generated files to the track-specific directory
    # -----------------------------------------------------------------
    files_2_move = [
        'precip_2d.nc',
        'wind_2d.nc',
        'sfincs.bnd',
        'sfincs.bzs',
        'sfincs.inp'
    ]

    for file in files_2_move:
        file_src = os.path.join(mod.root, file)
        file_dst = os.path.join(track_dir, file)
        _ = shutil.copyfile(file_src, file_dst)
        print(f'Created {file_dst}')

    # -----------------------------------------------------------------
    # Generate tide-only boundary conditions (runoff-only scenario)
    # -----------------------------------------------------------------
    tides_only = track.reanalysis_data.sel(
        index=track.coastal_locations_mapped['Point'].tolist()
    )

    tide_stations_int = [
        int(s.split('P')[1])
        for s in track.coastal_locations_mapped['Point'].tolist()
    ]

    tides_only2 = tides_only.assign_coords(index=tide_stations_int)

    mod.setup_waterlevel_forcing(
        geodataset=tides_only2,
        offset='lmsl_to_navd88',
        buffer=100000,
        merge=False
    )

    mod.write_forcing()
    print('Tides only bnd and bzs files created.')

    files_2_move = ['sfincs.bnd', 'sfincs.bzs']
    for file in files_2_move:
        file_src = os.path.join(mod.root, file)
        file_dst = os.path.join(track_dir, f'tides_{file}')
        _ = shutil.copyfile(file_src, file_dst)
        print(f'Created {file_dst}')

except:
    # Catch-all to allow batch processing to continue
    print(f'Issue with TC: {tc_index}')
