"""
================================================================================
Description:
    This script processes synthetic tropical cyclone (TC) rainfall and wind data
    for the Carolinas region under a given climate scenario. It performs the 
    following key tasks:
      1. Sets up the working environment and imports necessary libraries.
      2. Reads ADCIRC modeled storm tide peaks to select relevant TCs.
      3. Computes storm total precipitation and maximum rain rate for selected TCs.
      4. Computes wind speed statistics (max, mean, thresholded) for selected TCs.
      5. Loads the study domain from a Hydromt data catalog.
      6. Creates raster masks of basins for rainfall and wind grids with spatial reference.
    
Inputs:
    - ADCIRC gage peak CSV
    - Gridded TCR rainfall NetCDF files (hourly)
    - Regridded wind NetCDF files
    - Hydromt data catalog YAML file: data_catalog_SFINCS_Carolinas.yml
Outputs:
    - NetCDF of TC precipitation stats
    - NetCDF of TC wind speed stats
    - Basin raster masks for precipitation and wind grids

Dependencies:
    - Local package: src.utils (TCR_precip_stats2netcdf, TC_windspd_stats2netcdf)
================================================================================
"""

import sys
import os
import pandas as pd
from pathlib import Path

# Add custom module path to sys.path for importing utils functions
sys.path.append(r'C:\Users\lelise\Documents\GitHub\flood_model_carolinas\synthetic_tc_cmpdfld')
from src.utils import TCR_precip_stats2netcdf, TC_windspd_stats2netcdf

import hydromt  # hydrometeorology package for data catalog and geospatial tools


# ------------------------
# Set working directory
# ------------------------
os.chdir(r'Z:\Data-Expansion\users\lelise\projects\Carolinas_SFINCS\Chapter3_SyntheticTCs\02_DATA\CMIP6_585')


# ------------------------
# Load ADCIRC storm tide peaks and select TCs
# ------------------------
gage_peaks = pd.read_csv('.\stormTide\gage_peaks_ZerosRemoved_canesm_ssp585.csv', index_col=0)
selected_tcs = gage_peaks.index.tolist()  # List of TC IDs to process


# ------------------------
# Process TCR rainfall for selected TCs
# ------------------------
# Calculate cumulative precipitation, max and mean rain rates per TC
ds = TCR_precip_stats2netcdf(
    tc_ids=selected_tcs,
    rr_threshold=5,
    inputdir=Path(r'.\rain\03_TCR_RainOutput_Gridded_hourly'),
    outputdir=Path(r'.\rain')
)


# ------------------------
# Process TC wind statistics for selected TCs
# ------------------------
ds2 = TC_windspd_stats2netcdf(
    tc_ids=selected_tcs,
    inputdir=Path(r'.\wind\CLE15_ReGridded_canesm_ssp585cal'),
    outputdir=Path(r'.\wind')
)


# ------------------------
# Load study domain from Hydromt data catalog
# ------------------------
data_catalog_yml = r'Z:\Data-Expansion\users\lelise\data\data_catalog_SFINCS_Carolinas.yml'
cat = hydromt.DataCatalog(data_libs=[data_catalog_yml])
mod_domain = cat.get_geodataframe(data_like='enc_domain_HUC6_clipped').to_crs(4326)


# ------------------------
# Create basin mask for rainfall grids
# ------------------------
grid = ds.sel(tc_id=selected_tcs[0])['total_precip'].rio.write_crs(4326)
mask = grid.raster.rasterize(mod_domain, 'index', nodata=-9999.0, all_touched=False)
mask.raster.to_raster(r'.\rain\basin_mask_TCR.tif', nodata=-9999.0)


# ------------------------
# Create basin mask for wind grids
# ------------------------
grid = ds2.sel(tc_id=selected_tcs[0]).rio.write_crs(4326)
mask = grid.raster.rasterize(mod_domain, 'index', nodata=-9999.0, all_touched=False)
mask.raster.to_raster(r'.\wind\basin_mask_wind.tif', nodata=-9999.0)
