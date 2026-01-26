"""
Description:
    This script generates a YAML data catalog for synthetic tropical cyclone
    forcing datasets used in SFINCS modeling.
    The catalog includes references to gridded rainfall, wind fields, storm tide
    water levels, and a tidal reanalysis dataset.

    For each dataset type, the script:
    - Iterates over NetCDF files in predefined directories
    - Extracts a tropical cyclone ID from the filename
    - Writes a structured YAML entry describing the dataset

Notes:
    - The output YAML file is overwritten initially and then appended to
    - Directory structure and file naming conventions are assumed to be consistent
    - Paths are written exactly as provided for downstream ingestion
"""

import os

# Change working directory to the data location
os.chdir(r'Z:\Data-Expansion\users\lelise\projects\Carolinas_SFINCS\Chapter3_SyntheticTCs\02_DATA\CMIP6_585')

# Name of the output YAML data catalog
data_lib = 'data_catalog_canesm_ssp585.yml'

# ---------------------------------------------------------------------
# Rainfall datasets: hourly gridded TCR precipitation
# ---------------------------------------------------------------------
data_dir = r'.\rain\TCR_Gridded_canesm_ssp585cal_hourly'
with open(data_lib, mode="w+") as fcat:
    for file in os.listdir(data_dir):
        if file.endswith('.nc'):
            # Extract tropical cyclone ID from filename
            tc_id = file.split('.')[0]

            # YAML entry for rainfall dataset
            yml_str = f"""
precip_{tc_id}:
  path: {os.path.join(data_dir, file)}
  data_type: RasterDataset
  driver: netcdf
  crs: 4326
  meta:
    category: meteo
    units: mm/hr
    description: hourly gridded rain rates for synthetic tropical cyclone tracks
    title: Synthetic TCR rainfall
    model_ref: https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2013GL058284
    publication_date: 2022
    model: canesm ssp585cal
    data_ref: https://www.nature.com/articles/s41558-021-01272-7"""

            # Write entry to catalog
            fcat.write(yml_str)

# ---------------------------------------------------------------------
# Wind datasets: parametric wind fields
# ---------------------------------------------------------------------
data_dir = r'.\wind\CLE15_ReGridded_canesm_ssp585cal'
with open(data_lib, mode="a") as fcat:
    for file in os.listdir(data_dir):
        if file.endswith('.nc'):
            # Extract tropical cyclone ID from filename
            tc_id = file.split('.')[0]

            # YAML entry for wind dataset
            yml_str = f"""
wind_{tc_id}:
  path: {os.path.join(data_dir, file)}
  data_type: RasterDataset
  driver: netcdf
  crs: 4326
  meta:
    category: meteo
    units: m/s
    description: hourly gridded wind speeds for synthetic tropical cyclone tracks
    title: TC parametric wind fields
    model_ref: https://journals.ametsoc.org/view/journals/mwre/108/8/1520-0493_1980_108_1212_aamotw_2_0_co_2.xml
    publication_date: 2022
    data_ref: https://www.nature.com/articles/s41558-021-01272-7"""

            # Append entry to catalog
            fcat.write(yml_str)

# ---------------------------------------------------------------------
# Storm tide datasets: ADCIRC water levels
# ---------------------------------------------------------------------
data_dir = r'.\stormTide\adcirc_waterlevel_netcdf_canesm_ssp585'
with open(data_lib, mode="a") as fcat:
    for file in os.listdir(data_dir):
        if file.endswith('.nc'):
            # Extract tropical cyclone ID from filename
            tc_id = file.split('.')[0]

            # YAML entry for storm tide dataset
            yml_str = f"""
stormTide_{tc_id}:
  path: {os.path.join(data_dir, file)}
  data_type: GeoDataset
  driver: netcdf
  crs: 4326
  meta:
    category: waterlevel
    units: m+MSL
    description: hourly modeled water levels for synthetic tropical cyclone tracks
    title: Synthetic TC ADCIRC StormTide
    model_ref: ADCIRC
    publication_date: 2022
    data_ref: https://www.nature.com/articles/s41558-021-01272-7"""

            # Append entry to catalog
            fcat.write(yml_str)

# ---------------------------------------------------------------------
# Tidal reanalysis dataset (single file entry)
# ---------------------------------------------------------------------
data_dir = r'Z:\Data-Expansion\users\lelise\projects\Carolinas_SFINCS\Chapter3_SyntheticTCs\EDSReanalysis_data'
file = 'EDSReanalysis_V2_1992_2022_with90yrOffset.nc'
with open(data_lib, mode="a") as fcat:
    # YAML entry for tide reanalysis dataset
    yml_str = f"""
tide_reanalysis:
  path: {os.path.join(data_dir, file)}
  data_type: GeoDataset
  driver: netcdf
  crs: 4326
  meta:
    category: waterlevel
    units: m+MSL
    description: https://github.com/RENCI/EDSReanalysis
    title: ADCIRC tide reanalysis V2 with +90yr offset
    model_ref: ADCIRC
    publication_date: 2024
    data_ref: https://renci.github.io/edsreanalysisdoc/"""

    # Append entry to catalog
    fcat.write(yml_str)
