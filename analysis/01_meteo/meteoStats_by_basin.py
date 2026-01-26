"""
Basin-Level TC Rainfall and Wind Statistics Processing

Description:
This script processes tropical cyclone (TC) precipitation and wind data at the basin level
for the Carolinas using TCR (Tropical Cyclone Rainfall) and NCEP reanalysis wind data. 

It performs the following steps:
1. Loads a hydrometeorological data catalog and domain shapefiles.
2. Loads precomputed TC rainfall and wind datasets.
3. Applies a basin mask to extract values corresponding to individual basins.
4. Computes basin-level statistics for rainfall:
   - Cumulative precipitation (km³)
   - Mean, max, and thresholded mean rain rates (mm/hr)
5. Computes basin-level statistics for wind:
   - Maximum wind speed
   - Mean maximum wind speed
   - Mean wind speed (overall and thresholded)
   - Mean wind direction
6. Aggregates basin and domain-wide statistics and saves results to CSV files.
"""

import sys
import os
import pandas as pd
sys.path.append(r'C:\Users\lelise\Documents\GitHub\flood_model_carolinas\synthetic_tc_cmpdfld')
import hydromt
import xarray as xr

# -----------------------------
# Load hydrometeorological data catalog
# -----------------------------
data_catalog_yml = r'Z:\Data-Expansion\users\lelise\data\data_catalog_SFINCS_Carolinas.yml'
cat = hydromt.DataCatalog(data_libs=[data_catalog_yml])
# Load the modeling domain as a GeoDataFrame in WGS84 (lat/lon)
mod_domain = cat.get_geodataframe(data_like='enc_domain_HUC6_clipped').to_crs(4326)

# -----------------------------
# Process TCR rainfall by basin
# -----------------------------
os.chdir(r'Z:\Data-Expansion\users\lelise\projects\Carolinas_SFINCS\Chapter3_SyntheticTCs\02_DATA\CMIP6_585')
# Load precomputed TC precipitation dataset
ds = xr.open_dataset(r'.\rain\tc_precipitation_stats.nc').compute()
# Load basin mask raster where values correspond to mod_domain index
mask = xr.open_dataarray(r'.\rain\basin_mask_TCR.tif')
selected_tcs = ds['tc_id']  # List of TC IDs

# Initialize list to store per-basin statistics
df_list = []

# Loop through each basin to compute precipitation statistics
for i in range(len(mod_domain.index)):
    basin = mod_domain['Name'][i].replace(" ","")  # Basin name without spaces

    # Cumulative precipitation (adjusted for basin area in km³)
    df = ds['total_precip'].where(mask == i).sum(dim=['x', 'y']).to_dataframe()
    df.drop(columns='spatial_ref', inplace=True)
    basin_area = mod_domain['area_sqkm'][i]
    df = df * 1e-6 * basin_area
    df.columns = [f'{basin}_CumPrecipKM3']
    df_list.append(df)

    # Mean total precipitation per cell (mm)
    df = ds['total_precip'].where(mask == i).mean(dim=['x', 'y']).to_dataframe()
    df.drop(columns='spatial_ref', inplace=True)
    df.columns = [f'{basin}_MeanTotPrecipMM']
    df_list.append(df)

    # Max total precipitation per cell (mm)
    df = ds['total_precip'].where(mask == i).max(dim=['x', 'y']).to_dataframe()
    df.drop(columns='spatial_ref', inplace=True)
    df.columns = [f'{basin}_MaxTotPrecipMM']
    df_list.append(df)

    # Max rain rate (mm/hr)
    df = ds['max_RR'].where(mask == i).max(dim=['x', 'y']).to_dataframe()
    df.drop(columns='spatial_ref', inplace=True)
    df.columns = [f'{basin}_maxRR']
    df_list.append(df)

    # Average max rain rate (mm/hr)
    df = ds['max_RR'].where(mask == i).mean(dim=['x', 'y']).to_dataframe()
    df.drop(columns='spatial_ref', inplace=True)
    df.columns = [f'{basin}_AvgmaxRR']
    df_list.append(df)

    # Mean rain rate without threshold
    df = ds['mean_RR'].where(mask == i).mean(dim=['x', 'y']).to_dataframe()
    df.drop(columns='spatial_ref', inplace=True)
    df.columns = [f'{basin}_meanRR']
    df_list.append(df)

    # Mean rain rate with threshold >= 5 mm/hr
    df = ds['mean_RR_thresh'].where(mask == i).mean(dim=['x', 'y']).to_dataframe()
    df.drop(columns='spatial_ref', inplace=True)
    df.columns = [f'{basin}_meanRRthresh']
    df_list.append(df)

    print(f'Done processing {basin}')

# Combine all per-basin precipitation statistics
basin_precip_stats = pd.concat(objs=df_list, axis=1, ignore_index=False)
basin_precip_stats = basin_precip_stats.droplevel(1)

# -----------------------------
# Domain-wide precipitation statistics
# -----------------------------
basin = 'Domain'
df_list = []

# Repeat the same calculations for the entire domain
df = ds['total_precip'].sum(dim=['x', 'y']).to_dataframe()
df.drop(columns='spatial_ref', inplace=True)
basin_area = mod_domain['area_sqkm'][i]
df = df * 1e-6 * basin_area
df.columns = [f'{basin}_CumPrecipKM3']
df_list.append(df)

df = ds['total_precip'].mean(dim=['x', 'y']).to_dataframe()
df.drop(columns='spatial_ref', inplace=True)
df.columns = [f'{basin}_MeanTotPrecipMM']
df_list.append(df)

df = ds['total_precip'].max(dim=['x', 'y']).to_dataframe()
df.drop(columns='spatial_ref', inplace=True)
df.columns = [f'{basin}_MaxTotPrecipMM']
df_list.append(df)

df = ds['max_RR'].max(dim=['x', 'y']).to_dataframe()
df.drop(columns='spatial_ref', inplace=True)
df.columns = [f'{basin}_maxRR']
df_list.append(df)

df = ds['max_RR'].mean(dim=['x', 'y']).to_dataframe()
df.drop(columns='spatial_ref', inplace=True)
df.columns = [f'{basin}_AvgmaxRR']
df_list.append(df)

df = ds['mean_RR'].mean(dim=['x', 'y']).to_dataframe()
df.drop(columns='spatial_ref', inplace=True)
df.columns = [f'{basin}_meanRR']
df_list.append(df)

df = ds['mean_RR_thresh'].mean(dim=['x', 'y']).to_dataframe()
df.drop(columns='spatial_ref', inplace=True)
df.columns = [f'{basin}_meanRRthresh']
df_list.append(df)

rain_stats = pd.concat(objs=df_list, axis=1, ignore_index=False)
print(f'Done processing {basin}')

# Combine per-basin and domain-wide precipitation statistics
final = pd.concat(objs=[basin_precip_stats, rain_stats], axis=1, ignore_index=False).round(2)
final.to_csv(r'.\rain\basin_tc_precipitation_stats_v2.csv')

# -----------------------------
# Process wind by basin
# -----------------------------
os.chdir(r'Z:\Data-Expansion\users\lelise\projects\Carolinas_SFINCS\Chapter3_SyntheticTCs\02_DATA\NCEP_Reanalysis')
ds = xr.open_dataset(r'.\wind\tc_windspeed_stats.nc').compute()
mask = xr.open_dataarray(r'.\wind\basin_mask_wind.tif')
selected_tcs = ds['tc_id']

# Loop through each basin to compute wind statistics
df_list = []
for i in range(len(mod_domain.index)):
    basin = mod_domain['Name'][i].replace(" ","")

    # Max wind speed across basin
    df = ds['max_windspd'].where(mask == i).max(dim=['x', 'y']).to_dataframe()
    df.drop(columns=['spatial_ref'], inplace=True)
    df.columns = [f'{basin}_maxWS']
    df_list.append(df)
    print('1')

    # Mean of max wind speeds
    df = ds['max_windspd'].where(mask == i).mean(dim=['x', 'y']).to_dataframe()
    df.drop(columns=['spatial_ref'], inplace=True)
    df.columns = [f'{basin}_meanMaxWS']
    df_list.append(df)
    print('2')

    # Mean wind speed
    df = ds['mean_windspd'].where(mask == i).mean(dim=['x', 'y']).to_dataframe()
    df.drop(columns=['spatial_ref'], inplace=True)
    df.columns = [f'{basin}_meanWS']
    df_list.append(df)
    print('3')

    # Mean wind speed with threshold
    df = ds['mean_windspd_threshold'].where(mask == i).mean(dim=['x', 'y']).to_dataframe()
    df.drop(columns=['spatial_ref'], inplace=True)
    df.columns = [f'{basin}_meanWSthresh']
    df_list.append(df)
    print('4')

    # Mean wind direction (degrees)
    # Wind conventions: 0° = from north, 90° = from east, etc.
    df = ds['mean_direction_deg'].where(mask == i).mean(dim=['x', 'y']).to_dataframe()
    df.drop(columns=['spatial_ref'], inplace=True)
    df.columns = [f'{basin}_meanDirection']
    df_list.append(df)
    print('5')
    print(f'Done processing {basin}')

# Combine per-basin wind statistics
basin_windspd_stats = pd.concat(objs=df_list, axis=1, ignore_index=False)
basin_windspd_stats = basin_windspd_stats.droplevel(1)

# -----------------------------
# Domain-wide wind statistics
# -----------------------------
df_list = []
basin='Domain'

df = ds['max_windspd'].max(dim=['x', 'y']).to_dataframe()
df.columns = [f'{basin}_maxWS']
df_list.append(df)

df = ds['max_windspd'].mean(dim=['x', 'y']).to_dataframe()
df.columns = [f'{basin}_meanMaxWS']
df_list.append(df)

df = ds['mean_windspd'].mean(dim=['x', 'y']).to_dataframe()
df.columns = [f'{basin}_meanWS']
df_list.append(df)

df = ds['mean_windspd_threshold'].mean(dim=['x', 'y']).to_dataframe()
df.columns = [f'{basin}_meanWSthresh']
df_list.append(df)

df = ds['mean_direction_deg'].mean(dim=['x', 'y']).to_dataframe()
df.columns = [f'{basin}_meanDirection']
df_list.append(df)

windspd_stats = pd.concat(objs=df_list, axis=1, ignore_index=False)
print(f'Done processing {basin}')

# Combine per-basin and domain-wide wind statistics
final = pd.concat(objs=[basin_windspd_stats, windspd_stats], axis=1, ignore_index=False).round(2)
final.to_csv(r'.\wind\basin_tc_windspd_stats.csv')
