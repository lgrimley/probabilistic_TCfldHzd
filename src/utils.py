"""
TC and Flood Analysis Utilities

This module provides a set of tools for processing tropical cyclone (TC) track data,
precipitation, wind, tidal, and flood impact datasets. It includes functions to:

- Extract and process TC tracks from raw data files
- Convert track points to GeoDataFrames and LineStrings
- Select TCs by intersection with a study area
- Calculate TC precipitation and wind statistics, including threshold-based metrics
- Process station-based tidal data for high tide timing
- Calculate flood extents and classify inundation by driver (coastal/runoff/compound)
- Identify landfall locations relative to gauge networks

Dependencies:
- pandas, numpy, xarray, geopandas, shapely, matplotlib, cartopy, pathlib, os, datetime

"""

from shapely.geometry import LineString
import datetime
import os
import xarray as xr
import numpy as np
import geopandas as gpd
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from pathlib import Path


# ------------------------
# Processing TC tracks
# ------------------------

def get_track_datetime(df):
    """
    Converts separate year/month/day/hour columns to a single datetime column.
    Non-convertible entries are replaced with 0.
    
    Args:
        df (pd.DataFrame): DataFrame containing columns 'year100', 'month100', 'day100', 'hour100'.
        
    Returns:
        pd.DataFrame: Original DataFrame with a new 'datetime' column.
    """
    datetime_list = []
    for i in range(len(df)):
        year = df['year100'][1]  # Always take row 1 for the year
        month = df['month100'][i]
        day = df['day100'][i]
        hour = df['hour100'][i]
        try:
            x = pd.to_datetime(datetime.datetime(year, month, day, hour))
        except:
            x = 0  # fallback if date creation fails
        datetime_list.append(x)
    df['datetime'] = datetime_list
    return df


def get_track_info_in_df(tc_id, tc_tracks):
    """
    Extracts TC track information for a given TC ID from a MATLAB-like structure.
    
    Args:
        tc_id (int): Tropical cyclone ID.
        tc_tracks (dict): Dictionary containing TC track arrays.
        
    Returns:
        pd.DataFrame: Track information with datetime column added.
    """
    df = pd.DataFrame()
    py_index = tc_id - 1  # Python index correction
    for v in tc_tracks.keys():
        if v in ['__header__', '__version__', '__globals__', 'freq']:
            continue
        else:
            if len(tc_tracks[v][py_index]) == 1:
                df[v] = tc_tracks[v][py_index].item()
            else:
                df[v] = tc_tracks[v][py_index, :]

    df = get_track_datetime(df=df)
    return df


def track_as_gdf(df):
    """
    Converts a TC track DataFrame into a GeoDataFrame with Point geometries.
    
    Args:
        df (pd.DataFrame): DataFrame containing 'lon100' and 'lat100'.
    
    Returns:
        gpd.GeoDataFrame: GeoDataFrame with Point geometry.
    """
    df = df[df['datetime'] != 0]  # Remove invalid datetime entries
    track_gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(x=df['lon100'], y=df['lat100'], crs=4326))
    return track_gdf


def track_points_to_linestring(tc_tracks, output_shpfile=None):
    """
    Converts TC track points into LineString geometries and optionally saves as a shapefile.
    
    Args:
        tc_tracks (dict): TC track dictionary containing lat/lon arrays.
        output_shpfile (str, optional): Filepath to save output shapefile.
    
    Returns:
        gpd.GeoDataFrame: TC tracks as LineString geometries with 'tc_id' column.
    """
    lat = tc_tracks['lat100']
    lon = tc_tracks['lon100']
    tc_ids = np.arange(start=1, stop=len(lat) + 1, step=1)

    tc_tracks_gdf = pd.DataFrame()
    for tc_id in tc_ids:
        py_index_loc = tc_id - 1
        stm1_lat = lat[py_index_loc, :]
        stm1_lon = lon[py_index_loc, :]

        # Remove points with zero values
        stm1_lat = stm1_lat[stm1_lat != 0]
        stm1_lon = stm1_lon[stm1_lon != 0]

        # Convert points to GeoDataFrame
        track = gpd.GeoDataFrame(geometry=gpd.points_from_xy(x=stm1_lon, y=stm1_lat, crs=4326))
        track['tc_id'] = tc_id

        # Convert points to a LineString per TC
        track_line = track.groupby('tc_id')['geometry'].apply(lambda x: LineString(x.tolist()))
        track_line = gpd.GeoDataFrame(track_line, geometry='geometry', crs=4326)

        tc_tracks_gdf = pd.concat(objs=[tc_tracks_gdf, track_line], axis=0)

    tc_tracks_gdf['tc_id'] = tc_tracks_gdf.index
    if output_shpfile is not None:
        tc_tracks_gdf.to_file(output_shpfile)

    return tc_tracks_gdf


def select_TCs_by_AOI_intersect(aoi, tc_tracks_gdf, buffer_km=0):
    """
    Selects TC tracks that intersect with an area of interest (AOI) polygon.
    
    Args:
        aoi (gpd.GeoDataFrame or GeoSeries): Area of interest polygon.
        tc_tracks_gdf (gpd.GeoDataFrame): TC tracks as LineStrings.
        buffer_km (float): Buffer around AOI in km.
    
    Returns:
        tuple: Buffered AOI and TC tracks GeoDataFrame with intersection flag.
    """
    # Buffer AOI polygon
    aoi_buff = aoi.buffer(distance=buffer_km * 1000).to_crs(crs=4326)

    # Check for intersection of each TC track with buffered AOI
    select_tc = tc_tracks_gdf.geometry.apply(lambda g: aoi_buff.intersects(g))
    values, counts = np.unique(select_tc, return_counts=True)
    print(f'{counts[1]} of {len(tc_tracks_gdf)} TC tracks intersect with the study area')

    tc_tracks_gdf['aoi_sel'] = select_tc
    return aoi_buff, tc_tracks_gdf


# ------------------------
# TC Precipitation Statistics
# ------------------------

def TCR_precip_stats2netcdf(tc_ids: list, inputdir: Path, outputdir: Path,
                            rr_threshold: int=5) -> xr.Dataset:
    """
    Calculates cumulative, max, mean, and thresholded precipitation for selected TCs and saves to NetCDF.
    
    Args:
        tc_ids (list): List of TC IDs.
        inputdir (Path): Directory containing input NetCDF files.
        outputdir (Path): Directory to save output NetCDF.
        rr_threshold (int): Rainfall threshold for calculating thresholded mean.
    
    Returns:
        xr.Dataset: Dataset containing precipitation statistics for each TC.
    """
    cumP_list = []
    maxRR_list = []
    meanRR_list = []
    meanRR_thresh_list = []
    counter = 0
    for tc_id in tc_ids:
        filename = f'{str(tc_id).zfill(4)}.nc'
        d = xr.open_dataset(os.path.join(inputdir, filename))

        # Cumulative precipitation
        precip_sum = d.sum(dim='time')
        cumP_list.append(precip_sum)

        # Max rain rate across grid
        max_rain_rate = d.max(dim='time')
        maxRR_list.append(max_rain_rate)

        # Mean rain rate
        mean_rain_rate = d.mean(dim='time')
        meanRR_list.append(mean_rain_rate)

        # Mean rain rate above threshold
        mask = d > rr_threshold
        mean_rain_rate_thresh = d.where(mask).mean(dim='time')
        meanRR_thresh_list.append(mean_rain_rate_thresh)

        print(f'Processed {counter} out of {len(tc_ids)}')
        counter += 1

    # Concatenate along TC dimension and rename variables
    ds_sum = xr.concat(cumP_list, dim='tc_id')
    ds_sum['tc_id'] = xr.IndexVariable('tc_id', tc_ids)
    ds_sum = ds_sum.rename({'precip':'total_precip'})

    ds_max = xr.concat(maxRR_list, dim='tc_id')
    ds_max['tc_id'] = xr.IndexVariable('tc_id', tc_ids)
    ds_max = ds_max.rename({'precip':'max_RR'})

    ds_mean = xr.concat(meanRR_list, dim='tc_id')
    ds_mean['tc_id'] = xr.IndexVariable('tc_id', tc_ids)
    ds_mean = ds_mean.rename({'precip':'mean_RR'})

    ds_mean_thresh = xr.concat(meanRR_thresh_list, dim='tc_id')
    ds_mean_thresh['tc_id'] = xr.IndexVariable('tc_id', tc_ids)
    ds_mean_thresh = ds_mean_thresh.rename({'precip':'mean_RR_thresh'})

    ds = xr.merge([ds_sum, ds_max, ds_mean, ds_mean_thresh])
    ds.assign_attrs(RR_threshold=rr_threshold)
    outfile = os.path.join(outputdir, 'tc_precipitation_stats.nc')
    ds.to_netcdf(outfile)
    print(f'Created {outfile}')

    return ds


# ------------------------
# TC Wind Statistics
# ------------------------

def TC_windspd_stats2netcdf(tc_ids: list, inputdir: Path,
                            outputdir: Path) -> xr.Dataset:
    """
    Calculates max, mean, thresholded wind speed and mean wind direction for selected TCs.
    
    Args:
        tc_ids (list): List of TC IDs.
        inputdir (Path): Directory containing wind NetCDF files.
        outputdir (Path): Directory to save output NetCDF.
    
    Returns:
        xr.Dataset: Dataset containing wind speed and direction statistics for each TC.
    """
    vmax_list = []
    vmean_list = []
    wnd_direction_list = []
    mean_wndspd_thresh_list = []
    counter = 0
    for tc_id in tc_ids:
        filename = f'{str(tc_id).zfill(4)}.nc'
        d = xr.open_dataset(os.path.join(inputdir, filename))

        # Calculate wind direction (0-360 deg)
        wind_direction_rad = np.arctan2(d['wind10_v'], d['wind10_u'])
        wind_direction_deg = np.degrees(wind_direction_rad)
        normalized_direction = (wind_direction_deg + 360) % 360
        wnd_direction = normalized_direction.mean(dim='time')
        wnd_direction_list.append(wnd_direction)

        # Max wind speed
        max_wndspd = d['wind_speed'].max(dim='time')
        vmax_list.append(max_wndspd)

        # Mean wind speed
        mean_wndspd = d['wind_speed'].mean(dim='time')
        vmean_list.append(mean_wndspd)

        # Mean wind speed above threshold
        mask = d['wind_speed'] > 5
        mean_wndspd_thresh = d['wind_speed'].where(mask).mean(dim='time')
        mean_wndspd_thresh_list.append(mean_wndspd_thresh)

        print(f'Processed {counter} out of {len(tc_ids)}')
        counter += 1

    # Concatenate and merge all datasets
    ds1 = xr.concat(vmax_list, dim='tc_id').to_dataset().rename({'wind_speed':'max_windspd'})
    ds1['tc_id'] = xr.IndexVariable('tc_id', tc_ids)

    ds2 = xr.concat(vmean_list, dim='tc_id').to_dataset().rename({'wind_speed':'mean_windspd'})
    ds2['tc_id'] = xr.IndexVariable('tc_id', tc_ids)

    ds4 = xr.concat(mean_wndspd_thresh_list, dim='tc_id').to_dataset().rename({'wind_speed':'mean_windspd_threshold'})
    ds4['tc_id'] = xr.IndexVariable('tc_id', tc_ids)

    ds3 = xr.concat(wnd_direction_list, dim='tc_id').to_dataset().rename({'wind10_v':'mean_direction_deg'})
    ds3['tc_id'] = xr.IndexVariable('tc_id', tc_ids)

    ds = xr.merge([ds1, ds2, ds4, ds3])
    outfile = os.path.join(outputdir, 'tc_windspeed_stats.nc')
    ds.to_netcdf(outfile)
    print(f'Created {outfile}')

    return ds


# ------------------------
# Miscellaneous TC / Flood Utilities
# ------------------------

def find_nearest_non_nan_index(df, col_name, index):
    """
    Finds the nearest non-NaN index in a DataFrame column.
    
    Args:
        df (pd.DataFrame): Input DataFrame.
        col_name (str): Column to search.
        index (int): Reference index.
    
    Returns:
        int or None: Index of nearest non-NaN value or None if not found.
    """
    non_nan_indices = df[col_name].dropna().index
    if non_nan_indices.empty:
        return None
    else:
        return non_nan_indices[np.abs(non_nan_indices - index).argmin()]


def calculate_station_highTide_time(da: xr.Dataset,
                                    tref: pd.DatetimeIndex,
                                    station_names: list = None) -> pd.DataFrame:
    """
    Determines the high tide time for each station in the dataset within a 12-hour window.
    
    Args:
        da (xr.Dataset): Tidal water level dataset.
        tref (pd.DatetimeIndex): Reference start time.
        station_names (list): List of stations to process.
    
    Returns:
        pd.DataFrame: High tide times per station.
    """
    if station_names is None:
        station_names = da.index.values

    highTide_times = []
    sta_keep = []
    for station in station_names:
        data = da.sel(index=station).sel(time=slice(tref, tref + pd.to_timedelta('12h')))
        df = data.to_dataframe()
        df['delta_waterlevel'] = df['waterlevel'].diff()
        increasing_trend = df[df['delta_waterlevel'] > 0]
        try:
            if len(increasing_trend) > 0:
                ind = np.argmax(increasing_trend['waterlevel'])
                station_highTide_time = increasing_trend.index[ind]
            else:
                ind = np.argmax(df['waterlevel'])
                station_highTide_time = df.index[ind]
            highTide_times.append(station_highTide_time)
            sta_keep.append(station)
        except:
            print(f"Issue with {station}")
            pass

    highTide_df = pd.DataFrame()
    highTide_df['station'] = sta_keep
    highTide_df['highTide_time'] = highTide_times
    highTide_df.set_index(keys='station', drop=True, inplace=True)
    return highTide_df

def process_tmax_in_hours(tmax: xr.DataArray, time_min_hr: int = 0) -> xr.DataArray:
    """
    Converts the 'tmax' DataArray (time of maximum inundation) to hours,
    masks invalid values, and applies a minimum threshold for interest.
    
    Args:
        tmax (xr.DataArray): Maximum inundation time array (likely in ns or timedelta).
        time_min_hr (int): Minimum number of hours to consider.
    
    Returns:
        xr.DataArray: Processed tmax in hours with NaNs where below threshold or invalid.
    """
    # Identify invalid (NaT) values
    mask = np.isnan(tmax.values)

    # Convert to float for calculation
    tmax = tmax.astype(float)

    # Mask out invalid values
    tmax = tmax.where(~mask, np.nan)

    # Convert from nanoseconds to hours
    tmax = tmax / (3.6 * 10 ** 12)

    # Apply minimum threshold
    tmax = xr.where(tmax >= time_min_hr, tmax, np.nan)

    return tmax


def calculate_flooded_area_by_process(da: xr.DataArray, tc_index: int=0):
    """
    Calculates flooded area per process class from a process-coded DataArray.
    
    Args:
        da (xr.DataArray): DataArray with process codes (e.g., 0=NoFlood, 1=Coastal, etc.).
        tc_index (int): Index of the TC (used as DataFrame index).
    
    Returns:
        pd.DataFrame: Flooded area per process in square km, total flooded, compound area.
    """
    unique_codes, cell_counts = np.unique(da.data, return_counts=True)
    fld_area1 = cell_counts.copy()
    res = 200  # grid resolution in meters
    fld_area1 = fld_area1 * (res * res) / (1000 ** 2)  # convert to km^2

    # Cleanup and align columns
    fld_area1 = pd.DataFrame(fld_area1).T
    fld_area1.columns = unique_codes
    expected = pd.DataFrame(columns=[0.0, 1.0, 2.0, 3.0, 4.0, np.NAN])
    fld_area = pd.concat(objs=[expected, fld_area1], axis=0, ignore_index=True)

    # Rename columns for clarity
    fld_area.columns = ['NoFlood', 'Coastal', 'Coastal-Compound', 'Runoff', 'Runoff-Compound', 'Inactive']
    fld_area['Total_Area'] = fld_area.sum(axis=1)
    fld_area['Compound'] = fld_area['Coastal-Compound'] + fld_area['Runoff-Compound']
    fld_area['Total_Flooded'] = fld_area[['Coastal', 'Runoff', 'Compound']].sum(axis=1)
    fld_area['tc_index'] = tc_index
    fld_area.set_index('tc_index', inplace=True, drop=True)

    return fld_area


def get_landfall_info(tc_df, gage_locs, clip_gdf):
    """
    Determines the closest landfall location of a TC track relative to gauge points.
    
    Args:
        tc_df (pd.DataFrame): TC track points with 'lon100' and 'lat100'.
        gage_locs (pd.DataFrame): Gauge locations with 'x', 'y', and 'gage_id'.
        clip_gdf (gpd.GeoDataFrame): Polygon to clip the TC track (study area).
    
    Returns:
        gpd.GeoDataFrame: Single TC point of landfall with distance to nearest gauge.
    """
    # Convert TC track and gauge points to GeoDataFrames
    tc_gdf = gpd.GeoDataFrame(tc_df, geometry=gpd.points_from_xy(x=tc_df['lon100'], y=tc_df['lat100'], crs=4326))
    tc_gdf = tc_gdf.clip(clip_gdf)
    gage_locs_gdf = gpd.GeoDataFrame(gage_locs,
                                     geometry=gpd.points_from_xy(x=gage_locs['x'], y=gage_locs['y'], crs=4326))

    # Project to a metric CRS for accurate distance calculation
    tc_gdf = tc_gdf.to_crs(32617)
    gage_locs_gdf = gage_locs_gdf.to_crs(32617)

    min_distances = []

    # Loop through TC points to find closest gauge
    for i, point1 in tc_gdf.iterrows():
        distances = gage_locs_gdf.geometry.apply(lambda point2: point1.geometry.distance(point2))
        min_distance = np.round(distances.min(), 5)
        closest_point = gage_locs_gdf.loc[distances.idxmin()]
        min_distances.append((i, min_distance, closest_point['gage_id']))

    # Create a DataFrame with results
    result_df = gpd.GeoDataFrame(min_distances, columns=["Index_in_df1", "Min_Distance", "Closest_Point_in_df2"])

    # Select the point with minimum distance as landfall
    lf_idx = result_df.loc[result_df['Min_Distance'].idxmin()]
    landfall_info = tc_gdf[tc_gdf.index == lf_idx['Index_in_df1']]
    landfall_info['Min_Distance'] = lf_idx['Min_Distance']
    landfall_info['Closest_Point_in_df2'] = lf_idx['Closest_Point_in_df2']

    return landfall_info


def classify_zsmax_by_process(da: xr.DataArray, da_diff: xr.DataArray, hmin: float=0.1) -> xr.Dataset:
    """
    Classifies maximum water level (zsmax) by the driver process: coastal, runoff, or compound.
    
    Args:
        da (xr.DataArray): Water level DataArray with scenario dimension ['coastal','runoff','compound'].
        da_diff (xr.DataArray): Difference between compound and max single driver.
        hmin (float): Minimum water level threshold to assign compound effect.
    
    Returns:
        xr.Dataset: Dataset with classified zsmax per process.
    """
    compound_mask = da_diff > hmin
    coastal_mask = da.sel(scenario='coastal').fillna(0) > da.sel(scenario=['runoff']).fillna(0).max('scenario')
    runoff_mask = da.sel(scenario='runoff').fillna(0) > da.sel(scenario=['coastal']).fillna(0).max('scenario')
    assert ~np.logical_and(runoff_mask, coastal_mask).any()  # Ensure no overlap

    # Apply classification rules
    da_classified = (xr.where(coastal_mask, x=compound_mask + 1, y=0)
                     + xr.where(runoff_mask, x=compound_mask + 3, y=0)).compute()
    da_classified.name = 'zsmax_classified'
    da_classified = da_classified.assign_attrs(hmin=hmin,
                                               no_class=0, coast_class=1, coast_compound_class=2,
                                               runoff_class=3, runoff_compound_class=4)
    da_classified = da_classified.astype(int)

    ds = xr.Dataset({'zsmax_diff': da_diff,
                     'zsmax_attr': da_classified})
    return ds


def calc_diff_in_zsmax_compound_minus_max_individual(da: xr.DataArray) -> xr.DataArray:
    """
    Computes the difference between the compound scenario zsmax and the maximum of individual drivers.
    
    Args:
        da (xr.DataArray): Water level DataArray with scenario dimension ['coastal','runoff','compound'].
    
    Returns:
        xr.DataArray: Difference array (compound - max individual driver) per grid cell.
    """
    da_single_max = da.sel(scenario=['runoff', 'coastal']).max('scenario')
    da_diff = (da.sel(scenario='compound') - da_single_max).compute()
    da_diff.name = 'zsmax_diff'
    da_diff.attrs = da.attrs
    da_diff.assign_attrs(units='m')
    da_diff.assign_attrs(description='diff in waterlevel compound minus max. single driver')

    return da_diff
