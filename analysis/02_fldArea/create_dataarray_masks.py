"""
================================================================================
Script Overview: Create Basin and Water Body Masks for SFINCS Subgrid

Description:
------------
This script generates spatial masks used in processing SFINCS outputs
at a specified subgrid resolution.

1. Loads HydroMT data catalogs and required geospatial datasets
2. Reads basin geometries (HUC6) and reprojects them to the model CRS
3. Rasterizes basin boundaries to create a basin mask aligned with the model grid
4. Rasterizes coastal water bodies and large rivers to create a water body mask
5. Writes both basin and water masks to GeoTIFF files for use in modeling

These masks are used to:
- Identify basin membership of grid cells
- Exclude coastal and river water bodies from certain analyses

Assumptions:
------------
- All referenced datasets exist and are correctly registered in the HydroMT catalog
- The elevation subgrid raster defines the target grid and CRS
- EPSG:32617 (UTM Zone 17N) is the working projection

================================================================================
"""

import sys
import os
import hydromt
import xarray as xr

# Add local repository path for custom utilities (if needed downstream)
sys.path.append(
    r'C:\Users\lelise\Documents\GitHub\flood_model_carolinas\syntheticTCs_cmpdfld'
)

# =========================
# Set working directory
# =========================
wdir = r'Z:\Data-Expansion\users\lelise\projects\Carolinas_SFINCS\Chapter3_SyntheticTCs\04_MODEL_OUTPUTS'
os.chdir(wdir)

# ==========================================
# Load HydroMT data catalogs
# ==========================================
data_catalog_yml = r'Z:\Data-Expansion\users\lelise\data\data_catalog_SFINCS_Carolinas.yml'
yml_base = r'Z:\Data-Expansion\users\lelise\data\data_catalog_BASE_Carolinas.yml'
cat = hydromt.DataCatalog(data_libs=[data_catalog_yml, yml_base])

# ==========================================
# Load ADCIRC water level dataset (for testing / reference)
# ==========================================
file = r'Z:\Data-Expansion\users\lelise\adcirc_wl_sept2024.nc'
wl = cat.get_geodataset(file)
test = xr.open_dataset(file)

# ==========================================
# Load basin geometries (HUC6) and reproject
# ==========================================
basins = cat.get_geodataframe(
    data_like=r'Z:\Data-Expansion\users\lelise\projects\Carolinas_SFINCS\Chapter2_PGW'
              r'\sfincs\03_OBS\analysis_final\downscale_test\masks\huc6_basins.shp'
)
basins = basins.to_crs(epsg=32617)

# =========================
# Define subgrid resolution
# =========================
res = 5

'''' Create basin and water body masks  '''

# ==========================================
# Define input and output file paths
# ==========================================
elevation_file = rf'..\03_MODEL_RUNS\subgrid\dep_subgrid_{res}m.tif'
water_mask = rf'.\masks\water_mask_sbgRes{res}m.tif'
basin_mask = rf'.\masks\basin_mask_sbgRes{res}m.tif'

# ==========================================
# Load elevation raster (defines target grid)
# ==========================================
dep = cat.get_rasterdataset(elevation_file)

# ==========================================
# Create basin mask raster
# ==========================================
# Rasterize HUC6 basin polygons onto the model grid
mask = dep.raster.rasterize(
    basins,
    'index',
    nodata=-128,
    all_touched=False
)

# Assign CRS and write basin mask to disk
mask.rio.write_crs('EPSG:32617', inplace=True)
mask.raster.to_raster(basin_mask, nodata=-128, dtype='int8')

# ==========================================
# Create water body mask (coastal + large rivers)
# ==========================================

# Load and clip coastal water bodies
coastal_wb = cat.get_geodataframe('carolinas_coastal_wb', geom=basins)
coastal_wb = coastal_wb.to_crs(32617)
coastal_wb_clip = coastal_wb.clip(basins)
coastal_wb_clip['mask'] = 1

# Rasterize coastal water bodies
mask1 = dep.raster.rasterize(
    coastal_wb_clip,
    "mask",
    nodata=0,
    all_touched=True
)

# Load and rasterize NHD area rivers
carolinas_nhd_area_rivers = cat.get_geodataframe(
    'carolinas_nhd_area_rivers',
    geom=basins
)
carolinas_nhd_area_rivers = carolinas_nhd_area_rivers.to_crs(32617)
carolinas_nhd_area_rivers['mask'] = 1
mask2 = dep.raster.rasterize(
    carolinas_nhd_area_rivers,
    "mask",
    nodata=0,
    all_touched=True
)

# ==========================================
# Combine water body masks and finalize
# ==========================================
# Merge coastal and river masks
mask = (mask1 + mask2).compute()

# Convert to binary mask (1 = water, 0 = land)
mask = xr.where(cond=mask > 0, x=1, y=0)
mask = mask.astype('int8')

# Assign CRS and write water mask to disk
mask.rio.write_crs('EPSG:32617', inplace=True)
mask.raster.to_raster(water_mask, nodata=0)
