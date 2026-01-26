
## Formating TCR Rainfall for input to SFINCS

1. `01_rain_TCR_to_staticGrid.m`
- Overview: The TCR Model rainfall grid moves with the storm. This script embeds it into a very large fixed grid (Lgrid). The placement at each time step is done by matching the lower-left corner of the dynamic grid to the nearest static grid cell. This means small alignment offsets may occur if grid spacing or origin differs slightly. (see pseudo-code below)
- Input: 2-hr TCR MAT-files from Gori et al. 2022
- Output: 2-hr TCR rainfall NetCDF files on a static grid

2. `02_rain_process_grids.py`
- Overview: Converts gridded TCR rainfall NetCDF files from 2-hr resolution to hourly resolution and aligns rainfall time coordinates with storm track datetimes. Outputs are written as hourly NetCDF files suitable for use as rainfall forcing in SFINCS.
- Input: 2-hr TCR rainfall NetCDF files on a static grid
- Output: 1-hr TCR rainfall NetCDF files on a static grid w/ datetime assigned

3. `diagnosticCheck_plotTCRoutput.m` plotting TCR rainfall

4. `rainfallStats.py`
- Overview: This script computes storm-level rainfall summary statistics and total precipitation grids from gridded TCR rainfall NetCDF files. The outputs are organized by storm selection category and include both tabular statistics and spatial NetCDF products. This information is used in downstream analysis to understand how TC characteristics are changing and relate to flooding.
- Input: 1-hr TCR rainfall NetCDF files on a static grid w/ datetime assigned
- Output: NetCDF files of rainfall metrics

<!-- Pseudo-code -->

```text
FOR each GCM
    Read list of selected storms
    Create output directory if needed

    Define static latitude/longitude grid
    Define bounding box limits

    FOR each storm
        IF output NetCDF already exists
            Skip storm
        END

        Load TCR rainfall MAT file
        Convert longitude to −180 to 180

        Initialize large static rainfall grid

        FOR each timestep
            Extract dynamic TCR grid
            IF rainfall is all zero
                Stop processing this storm
            END

            Locate dynamic grid position in static grid
            Insert rainfall into static grid
        END

        Clip static grid to bounding box
        Write clipped rainfall to NetCDF
    END
END
