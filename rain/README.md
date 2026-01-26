
Formating TCR Rainfall for input to SFINCS


## Mapping rainfall to a static grid 
Script: `01_rain_TCR_to_staticGrid.m`
Overview: The TCR Model rainfall grid moves with the storm. This script embeds it into a very large fixed grid (Lgrid). The placement at each time step is done by matching the lower-left corner of the dynamic grid to the nearest static grid cell. This means small alignment offsets may occur if grid spacing or origin differs slightly.
Input: TCR MAT-files from Gori et al. 2022
Output: Rainfall netcdf-files at a static grid across the study area

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
