### Coastal Water Levels

The scripts in this folder are used for processing coastal water level datasets including:
1. ADCIRC modeled storm tides (Gori et al. 2022) as MAT-files
2. RENCI ADCIRC Reanalysis V2 as a CSV-file.
- This data was downloaded from https://renci.github.io/edsreanalysisdoc/. We downloaded the reanalysis at the same x,y locations as the storm tide gages.

Scripts and their usage:

1. `stormTide_adcirc_wl_output_locs.ipynb` is used to plot the spatial location of the ADCIRC storm tide gages available in the TC Dataset
2. `stormTide_plot_hydrograph.py` is used for plotting the hydrographs for the synthetic TCs at each ADCIRC storm tide gage 
3. `stormTide_process_adcircWL.py` is used for writing the water levels (stored in MAT-files) to netcdf for intput to SFINCS
4. `EDSReanlysis_csv2netcdf.py` is used to convert the RENCI ADCIRC Reanalysis V2 data stored as a CSV to a netcdf.
