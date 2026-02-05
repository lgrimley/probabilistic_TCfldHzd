The scripts in this folder are used for processing model output to analyze the TC track characteristics and meteorology (i.e., wind, rain)
at the basin and domain scale. A general description of the scripts and the order in which to run them is described below:

1. `meteoStats_to_netcdf.py`
This script processes the synthetic tropical cyclone (TC) rainfall and wind data to calculate storm-relevant information 
(i.e., peak rain rate/wind speed, total precipitation, etc.) using functions in `src\utils.py`. 
It generates NetCDF files combining all of these for each TC ID. These NetCDFs are used to calculate storm-intensity information
in the script `meteoStats_by_basin.py`. In this script, we also create rasters that can be used to mask the data by basis (HUC6).

2. `meteoStats_by_basin.py`
This script reads the NetCDFs and masks generated using the previous script and aggregates the meteo information to get basin and domain-wide statistics. 
The results are output to CSV files.

3. `meteo_violin_plot.py`
This script reads the CSV files generated using the script above to create violin plots of the storm-related intensity metrics.
