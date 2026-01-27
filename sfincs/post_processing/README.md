
## Processing SFINCS Outputs for Flood Hazard Analysis
**process_SFINCS_outputs.py** and ** process_SFINCs_outputs_SLR.py**
This script does the following:
1. Loads a base SFINCS model configuration using HydroMT
2. Identifies synthetic storms that:
   - Have landfall wind speed metadata 
   - Have valid ADCIRC storm tide outputs
3. Sorts storms by landfall maximum wind speed
4. Groups storms into fixed-size batches (default: 20 storms per group)
5. For each storm in the group:
   - Reads SFINCS outputs 
   - Extracts flood hazard results using `TCFloodHazard` class
   - Writes grouped outputs to NetCDF files

**plot_SFINCS_output_hydrographs.py** is a diagnostic script for plotting SFINCS modeled water level time series at points of interest for specified TC IDs