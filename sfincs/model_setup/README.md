**The scripts in this folder are used to build and update the SFINCS base model.***

## Build/update the base SFINCS model:
- `build_sfincs_ENC_200m_sbg5m_Ch3.py` build model files
- `run_build_sfincs_py.sh` bash script to submit the python script to HPC
- `add_SLR_to_stormtide.py` adds a given sea level rise amount (m) to an existing SFINCS water level file (.bnd and .bzs)
- `create_obsfile.py` adding x,y locations to the sfincs.obs file that SFINCS writes time series data in the sfincs_his.nc output
  
## Write SFINCS input files for individual TCs
- `create_BC_dataCatalog.py` Creating a data catalog and model inputs for the synthetic TCs
- `write_sfincs_track_inputs.py` This script generates SFINCS boundary condition files for a single synthetic TC track. It run from the command line with arguments specifying the TC index and output directory.
