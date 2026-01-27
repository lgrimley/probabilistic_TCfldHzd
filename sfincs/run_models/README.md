
## Batch running SFINCS models

Use **submit_runs.sh** to submit **run_sfincs_syntheticTCs.py** using slurm job arrays.
Each Slurm task processes one tropical cyclone and runs the three scenarios. 
**run.sh** can be used to submit the python script for a single SFINCS run.

This script executes SFINCS hydrodynamic simulations for synthetic tropical cyclones using a Slurm job array and a Singularity container. For each storm, it:
- Ranks storms by landfall maximum wind speed
- Selects a storm based on the Slurm array task ID 
- Generates scenario-specific sfincs.inp configuration files 
  - compound
  - coastal
  - runoff
- Runs SFINCS for multiple flooding scenarios 
- Logs model output and resource usage (i.e., see **python_output_file.txt**)
- Moves completed runs to a centralized archive directory




