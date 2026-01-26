"""
Script Name: run_writeScript.py

Description:
    This script batch-generates SFINCS boundary condition inputs for a selected
    list of synthetic tropical cyclones. For each storm, it calls an external Python
    script (`write_sfincs_track_inputs.py`) that writes SFINCS forcing files, 
    manages output directories, and logs stdout to a file.

    The script: 
    - Defines a list of synthetic TC indices to process
    - Creates and manages run directories
    - Skips storms that have already completed
    - Calls the SFINCS input-generation script via subprocess
    - Applies a timeout to prevent stalled runs
    - Tracks failed or timed-out storms

Expected Outputs:
    For each successfully processed tropical cyclone (TC_####), the following
    files are created in:

        <root>/<gcm>_<ssp>_runs/runs/TC_####/sfincs_bc_inputs/

    - precip_2d.nc     : Gridded precipitation forcing for SFINCS
    - wind_2d.nc       : Gridded wind forcing for SFINCS
    - sfincs.bnd       : SFINCS water level boundary condition file (x,y locations)
    - sfincs.bzs       : SFINCS water level boundary condition file (timeseries)
    - sfincs.inp       : SFINCS configuration file

    Additionally:
    - Standard output from all subprocess calls is appended to:
        output.txt
    - Successfully completed storms may be moved or flagged under:
        <root>/<gcm>_<ssp>_runs/completed_runs/
    - Storms that fail or exceed the timeout duration are tracked in-memory
      via the `failed` list during execution.

Notes:
    - Intended for long batch runs on shared or networked filesystems
    - Output is appended to a persistent log file (output.txt)
    - Assumes downstream script handles all SFINCS-specific logic
"""

import subprocess
import pandas as pd
import os
import shutil


def execute(cmd):
    popen = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True)
    for stdout_line in iter(popen.stdout.readline, ""):
        yield stdout_line
    popen.stdout.close()
    return_code = popen.wait()
    if return_code:
        raise subprocess.CalledProcessError(return_code, cmd)


module_path = r'C:\Users\lelise\Documents\GitHub\flood_model_carolinas\syntheticTCs_cmpdfld'
root = r'Z:\Data-Expansion\users\lelise\projects\Carolinas_SFINCS\Chapter3_SyntheticTCs\03_MODEL'
gcm = 'canesm'
ssp = 'ssp585'
dd = os.path.join(root, f'{gcm}_{ssp}_runs')
if os.path.exists(dd) is False:
    os.makedirs(dd)
os.chdir(dd)

storm_list = [1252, 1583, 2723, 4803, 5510, 6034, 6040]
timeout_duration = 1200
failed = []

with open(fr"output.txt", "a") as f:
    for tc_index in storm_list:
        if os.path.exists(os.path.join(dd, 'completed_runs', f'TC_{str(tc_index).zfill(4)}')):
            print(f'{tc_index} already run in SFINCS, skipping.')
            continue
        else:
            track_dir = os.path.join(dd, 'runs', f'TC_{str(tc_index).zfill(4)}', 'sfincs_bc_inputs')
            if os.path.exists(track_dir) is False:
                os.makedirs(track_dir)
            files = os.listdir(track_dir)
            if len(files) < 6:
                print(f'Working on {tc_index}')
                try:
                    subprocess.run(
                        [
                            'python',
                            r'C:\Users\lelise\Documents\GitHub\flood_model_carolinas\syntheticTCs_cmpdfld\sfincs\model_setup\write_sfincs_track_inputs.py',
                            str(tc_index),
                            track_dir
                        ],
                        env={**os.environ, 'PYTHONPATH': module_path},
                        stdout=f,
                        text=True,
                        timeout=timeout_duration
                    )
                except subprocess.TimeoutExpired:
                    print(f"Subprocess timed out after {timeout_duration} seconds for {tc_index}")
                    failed.append(tc_index)
                except Exception as e:
                    print(f"An error occurred: {e}")
                    failed.append(tc_index)
            else:
                print(f'Already processed {tc_index}')
