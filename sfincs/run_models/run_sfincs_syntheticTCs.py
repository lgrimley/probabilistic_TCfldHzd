"""
===============================================================================
Description:
------------
This script runs SFINCS hydrodynamic simulations for synthetic tropical cyclones
using a Slurm job array and a Singularity container. Each Slurm array task
processes one storm, selected based on landfall maximum wind speed ranking.

For each storm, the script:
- Selects a TC using the Slurm task ID
- Generates scenario-specific SFINCS configuration files
- Runs SFINCS inside a Singularity container for multiple scenarios
- Captures stdout/stderr to log files
- Moves completed runs to a finalized directory

Scenarios Run:
--------------
- runoff    : River discharge + precipitation
- coastal   : Storm surge + tides + wind
- compound  : Combined runoff and coastal forcing

Execution:
----------
Designed to be executed via Slurm job arrays.

Author: Lauren Grimley
===============================================================================
"""

import sys
import subprocess
import os
import traceback
import psutil
import pandas as pd
import shutil

# ------------------------------------------------------------------------------
# Read Slurm job array task ID (passed as command-line argument)
# ------------------------------------------------------------------------------

task_id = int(sys.argv[1])

# ------------------------------------------------------------------------------
# Load tropical cyclone metadata and rank storms by landfall Vmax
# ------------------------------------------------------------------------------

# CSV containing TC IDs and landfall maximum wind speed (vstore100)
tc_ids = pd.read_csv(
    r'/projects/sfincs/syntheticTCs_cmpdfld/SFINCS_model/ncep_TCs_vstore100.csv',
    index_col=0
)

# Rank storms from strongest to weakest landfall wind speed
tcs_ranked = tc_ids.sort_values(by='vstore100', ascending=False)

# Ordered list of TC IDs used to select storm by task ID
ordered_list_track = tcs_ranked['tc_id'].tolist()

# ------------------------------------------------------------------------------
# Utility functions
# ------------------------------------------------------------------------------

def get_usage(proc):
    """
    Return CPU and memory usage for a running process.

    Parameters
    ----------
    proc : psutil.Process
        Process object to monitor

    Returns
    -------
    cpu_percent : float
        CPU usage percentage
    mem_usage : float
        Memory usage in MB
    """
    mem_info = proc.memory_info()
    cpu_percent = proc.cpu_percent(interval=1)
    mem_usage = mem_info.rss / (1024 * 1024)
    return cpu_percent, mem_usage


def run_sfincs_models(scenario_dir, input_dir, base_model_dir, task_id):
    """
    Run SFINCS inside a Singularity container for a single scenario.

    Parameters
    ----------
    scenario_dir : str
        Directory where SFINCS is executed
    input_dir : str
        Directory containing boundary condition inputs
    base_model_dir : str
        Directory containing base model static files
    task_id : int
        Slurm job array task ID
    """
    # Prepare bind mounts for the Singularity container
    bind_str = f'-B {scenario_dir},{input_dir},{base_model_dir}'
    singularity_command = [
        "singularity", "run",
        bind_str,
        "/projects/sfincs/docker/sfincs.sif",
    ]
    print(singularity_command)

    try:
        # Change working directory to the scenario folder
        os.chdir(scenario_dir)

        # Open log file for this Slurm task
        with open(f'sfincs_log_{task_id}.log', 'w') as log_file:

            # Log available CPU resources
            logical_cores = psutil.cpu_count(logical=True)
            physical_cores = psutil.cpu_count(logical=False)
            log_file.write(f"Logical CPU cores: {logical_cores}\n")
            log_file.write(f"Physical CPU cores: {physical_cores}\n")
            log_file.flush()

            # Launch the Singularity container
            process = subprocess.Popen(
                singularity_command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )

            # Stream stdout to log file
            for line in process.stdout:
                log_file.write(line)
                log_file.flush()

            # Stream stderr to log file
            for line in process.stderr:
                log_file.write(line)
                log_file.flush()

            # Wait for SFINCS to finish
            process.wait()

            # Report success or failure
            if process.returncode == 0:
                log_file.write(f"Successfully completed SFINCS run: TaskID {task_id}")
                log_file.flush()
            else:
                log_file.write(
                    f"SFINCS {task_id} failed with exit code {process.returncode}"
                )
                sys.exit(1)

    except Exception as e:
        # Print detailed error diagnostics
        print(f"Error running {scenario_dir}: {e}")
        exc_type, exc_value, exc_traceback = sys.exc_info()
        tb = traceback.extract_tb(exc_traceback)[-1]
        line_number = tb.lineno

        print(f"Error Type: {exc_type.__name__}")
        print(f"Error Message: {exc_value}")
        print(f"Occurred at line: {line_number}")
        return None
        pass


def modify_config_file(template_file, output_file, updates_dict):
    """
    Modify a SFINCS configuration file (sfincs.inp).

    - Updates values on the right-hand side of '='
    - Adds new parameters if not present
    - Skips writing parameters with empty values

    Parameters
    ----------
    template_file : str
        Path to the template sfincs.inp file
    output_file : str
        Path where the modified file will be written
    updates_dict : dict
        Dictionary of parameter updates
    """
    try:
        # Read template file
        with open(template_file, 'r') as infile:
            lines = infile.readlines()

        existing_keys = set()
        updated_lines = []

        # Update existing parameters
        for line in lines:
            parts = line.split("=")

            if len(parts) == 2:
                key = parts[0].strip()
                value = parts[1].strip()

                if key in updates_dict:
                    value = updates_dict[key]

                existing_keys.add(key)

                if value:
                    k = "{:<20}".format(key)
                    updated_lines.append(f"{k} = {value}\n")
            else:
                updated_lines.append(line)

        # Add new parameters not present in template
        for key, value in updates_dict.items():
            if key not in existing_keys and value:
                k = "{:<20}".format(key)
                updated_lines.append(f"{k} = {value}\n")

        # Write updated configuration file
        with open(output_file, 'w') as outfile:
            outfile.writelines(updated_lines)

        print(f"File modified and saved to {output_file}")

    except Exception as e:
        print(f"An error occurred: {e}")

# ------------------------------------------------------------------------------
# Main execution block
# ------------------------------------------------------------------------------

with open(f'python_output_file.txt', 'a') as python_output_file:

    # Select storm corresponding to this Slurm task
    track_id = ordered_list_track[task_id]

    # Log storm assignment
    line_entry = (
        f'Running {track_id} with Slurm Job Array Task ID {task_id}\n'
    )
    python_output_file.write(line_entry)
    python_output_file.flush()

    # Define key directories
    base_model_dir = (
        r'/projects/sfincs/syntheticTCs_cmpdfld/SFINCS_model/base_model'
    )
    track_dir = (
        rf'/projects/sfincs/syntheticTCs_cmpdfld/SFINCS_model/runs/TC_{str(track_id).zfill(4)}'
    )

    # Exit if storm directory does not exist
    if not os.path.exists(track_dir):
        print(f"The folder does not exist for {track_id}. Exiting the python script.")
        sys.exit()

    # Verify boundary condition inputs exist
    bc_inputs_dir = os.path.join(track_dir, 'sfincs_bc_inputs')
    if len(os.listdir(bc_inputs_dir)) < 2:
        print(
            f"Boundary condition inputs not created for {track_id}. Exiting the python script."
        )
        sys.exit()

    # Define flooding scenarios
    scenarios = ['runoff', 'coastal', 'compound']

    # Base SFINCS configuration updates (applied to all scenarios)
    basemodel_config_updates = {
        "depfile": os.path.join(base_model_dir, 'sfincs.dep'),
        "mskfile": os.path.join(base_model_dir, 'sfincs.msk'),
        "indexfile": os.path.join(base_model_dir, 'sfincs.ind'),
        "sbgfile": os.path.join(base_model_dir, 'sfincs_subgrid.nc'),
        "obsfile": os.path.join(base_model_dir, 'sfincs_obsLocs_gauge_riv_v2.obs'),
        "smaxfile": os.path.join(base_model_dir, 'sfincs.smax'),
        "sefffile": os.path.join(base_model_dir, 'sfincs.seff'),
        "ksfile": os.path.join(base_model_dir, 'sfincs.ks'),
        "weirfile": os.path.join(base_model_dir, 'sfincs.weir'),
        "rstfile": os.path.join(base_model_dir, 'sfincs.19990119.000000.rst'),
        "tspinup": '60',
        "dtout": '432000',
        "dtrstout": '',
        "twet_threshold": '0.1',
        "storetwet": '1',
        "storevelmax": '1',
        "storefluxmax": '1',
        'zsini': '',
        'dthisout': '1800',
        'dtmaxout': '99999999',
    }

    # Loop over scenarios and execute SFINCS
    for scenario in scenarios:

        scenario_dir = os.path.join(track_dir, scenario)
        if os.path.exists(scenario_dir) is False:
            os.makedirs(scenario_dir)

        # Template and output configuration files
        template_inp = os.path.join(bc_inputs_dir, 'sfincs.inp')
        output_inp = os.path.join(scenario_dir, 'sfincs.inp')

        # Scenario-specific configuration updates
        if scenario == 'runoff':
            config_updates = {
                "bndfile": os.path.join(bc_inputs_dir, 'tides_sfincs.bnd'),
                "bzsfile": os.path.join(bc_inputs_dir, 'tides_sfincs.bzs'),
                "srcfile": os.path.join(base_model_dir, 'sfincs.src'),
                "disfile": os.path.join(base_model_dir, 'sfincs.dis'),
                "netamprfile": os.path.join(bc_inputs_dir, 'precip_2d.nc'),
                "netamuamvfile": ''
            }
        elif scenario == 'coastal':
            config_updates = {
                "bndfile": os.path.join(bc_inputs_dir, 'sfincs.bnd'),
                "bzsfile": os.path.join(bc_inputs_dir, 'sfincs.bzs'),
                "srcfile": '',
                "disfile": '',
                "netamprfile": '',
                "netamuamvfile": os.path.join(bc_inputs_dir, 'wind_2d.nc')
            }
        else:
            config_updates = {
                "bndfile": os.path.join(bc_inputs_dir, 'sfincs.bnd'),
                "bzsfile": os.path.join(bc_inputs_dir, 'sfincs.bzs'),
                "srcfile": os.path.join(base_model_dir, 'sfincs.src'),
                "disfile": os.path.join(base_model_dir, 'sfincs.dis'),
                "netamprfile": os.path.join(bc_inputs_dir, 'precip_2d.nc'),
                "netamuamvfile": os.path.join(bc_inputs_dir, 'wind_2d.nc')
            }

        # Write updated sfincs.inp file
        modify_config_file(
            template_file=template_inp,
            output_file=output_inp,
            updates_dict={**basemodel_config_updates, **config_updates}
        )

        print(f'SFINCS config file updated for {track_id} - {scenario}')

        # Log execution details
        print(
            f"Running {scenario} for track {track_id} "
            f"(Slurm task {os.getenv('SLURM_ARRAY_TASK_ID')}) in\n {scenario_dir}"
        )

        # Run SFINCS for the scenario
        run_sfincs_models(scenario_dir, bc_inputs_dir, base_model_dir, task_id)

    # Finalize storm processing
    print(
        f"All scenarios for storm {track_id} "
        f"(Slurm task {os.getenv('SLURM_ARRAY_TASK_ID')}) completed."
    )

    # Move completed storm runs to archive directory
    destination_folder = (
        rf'/projects/sfincs/syntheticTCs_cmpdfld/SFINCS_model/completed_runs/TC_{str(track_id).zfill(4)}'
    )
    shutil.move(track_dir, destination_folder)

    python_output_file.write(
        f'Completed run moved to: {destination_folder}\n'
    )
    python_output_file.flush()
