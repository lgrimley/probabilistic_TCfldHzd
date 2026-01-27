"""
Description:
    This script applies a uniform sea-level rise (SLR) offset to an existing
    SFINCS tidal boundary condition file (.bzs). The modified water levels are
    written to a new file using the same fixed-width formatting required by SFINCS.

    The script:
    - Reads an existing tidal boundary condition file
    - Applies a specified SLR value (in meters)
    - Writes a new .bzs file with the updated water levels

Notes:
    - SLR value represents mean estimate with uncertainty noted in comments
"""

import os
import pandas as pd

# Historical tropical cyclone IDs (not used directly in this script)
tc_ids_hist = [4347, 4796]

# Future/synthetic tropical cyclone IDs (not used directly in this script)
tc_ids_fut = [4550, 671, 5547, 5719, 4133, 4970, 1831, 1545, 3214, 4823, 4762, 5701, 2041, 5902, 2818]

# Directory containing SFINCS boundary condition files
bc_dir = r'Z:\Data-Expansion\users\lelise\projects\Carolinas_SFINCS\Chapter3_SyntheticTCs\03_MODEL\sfincs_initcond_mod'

# Sea-level rise value in meters
# Mean SLR = 1.12 m, 95% confidence interval: 0.68–1.57 m
slr = 1.12  

# Read tidal water level boundary conditions
waterlevel = pd.read_csv(
    os.path.join(bc_dir, 'tides_sfincs.bzs'),
    header=None,
    sep="\s+",
    index_col=0
)

# Ensure all values are floating-point
waterlevel = waterlevel.astype(float)

# Apply sea-level rise offset
waterlevel_slr = waterlevel + slr

# Change working directory to boundary condition directory
os.chdir(bc_dir)

# Define output filename including SLR value in centimeters
filename = f'tides_sfincs_SLR{int(slr*100)}.bzs'

# Open the output file in write mode
with open(filename, 'w') as f:
    # Iterate through each row of the modified DataFrame
    for index, row in waterlevel_slr.iterrows():
        # Format the index with fixed width and one decimal place
        # Format each water level value with fixed width and two decimal places
        formatted_row = f"{index:9.1f}" + "".join([f"{val:8.2f}" for val in row]) + "\n"
        
        # Write the formatted row to the output file
        f.write(formatted_row)
