# Synthetic TC Wind Fields

This codes in this folder are used for generating wind fields for TC tracks that are formatted for reading into SFINCS.
1. Select the TC IDs from TC Track files
   - Option by point: `wind/select_TCs_by_buffer.m`
   - Option by area: `tracks/tracks_select_by_polygon.ipynb`
3. Generate wind fields for selected TC IDs using the CLE15 Model (`CLE15_model`)
4. Reformat the wind fields for SFINCS (`wind/wind_process_grids.py`)

