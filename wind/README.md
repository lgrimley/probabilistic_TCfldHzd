# Synthetic TC Wind Fields

Description

Data Requirements:

---
### Workflow
- Get a list of TC IDs to processes using **wind/select_TCs_by_buffer.m** or the scripts in **tracks/tracks_select_by_polygon.ipynb**
- Apply wind model (**CLE15_model** folder) to generate wind fields for the select TC IDs using **wind/CLE15_model/getmaxwind_county_SynStorms_grid_LG.m**. 
- Reformat the gridded wind fields for application to SFINCS models using **wind/wind_process_grids.py**
---



