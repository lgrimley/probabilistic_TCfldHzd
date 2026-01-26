# Probablistic TC Flood Hazards

This repository contains the codes developed for probabilistic simulations of TC flooding leveraging multiple physics-based models. The engine of this framework is SFINCS, an open-access 2D hydrodynamic model (https://sfincs.readthedocs.io/en/latest/index.html). When possible, we focused on using open-source python packages including HydroMT (https://deltares.github.io/hydromt/latest/index.html). The codes used to process inputs and outputs from the SFINCS model can be adapted for other case studies. The codes in this repository are configured for a case study in North and South Carolina (Grimley et al. 2026). 

Key workflows:
1. Processing TC Dataset (Gori et al. 2022)
2. Setting up the SFINCS model runs
3. Submitting thousands of SFINCS model runs in parallel
4. Post-processing SFINCS outputs

---

## Repository Structure

```
.
├── .idea/
├── analysis/
├── rain/
├── sfincs/
├── src/
├── stormTide/
├── tracks/
├── wind/
├── LICENSE
├── README.md
├── __init__.py

```

---

## Dependencies

### MATLAB
- MATLAB (recent versions)

### Python
- Python ≥ 3.8
- `xarray`
- `scipy`
- `matplotlib`
- `rioxarray`
  
- Custom modules:
  - `src.utils`
  - `src.core`

---

## Data and Model Availability

Note: input datasets are not included in this repository but can be requested or accessed using the information below:

**TC Dataset**  
- Source: Gori et al. 2022 (https://www.nature.com/articles/s41558-021-01272-7)
- Data as MATLAB files:
  - TC Tracks (MIT Model) 
  - Storm Tide (ADCIRC Model)
  - Rainfall Fields (TCR Model)
  - Wind Fields (CLE15 Model)

**SFINCS Model**
- Software: https://sfincs.readthedocs.io/en/latest/index.html
- North and South Carolina Model
  - Grimley, L., A. Sebastian (2025). Topobathymetric Digital Elevation Models (DEM) for Flood Modeling in the Carolinas [Version 2]. DesignSafe-CI. https://doi.org/10.17603/ds2-mzc8-s589
  - https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2023WR036727

**SFINCS Outputs**
Available to download: Grimley, L., A. Sebastian, A. Gori (2025). Historic and Future Probabilistic Tropical Cyclone Flood Hazards for the Carolinas. DesignSafe-CI. https://doi.org/10.17603/ds2-kbhj-k266

---
