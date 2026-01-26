# Probablistic TC Flood Hazards

This repository contains the codes developed for probabilistic simulations of TC flooding leverage multiple physics-based models. The engine of this workflow is SFINCS which is an open-access 2D hydrodynamic model (link). These codes are configured for the study Grimley, Sebastian, and Gori (2026). The codes and workflows can be adapted for probablistic modeling and efficient large-scale data analysis. 

Key workflows:
1. Processing a TC Dataset forom Gori et al. (2022)
2. Setting up the SFINCS model inputs
3. Submitting thousands of SFINCS model runs in parallel
4. Post-processing of the SFINCS outputs

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

## Data Availability

Note: input datasets are not included in this repository.

**TC Dataset**  
- Source: Gori et al. 2022
- Data included (as MATLAB files):
  - TC Tracks (MIT Model) 
  - Storm Tide (ADCIRC Model)
  - Rainfall Fields (TCR Model)
  - Wind Fields (CLE15 Model)

- **SFINCS Model**
- Available for download: Grimley, L., A. Sebastian (2025). Topobathymetric Digital Elevation Models (DEM) for Flood Modeling in the Carolinas [Version 2]. DesignSafe-CI. https://doi.org/10.17603/ds2-mzc8-s589

- **SFINCS Model Outputs**
- Available to download: Grimley, L., A. Sebastian, A. Gori (2025). Historic and Future Probabilistic Tropical Cyclone Flood Hazards for the Carolinas. DesignSafe-CI. https://doi.org/10.17603/ds2-kbhj-k266

---
