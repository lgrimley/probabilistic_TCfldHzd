import os
import xarray as xr
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('TkAgg')
plt.ion()

# ----------------------------------------------
# HISTORIC Period
# ------------------------------------------------
#
# def get_sorted_aep_waterlevels(
#     sorted_data: xr.DataArray,
#     return_periods: list,
#     outputfile: str,
#     basin_wide_arrival_rate: float = 3.38,
#     basin_n_storms: int = 5018
# ) -> xr.DataArray:
#     """
#     Compute AEP water levels for historical data using Poisson storm arrival model.
#     sorted_data must have a 'probability' coordinate.
#     """
#     n_storms = len(sorted_data)
#     lambda_param = basin_wide_arrival_rate * (n_storms / basin_n_storms)
#     aep_probabilities = np.array([1 - np.exp(-lambda_param / T) for T in return_periods])
#
#     # Clean up the information in the sorted data
#     sorted_data = sorted_data.rename(
#         {"probability": "storm_exceedance_probability"}
#     )
#
#     sorted_data["storm_exceedance_probability"].attrs = {
#         "long_name": "Storm exceedance probability",
#         "standard_name": "exceedance_probability",
#         "units": "1",
#         "description": (
#             "Exceedance probability of a water level for an individual storm event, "
#             "derived from the empirical distribution of simulated storms. "
#             "This probability does not account for storm arrival frequency."
#         )
#     }
#
#     water_levels = []
#     for T, aep in zip(return_periods, aep_probabilities):
#         closest_idx = np.abs(
#             sorted_data["storm_exceedance_probability"] - aep
#         ).argmin(dim="storm_exceedance_probability")
#         rp_data = sorted_data.isel(
#             storm_exceedance_probability=closest_idx
#         ).astype('float32')
#         water_levels.append(rp_data)
#
#     AEP_waterlevel = xr.concat(water_levels, dim="return_period")
#     AEP_waterlevel.name = "AEP_waterlevel"
#     AEP_waterlevel = AEP_waterlevel.assign_coords(
#         return_period=("return_period", return_periods),
#         aep=("return_period", aep_probabilities)
#     )
#
#     AEP_waterlevel.attrs = {
#         "long_name": "SFINCS Modeled AEP Water Levels",
#         "standard_name": "water_surface_height",
#         "units": sorted_data.attrs.get("units", "m"),
#         "description": (
#             "Synthetic TC Water levels corresponding to specified return periods, "
#             "derived from a Poisson storm arrival model."
#         ),
#         "basin_wide_arrival_rate": basin_wide_arrival_rate,
#         "basin_number_of_storms": basin_n_storms,
#         "method": "Closest exceedance probability matching"
#     }
#     AEP_waterlevel["aep"].attrs = {
#         "long_name": "Annual exceedance probability",
#         "standard_name": "annual_exceedance_probability",
#         "units": "1",
#         "description": (
#             "Probability that the specified water level is exceeded at least once "
#             "in a given year, computed assuming storm arrivals follow a Poisson process."
#         ),
#         "arrival_rate_assumption": "Poisson"
#     }
#     AEP_waterlevel["spatial_ref"] = sorted_data["spatial_ref"]
#     print(f"Creating {outputfile}")
#     AEP_waterlevel.to_netcdf(outputfile)
#
#     return AEP_waterlevel
#
# base_hist = r'Z:\Data-Expansion\users\lelise\projects\Carolinas_SFINCS\Chapter3_SyntheticTCs\04_MODEL_OUTPUTS\ncep\aep'
# os.makedirs(os.path.join(base_hist, 'probabilistic_WSE'), exist_ok=True)
#
# # Return periods (years) for which AEP water levels are computed
# rps = [5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200, 250, 300, 350, 400, 500]
# for scenario in ['compound', 'coastal', 'runoff']:
#     # Historical
#     p = os.path.join(base_hist, 'sorted_peakWSE', f'ncep_sorted_peakWSE_{scenario}.nc')
#     sorted_data = xr.open_dataarray(p)
#     rp_hist = get_sorted_aep_waterlevels(
#         sorted_data=sorted_data,
#         return_periods=rps,
#         outputfile=os.path.join(base_hist, 'probabilistic_WSE', f'ncep_AEPs_WSE_{scenario}.nc'),
#         basin_wide_arrival_rate=3.38,
#         basin_n_storms=5018
#     )
#     print(scenario)
#
# rp100 = rp_hist.sel(return_period=100)
# rp100.plot()

# ----------------------------------------------
# FUTURE Period - need to consider storm weights
# ------------------------------------------------

base_future = r'Z:\Data-Expansion\users\lelise\projects\Carolinas_SFINCS\Chapter3_SyntheticTCs\04_MODEL_OUTPUTS\canesm_ssp585\aep'
os.makedirs(os.path.join(base_future, 'probabilistic_WSE'), exist_ok=True)
basin_wide_arrival_rate = 3.38
basin_n_storms = 6200
return_periods = [5, 10, 25, 30, 50, 100, 250, 500]

for scenario in ['runoff','compound', 'coastal']:
    if scenario == 'compound':
        outscenario = 'total'
    else:
        outscenario = scenario

    outputfile = os.path.join(base_future, 'probabilistic_WSE', f'canesm_ssp585_AEPs_WSE_{outscenario}.nc')
    sorted_data_future  = xr.open_dataarray(os.path.join(base_future, 'sorted_data',f'canesm_sorted_peakWSE_{scenario}.nc'))
    sorted_weights = xr.open_dataarray(os.path.join(base_future,'sorted_data', f'canesm_sorted_weights_{scenario}.nc'))
    cumulative_weights = sorted_weights.cumsum(dim='rank', skipna=False) # Calculate the cumulative sum of the sorted weights to form the CDF
    print('Finished calculating cumulative probability...')

    n_storms = len(sorted_data_future)
    lambda_param = basin_wide_arrival_rate * (n_storms / basin_n_storms)
    aep_probabilities = np.array([1 - np.exp(-lambda_param / T) for T in return_periods])

    return_period_water_levels = []
    for T in return_periods:
        exceedance_probability = np.round(np.exp(-lambda_param / T), decimals=5)
        probability = np.round(1 - exceedance_probability, decimals=5)

        # Calculate the difference between the cdf and the chosen probability
        diff = np.abs(cumulative_weights - exceedance_probability)

        # Find the closest rank
        closest_rank = diff.argmin(dim='rank').compute()

        # Select these water levels
        return_period_data = sorted_data_future.isel(rank=closest_rank)
        return_period_water_levels.append(return_period_data)
        print(f'Finished processing return period {T} for {outscenario} scenario')

    # Stack return periods into a single DataArray
    AEP_waterlevel = xr.concat(return_period_water_levels, dim="return_period")
    AEP_waterlevel = AEP_waterlevel.astype("float32")
    AEP_waterlevel.name = "AEP_waterlevel"

    # Assign coordinates
    AEP_waterlevel = AEP_waterlevel.assign_coords(
        return_period=("return_period", return_periods),
        aep=("return_period", aep_probabilities),
    )

    # Metadata
    AEP_waterlevel.attrs = {
        "long_name": "SFINCS Modeled Annual Exceedance Probability (AEP) TC Water Levels",
        "standard_name": "water_surface_height",
        "units": sorted_data_future.attrs.get("units", "m"),
        "description": (
            "Synthetic tropical cyclone water levels corresponding to specified return periods, "
            "derived using a Poisson storm arrival model."
        ),
        "basin_wide_arrival_rate": basin_wide_arrival_rate,
        "basin_number_of_storms": basin_n_storms,
        "method": "Closest exceedance probability matching",
        "climate_scenario":"High-emissions scenario SSP585 CMIP6 CanESM",
        "SFINCS_scenario":outscenario
    }

    AEP_waterlevel["aep"].attrs = {
        "long_name": "Annual exceedance probability",
        "standard_name": "annual_exceedance_probability",
        "units": "1",
        "description": (
            "Probability that the specified water level is exceeded at least once "
            "in a given year, assuming storm arrivals follow a Poisson process."
        ),
        "arrival_rate_assumption": "Poisson",
    }

    # Preserve spatial reference if available
    if "spatial_ref" in sorted_data_future.coords:
        AEP_waterlevel["spatial_ref"] = sorted_data_future["spatial_ref"]

    AEP_waterlevel2 = AEP_waterlevel.drop_vars('rank')

    print(f"Creating NetCDF: {outputfile}")
    AEP_waterlevel2.to_netcdf(outputfile)