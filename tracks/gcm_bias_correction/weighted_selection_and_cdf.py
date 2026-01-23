import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import random
import sys

# Add root to sys.path so custom modules could be imported (if needed)
sys.path.append(r'/')

import matplotlib as mpl
# Use TkAgg backend for interactive plotting
mpl.use('TkAgg')
plt.ion()  # interactive mode

# Set font parameters for plots
font = {'family': 'Arial', 'size': 10}
mpl.rc('font', **font)
mpl.rcParams.update({'axes.titlesize': 10})
mpl.rcParams["figure.autolayout"] = True

# -------------------------------------------------------------------------------------------------
# Define function for weighted sampling without replacement
# -------------------------------------------------------------------------------------------------
def weighted_sampling_without_replacement(data, weights, num_samples):
    """
    Sample `num_samples` points from `data` according to `weights` without replacement.
    
    Parameters:
        data (array-like): Data points to sample from
        weights (array-like): Weights corresponding to each data point
        num_samples (int): Number of points to sample

    Returns:
        sampled_indices (list): Indices of sampled data points
    """
    if num_samples > len(data):
        raise ValueError("Number of samples cannot be greater than the number of data points.")

    # Pair each data point with its weight and index
    population = list(enumerate(zip(data, weights)))  # (index, (data, weight))
    sampled_indices = []

    # Loop to sample without replacement
    for _ in range(num_samples):
        # sum of weights for remaining population
        total_weight = sum(w for _, (_, w) in population)
        # Select an index based on weights
        selected_index = random.choices(
            range(len(population)), weights=[w for _, (_, w) in population], k=1
        )[0]

        # Remove the selected item from population to prevent replacement
        sampled_item, _ = population.pop(selected_index)
        sampled_indices.append(sampled_item)  # store original index

    return sampled_indices

# -------------------------------------------------------------------------------------------------
# Set working directory for bias correction outputs
# -------------------------------------------------------------------------------------------------
os.chdir(r'Z:\Data-Expansion\users\lelise\projects\Carolinas_SFINCS\Chapter3_SyntheticTCs\02_DATA\BiasCorrection')

# -------------------------------------------------------------------------------------------------
# Load observed NCEP landfall wind speeds
# -------------------------------------------------------------------------------------------------
ncep = pd.read_csv(
    r'Z:\Data-Expansion\users\lelise\projects\Carolinas_SFINCS\Chapter3_SyntheticTCs\02_DATA\NCEP_Reanalysis\tracks\ncep_landfall_vmax_ZerosRemoved.csv',
    index_col=0
)
ncep['vmax'] = ncep['vstore100'].astype(float)  # create numeric vmax column
ncep_sorted = ncep.sort_values(by='vmax', axis=0, ascending=True)
# Empirical CDF for NCEP
ncep_cdf = np.arange(1, len(ncep_sorted) + 1) / len(ncep_sorted)

# -------------------------------------------------------------------------------------------------
# Load GCM bias-corrected tracks for CANESM
# -------------------------------------------------------------------------------------------------
gcm = 'canesm'
df = pd.read_csv(f'{gcm}_ssp585_weighted.csv', header=None)
df.columns = ['tc_id','vmax', 'gcm_weight']  # assign column names
df['vmax'] = df['vmax'].astype(float)
df['gcm_weight'] = df['gcm_weight'].astype(float)
df['gcm_tcid'] = [f'{gcm}_{x}' for x in df['tc_id']]  # unique TC identifier
sorted_data = df.sort_values(by='vmax', axis=0, ascending=True)
cdf = np.arange(1, len(sorted_data) + 1) / len(sorted_data)

# -------------------------------------------------------------------------------------------------
# Plot empirical CDF comparing NCEP and CANESM bias-corrected wind speeds
# -------------------------------------------------------------------------------------------------
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5, 3.5), sharex=True, sharey=True)
ax.plot(ncep_sorted['vmax'], ncep_cdf,
        linestyle='-', label='NCEP', linewidth=3, alpha=1, color='grey')
# Raw CANESM CDF plot commented out
# ax.plot(sorted_data['vmax'], cdf,
#         linestyle='-', label='CANESM', alpha=0.7, color='red')
# Use cumulative sum of bias-corrected weights for CDF
ax.plot(sorted_data['vmax'], sorted_data['gcm_weight'].cumsum(),
        linestyle='-', label=f'CANESM Bias Corrected', linewidth=3, alpha=1, color='black')
ax.set_xlabel('Maximum Sustained Wind Speeds (knots)', fontsize=12)
ax.set_ylabel('CDF', fontsize=12)
ax.legend()
ax.grid(False)
plt.subplots_adjust(wspace=0, hspace=0)
plt.margins(0,0)
plt.savefig(
    fr'Z:\Data-Expansion\users\lelise\projects\Carolinas_SFINCS\Chapter3_SyntheticTCs\02_DATA\BiasCorrection\tc_event_sets_vmax_cdf.png',
    dpi=300, bbox_inches="tight"
)
plt.close()

# -------------------------------------------------------------------------------------------------
# Combine multiple GCMS (currently only CANESM)
# -------------------------------------------------------------------------------------------------
gcms = ['canesm']  # can add more: 'cnrm6', 'ecearth6', etc.
full_set = pd.DataFrame()

for gcm in gcms:
    # Load GCM weighted tracks
    df = pd.read_csv(f'{gcm}_ssp585_weighted.csv', header=None)
    df.columns = ['tc_id','vmax', 'gcm_weight']
    df['vmax'] = df['vmax'].astype(float)
    df['gcm_weight'] = df['gcm_weight'].astype(float)
    df['gcm_tcid'] = [f'{gcm}_{x}' for x in df['tc_id']]

    # Append to full set
    full_set = pd.concat(objs=[full_set,df], axis = 0, ignore_index=True)

    # Plot raw vs bias-corrected CDF for each GCM
    sorted_data = df.sort_values(by='vmax', axis=0, ascending=True)
    cdf = np.arange(1, len(sorted_data) + 1) / len(sorted_data)
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6, 4), sharex=True, sharey=True)
    ax.plot(sorted_data['vmax'], cdf,
            linestyle='-', label='Raw Full Set', alpha=0.7, color='red')
    ax.plot(sorted_data['vmax'], sorted_data['gcm_weight'].cumsum(),
            linestyle='-', label=f'Bias Corrected', alpha=0.7, color='blue')
    ax.set_xlabel('Vmax', fontsize=12)
    ax.set_ylabel('CDF', fontsize=12)
    ax.legend()
    ax.set_title(f'{gcm}')
    ax.grid(False)
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.margins(0,0)
    #plt.savefig(fr'track_vmax_cdf_{gcm}.png', dpi=300, bbox_inches="tight")
    #plt.close()

# -------------------------------------------------------------------------------------------------
# Normalize weights across full set
# -------------------------------------------------------------------------------------------------
full_set['full_set_weight'] = full_set['gcm_weight'].to_numpy() / full_set['gcm_weight'].cumsum().values[-1]
full_set_ordered = full_set.sort_values(by='vmax', axis=0, ascending=True)
full_set_ordered.reset_index(inplace=True)
full_set_ordered.drop(columns='index', inplace=True)
#full_set_ordered.to_csv('full_set.csv', index=False)
print(full_set_ordered['full_set_weight'].sum())  # should be 1
print(full_set_ordered['full_set_weight'].cumsum().values[-1])  # should also be 1

# -------------------------------------------------------------------------------------------------
# Weighted sampling without replacement for different sample sizes
# -------------------------------------------------------------------------------------------------
for num_samples in np.arange(5000, 9500, 500):
    # Use custom weighted sampling function (without replacement)
    selected3 = weighted_sampling_without_replacement(
        data=full_set_ordered['vmax'].to_numpy(),
        weights=full_set_ordered['full_set_weight'].to_numpy(),
        num_samples=num_samples
    )

    # Subset the sampled storms
    subset = full_set_ordered[full_set_ordered.index.isin(selected3)].copy()
    subset_ordered = subset.sort_values(by='vmax', axis=0, ascending=True)

    # Normalize the weights for the sampled subset
    subset_ordered['sample_set_weight'] = subset_ordered['full_set_weight'].to_numpy() / subset_ordered['full_set_weight'].cumsum().values[-1]
    subset_ordered.to_csv(f'samples_set_{num_samples}.csv', index=False)
    print(subset_ordered['sample_set_weight'].sum())
    print(subset_ordered['sample_set_weight'].cumsum().values[-1])

    # Kolmogorov-Smirnov test between full set and sample
    ks_statistic, p_value = stats.ks_2samp(full_set_ordered['vmax'], subset_ordered['vmax'])

    # CDF for full set
    cdf = np.arange(1, len(full_set_ordered) + 1) / len(full_set_ordered)

    # Plot CDF comparison
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6, 4), sharex=True, sharey=True)
    ax.plot(full_set_ordered['vmax'].to_numpy(), cdf,
            linestyle='-', label=f'Raw Full Set ({len(full_set_ordered)})', alpha=0.7, color='red')
    ax.plot(full_set_ordered['vmax'].to_numpy(), full_set_ordered['full_set_weight'].cumsum().to_numpy(),
            linestyle='-', label=f'Bias Corrected Full Set ({len(full_set_ordered)})', alpha=0.7, color='green')
    ax.plot(subset_ordered['vmax'].to_numpy(), subset_ordered['sample_set_weight'].cumsum().to_numpy(),
            linestyle='-', label=f'Bias Corrected Sample ({num_samples})', alpha=0.7, color='blue')
    ax.set_xlabel('Vmax', fontsize=12)
    ax.set_ylabel('CDF', fontsize=12)
    ax.legend()
    ax.set_title(f'ks_stat:{np.round(ks_statistic,3)}\np-value:{np.format_float_scientific(p_value, 3)}')
    ax.grid(False)
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.margins(0,0)
    plt.savefig(fr'track_vmax_cdf_{num_samples}_v3.png', dpi=300, bbox_inches="tight")
    plt.close()
