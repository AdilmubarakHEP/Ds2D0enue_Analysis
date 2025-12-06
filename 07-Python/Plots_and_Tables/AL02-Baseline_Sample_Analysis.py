# %%
import uproot
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm.auto import tqdm

# %%
plt.rcParams.update({
    "axes.labelsize": 14,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 12,
    "figure.titlesize": 16
})

# %% [markdown]
#  # Configuration

# %%
# ========================================
# BACKGROUND CATEGORY TOGGLE
# ========================================
USE_4_CATEGORIES = False  # Set to True for 4 categories (split NaN), False for 3 categories (NaN in Comb)

# %%
# Define decay modes, trees, and cuts
decay_config = {
    'kmpip': {
        'ds_tree': 'DstreeCh1',
        'cut': "(-0.014291 <= D0_dM) & (D0_dM <= 0.014152)"
    },
    'km3pi': {
        'ds_tree': 'DstreeCh3',
        'cut': "(-0.013093 <= D0_dM) & (D0_dM <= 0.012520)"
    },
    'kmpippi0_eff20_May2020': {
        'ds_tree': 'DstreeCh2',
        'cut': "(-0.052152 <= D0_dM) & (D0_dM <= 0.024237)"
    }
}

# Background samples
background_samples = ["ccbar", "charged", "ddbar", "mixed", "ssbar", "uubar"]

# %% [markdown]
#  # Prep-Work

# %% [markdown]
#  ### Import Data

# %%
DataFrames = {}  # define empty dictionary to hold dataframes

# Load Signal files
print("Loading Signal files...")
for mode, config in tqdm(list(decay_config.items()), desc="Signal modes"):
    signal_file = f"/home/belle2/amubarak/C01-Simulated_Events/Signal/output_test_{mode}.root"
    df = uproot.concatenate(f"{signal_file}:{config['ds_tree']}", library='pd')
    df = df.query(config['cut'])
    DataFrames[f"Signal_{mode}"] = df

# Load Background files
print("\nLoading Background files...")
for sample in tqdm(background_samples, desc="Background samples"):
    for mode, config in decay_config.items():
        generic_file = f"/group/belle/users/amubarak/03-KEKCC/Ds2D0e-Generic_Ds_120325_1_{sample}_{mode}.root"
        df = uproot.concatenate(f"{generic_file}:{config['ds_tree']}", library='pd')
        df = df.query(config['cut'])
        DataFrames[f"{sample}_{mode}"] = df

# Combine all backgrounds for each mode
print("\nCombining background samples by mode...")
for mode in decay_config.keys():
    dfs_list = [DataFrames[f"{sample}_{mode}"] for sample in background_samples]
    DataFrames[f"All_{mode}"] = pd.concat(dfs_list, ignore_index=True)

print("\nData loading complete!")
print(f"Successfully loaded {len(DataFrames)} dataframes")
print(f"\nUsing {'4' if USE_4_CATEGORIES else '3'} background categories")

# %%
pd.set_option('display.max_rows', 200000)
pd.set_option('display.max_columns', 200000)

# %% [markdown]
# # PRELIMINARY RESULTS
# ---

# %% [markdown]
#  ## Value Counts - Ds_mcPDG

# %% [markdown]
#  ### Signal Samples

# %%
for mode in decay_config.keys():
    print(f"\n{'='*80}")
    print(f"Signal_{mode}")
    print('='*80)
    print(abs(DataFrames[f"Signal_{mode}"]['Ds_mcPDG']).value_counts(normalize=True, dropna=False).apply(lambda x: f"{x:.6f}"))

# %% [markdown]
#  ### Background Samples (Combined by Mode)

# %%
for mode in decay_config.keys():
    print(f"\n{'='*80}")
    print(f"All_{mode} (ccbar + charged + ddbar + mixed + ssbar + uubar)")
    print('='*80)
    print(abs(DataFrames[f"All_{mode}"]['Ds_mcPDG']).value_counts(normalize=True, dropna=False).apply(lambda x: f"{x:.6f}"))

# %% [markdown]
#  ### Individual Background Samples

# %%
for sample in background_samples:
    print(f"\n{'='*80}")
    print(f"{sample.upper()}")
    print('='*80)
    for mode in decay_config.keys():
        key = f"{sample}_{mode}"
        print(f"\n  {key}:")
        print(f"  {'-'*76}")
        print(abs(DataFrames[key]['Ds_mcPDG']).value_counts(normalize=True, dropna=False).head(10).apply(lambda x: f"{x:.6f}"))

# %% [markdown]
#  ## Basic Distribution Plots

# %% [markdown]
#  ### Signal Plots

# %%
# Signal plotting parameters
dM = 0.1  # gammaveto_M_Correction cut value

# %%
# Plot Signal - Ds_diff_D0pi
Bins = 50
Density = False
Stacked = False
Range = [0.1, 0.7]
var = 'Ds_diff_D0pi'

for mode in decay_config.keys():
    perBin = ((Range[1] - Range[0])/Bins)*1000
    print(f"Mode: {mode}, Width Per Bin: {perBin:.2f} MeV")
    
    label1 = r'$isSignal(D_s^{+})=1$'
    label2 = r'$isSignal(D_s^{+})=0$'
    label3 = r'$NaN$'
    
    labels = [label1, label2, label3]
    colors = ['#fd7f6f', '#7eb0d5', 'purple']
    
    df = DataFrames[f"Signal_{mode}"]
    df_cut = df[df['gammaveto_M_Correction'] >= dM]
    
    data = [
        df_cut[df_cut['Ds_isSignal'] == 1][var],
        df_cut[df_cut['Ds_isSignal'] == 0][var],
        df_cut[df_cut['Ds_isSignal'].isna()][var]
    ]
    
    plt.hist(data[::-1], color=colors[::-1], label=labels[::-1], alpha=1, 
             range=Range, stacked=Stacked, density=Density, bins=Bins, 
             histtype='step', linewidth=1.5)
    
    # Title
    plt.title(f'Signal - {mode}', loc="left")
    plt.title(r'$\bf Signal\;Events$', loc="right")
    # Label
    plt.ylabel(r'$Entries/(\; {width:.2f}\;MeV/c^2)$'.format(width=perBin))
    plt.xlabel(r'$\Delta m(D_s^{+} - D^{0})\;[GeV/c^{2}]$')
    # plt.yscale("log")
    plt.legend()
    plt.show()

# %%
# Plot Signal - Ds_massDifference_0
Bins = 50
Density = False
Stacked = False
Range = [0.0, 0.25]
var = 'Ds_massDifference_0'

for mode in decay_config.keys():
    perBin = ((Range[1] - Range[0])/Bins)*1000
    print(f"Mode: {mode}, Width Per Bin: {perBin:.2f} MeV")
    
    label1 = r'$isSignal(D_s^{+})=1$'
    label2 = r'$isSignal(D_s^{+})=0$'
    label3 = r'$NaN$'
    
    labels = [label1, label2, label3]
    colors = ['#fd7f6f', '#7eb0d5', 'purple']
    
    df = DataFrames[f"Signal_{mode}"]
    df_cut = df[df['gammaveto_M_Correction'] >= dM]
    
    data = [
        df_cut[df_cut['Ds_isSignal'] == 1][var],
        df_cut[df_cut['Ds_isSignal'] == 0][var],
        df_cut[df_cut['Ds_isSignal'].isna()][var]
    ]
    
    plt.hist(data[::-1], color=colors[::-1], label=labels[::-1], alpha=1, 
             range=Range, stacked=Stacked, density=Density, bins=Bins, 
             histtype='step', linewidth=1.5)
    
    # Title
    plt.title(f'Signal - {mode}', loc="left")
    plt.title(r'$\bf Signal\;Events$', loc="right")
    # Label
    plt.ylabel(r'$Entries/(\; {width:.2f}\;MeV/c^2)$'.format(width=perBin))
    plt.xlabel(r'$\Delta m(D_s^{+} - D^{0})\;[GeV/c^{2}]$')
    # plt.yscale("log")
    plt.legend()
    plt.show()

# %% [markdown]
#  ### Background Plots

# %%
# Background plotting parameters
dM = 0.1  # gammaveto_M_Correction cut value

# %%
# Plot Background - Ds_diff_D0pi
Bins = 50
Density = False
Stacked = True
Range = [0.1, 0.7]
var = 'Ds_diff_D0pi'

for mode in decay_config.keys():
    perBin = ((Range[1] - Range[0])/Bins)*1000
    print(f"Mode: {mode}, Width Per Bin: {perBin:.2f} MeV")
    
    df = DataFrames[f"All_{mode}"]
    df_cut = df[df['gammaveto_M_Correction'] >= dM]
    
    if USE_4_CATEGORIES:
        label1 = r'$Comb.$'
        label2 = r'$NaN$'
        label3 = r'$D^{*0}$'
        label4 = r'$D^{*+} \rightarrow D^{0} \pi^{+}$'
        
        labels = [label1, label2, label3, label4]
        colors = ["#DD8452", "#C44E52", "#55A868", "#4C72B0"]
        
        data = [
            df_cut[(abs(df_cut["Ds_mcPDG"]) != 413) & 
                   (abs(df_cut["Ds_mcPDG"]) != 423) & 
                   (~df_cut["Ds_mcPDG"].isna())][var],
            df_cut[df_cut["Ds_mcPDG"].isna()][var],
            df_cut[abs(df_cut["Ds_mcPDG"]) == 423][var],
            df_cut[abs(df_cut["Ds_mcPDG"]) == 413][var]
        ]
    else:  # 3 categories
        label1 = r'$Comb.$'
        label2 = r'$D^{*0}$'
        label3 = r'$D^{*+} \rightarrow D^{0} \pi^{+}$'
        
        labels = [label1, label2, label3]
        colors = ["#DD8452", "#55A868", "#4C72B0"]
        
        data = [
            df_cut[(abs(df_cut["Ds_mcPDG"]) != 413) & 
                   (abs(df_cut["Ds_mcPDG"]) != 423)][var],
            df_cut[abs(df_cut["Ds_mcPDG"]) == 423][var],
            df_cut[abs(df_cut["Ds_mcPDG"]) == 413][var]
        ]
    
    plt.hist(data, color=colors, label=labels, density=Density, 
             stacked=Stacked, bins=Bins, alpha=1, histtype='step', 
             linewidth=1.5, range=Range)
    
    # Title
    plt.title(f'Background - {mode}', loc="left")
    plt.title(r'$\int\mathcal{L}dt\approx\;1444$ fb$^{-1}$', loc="right")
    # Label
    plt.ylabel(r'$Entries/(\; {width:.2f}\;MeV/c^2)$'.format(width=perBin))
    plt.xlabel(r'$\Delta m(D_s^{+} - D^{0})\;[GeV/c^{2}]$')
    plt.legend(loc='upper right')
    plt.show()

# %%
# Plot Background - Ds_massDifference_0
Bins = 50
Density = False
Stacked = True
Range = [0.0, 0.25]
var = 'Ds_massDifference_0'

for mode in decay_config.keys():
    perBin = ((Range[1] - Range[0])/Bins)*1000
    print(f"Mode: {mode}, Width Per Bin: {perBin:.2f} MeV")
    
    df = DataFrames[f"All_{mode}"]
    df_cut = df[df['gammaveto_M_Correction'] >= dM]
    
    if USE_4_CATEGORIES:
        label1 = r'$Comb.$'
        label2 = r'$NaN$'
        label3 = r'$D^{*0}$'
        label4 = r'$D^{*+} \rightarrow D^{0} \pi^{+}$'
        
        labels = [label1, label2, label3, label4]
        colors = ["#DD8452", "#C44E52", "#55A868", "#4C72B0"]
        
        data = [
            df_cut[(abs(df_cut["Ds_mcPDG"]) != 413) & 
                   (abs(df_cut["Ds_mcPDG"]) != 423) & 
                   (~df_cut["Ds_mcPDG"].isna())][var],
            df_cut[df_cut["Ds_mcPDG"].isna()][var],
            df_cut[abs(df_cut["Ds_mcPDG"]) == 423][var],
            df_cut[abs(df_cut["Ds_mcPDG"]) == 413][var]
        ]
    else:  # 3 categories
        label1 = r'$Comb.$'
        label2 = r'$D^{*0}$'
        label3 = r'$D^{*+} \rightarrow D^{0} \pi^{+}$'
        
        labels = [label1, label2, label3]
        colors = ["#DD8452", "#55A868", "#4C72B0"]
        
        data = [
            df_cut[(abs(df_cut["Ds_mcPDG"]) != 413) & 
                   (abs(df_cut["Ds_mcPDG"]) != 423)][var],
            df_cut[abs(df_cut["Ds_mcPDG"]) == 423][var],
            df_cut[abs(df_cut["Ds_mcPDG"]) == 413][var]
        ]
    
    plt.hist(data, color=colors, label=labels, density=Density, 
             stacked=Stacked, bins=Bins, alpha=1, histtype='step', 
             linewidth=1.5, range=Range)
    
    # Title
    plt.title(f'Background - {mode}', loc="left")
    plt.title(r'$\int\mathcal{L}dt\approx\;1444$ fb$^{-1}$', loc="right")
    # Label
    plt.ylabel(r'$Entries/(\; {width:.2f}\;MeV/c^2)$'.format(width=perBin))
    plt.xlabel(r'$\Delta m(D_s^{+} - D^{0})\;[GeV/c^{2}]$')
    plt.legend(loc='upper right')
    plt.show()

# %% [markdown]
#  # DEEP ANALYSIS: Δm_e vs Δm_π SANITY CHECKS
#  ---
# 
#  **Goal**: Verify whether a 2D fit in (Δm_e, Δm_π) is justified before investing in BDT and fit development.
# 
#  ## Definitions
# 
#  * **Δm_e** = `Ds_massDifference_0` (electron mass on lepton track)
#  * **Δm_π** = `Ds_diff_D0pi` (pion mass on same track)
#  * **Δm_shift** = Δm_π - Δm_e
#  * **D*+ veto**: 0.142 < Δm_π ≤ 0.15 GeV/c²
# 
#  **Truth categories**:
#  - Signal: isSignal(D_s^+) = 1
#  - Bkg1: D*+ → D0 π+ (mcPDG = 413)
#  - Bkg2: D*0 → D0 π0/γ (mcPDG = 423)
#  - Bkg3: Combinatorial (everything else)

# %% [markdown]
#  ## Stage 1: Global Δm_e Efficiency of the Δm_π Veto
# 
#  Check if the D*+ veto (0.142 < Δm_π ≤ 0.15) sculpts the Δm_e distribution.

# %%
# Define D*+ veto
dstar_veto_min = 0.142
dstar_veto_max = 0.15
dM = 0.1  # gammaveto_M_Correction cut

# %%
# Stage 1: Signal - Δm_e efficiency plots
for mode in decay_config.keys():
    print(f"\n{'='*80}")
    print(f"Stage 1 - Signal - {mode}")
    print('='*80)
    
    df = DataFrames[f"Signal_{mode}"]
    df = df[df['gammaveto_M_Correction'] >= dM]
    df_signal = df[df['Ds_isSignal'] == 1]
    
    # Before veto
    df_before = df_signal
    
    # After veto (OUTSIDE the D*+ window)
    df_after = df_signal[~((df_signal['Ds_diff_D0pi'] > dstar_veto_min) & 
                           (df_signal['Ds_diff_D0pi'] <= dstar_veto_max))]
    
    # Make histograms
    bins = 50
    range_dm = [0.0, 0.25]
    
    hist_before, bin_edges = np.histogram(df_before['Ds_massDifference_0'], 
                                          bins=bins, range=range_dm)
    hist_after, _ = np.histogram(df_after['Ds_massDifference_0'], 
                                 bins=bins, range=range_dm)
    
    # Calculate efficiency
    with np.errstate(divide='ignore', invalid='ignore'):
        efficiency = np.where(hist_before > 0, hist_after / hist_before, 0)
        # Binomial errors
        eff_err = np.where(hist_before > 0, 
                          np.sqrt(efficiency * (1 - efficiency) / hist_before), 0)
    
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    
    # Plot before and after
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))
    
    # Top: distributions
    ax1.hist(df_before['Ds_massDifference_0'], bins=bins, range=range_dm, 
             histtype='step', linewidth=1.5, label='Before D*+ veto', color='blue')
    ax1.hist(df_after['Ds_massDifference_0'], bins=bins, range=range_dm, 
             histtype='step', linewidth=1.5, label='After D*+ veto', color='red')
    ax1.set_xlabel(r'$\Delta m_e(D_s^{+} - D^{0})\;[GeV/c^{2}]$')
    ax1.set_ylabel('Entries')
    ax1.set_title(f'Signal - {mode}', loc='left')
    ax1.legend()
    ax1.grid(alpha=0.3)
    
    # Bottom: efficiency
    ax2.errorbar(bin_centers, efficiency, yerr=eff_err, fmt='o', 
                markersize=4, capsize=3, color='black')
    ax2.axhline(y=1.0, color='gray', linestyle='--', alpha=0.5)
    ax2.set_xlabel(r'$\Delta m_e(D_s^{+} - D^{0})\;[GeV/c^{2}]$')
    ax2.set_ylabel(r'Efficiency $\varepsilon(\Delta m_e)$')
    ax2.set_ylim([0, 1.1])
    ax2.grid(alpha=0.3)
    
    plt.tight_layout()
    plt.show()
    
    # Print summary statistics
    overall_eff = len(df_after) / len(df_before) if len(df_before) > 0 else 0
    print(f"Overall efficiency: {overall_eff:.4f}")
    print(f"Events before: {len(df_before)}")
    print(f"Events after: {len(df_after)}")

# %%
# Stage 1: Background categories - Δm_e efficiency plots
for mode in decay_config.keys():
    print(f"\n{'='*80}")
    print(f"Stage 1 - Background - {mode}")
    print('='*80)
    
    df = DataFrames[f"All_{mode}"]
    df = df[df['gammaveto_M_Correction'] >= dM]
    
    # Define background categories based on toggle
    if USE_4_CATEGORIES:
        bkg_categories = {
            'D*+ → D0π+': df[abs(df['Ds_mcPDG']) == 413],
            'D*0 → D0π0/γ': df[abs(df['Ds_mcPDG']) == 423],
            'Comb': df[(abs(df['Ds_mcPDG']) != 413) & 
                       (abs(df['Ds_mcPDG']) != 423) & 
                       (~df['Ds_mcPDG'].isna())],
            'NaN': df[df['Ds_mcPDG'].isna()]
        }
        fig, axes = plt.subplots(4, 2, figsize=(15, 16))
    else:
        bkg_categories = {
            'D*+ → D0π+': df[abs(df['Ds_mcPDG']) == 413],
            'D*0 → D0π0/γ': df[abs(df['Ds_mcPDG']) == 423],
            'Comb': df[(abs(df['Ds_mcPDG']) != 413) & 
                       (abs(df['Ds_mcPDG']) != 423)]
        }
        fig, axes = plt.subplots(3, 2, figsize=(15, 12))
    
    for idx, (cat_name, df_cat) in enumerate(bkg_categories.items()):
        # Before veto
        df_before = df_cat
        
        # After veto
        df_after = df_cat[~((df_cat['Ds_diff_D0pi'] > dstar_veto_min) & 
                            (df_cat['Ds_diff_D0pi'] <= dstar_veto_max))]
        
        # Make histograms
        bins = 50
        range_dm = [0.0, 0.25]
        
        hist_before, bin_edges = np.histogram(df_before['Ds_massDifference_0'], 
                                              bins=bins, range=range_dm)
        hist_after, _ = np.histogram(df_after['Ds_massDifference_0'], 
                                     bins=bins, range=range_dm)
        
        # Calculate efficiency
        with np.errstate(divide='ignore', invalid='ignore'):
            efficiency = np.where(hist_before > 0, hist_after / hist_before, 0)
            eff_err = np.where(hist_before > 0, 
                              np.sqrt(efficiency * (1 - efficiency) / hist_before), 0)
        
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        
        # Plot distributions
        ax1 = axes[idx, 0]
        ax1.hist(df_before['Ds_massDifference_0'], bins=bins, range=range_dm, 
                histtype='step', linewidth=1.5, label='Before', color='blue')
        ax1.hist(df_after['Ds_massDifference_0'], bins=bins, range=range_dm, 
                histtype='step', linewidth=1.5, label='After', color='red')
        ax1.set_ylabel('Entries')
        ax1.set_title(f'{cat_name}', loc='left')
        ax1.legend()
        ax1.grid(alpha=0.3)
        
        # Plot efficiency
        ax2 = axes[idx, 1]
        ax2.errorbar(bin_centers, efficiency, yerr=eff_err, fmt='o', 
                    markersize=4, capsize=3, color='black')
        ax2.axhline(y=1.0, color='gray', linestyle='--', alpha=0.5)
        ax2.set_ylabel(r'$\varepsilon(\Delta m_e)$')
        ax2.set_ylim([0, 1.1])
        ax2.grid(alpha=0.3)
        
        if idx == len(bkg_categories) - 1:
            ax1.set_xlabel(r'$\Delta m_e\;[GeV/c^{2}]$')
            ax2.set_xlabel(r'$\Delta m_e\;[GeV/c^{2}]$')
        
        overall_eff = len(df_after) / len(df_before) if len(df_before) > 0 else 0
        print(f"{cat_name}: Overall eff = {overall_eff:.4f}, N_before = {len(df_before)}, N_after = {len(df_after)}")
    
    plt.suptitle(f'Background - {mode}', fontsize=16)
    plt.tight_layout()
    plt.show()

# %% [markdown]
#  ## Stage 2: Momentum-Sliced Efficiency
# 
#  Check if the veto behaves differently at low and high p_e.

# %%
# Define momentum bins
p_bins = [0.05, 0.15, 0.25, 0.35, 0.6]  # GeV
p_labels = ['0.05-0.15', '0.15-0.25', '0.25-0.35', '>0.35']

# %%
# Stage 2: Signal - momentum sliced efficiency
for mode in decay_config.keys():
    print(f"\n{'='*80}")
    print(f"Stage 2 - Signal - {mode}")
    print('='*80)
    
    df = DataFrames[f"Signal_{mode}"]
    df = df[df['gammaveto_M_Correction'] >= dM]
    df_signal = df[df['Ds_isSignal'] == 1]
    
    fig, axes = plt.subplots(len(p_bins)-1, 2, figsize=(15, 12))
    
    for idx in range(len(p_bins)-1):
        p_min = p_bins[idx]
        p_max = p_bins[idx+1]
        
        # Select momentum slice
        df_slice = df_signal[(df_signal['e_p'] >= p_min) & (df_signal['e_p'] < p_max)]
        
        # Before and after veto
        df_before = df_slice
        df_after = df_slice[~((df_slice['Ds_diff_D0pi'] > dstar_veto_min) & 
                              (df_slice['Ds_diff_D0pi'] <= dstar_veto_max))]
        
        # Make histograms
        bins = 30
        range_dm = [0.0, 0.25]
        
        hist_before, bin_edges = np.histogram(df_before['Ds_massDifference_0'], 
                                              bins=bins, range=range_dm)
        hist_after, _ = np.histogram(df_after['Ds_massDifference_0'], 
                                     bins=bins, range=range_dm)
        
        # Calculate efficiency
        with np.errstate(divide='ignore', invalid='ignore'):
            efficiency = np.where(hist_before > 0, hist_after / hist_before, 0)
            eff_err = np.where(hist_before > 0, 
                              np.sqrt(efficiency * (1 - efficiency) / hist_before), 0)
        
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        
        # Plot distributions
        ax1 = axes[idx, 0]
        ax1.hist(df_before['Ds_massDifference_0'], bins=bins, range=range_dm, 
                histtype='step', linewidth=1.5, label='Before', color='blue')
        ax1.hist(df_after['Ds_massDifference_0'], bins=bins, range=range_dm, 
                histtype='step', linewidth=1.5, label='After', color='red')
        ax1.set_ylabel('Entries')
        ax1.set_title(f'$p_e$ = {p_labels[idx]} GeV', loc='left')
        ax1.legend()
        ax1.grid(alpha=0.3)
        
        # Plot efficiency
        ax2 = axes[idx, 1]
        ax2.errorbar(bin_centers, efficiency, yerr=eff_err, fmt='o', 
                    markersize=4, capsize=3, color='black')
        ax2.axhline(y=1.0, color='gray', linestyle='--', alpha=0.5)
        ax2.set_ylabel(r'$\varepsilon(\Delta m_e | p_e)$')
        ax2.set_ylim([0, 1.1])
        ax2.grid(alpha=0.3)
        
        if idx == len(p_bins)-2:
            ax1.set_xlabel(r'$\Delta m_e\;[GeV/c^{2}]$')
            ax2.set_xlabel(r'$\Delta m_e\;[GeV/c^{2}]$')
        
        overall_eff = len(df_after) / len(df_before) if len(df_before) > 0 else 0
        print(f"p_e {p_labels[idx]} GeV: eff = {overall_eff:.4f}, N = {len(df_before)}")
    
    plt.suptitle(f'Signal Momentum Slices - {mode}', fontsize=16)
    plt.tight_layout()
    plt.show()

# %% [markdown]
#  ## Stage 3: 2D Structure Study with Δm_shift
# 
#  Check if Δm_π gives independent information once Δm_e is known.

# %%
# Stage 3: Create Δm_shift variable and make 2D plots
for mode in decay_config.keys():
    print(f"\n{'='*80}")
    print(f"Stage 3 - 2D Structure - {mode}")
    print('='*80)
    
    df_signal = DataFrames[f"Signal_{mode}"]
    df_signal = df_signal[df_signal['gammaveto_M_Correction'] >= dM]
    df_signal = df_signal[df_signal['Ds_isSignal'] == 1].copy()
    
    df_bkg = DataFrames[f"All_{mode}"]
    df_bkg = df_bkg[df_bkg['gammaveto_M_Correction'] >= dM].copy()
    
    # Calculate Δm_shift
    df_signal['Dm_shift'] = df_signal['Ds_diff_D0pi'] - df_signal['Ds_massDifference_0']
    df_bkg['Dm_shift'] = df_bkg['Ds_diff_D0pi'] - df_bkg['Ds_massDifference_0']
    
    # Define background categories based on toggle
    df_bkg1 = df_bkg[abs(df_bkg['Ds_mcPDG']) == 413].copy()
    df_bkg2 = df_bkg[abs(df_bkg['Ds_mcPDG']) == 423].copy()
    df_bkg3 = df_bkg[(abs(df_bkg['Ds_mcPDG']) != 413) & 
                     (abs(df_bkg['Ds_mcPDG']) != 423)].copy()
    
    # 2D histograms
    bins_2d = [40, 40]
    range_2d = [[0.0, 0.25], [0.0, 0.6]]
    
    datasets = {
        'Signal': df_signal,
        'Bkg1: D*+ → D0π+': df_bkg1,
        'Bkg2: D*0 → D0π0/γ': df_bkg2,
        'Bkg3: Comb': df_bkg3
    }
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 14))
    axes = axes.flatten()
    
    for idx, (name, df) in enumerate(datasets.items()):
        if len(df) == 0:
            continue
            
        H, xedges, yedges = np.histogram2d(
            df['Ds_massDifference_0'], 
            df['Dm_shift'],
            bins=bins_2d, 
            range=range_2d
        )
        
        # Normalize
        H = H / np.sum(H) if np.sum(H) > 0 else H
        
        im = axes[idx].imshow(H.T, origin='lower', aspect='auto', 
                             extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
                             cmap='viridis', interpolation='nearest')
        
        # Mark D*+ veto region
        axes[idx].axhline(y=dstar_veto_min, color='red', linestyle='--', 
                         linewidth=2, alpha=0.7, label='D*+ veto')
        axes[idx].axhline(y=dstar_veto_max, color='red', linestyle='--', 
                         linewidth=2, alpha=0.7)
        
        axes[idx].set_xlabel(r'$\Delta m_e\;[GeV/c^{2}]$')
        axes[idx].set_ylabel(r'$\Delta m_\pi\;[GeV/c^{2}]$')
        axes[idx].set_title(f'{name} (N={len(df)})', loc='left')
        axes[idx].legend(loc='upper right')
        plt.colorbar(im, ax=axes[idx], label='Normalized density')
    
    plt.suptitle(f'2D Structure: ($\\Delta m_e$, $\\Delta m_\\pi$) - {mode}', fontsize=16)
    plt.tight_layout()
    plt.show()

# %%
# Stage 3: Δm_shift projections at fixed Δm_e slices
for mode in decay_config.keys():
    print(f"\n{'='*80}")
    print(f"Stage 3 - Δm_shift Projections - {mode}")
    print('='*80)
    
    df_signal = DataFrames[f"Signal_{mode}"]
    df_signal = df_signal[df_signal['gammaveto_M_Correction'] >= dM]
    df_signal = df_signal[df_signal['Ds_isSignal'] == 1].copy()
    
    df_bkg = DataFrames[f"All_{mode}"]
    df_bkg = df_bkg[df_bkg['gammaveto_M_Correction'] >= dM].copy()
    
    # Calculate Δm_shift
    df_signal['Dm_shift'] = df_signal['Ds_diff_D0pi'] - df_signal['Ds_massDifference_0']
    df_bkg['Dm_shift'] = df_bkg['Ds_diff_D0pi'] - df_bkg['Ds_massDifference_0']
    
    # Define background categories
    df_bkg1 = df_bkg[abs(df_bkg['Ds_mcPDG']) == 413].copy()
    df_bkg2 = df_bkg[abs(df_bkg['Ds_mcPDG']) == 423].copy()
    
    # Define Δm_e slices
    dm_e_slices = [(0.00, 0.05), (0.05, 0.10), (0.10, 0.15), (0.15, 0.20)]
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    axes = axes.flatten()
    
    for idx, (dm_min, dm_max) in enumerate(dm_e_slices):
        ax = axes[idx]
        
        # Select Δm_e slice
        sig_slice = df_signal[(df_signal['Ds_massDifference_0'] >= dm_min) & 
                              (df_signal['Ds_massDifference_0'] < dm_max)]
        bkg1_slice = df_bkg1[(df_bkg1['Ds_massDifference_0'] >= dm_min) & 
                             (df_bkg1['Ds_massDifference_0'] < dm_max)]
        bkg2_slice = df_bkg2[(df_bkg2['Ds_massDifference_0'] >= dm_min) & 
                             (df_bkg2['Ds_massDifference_0'] < dm_max)]
        
        # Plot Δm_shift distributions
        bins = 40
        range_shift = [0.0, 0.6]
        
        ax.hist(sig_slice['Dm_shift'], bins=bins, range=range_shift, 
               histtype='step', linewidth=2, label='Signal', 
               color='blue', density=True)
        ax.hist(bkg1_slice['Dm_shift'], bins=bins, range=range_shift, 
               histtype='step', linewidth=2, label='Bkg1: D*+', 
               color='red', density=True)
        ax.hist(bkg2_slice['Dm_shift'], bins=bins, range=range_shift, 
               histtype='step', linewidth=2, label='Bkg2: D*0', 
               color='green', density=True)
        
        # Mark D*+ veto
        ax.axvspan(dstar_veto_min, dstar_veto_max, alpha=0.2, 
                  color='gray', label='D*+ veto')
        
        ax.set_xlabel(r'$\Delta m_\pi\;[GeV/c^{2}]$')
        ax.set_ylabel('Normalized entries')
        ax.set_title(f'$\\Delta m_e$ = [{dm_min:.2f}, {dm_max:.2f}] GeV/c²', 
                    loc='left')
        ax.legend()
        ax.grid(alpha=0.3)
        
        print(f"Δm_e [{dm_min:.2f}, {dm_max:.2f}]: Signal N={len(sig_slice)}, "
              f"Bkg1 N={len(bkg1_slice)}, Bkg2 N={len(bkg2_slice)}")
    
    plt.suptitle(f'$\\Delta m_\\pi$ Projections at Fixed $\\Delta m_e$ - {mode}', 
                fontsize=16)
    plt.tight_layout()
    plt.show()

# %% [markdown]
#  ## Stage 3 Summary: Correlation Analysis

# %%
# Calculate correlation between Δm_e and Δm_π for each category
for mode in decay_config.keys():
    print(f"\n{'='*80}")
    print(f"Stage 3 - Correlation Analysis - {mode}")
    print('='*80)
    
    df_signal = DataFrames[f"Signal_{mode}"]
    df_signal = df_signal[df_signal['gammaveto_M_Correction'] >= dM]
    df_signal = df_signal[df_signal['Ds_isSignal'] == 1]
    
    df_bkg = DataFrames[f"All_{mode}"]
    df_bkg = df_bkg[df_bkg['gammaveto_M_Correction'] >= dM]
    
    df_bkg1 = df_bkg[abs(df_bkg['Ds_mcPDG']) == 413]
    df_bkg2 = df_bkg[abs(df_bkg['Ds_mcPDG']) == 423]
    df_bkg3 = df_bkg[(abs(df_bkg['Ds_mcPDG']) != 413) & 
                     (abs(df_bkg['Ds_mcPDG']) != 423)]
    
    datasets = {
        'Signal': df_signal,
        'Bkg1: D*+ → D0π+': df_bkg1,
        'Bkg2: D*0 → D0π0/γ': df_bkg2,
        'Bkg3: Comb': df_bkg3
    }
    
    for name, df in datasets.items():
        if len(df) > 0:
            corr = df[['Ds_massDifference_0', 'Ds_diff_D0pi']].corr().iloc[0, 1]
            print(f"{name:25s}: correlation(Δm_e, Δm_π) = {corr:.4f}")

# %% [markdown]
# ## MC-truth based veto test in Δm_e

signal_dm_min = 0.05
signal_dm_max = 0.20

dstar_veto_min = 0.142
dstar_veto_max = 0.15
dM = 0.1  # gammaveto_M_Correction cut

def apply_base_cuts(df):
    return df[df["gammaveto_M_Correction"] >= dM]

def in_dm_window(df):
    return df[(df["Ds_massDifference_0"] >= signal_dm_min) &
              (df["Ds_massDifference_0"] <= signal_dm_max)]

for mode in decay_config.keys():
    print("\n" + "="*80)
    print(f"MC-truth veto test - mode: {mode}")
    print("="*80)

    # 1. Signal: Ds_isSignal == 1
    df_sig = DataFrames[f"Signal_{mode}"]
    df_sig = apply_base_cuts(df_sig)
    df_sig = df_sig[df_sig["Ds_isSignal"] == 1]

    # no-veto selection in dm window
    sig_before = in_dm_window(df_sig)
    # veto selection in dm window
    sig_after = sig_before[~((sig_before["Ds_diff_D0pi"] > dstar_veto_min) &
                             (sig_before["Ds_diff_D0pi"] <= dstar_veto_max))]

    N_sig_before = len(sig_before)
    N_sig_after = len(sig_after)
    f_eps = N_sig_after / N_sig_before if N_sig_before > 0 else 0.0

    print(f"Signal: N_before = {N_sig_before}, N_after = {N_sig_after}, f_eps = {f_eps:.3f}")

    # 2. Background: split MC truth categories in All_mode
    df_bkg_all = DataFrames[f"All_{mode}"]
    df_bkg_all = apply_base_cuts(df_bkg_all)

    # MC truth categories
    df_bkg_Dsp = df_bkg_all[abs(df_bkg_all["Ds_mcPDG"]) == 413]  # D*+
    df_bkg_Ds0 = df_bkg_all[abs(df_bkg_all["Ds_mcPDG"]) == 423]  # D*0
    df_bkg_comb = df_bkg_all[(abs(df_bkg_all["Ds_mcPDG"]) != 413) &
                             (abs(df_bkg_all["Ds_mcPDG"]) != 423)]

    cats = {
        "D*+": df_bkg_Dsp,
        "D*0": df_bkg_Ds0,
        "Comb": df_bkg_comb,
    }

    B_before_total = 0
    B_after_total = 0

    for name, df_cat in cats.items():
        cat_before = in_dm_window(df_cat)
        cat_after = cat_before[~((cat_before["Ds_diff_D0pi"] > dstar_veto_min) &
                                 (cat_before["Ds_diff_D0pi"] <= dstar_veto_max))]

        N_b_before = len(cat_before)
        N_b_after = len(cat_after)
        f_b_cat = N_b_after / N_b_before if N_b_before > 0 else 0.0

        B_before_total += N_b_before
        B_after_total += N_b_after

        print(f"{name:4s}: N_before = {N_b_before:7d}, N_after = {N_b_after:7d}, f_b_cat = {f_b_cat:.3f}")

    f_b_tot = B_after_total / B_before_total if B_before_total > 0 else 0.0
    R = (np.sqrt(f_b_tot) / f_eps) if (f_eps > 0 and f_b_tot > 0) else np.inf

    print("-"*80)
    print(f"Total background: B_before = {B_before_total}, B_after = {B_after_total}, f_b_tot = {f_b_tot:.3f}")
    print(f"Veto FoM R = sqrt(f_b_tot) / f_eps = {R:.3f}")
    if R < 1.0:
        print("=> In this MC truth scenario, veto improves expected UL (statistically).")
    else:
        print("=> In this MC truth scenario, veto worsens or does not improve expected UL.")

# %% [markdown]
# ## Build pyhf toys: with and without veto

import pyhf
pyhf.set_backend("numpy")

# binning for Δm_e
nbins = 25
dm_range = (0.0, 0.25)

def hist_dm(series):
    counts, edges = np.histogram(series, bins=nbins, range=dm_range)
    return counts.astype(float), edges

def make_templates_for_mode(mode, apply_veto):
    """Return dict of per-process Δm_e histograms for a given mode and selection."""
    # Signal
    df_sig = DataFrames[f"Signal_{mode}"]
    df_sig = apply_base_cuts(df_sig)
    df_sig = df_sig[df_sig["Ds_isSignal"] == 1]
    if apply_veto:
        df_sig = df_sig[~((df_sig["Ds_diff_D0pi"] > dstar_veto_min) &
                          (df_sig["Ds_diff_D0pi"] <= dstar_veto_max))]
    sig_counts, edges = hist_dm(df_sig["Ds_massDifference_0"])

    # Backgrounds
    df_bkg_all = DataFrames[f"All_{mode}"]
    df_bkg_all = apply_base_cuts(df_bkg_all)
    if apply_veto:
        df_bkg_all = df_bkg_all[~((df_bkg_all["Ds_diff_D0pi"] > dstar_veto_min) &
                                  (df_bkg_all["Ds_diff_D0pi"] <= dstar_veto_max))]

    df_Dsp  = df_bkg_all[abs(df_bkg_all["Ds_mcPDG"]) == 413]
    df_Ds0  = df_bkg_all[abs(df_bkg_all["Ds_mcPDG"]) == 423]
    df_comb = df_bkg_all[(abs(df_bkg_all["Ds_mcPDG"]) != 413) &
                         (abs(df_bkg_all["Ds_mcPDG"]) != 423)]

    Dsp_counts, _  = hist_dm(df_Dsp["Ds_massDifference_0"])
    Ds0_counts, _  = hist_dm(df_Ds0["Ds_massDifference_0"])
    comb_counts, _ = hist_dm(df_comb["Ds_massDifference_0"])

    return {
        "edges": edges,
        "signal": sig_counts,
        "DstarPlus": Dsp_counts,
        "Dstar0": Ds0_counts,
        "comb": comb_counts,
    }

# Example: combine all modes into one channel by summing histograms
def build_combined_templates(apply_veto):
    templates = None
    for mode in decay_config.keys():
        t = make_templates_for_mode(mode, apply_veto)
        if templates is None:
            templates = {k: v.copy() for k, v in t.items()}
        else:
            for k in ["signal", "DstarPlus", "Dstar0", "comb"]:
                templates[k] += t[k]
    return templates

templates_no_veto = build_combined_templates(apply_veto=False)
templates_with_veto = build_combined_templates(apply_veto=True)

# Build pyhf models

def build_pyhf_model(templates):
    # Convert to numpy arrays
    sig      = np.asarray(templates["signal"],     dtype=float)
    bkg_Dsp  = np.asarray(templates["DstarPlus"],  dtype=float)
    bkg_Ds0  = np.asarray(templates["Dstar0"],     dtype=float)
    bkg_comb = np.asarray(templates["comb"],       dtype=float)

    # Single background template = sum of all three
    bkg_total = bkg_Dsp + bkg_Ds0 + bkg_comb

    # Convert to plain lists of numbers for pyhf spec
    signal = sig.tolist()
    bkg    = bkg_total.tolist()
    # 10% relative uncertainty per bin
    bkg_uncerts = (0.1 * np.ones_like(bkg_total, dtype=float)).tolist()

    model = pyhf.simplemodels.uncorrelated_background(
        signal=signal,
        bkg=bkg,
        bkg_uncertainty=bkg_uncerts,
        # signal_modifier="mu",  # only needed if your pyhf version supports it and you care about the name
    )
    return model

model_no_veto = build_pyhf_model(templates_no_veto)
model_with_veto = build_pyhf_model(templates_with_veto)

# Asimov background-only data
asimov_no_veto = model_no_veto.expected_data(0.0)
asimov_with_veto = model_with_veto.expected_data(0.0)

# Compute expected 95% CL upper limit on mu using CLs
def expected_mu_up(model, asimov_data):
    scan = np.linspace(0.0, 10.0, 101)
    best_mu = None
    for mu in scan:
        cls = pyhf.infer.hypotest(
            mu,
            asimov_data,
            model,
            test_stat="qtilde",
            return_expected=True,
        )
        # hypotest can return a scalar or array; make it a float
        cls = float(cls[0][2]) if isinstance(cls, (list, tuple, np.ndarray)) else float(cls)
        if cls < 0.05:
            best_mu = mu
            break
    return best_mu

mu_up_no_veto = expected_mu_up(model_no_veto, asimov_no_veto)
mu_up_with_veto = expected_mu_up(model_with_veto, asimov_with_veto)

print("\npyhf Asimov expected UL comparison (mu is signal strength):")
print(f"  no veto   : mu_up ≈ {mu_up_no_veto}")
print(f"  with veto : mu_up ≈ {mu_up_with_veto}")
if mu_up_with_veto is not None and mu_up_no_veto is not None:
    if mu_up_with_veto < mu_up_no_veto:
        print("=> With current MC, veto improves expected UL.")
    else:
        print("=> With current MC, veto does not improve expected UL.")
else:
    print("Scan did not cross CLs = 0.05 in one of the cases. Extend the scan range.")
