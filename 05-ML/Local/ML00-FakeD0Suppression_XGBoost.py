# %%
import os

os.cpu_count()

# %%
import psutil

mem = psutil.virtual_memory()
print(f"Total RAM: {mem.total / 1e9:.2f} GB")
print(f"Available: {mem.available / 1e9:.2f} GB")

# %%
# !pip uninstall -y scikit-learn
# !pip install scikit-learn==1.3.1

# %%
# ! pip install --upgrade pip
# ! pip install --user xgboost seaborn
# ! pip install --user bayesian-optimization

# %%
# import mplhep
import sys

import seaborn as sns

import numpy as np
import pandas as pd
import uproot
from matplotlib import pyplot as plt

from sklearn.datasets import make_classification,make_regression
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.metrics import auc,roc_curve,confusion_matrix,classification_report,precision_recall_curve,mean_squared_error,accuracy_score,roc_auc_score
from sklearn.model_selection import GridSearchCV, cross_validate, validation_curve,train_test_split,KFold,learning_curve,cross_val_score
from sklearn.utils import compute_sample_weight
from scipy.stats import ks_2samp

import xgboost
from xgboost import XGBClassifier

from sklearn.model_selection import RandomizedSearchCV

from scipy.stats import randint, uniform

# %%
plt.rcParams.update({
    "axes.labelsize": 16,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 14,
    "figure.titlesize": 20
})

# %%
pd.set_option('display.max_rows', 200000)
pd.set_option('display.max_columns', 200000)

# %%
sys.path.append("/home/belle2/amubarak/Ds2D0enue_Analysis/07-Python_Functions/")

# %% [markdown]
# # Prep-Work

# %% [markdown]
# ### Import Data

# %% [markdown]
# Correct Charge

# %%
# # In this notebook we only process the main signal and the generic events,
# # for illustration purposes.
# # You can add other backgrounds after if you wish.
# samples = ["Signal","All","BB","ccbar","ddbar","ssbar","taupair","uubar","uds"]
# GenEvents = ["Signal","BB","ccbar","ddbar","ssbar","taupair","uubar"]

# DataFrames = {}  # define empty dictionary to hold dataframes

# # Signal:
# DataFrames[samples[0]] =  uproot.concatenate("/home/belle2/amubarak/C01-Simulated_Events/Ds2D0enu-Signal.root:Dstree",library='pd')
# # Background
# for s in samples[1:]: # loop over samples
#     DataFrames[s] =  uproot.concatenate("/group/belle2/users2022/amubarak/TopoAna/Completed_TopoAna/TopoAna_"+ s +".root:Dstree",library='pd')

# %%
import os
import pandas as pd
import uproot
from tqdm import tqdm

# === Load only selected branches ===
with open("/home/belle2/amubarak/Ds2D0enue_Analysis/03-Grid/Save_var.txt") as f:
    variables_to_load = [
        line.strip().strip(",").strip('"').strip("'")
        for line in f
        if line.strip() and not line.strip().startswith("#")
    ]

samples = ["Signal", "BB", "ccbar", "ddbar", "ssbar", "taupair", "uubar"]
GenEvents = ["Signal","BB","ccbar","ddbar","ssbar","taupair","uubar"]
Date = "0530"
Attempt = "0"

DataFrames = {}

# === Load each sample one by one with progress bar ===
for name in tqdm(samples, desc="Loading samples"):
    if name == "Signal":
        path = "/home/belle2/amubarak/C01-Simulated_Events/Ds2D0enu-Signal.root:Dstree"
    else:
        path = f"/group/belle/users/amubarak/02-Grid/Sample_Grid/Ds2D0e-Generic_Ds_{Date}25_{Attempt}_{name}.root:Dstree"

    try:
        df = uproot.concatenate(path, filter_name=variables_to_load, library='pd')
        print(f"✔️ Loaded {name} with {len(df):,} entries")
        DataFrames[name] = df
    except Exception as e:
        print(f"❌ Failed to load {name}: {e}")
        DataFrames[name] = pd.DataFrame()

# === Merge background categories ===
background_samples = ["BB", "ccbar", "ddbar", "ssbar", "taupair", "uubar"]
DataFrames["All"] = pd.concat([DataFrames[s] for s in background_samples], ignore_index=True)

# === Combine uds light-quark backgrounds ===
DataFrames["uds"] = pd.concat(
    [DataFrames["uubar"], DataFrames["ddbar"], DataFrames["ssbar"]],
    ignore_index=True
)

# %% [markdown]
# Incorrect Charge

# %%
import os
import pandas as pd
import uproot
from tqdm import tqdm

# === Load only selected branches ===
with open("/home/belle2/amubarak/Ds2D0enue_Analysis/03-Grid/Save_var.txt") as f:
    variables_to_load = [
        line.strip().strip(",").strip('"').strip("'")
        for line in f
        if line.strip() and not line.strip().startswith("#")
    ]

# === Configuration ===
samples_WCh = [
    "Signal_WCh", "BB_WCh", "ccbar_WCh", "ddbar_WCh",
    "ssbar_WCh", "taupair_WCh", "uubar_WCh", "Data_WCh"
]
Date_WCh = "0630"
Attempt_WCh = "0"

# === Load one sample at a time ===
for sample in tqdm(samples_WCh, desc="Loading WCh samples"):
    if sample == "Signal_WCh":
        path = "/home/belle2/amubarak/C01-Simulated_Events/Ds2D0enu-Signal_WCh.root:Dstree"
    else:
        path = f"/group/belle/users/amubarak/02-Grid/Sample_Grid_WCh/Ds2D0e-Generic_Ds_{Date_WCh}25_{Attempt_WCh}_{sample}.root:Dstree"

    try:
        df = uproot.concatenate(path, filter_name=variables_to_load, library='pd')
        DataFrames[sample] = df
        print(f"✔️ Loaded: {path} [{len(df):,} entries]")
    except Exception as e:
        print(f"❌ Failed: {sample} — {e}")

# %% [markdown]
# Reversed PID

# %%
import os
import glob
import uproot
import pandas as pd
from tqdm import tqdm

# === Configuration ===
Date_ReverseID = "0626"
Attempt_ReverseID = "0"

# === Prompt user for veto toggle ===
apply_veto = input("Apply veto cut on Ds_diff_D0pi? (y/n): ").strip().lower() == "y"

# === Define veto window ===
cut_low = 0.14543 - (3*0.00041124)
cut_high = 0.14543 + (3*0.00041124)

# === Variables to load ===
with open("/home/belle2/amubarak/Ds2D0enue_Analysis/03-Grid/Save_var.txt") as f:
    variables_to_load = [
        line.strip().strip(",").strip('"').strip("'")
        for line in f if line.strip() and not line.strip().startswith("#")
    ]

if apply_veto and "Ds_diff_D0pi" not in variables_to_load:
    variables_to_load.append("Ds_diff_D0pi")

# === Load merged background and signal samples ===
merged_samples = {
    "Signal_ReverseID": "/home/belle2/amubarak/C01-Simulated_Events/Ds2D0enu-Signal_ReverseID.root",
    "BB_ReverseID": f"/group/belle/users/amubarak/02-Grid/Sample_Grid_ReverseID/Ds2D0e-Generic_Ds_{Date_ReverseID}25_{Attempt_ReverseID}_BB_ReverseID.root",
    "ddbar_ReverseID": f"/group/belle/users/amubarak/02-Grid/Sample_Grid_ReverseID/Ds2D0e-Generic_Ds_{Date_ReverseID}25_{Attempt_ReverseID}_ddbar_ReverseID.root",
    "ssbar_ReverseID": f"/group/belle/users/amubarak/02-Grid/Sample_Grid_ReverseID/Ds2D0e-Generic_Ds_{Date_ReverseID}25_{Attempt_ReverseID}_ssbar_ReverseID.root",
    "taupair_ReverseID": f"/group/belle/users/amubarak/02-Grid/Sample_Grid_ReverseID/Ds2D0e-Generic_Ds_{Date_ReverseID}25_{Attempt_ReverseID}_taupair_ReverseID.root",
    "uubar_ReverseID": f"/group/belle/users/amubarak/02-Grid/Sample_Grid_ReverseID/Ds2D0e-Generic_Ds_{Date_ReverseID}25_{Attempt_ReverseID}_uubar_ReverseID.root",
    "Data_ReverseID": f"/group/belle/users/amubarak/02-Grid/Sample_Grid_ReverseID/Ds2D0e-Generic_Ds_{Date_ReverseID}25_{Attempt_ReverseID}_Data_ReverseID.root",
}

print("\n📦 Loading merged samples...")
for sample, path in tqdm(merged_samples.items(), desc="Merged Samples"):
    try:
        df = uproot.concatenate(f"{path}:Dstree", filter_name=variables_to_load, library="pd")
        if apply_veto:
            df = df[(df["Ds_diff_D0pi"] <= cut_low) | (df["Ds_diff_D0pi"] >= cut_high)]
        DataFrames[sample] = df
        print(f"✔️ {sample}: {len(df):,} entries")
    except Exception as e:
        print(f"❌ Failed to load {sample}: {e}")

# === Load ccbar_ReverseID chunks sequentially ===
chunk_dir = f"/group/belle/users/amubarak/02-Grid/Sample_Grid_ReverseID/Ds_{Date_ReverseID}25_{Attempt_ReverseID}_ccbar_ReverseID_Chunks"
chunk_paths = [os.path.join(chunk_dir, f"ccbar_chunk_{i:02d}.root") for i in range(40)]

DataFrames["ccbar_ReverseID"] = []

print("\n🧱 Loading ccbar_ReverseID chunks one by one...")
for path in tqdm(chunk_paths, desc="ccbar Chunks"):
    try:
        df = uproot.concatenate(f"{path}:Dstree", filter_name=variables_to_load, library="pd")
        if apply_veto:
            df = df[(df["Ds_diff_D0pi"] <= cut_low) | (df["Ds_diff_D0pi"] >= cut_high)]
        DataFrames["ccbar_ReverseID"].append(df)
    except Exception as e:
        print(f"❌ Error loading {path}: {e}")

# === Concatenate all loaded DataFrames ===
if DataFrames["ccbar_ReverseID"]:
    DataFrames["ccbar_ReverseID"] = pd.concat(DataFrames["ccbar_ReverseID"], ignore_index=True)
    print(f"✅ ccbar_ReverseID: {len(DataFrames['ccbar_ReverseID']):,} entries")
else:
    print("❌ ccbar_ReverseID failed to load any chunks.")

# %% [markdown]
# Reverse PID and Wrong Charge

# %%
import os
import glob
import uproot
import pandas as pd
from tqdm import tqdm

# === Configuration ===
Date_ReverseID = "0708"
Attempt_ReverseID = "0"

# === Prompt user for veto toggle ===
apply_veto = input("Apply veto cut on Ds_diff_D0pi? (y/n): ").strip().lower() == "y"

# === Define veto window ===
cut_low = 0.14543 - (3*0.00041124)
cut_high = 0.14543 + (3*0.00041124)

# === Variables to load ===
with open("/home/belle2/amubarak/Ds2D0enue_Analysis/03-Grid/Save_var.txt") as f:
    variables_to_load = [
        line.strip().strip(",").strip('"').strip("'")
        for line in f if line.strip() and not line.strip().startswith("#")
    ]

if apply_veto and "Ds_diff_D0pi" not in variables_to_load:
    variables_to_load.append("Ds_diff_D0pi")

# === Load merged background and signal samples ===
merged_samples = {
    "BB_ReverseID_WCh": f"/group/belle/users/amubarak/02-Grid/Sample_Grid_ReverseID_WCh/Ds2D0e-Generic_Ds_{Date_ReverseID}25_{Attempt_ReverseID}_BB_ReverseID_WCh.root",
    "ccbar_ReverseID_WCh": f"/group/belle/users/amubarak/02-Grid/Sample_Grid_ReverseID_WCh/Ds2D0e-Generic_Ds_{Date_ReverseID}25_{Attempt_ReverseID}_ccbar_ReverseID_WCh.root",
    "ddbar_ReverseID_WCh": f"/group/belle/users/amubarak/02-Grid/Sample_Grid_ReverseID_WCh/Ds2D0e-Generic_Ds_{Date_ReverseID}25_{Attempt_ReverseID}_ddbar_ReverseID_WCh.root",
    "ssbar_ReverseID_WCh": f"/group/belle/users/amubarak/02-Grid/Sample_Grid_ReverseID_WCh/Ds2D0e-Generic_Ds_{Date_ReverseID}25_{Attempt_ReverseID}_ssbar_ReverseID_WCh.root",
    "taupair_ReverseID_WCh": f"/group/belle/users/amubarak/02-Grid/Sample_Grid_ReverseID_WCh/Ds2D0e-Generic_Ds_{Date_ReverseID}25_{Attempt_ReverseID}_taupair_ReverseID_WCh.root",
    "uubar_ReverseID_WCh": f"/group/belle/users/amubarak/02-Grid/Sample_Grid_ReverseID_WCh/Ds2D0e-Generic_Ds_{Date_ReverseID}25_{Attempt_ReverseID}_uubar_ReverseID_WCh.root",
    "Data_ReverseID_WCh": f"/group/belle/users/amubarak/02-Grid/Sample_Grid_ReverseID_WCh/Ds2D0e-Generic_Ds_{Date_ReverseID}25_{Attempt_ReverseID}_Data_ReverseID_WCh.root",
}

print("\n📦 Loading merged samples...")
for sample, path in tqdm(merged_samples.items(), desc="Merged Samples"):
    try:
        df = uproot.concatenate(f"{path}:Dstree", filter_name=variables_to_load, library="pd")
        if apply_veto:
            df = df[(df["Ds_diff_D0pi"] <= cut_low) | (df["Ds_diff_D0pi"] >= cut_high)]
        DataFrames[sample] = df
        print(f"✔️ {sample}: {len(df):,} entries")
    except Exception as e:
        print(f"❌ Failed to load {sample}: {e}")

# # === Load ccbar_ReverseID_WCh chunks sequentially ===
# chunk_dir = f"/group/belle/users/amubarak/02-Grid/Sample_Grid_ReverseID_WCh/Ds_{Date_ReverseID}25_{Attempt_ReverseID}_ccbar_ReverseID_WCh_Chunks"
# chunk_paths = [os.path.join(chunk_dir, f"ccbar_chunk_{i:02d}.root") for i in range(40)]

# DataFrames["ccbar_ReverseID_WCh"] = []

# print("\n🧱 Loading ccbar_ReverseID chunks one by one...")
# for path in tqdm(chunk_paths, desc="ccbar Chunks"):
#     try:
#         df = uproot.concatenate(f"{path}:Dstree", filter_name=variables_to_load, library="pd")
#         if apply_veto:
#             df = df[(df["Ds_diff_D0pi"] <= cut_low) | (df["Ds_diff_D0pi"] >= cut_high)]
#         DataFrames["ccbar_ReverseID_WCh"].append(df)
#     except Exception as e:
#         print(f"❌ Error loading {path}: {e}")

# # === Concatenate all loaded DataFrames ===
# if DataFrames["ccbar_ReverseID_WCh"]:
#     DataFrames["ccbar_ReverseID_WCh"] = pd.concat(DataFrames["ccbar_ReverseID_WCh"], ignore_index=True)
#     print(f"✅ ccbar_ReverseID_WCh: {len(DataFrames['ccbar_ReverseID_WCh']):,} entries")
# else:
#     print("❌ ccbar_ReverseID_WCh failed to load any chunks.")

# %% [markdown]
# The line below is to look at the available variables.

# %%
print(DataFrames.keys())

# %%
DataFrames["All"].columns.tolist()

# %% [markdown]
# ### Setup
# The code below will be used to apply cuts to the data.  
# The range of the plots.

# %%
# Electron ID
#-------------------
# DataFrames["Signal"] = DataFrames["Signal"][DataFrames["Signal"]['e_electronID']>=0.95]
# DataFrames["ccbar"] = DataFrames["ccbar"][DataFrames["ccbar"]['e_electronID']>=0.95]
# DataFrames["Signal"] = DataFrames["Signal"][DataFrames["Signal"]['Ds_gammaveto_em_electronID']>=0.95]
# DataFrames["ccbar"] = DataFrames["ccbar"][DataFrames["ccbar"]['Ds_gammaveto_em_electronID']>=0.95]

# Photon Conversion
#-------------------
# DataFrames[samples[0]] = DataFrames[samples[0]][DataFrames[samples[0]]['Ds_gammaveto_M_Correction']>=0.1]
# DataFrames[samples[1]] = DataFrames[samples[1]][DataFrames[samples[1]]['Ds_gammaveto_M_Correction']>=0.1]

# Peaking Background Removal
#----------------------------
# DataFrames["ccbar"] = DataFrames["ccbar"][(DataFrames["ccbar"]['Ds_diff_D0pi']>=0.15)]
# DataFrames["Signal"] = DataFrames["Signal"][(DataFrames["Signal"]['Ds_diff_D0pi']>=0.15)]

# # Vertex Fitting
# #----------------
# DataFrames["Signal"] = DataFrames["Signal"][DataFrames["Signal"]['Ds_chiProb']>=0.01]
# DataFrames["ccbar"] = DataFrames["ccbar"][DataFrames["ccbar"]['Ds_chiProb']>=0.01]

# Dalitz Removal
#----------------------------
# DataFrames["ccbar"] = DataFrames["ccbar"][(DataFrames["ccbar"]['Ds_pi0veto_M_Correction']<=0.08) | (DataFrames["ccbar"]['Ds_pi0veto_M_Correction']>=0.16)]
# DataFrames["Signal"] = DataFrames["Signal"][(DataFrames["Signal"]['Ds_pi0veto_M_Correction']<=0.08) | (DataFrames["Signal"]['Ds_pi0veto_M_Correction']>=0.16)]

# Vertex Fit
#----------------
# DataFrames[samples[0]] = DataFrames[samples[0]][DataFrames[samples[0]]['Ds_chiProb_rank']==1]
# DataFrames[samples[1]] = DataFrames[samples[1]][DataFrames[samples[1]]['Ds_chiProb_rank']==1]

# D0 Invariant Mass
#-----------------------
# DataFrames[samples[0]] = DataFrames[samples[0]][(DataFrames[samples[0]]['Ds_D0_sideband']==1)]
# DataFrames[samples[1]] = DataFrames[samples[1]][(DataFrames[samples[1]]['Ds_D0_sideband']==1)]

# %% [markdown]
# ## Photon Conversion Veto

# %%
for key in DataFrames.keys():
    df = DataFrames[key]
    if "Ds_gammaveto_M_Correction" in df.columns:
        DataFrames[key] = df[
            (df["Ds_gammaveto_M_Correction"] >= 0.1)
        ]

# %% [markdown]
# ## $D^{*+}$ Veto

# %%
cut_low = 0.14541 - (3*0.00039706)
cut_high = 0.14541 + (3*0.00042495)

for key in DataFrames.keys():
    df = DataFrames[key]
    if "Ds_diff_D0pi" in df.columns:
        DataFrames[key] = df[
            (df["Ds_diff_D0pi"] <= cut_low) | (df["Ds_diff_D0pi"] >= cut_high)
        ]

# %% [markdown]
# # Fake $D^0$ Suppression

# %%
# # === Updated Variables ===
# Variables = [
#     'pi_dr',
# #     'pi_dz',
#     'K_dr',
# #     'K_dz',
#     'D0_dM',
#     'D0_chiProb',
#     'D0_flightDistance',
#     'D0_flightTime',
#     'D0_useCMSFrame_p',
#     'D0_cos_decayAngle_1',
# ]

# features = [
#     r'$dr(\pi^{+})\;[\mathrm{cm}]$',
# #     r'$dz(\pi^{+})\;[\mathrm{cm}]$',
#     r'$dr(K^{-})\;[\mathrm{cm}]$',
# #     r'$dz(K^{-})\;[\mathrm{cm}]$',
#     r'$m(D^{0}) - m_{\mathrm{PDG}}(D^{0})\;[\mathrm{GeV}/c^{2}]$',
#     r'p-value of $D^{0}$',
#     r'$Flight\;Distance(D^{0})\;[\mathrm{cm}]$',
#     r'$Flight\;Time(D^{0})\;[\mathrm{ns}]$',
#     r'$p^{*}(D^{0})\;[\mathrm{GeV}/c]$',
#     r'$\cos\theta^*_{\mathrm{daughter}_1}$',
# ]

# ranges = {
#     'pi_dr': [0, 0.1],
# #     'pi_dz': [-0.5, 0.5],
#     'K_dr': [0, 0.1],
# #     'K_dz': [-0.5, 0.5],
#     'D0_dM': [-0.02, 0.02],
#     'D0_chiProb': [0, 1],
#     'D0_flightDistance': [-0.4, 0.4],
#     'D0_flightTime': [-0.005, 0.005],
#     'D0_useCMSFrame_p': [2.5, 5.0],
#     'D0_cos_decayAngle_1': [-1, 1],
# }

# bins = 50
# density = True
# samples = "All"

# # === Fixed Color Scheme ===
# colors = {
#     'signal': '#007C91',   # Real signal
#     'other': '#2E2E2E',    # Everything else
# }

# bg_labels = [
#     r'$Other$',
#     r'$D^{0}$',
#     r'$D^{*0} \rightarrow D^{0} \; \pi^{0} / \gamma$',
#     r'$D^{*+} \rightarrow D^{0} \; \pi^{+}$'
# ]

# bg_masks = [
#     lambda df: df['D0_isSignal'].isna() | (df['D0_isSignal'] == 0),
#     lambda df: (df['Ds_D0_NoDstarplusDstar0'] == 1) & (df['D0_isSignal'] == 1),
#     lambda df: (df['Ds_D0_Dstar0'] == 1) & (df['D0_isSignal'] == 1),
#     lambda df: (df['Ds_D0_Dstarplus'] == 1) & (df['D0_isSignal'] == 1),
# ]

# # === Plotting ===
# if "All" in DataFrames and "Signal" in DataFrames:
#     for idx, (var, label) in enumerate(zip(Variables, features)):
#         var_range = ranges[var]
#         bin_width = (var_range[1] - var_range[0]) / bins

#         real_signal_data = DataFrames["Signal"][DataFrames["Signal"]['Ds_isSignal'] == 1][var]

#         for jdx, (mask, bg_label) in enumerate(zip(bg_masks, bg_labels)):
#             bg_data = DataFrames[samples][mask(DataFrames[samples])][var]

#             plt.hist(real_signal_data, label="Real Signal", histtype='step', density=density,
#                      bins=bins, alpha=1, range=var_range, linewidth=2, color=colors['signal'])

#             plt.hist(bg_data, label=bg_label, histtype='step', density=density,
#                      bins=bins, alpha=1, range=var_range, linewidth=2, color=colors['other'])

#             plt.xlabel(label)
#             if bin_width < 0.01:
#                 exponent = int(np.floor(np.log10(bin_width)))
#                 base = bin_width / (10**exponent)
#                 ylabel = r'$Norm.\;Entries/({:.2f} \times 10^{{{}}})$'.format(base, exponent)
#             else:
#                 ylabel = r'$Norm.\;Entries/({:.2f})$'.format(bin_width)

#             plt.ylabel(ylabel)
#             plt.legend(loc='upper right')
# #             plt.title(f"Real Signal vs {bg_label}", fontsize=15)
#             plt.show()

#         # Real vs Fake Signal
#         fake_signal = DataFrames["Signal"][DataFrames["Signal"]['Ds_isSignal'] == 0][var]

#         plt.hist(real_signal_data, label="Real Signal", histtype='step', density=density,
#                  bins=bins, alpha=1, range=var_range, linewidth=2, color=colors['signal'])

#         plt.hist(fake_signal, label="Fake Signal", histtype='step', density=density,
#                  bins=bins, alpha=1, range=var_range, linewidth=2, color=colors['other'])

#         plt.xlabel(label)
#         plt.ylabel(r'$Normalized\;Entries/({:.2f})$'.format(bin_width))
#         plt.legend(loc='upper right')
#         plt.title("Real vs Fake Signal", fontsize=15)
#         plt.show()
# else:
#     print("DataFrames['All'] and DataFrames['Signal'] must be defined.")

# %%
# === Input Variables for the BDT ===
Variables = [
    'K_dr',
    'pi_dr',
    'K_kaonID',
    'pi_pionID',
    'D0_dM',
    'D0_chiProb',
    'D0_flightDistance',
    'D0_useCMSFrame_p',
    'D0_cos_decayAngle_1',
]

features = [
    r'$dr(K^{-})\;[\mathrm{cm}]$',
    r'$dr(\pi^{+})\;[\mathrm{cm}]$',
    r'$kaonID(K^{-})$',  # unitless PID likelihood
    r'$pionID(\pi^{+})$',  # unitless PID likelihood
    r'$m(D^{0}) - m_{PDG}(D^{0})\;[\mathrm{GeV}/c^{2}]$',
    r'$p$-value$(D^{0})$',  # unitless
    r'$Flight \; Distance(D^{0})\;[\mathrm{cm}]$',
    r'$p^{*} (D^{0})\;[\mathrm{GeV}/c]$',
    r'$\cos\theta^*_{daughter_1}$',  # unitless angle cosine
]

# === Plot Ranges ===
ranges = {
    'K_dr': [0, 0.08],
    'pi_dr': [0, 0.08],
    'K_kaonID': [0.5, 1],
    'pi_pionID': [0.2, 1],
    'D0_dM': [-0.02, 0.02],
    'D0_chiProb': [0, 1],
    'D0_flightDistance': [-0.4, 0.4],
    'D0_useCMSFrame_p': [2.5, 5.0],
    'D0_cos_decayAngle_1': [-1, 1],
}

bins = 50
density = True

# === Scientific Colors ===
colors = {
    'signal': '#007C91',  # Blue for real signal
    'fake': '#C44E52',    # Red for fake D0
}

# === Extract Samples ===
df_signal = DataFrames["Signal"]
df_generic = DataFrames["All"]

df_true_signal = df_signal[
    (df_signal["Ds_isSignal"] == 1) & (df_signal["D0_isSignal"] == 1)
]
df_fake_d0 = df_generic[
    (df_generic["D0_isSignal"] == 0) | (df_generic["D0_isSignal"].isna())
]

# === Plotting ===
for var, label in zip(Variables, features):
    if var not in ranges:
        print(f"Skipping {var}: no range defined.")
        continue

    var_range = ranges[var]
    bin_width = (var_range[1] - var_range[0]) / bins

    signal_data = df_true_signal[var].dropna()
    fake_data = df_fake_d0[var].dropna()

    plt.hist(signal_data, label="Real $D^0$ (Signal MC)",
             histtype='step', density=density,
             bins=bins, range=var_range, linewidth=2, color=colors['signal'])

    plt.hist(fake_data, label="Fake $D^0$ (Generic MC)",
             histtype='step', density=density,
             bins=bins, range=var_range, linewidth=2, color=colors['fake'])

    plt.xlabel(label)

    # if bin_width < 0.01:
    #     exponent = int(np.floor(np.log10(bin_width)))
    #     base = bin_width / (10**exponent)
    #     ylabel = r'$Norm.\;Entries/({:.2f} \times 10^{{{}}})$'.format(base, exponent)
    # else:
    ylabel = r'$Norm.\;Entries/({:.3f})$'.format(bin_width)

    plt.ylabel(ylabel)
    plt.legend(loc='upper right')
#     plt.title("Real vs Fake $D^0$", fontsize=15)
    plt.show()


# %%
DataFrames["All"].isna().sum()

# %%
DataFrames["All"]["D0_isSignal"] = DataFrames["All"]["D0_isSignal"].replace(np.nan, 0)

for s in GenEvents[0:]: # loop over samples
    DataFrames[s]["D0_isSignal"] = DataFrames[s]["D0_isSignal"].replace(np.nan, 0)

# %%
DataFrames["All"]["Ds_isSignal"] = DataFrames["All"]["Ds_isSignal"].replace(np.nan, 0)

for s in GenEvents[0:]: # loop over samples
    DataFrames[s]["Ds_isSignal"] = DataFrames[s]["Ds_isSignal"].replace(np.nan, 0)

# %%
Variables = [
             'K_dr',
             'pi_dr',
             'D0_significanceOfDistance',
             'D0_chiProb',
             'D0_flightDistance',
             'D0_useCMSFrame_p',
             'D0_cos_decayAngle_1',
             ]

features = [
             r'$dr(K^{-})$',
             r'$dr(\pi^{+})$',
             'D0_significanceOfDistance',
             # r'$m(D^{0}) - m_{PDG}(D^{0})$',
             r'$p-value(D^{0})$',
             r'$Flight \; Distance(D^{0})$',
             r'$p^{*} (D^{0})$',
             r'$\cos\theta^*_{daughter_1}$',
             ]

# %%
plt.figure(figsize=(18, 15))

heatmap = sns.heatmap(DataFrames["Signal"][Variables].corr(), annot=True, cmap="coolwarm",vmin=-1, vmax=1)

heatmap.set_title('Signal Correlation Heatmap', fontdict={'fontsize':20}, pad=16)

# %%
plt.figure(figsize=(18, 15))

heatmap = sns.heatmap(DataFrames["All"][Variables].corr(), annot=True, cmap="coolwarm",vmin=-1, vmax=1)

heatmap.set_title('Background Correlation Heatmap', fontdict={'fontsize':20}, pad=16)

# %% [markdown]
# ## Data Preprocessing

# %%
# # Define your features and labels from the 'All' dataset
# X = DataFrames["All"][Variables].to_numpy(dtype=np.float32)
# y = DataFrames["All"]['D0_isSignal'].to_numpy(dtype=np.int64)

# # # Reference variable for decorrelation — this is what uBoost will try to flatten for background
# # ref_variable = DataFrames["All"]["D0_dM"]

# #splitting with  Holdout method for eval_set
# X_train, X_test, y_train, y_test = train_test_split(X, y,
#                                                     test_size=0.30,
#                                                     random_state=42,
#                                                     # stratify=y
#                                                     )

# %%
# Signal: Real D⁰ from signal MC
df_signal = DataFrames["Signal"]
real_signal_mask = (df_signal["Ds_isSignal"] == 1) & (df_signal["D0_isSignal"] == 1)
df_true_signal = df_signal[real_signal_mask]

# Background: Fake D⁰ from generic MC
df_generic = DataFrames["All"]
df_fake_d0 = df_generic[(df_generic["D0_isSignal"] == 0) | (df_generic["D0_isSignal"].isna())]

# Combine
df_train = pd.concat([df_true_signal, df_fake_d0], axis=0)

# Labels
labels = np.concatenate([
    np.ones(len(df_true_signal), dtype=np.int64),   # Signal = 1
    np.zeros(len(df_fake_d0), dtype=np.int64)       # Fake D⁰ = 0
])

# Features
X = df_train[Variables].to_numpy(dtype=np.float32)
y = labels

#splitting with  Holdout method for eval_set
X_train, X_test, y_train, y_test = train_test_split(X, y,
                                                    test_size=0.30,
                                                    random_state=42,
                                                    # stratify=y
                                                    )

# %% [markdown]
# ## Model Training

# %%
weights = compute_sample_weight('balanced', y_train)

# Create EarlyStopping callback
early_stop = xgboost.callback.EarlyStopping(
    rounds=10,
    metric_name='rmse',
    data_name="validation_0",
    save_best=True,
)

# %%
eval_set = [(X_test, y_test)]
bdt = XGBClassifier(objective="binary:logistic",
                    eval_metric="logloss",
                    # early_stopping_rounds=10,
                    # scale_pos_weight=pos_class_weight,
                #     scale_pos_weight=scale,
                    max_delta_step=1,
                    random_state=42,
                    n_estimators=100)

bdt.fit(X_train, y_train, 
        eval_set=[(X_train, y_train),(X_test, y_test)], 
        sample_weight=weights,
        verbose=0) 

# %%
# loss curve of xgboost
results = bdt.evals_result()

plt.figure(figsize=(10,7))
plt.plot(results["validation_0"]["logloss"], label="Training loss")
plt.plot(results["validation_1"]["logloss"], label="Validation loss")
plt.xlabel("Number of trees")
plt.ylabel("Loss")
plt.legend()

# %% [markdown]
# ## Parameter Optimization
# This optimization is pulling too much resources and ending the connection

# %%
param_dist = {
    "learning_rate": uniform(0.01, 0.2),        # e.g., 0.01 to 0.21
    "max_depth": randint(1, 5),                 # 1 to 4
    "n_estimators": randint(100, 201),          # 100 to 200
    "reg_lambda": randint(1, 5),
    "gamma": randint(0, 4),
    "subsample": uniform(0.5, 0.5),             # 0.5 to 1.0
    "min_child_weight": randint(1, 6),
    "colsample_bytree": uniform(0.3, 0.7)
}

random_search = RandomizedSearchCV(
    bdt,
    param_distributions=param_dist,
    n_iter=50,  # Try only 50 random combos (you can adjust)
    cv=5,
    n_jobs=-1,
    random_state=42,
    verbose=1
)

# After running RandomizedSearchCV:
random_search.fit(X_train, y_train, eval_set=[(X_test, y_test)], sample_weight=weights, verbose=0)

# Extract best parameters and apply them to the base model
xgbm_final = bdt.set_params(**random_search.best_params_, random_state=17).fit(X_train, y_train, sample_weight=weights)

# %%
# cv_results = cross_validate(xgbm_final, X, y, cv=10,
#                             scoring=["f1"],return_train_score=True)

# print(cv_results['train_f1'].mean())
# print(cv_results['test_f1'].mean())

# %%
print("Best parameters found by RandomizedSearchCV:")
for param, value in random_search.best_params_.items():
    print(f"{param:20s}: {value}")

# %% [markdown]
# ## Feature Importance

# %%
# Get feature importance scores
print(xgbm_final.feature_importances_)

# %%
def plot_importance(model, features, num=len(X), save=False):
    feature_imp = pd.DataFrame({'Value': model.feature_importances_, 'Feature': features})
    plt.figure(figsize=(16, 8))
    sns.set(font_scale=1)
    sns.barplot(x="Value", y="Feature", data=feature_imp.sort_values(by="Value",
                                                                     ascending=False)[0:num])
    plt.title('Feature Importance', fontsize=20)
    plt.xlabel('Importance Value', fontsize=16)
    plt.ylabel('Features', fontsize=20)
    # ax.tick_params(axis='both', labelsize=14)  # Tick labels
    if save:
        plt.savefig('importances.png')

# %%
plot_importance(xgbm_final, features)
# When the feature importance graph is observed, 
# it is seen that the variables other than a02 and a01 are important for the xgboost model.

# %% [markdown]
# ## Overfitting Check

# %%
from scipy import stats
def get_pulls(counts,errors,pdf):
    pull = (-pdf + counts) / errors
    return pull

# %%
def compare_train_test(clf, X_train, y_train, X_test, y_test):
    Density = True
    decisions = [] # list to hold decisions of classifier
    for X,y in ((X_train, y_train), (X_test, y_test)): # train and test
        if hasattr(clf, "predict_proba"): # if predict_proba function exists
            d1 = clf.predict_proba(X[y<0.5])[:, 1] # background
            d2 = clf.predict_proba(X[y>0.5])[:, 1] # signal
        else: # predict_proba function doesn't exist
            X_tensor = torch.as_tensor(X, dtype=torch.float) # make tensor from X_test_scaled
            y_tensor = torch.as_tensor(y, dtype=torch.long) # make tensor from y_test
            X_var, y_var = Variable(X_tensor), Variable(y_tensor) # make variables from tensors
            d1 = clf(X_var[y_var<0.5])[1][:, 1].cpu().detach().numpy() # background
            d2 = clf(X_var[y_var>0.5])[1][:, 1].cpu().detach().numpy() # signal
        decisions += [d1, d2] # add to list of classifier decision

    #pd.set_option('max_columns', None)
#     %config InlineBackend.figure_format = 'retina'
    # plt.style.use('belle2')
    lw=3

    fig,axs=plt.subplots(3,1,figsize=(10,10),gridspec_kw={'height_ratios':[1,0.2,0.2]})

    bins = 50
    bin_edges = np.linspace(0,1,bins)
    
    test_bkg_count_weight=bins/len(decisions[2])
    test_sig_count_weight=bins/len(decisions[3])
    test_bkg_counts,test_bkg_bins = np.histogram(decisions[2],bins=bins,range=(0,1))
    test_sig_counts,test_sig_bins = np.histogram(decisions[3],bins=bins,range=(0,1))

    train_bkg_counts,train_bkg_bins,_etc=axs[0].hist(decisions[0],color = 'tab:blue',
            histtype='step',bins=bins,density=Density,range=(0,1),linewidth=lw,label='Train Background')
    train_sig_counts,train_sig_bins,_etc=axs[0].hist(decisions[1],color = 'tab:red',
            histtype='step',bins=bins,density=Density,range=(0,1),linewidth=lw,label=r'Train Signal')
    axs[0].hist(decisions[0],color = 'tab:blue',
            histtype='stepfilled',alpha=0.4,bins=bins,density=Density,range=(0,1))
    axs[0].hist(decisions[1],color = 'tab:red',
            histtype='stepfilled',alpha=0.4,bins=bins,density=Density,range=(0,1))
    bin_width=test_bkg_bins[1]-test_bkg_bins[0]
    bin_centers=[el+(bin_width/2) for el in test_bkg_bins[:-1]]

    axs[0].errorbar(bin_centers,test_bkg_count_weight*test_bkg_counts,
                yerr=test_bkg_count_weight*np.sqrt(test_bkg_counts),label='Test Background',color='tab:blue',
                marker='o',linewidth=lw,ls='')
    axs[0].errorbar(bin_centers,test_sig_count_weight*test_sig_counts,
                yerr=test_sig_count_weight*np.sqrt(test_sig_counts),label='Test Signal',color='tab:red',
                marker='o',linewidth=lw,ls='')
    axs[0].set_title(r'$D_{s}^{+} \rightarrow D^{0} e^{+} \nu_{e}$',loc='left')
    axs[0].set_xlim(0,1)
    axs[0].set_ylim(0)
    axs[0].set_ylabel('Event Density')

    x= decisions[1]
    y=  decisions[3]
    ks_p_value_sig = ks_2samp(x, y)[1]

    x= decisions[0]
    y= decisions[2]
    ks_p_value_bkg = ks_2samp(x, y)[1]

    leg=axs[0].legend(loc='upper center',title=f"Sig K-S test score: {ks_p_value_sig:0.3f}"+
                      "\n"+f"Bkg K-S test score: {ks_p_value_bkg:0.3f}")
    leg._legend_box.align = "left"  

    pulls=get_pulls(test_bkg_count_weight*test_bkg_counts,test_bkg_count_weight*np.sqrt(test_bkg_counts),np.array(train_bkg_counts))
    axs[1].bar(bin_centers,pulls,width=bin_width)
    axs[1].set_xlim(0,1)
    axs[1].set_ylabel('Pulls')
    axs[1].set_ylim(-5,5)

    pulls=get_pulls(test_sig_count_weight*test_sig_counts,test_sig_count_weight*np.sqrt(test_sig_counts),np.array(train_sig_counts))
    axs[2].bar(bin_centers,pulls,width=bin_width,color='tab:red')
    axs[2].set_xlim(0,1)
    axs[2].set_ylabel('Pulls')
    axs[2].set_ylim(-5,5)
    axs[2].set_xlabel(r'BDT output')

    return decisions

# %%
decisions = compare_train_test(xgbm_final, X_train, y_train, X_test, y_test)

# %% [markdown]
# ## Model Check

# %% [markdown]
# ### Basf2 ROC

# %%
y_score_test = xgbm_final.predict_proba(X_test)[:, 1]
fpr_test, tpr_test, thresholds_test = roc_curve(y_test, y_score_test)
area_test = auc(fpr_test, tpr_test)

y_score_train = xgbm_final.predict_proba(X_train)[:, 1]
fpr_train, tpr_train, thresholds_train = roc_curve(y_train, y_score_train)
area_train = auc(fpr_train, tpr_train)

# Get classifier scores (probabilities for class 1)
train_scores = xgbm_final.predict_proba(X_train)[:, 1]
test_scores  = xgbm_final.predict_proba(X_test)[:, 1]

# Use y_train and y_test to separate signal/background
sig_train = train_scores[y_train == 1]
bkg_train = train_scores[y_train == 0]
sig_test  = test_scores[y_test == 1]
bkg_test  = test_scores[y_test == 0]

# Optionally, group them into one list like this:
decisions = [bkg_train, sig_train, bkg_test, sig_test]

bdt_cuts = np.linspace(0, 1, 100)

sig_eff_train = []
bkg_rej_train = []
sig_eff_test = []
bkg_rej_test = []
fom_vals = []

for cut in bdt_cuts:
    num_sig_train = np.sum(sig_train > cut)
    num_bkg_train = np.sum(bkg_train > cut)
    num_sig_test = np.sum(sig_test > cut)
    num_bkg_test = np.sum(bkg_test > cut)

    # FoM calculation
    fom = num_sig_test / np.sqrt(num_sig_test + num_bkg_test) if (num_sig_test + num_bkg_test) > 0 else 0
    fom_vals.append(fom)

    sig_eff_train.append(num_sig_train / len(sig_train))
    bkg_rej_train.append(1 - (num_bkg_train / len(bkg_train)))
    sig_eff_test.append(num_sig_test / len(sig_test))
    bkg_rej_test.append(1 - (num_bkg_test / len(bkg_test)))

# Find optimal FoM point
fom_vals = np.array(fom_vals)
best_idx = np.argmax(fom_vals)
best_cut = bdt_cuts[best_idx]

# Plot
fig, axs = plt.subplots(1, 1, figsize=(7, 6))
lw = 2

# axs.plot([0, 1], [0, 1], color='grey', linestyle='--', label='Random')
axs.plot(bkg_rej_train, sig_eff_train, color='tab:blue', linewidth=lw, label=f'Train (AUC = {area_train:.2f})')
axs.plot(bkg_rej_test, sig_eff_test, color='tab:red', linestyle='--', linewidth=lw, label=f'Test (AUC = {area_test:.2f})')

# ① Shade the overfit gap
axs.fill_between(bkg_rej_test,
                 sig_eff_train,
                 sig_eff_test,
                 where=(np.array(sig_eff_train) > np.array(sig_eff_test)),
                 color='gray', alpha=0.2, label='Overfit Gap')

# ② Mark the optimal cut point (from test curve)
axs.axhline(sig_eff_test[best_idx], color='black', ls='--', linewidth=1.6)
# axs.axhline(sig_eff_test[best_idx], color='black', ls='--', linewidth=1.6,
#             label=f'Best FoM Cut = {best_cut:.3f}')
axs.axvline(bkg_rej_test[best_idx], color='black', ls='--', linewidth=1.6)
axs.scatter(bkg_rej_test[best_idx], sig_eff_test[best_idx], color='green', s=50)

# Axis labels and formatting
axs.set_title(r'$D_{s}^{+} \rightarrow D^{0} e^{+} \nu_{e}$', loc='left')
axs.set_ylim(0, 1.05)
axs.set_xlim(0, 1.05)
axs.set_xlabel('Background rejection')
axs.set_ylabel('Signal efficiency')
axs.legend(loc='lower left')
axs.grid(True)
plt.tight_layout()
plt.show()

# %% [markdown]
# ### Machine Learing ROC

# %%
y_score_test = xgbm_final.predict_proba(X_test)[:, 1]
fpr_test, tpr_test, thresholds_test = roc_curve(y_test, y_score_test)
area_test = auc(fpr_test, tpr_test)

y_score_train = xgbm_final.predict_proba(X_train)[:, 1]
fpr_train, tpr_train, thresholds_train = roc_curve(y_train, y_score_train)
area_train = auc(fpr_train, tpr_train)

plt.plot([0, 1], [0, 1], color='grey', linestyle='--')
plt.plot(fpr_test, tpr_test, label=f'Test ROC curve (AUC = {area_test:.2f})')
plt.plot(fpr_train, tpr_train, label=f'Train ROC curve (AUC = {area_train:.2f})')
plt.xlim(0.0, 1.0)
plt.ylim(0.0, 1.0)
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.legend(loc='lower right')
# We can make the plot look nicer by forcing the grid to be square
plt.gca().set_aspect('equal', adjustable='box')

# %%
# Make predictions on the test set
y_pred_proba = xgbm_final.predict_proba(X_test)[:, 1]

# Calculate the ROC AUC score
roc_auc = roc_auc_score(y_test, y_pred_proba)

print(f"ROC AUC Score: {roc_auc:.2f}")

# %% [markdown]
# ### Other Checks

# %%
# Predict on training and validation sets
train_preds = xgbm_final.predict(X_train)
val_preds = xgbm_final.predict(X_test)

# Calculate accuracy scores
train_accuracy = accuracy_score(y_train, train_preds)
val_accuracy = accuracy_score(y_test, val_preds)

print(f"Training Accuracy: {train_accuracy:.4f}")
print(f"Validation Accuracy: {val_accuracy:.4f}")

# Check for large difference between train and validation accuracy
if train_accuracy - val_accuracy > 0.1:
    print("Warning: The model may be overfitting!")

# %% [markdown]
# Check if XGBoost Is Underfitting

# %%
# Predict on training and validation sets
train_preds = xgbm_final.predict(X_train)
val_preds = xgbm_final.predict(X_test)

# Calculate MSE for training and validation sets
train_mse = mean_squared_error(y_train, train_preds)
val_mse = mean_squared_error(y_test, val_preds)

print(f"Training MSE: {train_mse:.4f}")
print(f"Validation MSE: {val_mse:.4f}")

# Check if both training and validation MSE are high
if train_mse > 100 and val_mse > 100:
    print("Warning: The model may be underfitting!")
    print("Consider increasing model complexity by adding more estimators, reducing learning rate, or adjusting other hyperparameters.")

# %% [markdown]
# ## BDT Cut Optimization

# %%
# Apply BDT to all DataFrames that contain the required Variables
for key in DataFrames.keys():
    df = DataFrames[key]
    
    # Check: make sure all input BDT variables exist in this DataFrame
    if all(var in df.columns for var in Variables):
        # Apply BDT and store the result
        DataFrames[key]["Ds_FakeD0BDT"] = xgbm_final.predict_proba(df[Variables])[:, 1].astype(np.float32)

# %%
def compute_fom_curve(scores, labels, weights=None, n_thresholds=200):
    """
    Compute FoM (S / sqrt(S + B)) across multiple BDT score thresholds.
    
    Parameters:
        scores (np.array): BDT scores for the validation/test set
        labels (np.array): True labels (1 for real D0, 0 for fake)
        weights (np.array): Optional per-event weights
        n_thresholds (int): Number of thresholds to scan (default=200)

    Returns:
        thresholds (np.array), foms (np.array), best_threshold (float), best_fom (float)
    """
    thresholds = np.linspace(0, 1, n_thresholds)
    foms = []

    for t in thresholds:
        mask = scores > t

        if weights is not None:
            S = np.sum(weights[(labels == 1) & mask])
            B = np.sum(weights[(labels == 0) & mask])
        else:
            S = np.sum((labels == 1) & mask)
            B = np.sum((labels == 0) & mask)

        fom = S / np.sqrt(S + B) if (S + B) > 0 else 0
        foms.append(fom)

    foms = np.array(foms)
    best_idx = np.argmax(foms)
    return thresholds, foms, thresholds[best_idx], foms[best_idx]

# %%
# Predict scores from your trained model
scores = xgbm_final.predict_proba(X_test)[:, 1]

# Optionally define weights (or leave as None)
weights = np.ones_like(y_test)  # or from your MC truth if applicable

# Compute FoM curve
thresholds, foms, best_thresh, best_fom = compute_fom_curve(scores, y_test)

# Print results
print(f"Best threshold: {best_thresh:.3f}")
print(f"Best FoM: {best_fom:.3f}")

# Plot it
plt.plot(thresholds, foms)
plt.axvline(best_thresh, color='red', linestyle='--', label=f'Best = {best_thresh:.3f}')
plt.axvspan(0,best_thresh,color='gray',alpha=0.2)
plt.xlabel("BDT Threshold")
plt.ylabel("FoM = S / √(S + B)")
plt.title("FoM Scan vs BDT Threshold")
plt.legend()
plt.grid(True)
plt.show()

# %%
from Functions import optimize_cut, plot_save

cut = optimize_cut(
    df_sig=DataFrames["Signal"],                  # DataFrame used for signal plotting
    df_bkg=DataFrames["All"],                     # DataFrame used for background plotting
    Signal=DataFrames["Signal"],                  # DataFrame used for signal FoM calculation
    Background=DataFrames["All"],                 # DataFrame used for background FoM calculation
    var="Ds_FakeD0BDT",                           # Variable to plot
    FoM="Ds_FakeD0BDT",                           # Variable to optimize over (can be same as var)
    xlabel="Classifier Output",                   # X-axis label
    Bins=50,
    Range=[0, 1],
    varmin=0,
    varmax=0.99,
    select="right",                               # "right" for >= cut, "left" for <= cut
    Width=False,
    query_signal="Ds_isSignal == 1"               # Only consider true signal
)

print(f"Best cut is: {cut:.3f}")

# %% [markdown]
# ## Fake $D^0$ BDT Cut

# %%
# DataFrames["All"] = DataFrames["All"][(DataFrames["All"]["Ds_FakeD0BDT"]>=0.556)]

# for s in GenEvents[0:]: # loop over samples
#     DataFrames[s] = DataFrames[s][(DataFrames[s]["Ds_FakeD0BDT"]>=0.556)]

# %% [markdown]
# # Save BDT Output

# %% [markdown]
# Correct Charge

# %%
# import os
# import uproot

# # === Make sure samples is a list ===
# samples = ["Signal", "BB", "ccbar", "ddbar", "ssbar", "taupair", "uubar"]

# # === Output directory ===
# output_dir = "/group/belle/users/amubarak/03-ML/FakeD0/"
# os.makedirs(output_dir, exist_ok=True)

# # === Base input path for original files ===
# base_input_dir = "/group/belle/users/amubarak/02-Grid/Sample_Grid"
# Date = "0530"
# Attempt = "0"

# # === Save each DataFrame using original filename with _withBDT suffix ===
# for s in samples:
#     if s not in DataFrames:
#         print(f"Warning: {s} not in DataFrames — skipping.")
#         continue

#     # Convert Fake D⁰ BDT output to float32 if it exists
#     if "Ds_FakeD0BDT" in DataFrames[s].columns:
#         DataFrames[s]["Ds_FakeD0BDT"] = DataFrames[s]["Ds_FakeD0BDT"].astype(np.float32)

#     # Set original file path
#     if s == "Signal":
#         original_name = "Ds2D0enu-Signal.root"
#     else:
#         original_name = f"Ds2D0e-Generic_Ds_{Date}25_{Attempt}_{s}.root"

#     # Build output path with _withBDT suffix
#     output_name = original_name.replace(".root", "_withBDT.root")
#     out_path = os.path.join(output_dir, output_name)

#     # Save DataFrame to ROOT
#     with uproot.recreate(out_path) as f:
#         f["Dstree"] = DataFrames[s]

#     print(f"Saved: {out_path}")

# %% [markdown]
# Wrong Charge

# %%
# import os
# import uproot

# # === Make sure wrong-charge samples list is defined ===
# samples_WCh = ["Signal_WCh", "BB_WCh", "ccbar_WCh", "ddbar_WCh", "ssbar_WCh", "taupair_WCh", "uubar_WCh", "Data_WCh"]

# # === Output directory for wrong charge ===
# output_dir_WCh = "/group/belle/users/amubarak/03-ML/FakeD0_WCh/"
# os.makedirs(output_dir_WCh, exist_ok=True)

# # === Base input path for original files (wrong charge) ===
# base_input_dir_WCh = "/group/belle/users/amubarak/02-Grid/Sample_Grid_WCh"
# Date_WCh = "0630"
# Attempt_WCh = "0"

# # === Save each wrong-charge DataFrame using original filename with _withBDT suffix ===
# for s in samples_WCh:
#     if s not in DataFrames:
#         print(f"Warning: {s} not in DataFrames — skipping.")
#         continue

#     # Convert Fake D⁰ BDT output to float32 if it exists
#     if "Ds_FakeD0BDT" in DataFrames[s].columns:
#         DataFrames[s]["Ds_FakeD0BDT"] = DataFrames[s]["Ds_FakeD0BDT"].astype(np.float32)

#     # Set original file name
#     if s == "Signal_WCh":
#         original_name = "Ds2D0enu-Signal_WCh.root"
#     else:
#         tag = s.replace("_WCh", "")
#         original_name = f"Ds2D0e-Generic_Ds_{Date_WCh}25_{Attempt_WCh}_{tag}.root"

#     # Build output path with _withBDT suffix
#     output_name = original_name.replace(".root", "_withBDT.root")
#     out_path = os.path.join(output_dir_WCh, output_name)

#     # Save DataFrame to ROOT
#     with uproot.recreate(out_path) as f:
#         f["Dstree"] = DataFrames[s]

#     print(f"Saved: {out_path}")

# %% [markdown]
# Reverse PID

# %%
# import os
# import uproot

# # === Make sure ReverseID samples list is defined ===
# samples_ReverseID = ["Signal_ReverseID", "BB_ReverseID", "ccbar_ReverseID", "ddbar_ReverseID", "ssbar_ReverseID", "taupair_ReverseID", "uubar_ReverseID", "Data_ReverseID"]

# # === Output directory for ReverseID ===
# output_dir_ReverseID = "/group/belle/users/amubarak/03-ML/FakeD0_ReverseID/"
# os.makedirs(output_dir_ReverseID, exist_ok=True)

# # === Base input path for original files (ReverseID) ===
# base_input_dir_ReverseID = "/group/belle/users/amubarak/02-Grid/Sample_Grid_ReverseID"
# Date_ReverseID = "0626"
# Attempt_ReverseID = "0"

# # === Save each ReverseID DataFrame using original filename with _withBDT suffix ===
# for s in samples_ReverseID:
#     if s not in DataFrames:
#         print(f"Warning: {s} not in DataFrames — skipping.")
#         continue

#     # Convert Fake D⁰ BDT output to float32 if it exists
#     if "Ds_FakeD0BDT" in DataFrames[s].columns:
#         DataFrames[s]["Ds_FakeD0BDT"] = DataFrames[s]["Ds_FakeD0BDT"].astype(np.float32)

#     # Set original file name
#     if s == "Signal_ReverseID":
#         original_name = "Ds2D0enu-Signal_ReverseID.root"
#     else:
#         tag = s.replace("_ReverseID", "")
#         original_name = f"Ds2D0e-Generic_Ds_{Date_ReverseID}25_{Attempt_ReverseID}_{tag}.root"

#     # Build output path with _withBDT suffix
#     output_name = original_name.replace(".root", "_withBDT.root")
#     out_path = os.path.join(output_dir_ReverseID, output_name)

#     # Save DataFrame to ROOT
#     with uproot.recreate(out_path) as f:
#         f["Dstree"] = DataFrames[s]

#     print(f"Saved: {out_path}")

# %% [markdown]
# Reverse PID and Wrong Charge

# %%
import os
import uproot

# === Make sure ReverseID samples list is defined ===
samples_ReverseID_WCh = ["BB_ReverseID_WCh", "ccbar_ReverseID_WCh", "ddbar_ReverseID_WCh", "ssbar_ReverseID_WCh", 
                         "taupair_ReverseID_WCh", "uubar_ReverseID_WCh", "Data_ReverseID_WCh"]

# === Output directory for ReverseID ===
output_dir_ReverseID_WCh = "/group/belle/users/amubarak/03-ML/FakeD0_ReverseID_WCh/"
os.makedirs(output_dir_ReverseID_WCh, exist_ok=True)

# === Base input path for original files (ReverseID) ===
base_input_dir_ReverseID_WCh = "/group/belle/users/amubarak/02-Grid/Sample_Grid_ReverseID_WCh"
Date_ReverseID_WCh = "0708"
Attempt_ReverseID_WCh = "0"

# === Save each ReverseID DataFrame using original filename with _withBDT suffix ===
for s in samples_ReverseID_WCh:
    if s not in DataFrames:
        print(f"Warning: {s} not in DataFrames — skipping.")
        continue

    # Convert Fake D⁰ BDT output to float32 if it exists
    if "Ds_FakeD0BDT" in DataFrames[s].columns:
        DataFrames[s]["Ds_FakeD0BDT"] = DataFrames[s]["Ds_FakeD0BDT"].astype(np.float32)

    # Set original file name
    if s == "Signal_ReverseID_WCh":
        original_name = "Ds2D0enu-Signal_ReverseID_WCh.root"
    else:
        tag = s.replace("_ReverseID_WCh", "")
        original_name = f"Ds2D0e-Generic_Ds_{Date_ReverseID_WCh}25_{Attempt_ReverseID_WCh}_{tag}.root"

    # Build output path with _withBDT suffix
    output_name = original_name.replace(".root", "_withBDT.root")
    out_path = os.path.join(output_dir_ReverseID_WCh, output_name)

    # Save DataFrame to ROOT
    with uproot.recreate(out_path) as f:
        f["Dstree"] = DataFrames[s]

    print(f"Saved: {out_path}")


