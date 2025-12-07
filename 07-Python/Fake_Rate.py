import sys
import numpy as np
import matplotlib.pyplot as plt
from IPython.display import HTML 
import warnings
# Path to scripts on KEKCC
sys_path = '/group/belle2/dataprod/Systematics/systematic_corrections_framework/scripts'
# for NAF: 
# sys_path = '/nfs/dust/belle2/group/dataprod/Systematics/systematic_corrections_framework/scripts'
sys.path.insert(1, sys_path)

import weight_table as wm
from show_db_content import show_db_content
from show_variables import show_ntuple_variables
import show_collections as sc
import sysvar

import seaborn as sns

import numpy as np
import pandas as pd
import uproot
from matplotlib import pyplot as plt

plt.rcParams.update({
    "axes.labelsize": 16,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 16,
    "figure.titlesize": 20
})

pd.set_option('display.max_rows', 200000)
pd.set_option('display.max_columns', 200000)

import sys
sys.path.append("/home/belle2/amubarak/Ds2D0enue_Analysis/07-Python_Functions/")

from plotting_utils_pidopt import optimize_cut, plot_roc_curve

import os
import uproot
import pandas as pd
from tqdm import tqdm

# === Load only selected branches ===
with open("/home/belle2/amubarak/Ds2D0enue_Analysis/03-Grid/Save_var.txt") as f:
    variables_to_load = [
        line.strip().strip(",").strip('"').strip("'")
        for line in f
        if line.strip() and not line.strip().startswith("#")
    ]

# Ensure BDT variable is included
if "Ds_FakeD0BDT" not in variables_to_load:
    variables_to_load.append("Ds_FakeD0BDT")

# Make sure Ds_BkgBDT is included
if "Ds_BkgBDT" not in variables_to_load:
    variables_to_load.append("Ds_BkgBDT")

# === Sample list ===
samples = ["Signal", "BB", "ccbar", "ddbar", "ssbar", "taupair", "uubar"]
GenEvents = samples.copy()

# === Input configuration for BkgBDT ===
Date = "0530"
Attempt = "0"
input_dir = "/group/belle/users/amubarak/03-ML/BkgBDT/"

# === Load ROOT files into DataFrames ===
DataFrames = {}

for s in tqdm(samples, desc="Loading ROOT files"):
    if s == "Signal":
        file_path = os.path.join(input_dir, "Ds2D0enu-Signal_withBkgBDT.root")
    else:
        file_path = os.path.join(
            input_dir, f"Ds2D0e-Generic_Ds_{Date}25_{Attempt}_{s}_withBkgBDT.root"
        )

    df = uproot.concatenate(
        f"{file_path}:Dstree",
        filter_name=variables_to_load,
        library="pd"
    )

    print(f"✔️ Loaded {file_path} with {len(df):,} entries")
    DataFrames[s] = df

# === Define combined background ===
background_samples = ["BB", "ccbar", "ddbar", "ssbar", "taupair", "uubar"]
DataFrames["All"] = pd.concat([DataFrames[s] for s in background_samples], ignore_index=True)

# === Combine uds backgrounds for convenience ===
DataFrames["uds"] = pd.concat(
    [DataFrames["uubar"], DataFrames["ddbar"], DataFrames["ssbar"]],
    ignore_index=True
)

import os
import uproot
import numpy as np
import pandas as pd
from tqdm import tqdm

# === Load only selected branches ===
with open("/home/belle2/amubarak/Ds2D0enue_Analysis/03-Grid/Save_var.txt") as f:
    variables_to_load = [
        line.strip().strip(",").strip('"').strip("'")
        for line in f
        if line.strip() and not line.strip().startswith("#")
    ]

# Ensure critical vars included
for var in ["Ds_FakeD0BDT", "Ds_BkgBDT"]:
    if var not in variables_to_load:
        variables_to_load.append(var)

# === Ask user if they want to load every 4th entry
skip_by_4 = input("Load only every 4th entry of control sample? (y/n): ").strip().lower() == "y"

# === ReverseID sample list ===
samples_ReverseID = [
    "Signal_ReverseID", "BB_ReverseID", "ccbar_ReverseID",
    "ddbar_ReverseID", "ssbar_ReverseID", "taupair_ReverseID",
    "uubar_ReverseID", "Data_ReverseID"
]

# === Input configuration ===
Date = "0626"
Attempt = "0"
input_dir = "/group/belle/users/amubarak/03-ML/BkgBDT_ReverseID/"

# === Load ROOT files into DataFrames ===
for s in tqdm(samples_ReverseID, desc="Loading ReverseID"):
    if s == "Signal_ReverseID":
        file_path = os.path.join(input_dir, "Ds2D0enu-Signal_ReverseID_withBkgBDT.root")
    else:
        tag = s.replace("_ReverseID", "")
        file_path = os.path.join(
            input_dir, f"Ds2D0e-Generic_Ds_{Date}25_{Attempt}_{tag}_withBkgBDT.root"
        )

    # Load full tree first
    df = uproot.concatenate(
        f"{file_path}:Dstree",
        filter_name=variables_to_load,
        library="pd"
    )

    # Then skip every 4th entry if requested (not for signal)
    if skip_by_4 and s != "Signal_ReverseID":
        df = df.iloc[::4].reset_index(drop=True)

    print(f"✔️ Loaded {file_path} with {len(df):,} entries")
    DataFrames[s] = df

# === Combine reversed PID backgrounds
background_samples = [
    "BB_ReverseID", "ccbar_ReverseID", "ddbar_ReverseID",
    "ssbar_ReverseID", "taupair_ReverseID", "uubar_ReverseID"
]

DataFrames["All_ReverseID"] = pd.concat(
    [DataFrames[s] for s in background_samples],
    ignore_index=True
)

DataFrames["uds_ReverseID"] = pd.concat(
    [DataFrames["uubar_ReverseID"], DataFrames["ddbar_ReverseID"], DataFrames["ssbar_ReverseID"]],
    ignore_index=True
)

import os
import uproot
import numpy as np
import pandas as pd
from tqdm import tqdm

# === Load only selected branches ===
with open("/home/belle2/amubarak/Ds2D0enue_Analysis/03-Grid/Save_var.txt") as f:
    variables_to_load = [
        line.strip().strip(",").strip('"').strip("'")
        for line in f
        if line.strip() and not line.strip().startswith("#")
    ]

# Ensure critical vars included
for var in ["Ds_FakeD0BDT", "Ds_BkgBDT"]:
    if var not in variables_to_load:
        variables_to_load.append(var)

# === Ask user if they want to load every 4th entry
skip_by_4 = input("Load only every 4th entry of control sample? (y/n): ").strip().lower() == "y"

# === ReverseID_WCh sample list ===
samples_ReverseID_WCh = [
    "BB_ReverseID_WCh", "ccbar_ReverseID_WCh",
    "ddbar_ReverseID_WCh", "ssbar_ReverseID_WCh", "taupair_ReverseID_WCh",
    "uubar_ReverseID_WCh", "Data_ReverseID_WCh"
]

# === Input configuration ===
Date = "0708"
Attempt = "0"
input_dir = "/group/belle/users/amubarak/03-ML/BkgBDT_ReverseID_WCh/"

# === Load ROOT files into DataFrames ===
for s in tqdm(samples_ReverseID_WCh, desc="Loading ReverseID_WCh"):
    if s == "Signal_ReverseID_WCh":
        file_path = os.path.join(input_dir, "Ds2D0enu-Signal_ReverseID_WCh_withBkgBDT.root")
    else:
        tag = s.replace("_ReverseID_WCh", "")
        file_path = os.path.join(
            input_dir, f"Ds2D0e-Generic_Ds_{Date}25_{Attempt}_{tag}_withBkgBDT.root"
        )

    # Load full tree first
    df = uproot.concatenate(
        f"{file_path}:Dstree",
        filter_name=variables_to_load,
        library="pd"
    )

    # Then skip every 4th entry if requested (not for signal)
    if skip_by_4 and s != "Signal_ReverseID_WCh":
        df = df.iloc[::4].reset_index(drop=True)

    print(f"✔️ Loaded {file_path} with {len(df):,} entries")
    DataFrames[s] = df

# === Combine reversed PID backgrounds
background_samples = [
    "BB_ReverseID_WCh", "ccbar_ReverseID_WCh", "ddbar_ReverseID_WCh",
    "ssbar_ReverseID_WCh", "taupair_ReverseID_WCh", "uubar_ReverseID_WCh"
]

DataFrames["All_ReverseID_WCh"] = pd.concat(
    [DataFrames[s] for s in background_samples],
    ignore_index=True
)

DataFrames["uds_ReverseID_WCh"] = pd.concat(
    [DataFrames["uubar_ReverseID_WCh"], DataFrames["ddbar_ReverseID_WCh"], DataFrames["ssbar_ReverseID_WCh"]],
    ignore_index=True
)

print(DataFrames.keys())

# import numpy as np
# import pandas as pd
# import matplotlib.pyplot as plt
# import seaborn as sns
# from matplotlib.patches import Rectangle
# import matplotlib.colors as colors

# # === Use your input DataFrame ===
# df = DataFrames["All_ReverseID"]

# # === Configuration ===
# pid_cut = 0.5
# theta_vals = [0.22, 0.56, 1.13, 1.57, 1.88, 2.23, 2.71]
# cos_edges = np.round(np.cos(theta_vals), 4)[::-1]  # descending in cos(theta)
# charge_vals = [-2, 0]
# p_edges = np.linspace(0.0, 3.0, 11)
# x_bin_pairs = [(i, ch) for i in range(len(cos_edges) - 1) for ch in charge_vals]

# # === Compute Fake Rate Matrix ===
# fake_rate_matrix = np.full((len(p_edges) - 1, len(x_bin_pairs)), np.nan)
# for i in range(len(p_edges) - 1):
#     df_p = df[(df["e_p"] > p_edges[i]) & (df["e_p"] <= p_edges[i + 1])]
#     for j, (cos_idx, ch_val) in enumerate(x_bin_pairs):
#         cos_min = min(cos_edges[cos_idx], cos_edges[cos_idx+1])
#         cos_max = max(cos_edges[cos_idx], cos_edges[cos_idx+1])
#         cos_mask = (df_p["e_cos_theta"] > cos_min) & (df_p["e_cos_theta"] <= cos_max)
#         charge_mask = df_p["e_charge"] < 0 if ch_val == -2 else df_p["e_charge"] >= 0
#         bin_df = df_p[charge_mask & cos_mask]
#         n_pass = (bin_df["e_electronID"] > pid_cut).sum()
#         n_fail = (bin_df["e_electronID"] <= pid_cut).sum()
#         total = n_pass + n_fail
#         if total > 0:
#             fake_rate_matrix[i, j] = n_pass / total

# # === Helper for label formatting ===
# class cm:
#     @staticmethod
#     def compactify_numbered(x):
#         return x

# # === Final FakeRateTable class ===
# class FakeRateTable:
#     def __init__(self, values, binning, bin_variables, cut, particle_type, processing, title):
#         self.values = values
#         self.binning = binning
#         self.bin_variables = bin_variables
#         self.cut = cut
#         self.particle_type = particle_type
#         self.processing = processing
#         self.title = title

#     def compute_vmax(self, array):
#         return np.nanmax(array)

#     def plot(self, *args, **kw_args):
#         plot_args = kw_args.get('plot_args')
#         plot_percent_tuple = (100, ' [%]', 0.1) if kw_args.get('use_percent') else (1, '', 0.2)
#         if kw_args.get('display_plots') is None:
#             kw_args['display_plots'] = True

#         label = cm.compactify_numbered(self.processing)
#         if kw_args.get('label'):
#             label = kw_args.get('label')

#         object_name = self.title.replace('_', ' ')
#         plot_title = f'{self.particle_type} {object_name} {plot_percent_tuple[1]} for cut "{self.cut}"'

#         fig = kw_args.get('fig')
#         ax = kw_args.get('ax')
#         if fig is None or ax is None:
#             fig, ax = plt.subplots(1, 1, figsize=(9, 5), dpi=120)

#         if len(self.values.shape) != 2:
#             raise ValueError("Only 2D plotting is supported.")

#         # Annotations WITHOUT uncertainties
#         annot_labels = np.asarray([
#             fr'{value*plot_percent_tuple[0]:{plot_percent_tuple[2]}f}'
#             for value in self.values.flatten()
#         ]).reshape(self.values.shape)

#         if plot_args is None:
#             plot_args = {
#                 'linewidths': .5,
#                 'fmt': '',
#                 'cbar_kws': {'fraction': .05},
#                 'vmin': 0,
#                 'vmax': self.compute_vmax(self.values) * plot_percent_tuple[0],
#                 'cmap': 'plasma'
#             }

#         sns.heatmap(self.values * plot_percent_tuple[0], annot=annot_labels, ax=ax, **plot_args)

#         # X-axis: two-line labels
#         xtick_labels = []
#         for idx, (cos_idx, ch) in enumerate(x_bin_pairs):
#             cos_val = cos_edges[cos_idx]
#             if ch == -2:
#                 xtick_labels.append(f"{ch}\n{cos_val:.4f}")
#             else:
#                 xtick_labels.append(f"{ch}\n")
#         ax.set_xticks(np.arange(len(xtick_labels)))
#         ax.set_xticklabels(xtick_labels, rotation=0)

#         # Y-axis: momentum bin edges
#         ytick_labels = [f"{e:.1f}" for e in p_edges[:-1]]
#         ax.set_yticks(np.arange(len(ytick_labels)))
#         ax.set_yticklabels(ytick_labels, rotation=0, fontsize=10)

#         ax.set_xlabel(f'{self.bin_variables[1]} bins', fontsize=14)
#         ax.set_ylabel(f'{self.bin_variables[0]} bins', fontsize=14)
#         ax.set_title(plot_title)
#         fig.tight_layout()

#         if kw_args.get('save_plots'):
#             fig.savefig(kw_args['save_plots'])

#         if kw_args.get('display_plots'):
#             if not hasattr(__builtins__, '__IPYTHON__'):
#                 fig.show()
#                 input()
#         else:
#             plt.close()

# # === Create and Plot ===
# table = FakeRateTable(
#     values=fake_rate_matrix,
#     binning=[p_edges, cos_edges],
#     bin_variables=["e_p", "cosθ + charge"],
#     cut="e_electronID > 0.5",
#     particle_type="electron",
#     processing="proc13",
#     title="Fake Rate"
# )

# table.plot(display_plots=True)

import numpy as np
import pandas as pd

# === Use your input DataFrame ===
df = DataFrames["All_ReverseID"]

# === Configuration ===
pid_cut = 0.5
theta_vals = [0.22, 0.56, 1.13, 1.57, 1.88, 2.23, 2.71]
cos_edges = np.round(np.cos(theta_vals), 4)[::-1]  # descending in cos(theta)
p_edges = np.linspace(0.0, 3.0, 11)
charge_vals = [-2, 0]  # -2 = negative, 0 = positive

# === Build flat DataFrame ===
rows = []
for i in range(len(p_edges) - 1):
    p_min = p_edges[i]
    p_max = p_edges[i + 1]
    df_p = df[(df["e_p"] > p_min) & (df["e_p"] <= p_max)]

    for j in range(len(cos_edges) - 1):
        cos_min = min(cos_edges[j], cos_edges[j+1])
        cos_max = max(cos_edges[j], cos_edges[j+1])
        cos_mask = (df_p["e_cos_theta"] > cos_min) & (df_p["e_cos_theta"] <= cos_max)

        for ch_val in charge_vals:
            if ch_val == -2:
                charge_mask = df_p["e_charge"] < 0
            else:
                charge_mask = df_p["e_charge"] >= 0

            bin_df = df_p[cos_mask & charge_mask]
            n_pass = (bin_df["e_electronID"] > pid_cut).sum()
            n_fail = (bin_df["e_electronID"] <= pid_cut).sum()
            total = n_pass + n_fail

            fake_rate = n_pass / total if total > 0 else np.nan

            rows.append({
                "p_min": round(p_min, 4),
                "p_max": round(p_max, 4),
                "cosTheta_min": round(cos_min, 4),
                "cosTheta_max": round(cos_max, 4),
                "charge_min": ch_val,
                "charge_max": 0 if ch_val == -2 else 2,
                "fake_rate": fake_rate
            })

# === Final table ===
fake_rate_df = pd.DataFrame(rows)
# fake_rate_df

import pandas as pd
import numpy as np

def add_fake_rate_weights(df: pd.DataFrame,
                          fake_rate_table: pd.DataFrame,
                          p_col: str = "e_p",
                          cos_col: str = "e_cos_theta",
                          charge_col: str = "e_charge",
                          weight_col: str = "fake_weight",
                          fillna_value: float = 1.0,
                          epsilon: float = 1e-8) -> pd.DataFrame:
    """
    Assigns fake rate weights to a dataframe using a bin-by-bin match with a fake rate lookup table.

    Args:
        df (pd.DataFrame): Input DataFrame to be weighted.
        fake_rate_table (pd.DataFrame): Table of fake rates with bin edges and fake_rate column.
        p_col (str): Column in df for particle momentum.
        cos_col (str): Column in df for cos(theta).
        charge_col (str): Column in df for particle charge.
        weight_col (str): Name of the output weight column.
        fillna_value (float): Value to fill for unmatched entries (default: 1.0).
        epsilon (float): Small value to loosen bin edge comparisons (default: 1e-8).

    Returns:
        pd.DataFrame: The input DataFrame with an additional `fake_weight` column.
    """
    df = df.copy()
    df[weight_col] = np.nan  # Initialize weight column with NaN

    for _, row in fake_rate_table.iterrows():
        p_mask = (df[p_col] >= row["p_min"] - epsilon) & (df[p_col] < row["p_max"] + epsilon)
        cos_mask = (df[cos_col] >= row["cosTheta_min"] - epsilon) & (df[cos_col] < row["cosTheta_max"] + epsilon)

        if row["charge_min"] == -2:
            charge_mask = df[charge_col] < 0
        else:
            charge_mask = df[charge_col] >= 0

        total_mask = p_mask & cos_mask & charge_mask
        df.loc[total_mask, weight_col] = row["fake_rate"]

    df[weight_col] = df[weight_col].fillna(fillna_value)
    return df

# Input:
# - df: any DataFrame to apply weights to (signal, background, etc.)
# - fake_rate_df: already built using your code

DataFrames["All"] = add_fake_rate_weights(DataFrames["All"], fake_rate_df)
DataFrames["Signal_ReverseID"] = add_fake_rate_weights(DataFrames["Signal_ReverseID"], fake_rate_df)
DataFrames["All_ReverseID"] = add_fake_rate_weights(DataFrames["All_ReverseID"], fake_rate_df)
DataFrames["Data_ReverseID"] = add_fake_rate_weights(DataFrames["Data_ReverseID"], fake_rate_df)
DataFrames["All_ReverseID_WCh"] = add_fake_rate_weights(DataFrames["All_ReverseID_WCh"], fake_rate_df)
DataFrames["Data_ReverseID_WCh"] = add_fake_rate_weights(DataFrames["Data_ReverseID_WCh"], fake_rate_df)

# Optionally check result
# DataFrames["All_ReverseID"][["e_p", "e_cos_theta", "e_charge", "e_electronID", "fake_weight"]].head()

import os
import uproot
import numpy as np  # Ensure this is available for dtype conversion

# === Define ReverseID samples ===
samples_ReverseID = ["Signal_ReverseID", "All_ReverseID", "Data_ReverseID"]

# === Output directory for BkgBDT ReverseID files ===
output_dir_ReverseID = "/group/belle/users/amubarak/04-Reversed_PID_FakeRate/"
os.makedirs(output_dir_ReverseID, exist_ok=True)

# === Base input path for ReverseID ===
base_input_dir_ReverseID = "/group/belle/users/amubarak/02-Grid/Sample_Grid_ReverseID"
Date_ReverseID = "0626"
Attempt_ReverseID = "0"

# === Save each ReverseID DataFrame with Ds_BkgBDT ===
for s in samples_ReverseID:
    if s not in DataFrames:
        print(f"Warning: {s} not in DataFrames — skipping.")
        continue

    # Set original file name
    if s == "Signal_ReverseID":
        original_name = "Ds2D0enu-Signal_ReverseID.root"
    else:
        tag = s.replace("_ReverseID", "")
        original_name = f"Ds2D0e-Generic_Ds_{Date_ReverseID}25_{Attempt_ReverseID}_{tag}.root"

    # Build output file name with _withBkgBDT suffix
    output_name = original_name.replace(".root", "_withBkgBDT.root")
    out_path = os.path.join(output_dir_ReverseID, output_name)

    # Save DataFrame to ROOT
    with uproot.recreate(out_path) as f:
        f["Dstree"] = DataFrames[s]

    print(f"Saved: {out_path}")

import os
import uproot
import numpy as np  # Ensure this is available for dtype conversion

# === Define ReverseID_WCh samples ===
samples_ReverseID_WCh = ["All_ReverseID_WCh", "Data_ReverseID_WCh"]

# === Output directory for BkgBDT ReverseID files ===
output_dir_ReverseID_WCh = "/group/belle/users/amubarak/04-Reversed_PID_FakeRate/"
os.makedirs(output_dir_ReverseID_WCh, exist_ok=True)

# === Base input path for ReverseID ===
base_input_dir_ReverseID_WCh = "/group/belle/users/amubarak/02-Grid/Sample_Grid_ReverseID_WCh"
Date_ReverseID_WCh = "0708"
Attempt_ReverseID_WCh = "0"

# === Save each ReverseID_WCh DataFrame with Ds_BkgBDT ===
for s in samples_ReverseID_WCh:
    if s not in DataFrames:
        print(f"Warning: {s} not in DataFrames — skipping.")
        continue

    # Set original file name
    if s == "Signal_ReverseID_WCh":
        original_name = "Ds2D0enu-Signal_ReverseID_WCh.root"
    else:
        tag = s.replace("_ReverseID_WCh", "")
        original_name = f"Ds2D0e-Generic_Ds_{Date_ReverseID_WCh}25_{Attempt_ReverseID_WCh}_{tag}.root"

    # Build output file name with _withBkgBDT suffix
    output_name = original_name.replace(".root", "_withBkgBDT.root")
    out_path = os.path.join(output_dir_ReverseID_WCh, output_name)

    # Save DataFrame to ROOT
    with uproot.recreate(out_path) as f:
        f["Dstree"] = DataFrames[s]

    print(f"Saved: {out_path}")