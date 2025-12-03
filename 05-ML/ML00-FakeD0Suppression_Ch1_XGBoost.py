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

pd.set_option('display.max_rows', 200000)
pd.set_option('display.max_columns', 200000)

sys.path.append("/home/belle2/amubarak/Ds2D0enue_Analysis/07-Python_Functions/")

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
cut_low = 0.14541 - (3*0.00039706)
cut_high = 0.14541 + (3*0.00042495)

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
chunk_paths = [os.path.join(chunk_dir, f"ccbar_chunk_{i:02d}.root") for i in range(20)]

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

# import os
# import glob
# import uproot
# import pandas as pd
# from tqdm import tqdm

# # === Configuration ===
# Date_ReverseID = "0708"
# Attempt_ReverseID = "0"

# # === Prompt user for veto toggle ===
# apply_veto = input("Apply veto cut on Ds_diff_D0pi? (y/n): ").strip().lower() == "y"

# # === Define veto window ===
# cut_low = 0.14543 - (3*0.00041124)
# cut_high = 0.14543 + (3*0.00041124)

# # === Variables to load ===
# with open("/home/belle2/amubarak/Ds2D0enue_Analysis/03-Grid/Save_var.txt") as f:
#     variables_to_load = [
#         line.strip().strip(",").strip('"').strip("'")
#         for line in f if line.strip() and not line.strip().startswith("#")
#     ]

# if apply_veto and "Ds_diff_D0pi" not in variables_to_load:
#     variables_to_load.append("Ds_diff_D0pi")

# # === Load merged background and signal samples ===
# merged_samples = {
#     "BB_ReverseID_WCh": f"/group/belle/users/amubarak/02-Grid/Sample_Grid_ReverseID_WCh/Ds2D0e-Generic_Ds_{Date_ReverseID}25_{Attempt_ReverseID}_BB_ReverseID_WCh.root",
#     "ccbar_ReverseID_WCh": f"/group/belle/users/amubarak/02-Grid/Sample_Grid_ReverseID_WCh/Ds2D0e-Generic_Ds_{Date_ReverseID}25_{Attempt_ReverseID}_ccbar_ReverseID_WCh.root",
#     "ddbar_ReverseID_WCh": f"/group/belle/users/amubarak/02-Grid/Sample_Grid_ReverseID_WCh/Ds2D0e-Generic_Ds_{Date_ReverseID}25_{Attempt_ReverseID}_ddbar_ReverseID_WCh.root",
#     "ssbar_ReverseID_WCh": f"/group/belle/users/amubarak/02-Grid/Sample_Grid_ReverseID_WCh/Ds2D0e-Generic_Ds_{Date_ReverseID}25_{Attempt_ReverseID}_ssbar_ReverseID_WCh.root",
#     "taupair_ReverseID_WCh": f"/group/belle/users/amubarak/02-Grid/Sample_Grid_ReverseID_WCh/Ds2D0e-Generic_Ds_{Date_ReverseID}25_{Attempt_ReverseID}_taupair_ReverseID_WCh.root",
#     "uubar_ReverseID_WCh": f"/group/belle/users/amubarak/02-Grid/Sample_Grid_ReverseID_WCh/Ds2D0e-Generic_Ds_{Date_ReverseID}25_{Attempt_ReverseID}_uubar_ReverseID_WCh.root",
#     "Data_ReverseID_WCh": f"/group/belle/users/amubarak/02-Grid/Sample_Grid_ReverseID_WCh/Ds2D0e-Generic_Ds_{Date_ReverseID}25_{Attempt_ReverseID}_Data_ReverseID_WCh.root",
# }

# print("\n📦 Loading merged samples...")
# for sample, path in tqdm(merged_samples.items(), desc="Merged Samples"):
#     try:
#         df = uproot.concatenate(f"{path}:Dstree", filter_name=variables_to_load, library="pd")
#         if apply_veto:
#             df = df[(df["Ds_diff_D0pi"] <= cut_low) | (df["Ds_diff_D0pi"] >= cut_high)]
#         DataFrames[sample] = df
#         print(f"✔️ {sample}: {len(df):,} entries")
#     except Exception as e:
#         print(f"❌ Failed to load {sample}: {e}")

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

print(DataFrames.keys())

# Vertex Fit
#----------------
for key in DataFrames.keys():
    df = DataFrames[key]
    if "Ds_chiProb_Ds_rank" in df.columns:
        DataFrames[key] = df[
            (df['Ds_chiProb_Ds_rank'] == 1)
        ]

# Photon Conversion
#--------------------
for key in DataFrames.keys():
    df = DataFrames[key]
    if "Ds_gammaveto_M_Correction" in df.columns:
        DataFrames[key] = df[
            (df["Ds_gammaveto_M_Correction"] >= 0.1)
        ]

# cut_low = 0.14541 - (3*0.00039706)
# cut_high = 0.14541 + (3*0.00042495)

# for key in DataFrames.keys():
#     df = DataFrames[key]
#     if "Ds_diff_D0pi" in df.columns:
#         DataFrames[key] = df[
#             (df["Ds_diff_D0pi"] <= cut_low) | (df["Ds_diff_D0pi"] >= cut_high)
#         ]

DataFrames["All"]["D0_isSignal"] = DataFrames["All"]["D0_isSignal"].replace(np.nan, 0)

for s in GenEvents[0:]: # loop over samples
    DataFrames[s]["D0_isSignal"] = DataFrames[s]["D0_isSignal"].replace(np.nan, 0)

DataFrames["All"]["Ds_isSignal"] = DataFrames["All"]["Ds_isSignal"].replace(np.nan, 0)

for s in GenEvents[0:]: # loop over samples
    DataFrames[s]["Ds_isSignal"] = DataFrames[s]["Ds_isSignal"].replace(np.nan, 0)

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

weights = compute_sample_weight('balanced', y_train)

# Create EarlyStopping callback
early_stop = xgboost.callback.EarlyStopping(
    rounds=10,
    metric_name='rmse',
    data_name="validation_0",
    save_best=True,
)

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

print("Best parameters found by RandomizedSearchCV:")
for param, value in random_search.best_params_.items():
    print(f"{param:20s}: {value}")

# Make predictions on the test set
y_pred_proba = xgbm_final.predict_proba(X_test)[:, 1]

# Calculate the ROC AUC score
roc_auc = roc_auc_score(y_test, y_pred_proba)

print(f"ROC AUC Score: {roc_auc:.2f}")

# Apply BDT to all DataFrames that contain the required Variables
for key in DataFrames.keys():
    df = DataFrames[key]
    
    # Check: make sure all input BDT variables exist in this DataFrame
    if all(var in df.columns for var in Variables):
        # Apply BDT and store the result
        DataFrames[key]["Ds_FakeD0BDT"] = xgbm_final.predict_proba(df[Variables])[:, 1].astype(np.float32)

import os
import uproot

# === Make sure samples is a list ===
samples = ["Signal", "BB", "ccbar", "ddbar", "ssbar", "taupair", "uubar"]

# === Output directory ===
output_dir = "/group/belle/users/amubarak/03-ML/FakeD0/"
os.makedirs(output_dir, exist_ok=True)

# === Base input path for original files ===
base_input_dir = "/group/belle/users/amubarak/02-Grid/Sample_Grid"
Date = "0530"
Attempt = "0"

# === Save each DataFrame using original filename with _withBDT suffix ===
for s in samples:
    if s not in DataFrames:
        print(f"Warning: {s} not in DataFrames — skipping.")
        continue

    # Convert Fake D⁰ BDT output to float32 if it exists
    if "Ds_FakeD0BDT" in DataFrames[s].columns:
        DataFrames[s]["Ds_FakeD0BDT"] = DataFrames[s]["Ds_FakeD0BDT"].astype(np.float32)

    # Set original file path
    if s == "Signal":
        original_name = "Ds2D0enu-Signal.root"
    else:
        original_name = f"Ds2D0e-Generic_Ds_{Date}25_{Attempt}_{s}.root"

    # Build output path with _withBDT suffix
    output_name = original_name.replace(".root", "_withBDT.root")
    out_path = os.path.join(output_dir, output_name)

    # Save DataFrame to ROOT
    with uproot.recreate(out_path) as f:
        f["Dstree"] = DataFrames[s]

    print(f"Saved: {out_path}")

import os
import uproot

# === Make sure wrong-charge samples list is defined ===
samples_WCh = ["Signal_WCh", "BB_WCh", "ccbar_WCh", "ddbar_WCh", "ssbar_WCh", "taupair_WCh", "uubar_WCh", "Data_WCh"]

# === Output directory for wrong charge ===
output_dir_WCh = "/group/belle/users/amubarak/03-ML/FakeD0_WCh/"
os.makedirs(output_dir_WCh, exist_ok=True)

# === Base input path for original files (wrong charge) ===
base_input_dir_WCh = "/group/belle/users/amubarak/02-Grid/Sample_Grid_WCh"
Date_WCh = "0630"
Attempt_WCh = "0"

# === Save each wrong-charge DataFrame using original filename with _withBDT suffix ===
for s in samples_WCh:
    if s not in DataFrames:
        print(f"Warning: {s} not in DataFrames — skipping.")
        continue

    # Convert Fake D⁰ BDT output to float32 if it exists
    if "Ds_FakeD0BDT" in DataFrames[s].columns:
        DataFrames[s]["Ds_FakeD0BDT"] = DataFrames[s]["Ds_FakeD0BDT"].astype(np.float32)

    # Set original file name
    if s == "Signal_WCh":
        original_name = "Ds2D0enu-Signal_WCh.root"
    else:
        tag = s.replace("_WCh", "")
        original_name = f"Ds2D0e-Generic_Ds_{Date_WCh}25_{Attempt_WCh}_{tag}.root"

    # Build output path with _withBDT suffix
    output_name = original_name.replace(".root", "_withBDT.root")
    out_path = os.path.join(output_dir_WCh, output_name)

    # Save DataFrame to ROOT
    with uproot.recreate(out_path) as f:
        f["Dstree"] = DataFrames[s]

    print(f"Saved: {out_path}")

import os
import uproot

# === Make sure ReverseID samples list is defined ===
samples_ReverseID = ["Signal_ReverseID", "BB_ReverseID", "ccbar_ReverseID", "ddbar_ReverseID", "ssbar_ReverseID", "taupair_ReverseID", "uubar_ReverseID", "Data_ReverseID"]

# === Output directory for ReverseID ===
output_dir_ReverseID = "/group/belle/users/amubarak/03-ML/FakeD0_ReverseID/"
os.makedirs(output_dir_ReverseID, exist_ok=True)

# === Base input path for original files (ReverseID) ===
base_input_dir_ReverseID = "/group/belle/users/amubarak/02-Grid/Sample_Grid_ReverseID"
Date_ReverseID = "0626"
Attempt_ReverseID = "0"

# === Save each ReverseID DataFrame using original filename with _withBDT suffix ===
for s in samples_ReverseID:
    if s not in DataFrames:
        print(f"Warning: {s} not in DataFrames — skipping.")
        continue

    # Convert Fake D⁰ BDT output to float32 if it exists
    if "Ds_FakeD0BDT" in DataFrames[s].columns:
        DataFrames[s]["Ds_FakeD0BDT"] = DataFrames[s]["Ds_FakeD0BDT"].astype(np.float32)

    # Set original file name
    if s == "Signal_ReverseID":
        original_name = "Ds2D0enu-Signal_ReverseID.root"
    else:
        tag = s.replace("_ReverseID", "")
        original_name = f"Ds2D0e-Generic_Ds_{Date_ReverseID}25_{Attempt_ReverseID}_{tag}.root"

    # Build output path with _withBDT suffix
    output_name = original_name.replace(".root", "_withBDT.root")
    out_path = os.path.join(output_dir_ReverseID, output_name)

    # Save DataFrame to ROOT
    with uproot.recreate(out_path) as f:
        f["Dstree"] = DataFrames[s]

    print(f"Saved: {out_path}")

# import os
# import uproot

# # === Make sure ReverseID samples list is defined ===
# samples_ReverseID_WCh = ["BB_ReverseID_WCh", "ccbar_ReverseID_WCh", "ddbar_ReverseID_WCh", "ssbar_ReverseID_WCh", 
#                          "taupair_ReverseID_WCh", "uubar_ReverseID_WCh", "Data_ReverseID_WCh"]

# # === Output directory for ReverseID ===
# output_dir_ReverseID_WCh = "/group/belle/users/amubarak/03-ML/FakeD0_ReverseID_WCh/"
# os.makedirs(output_dir_ReverseID_WCh, exist_ok=True)

# # === Base input path for original files (ReverseID) ===
# base_input_dir_ReverseID_WCh = "/group/belle/users/amubarak/02-Grid/Sample_Grid_ReverseID_WCh"
# Date_ReverseID_WCh = "0708"
# Attempt_ReverseID_WCh = "0"

# # === Save each ReverseID DataFrame using original filename with _withBDT suffix ===
# for s in samples_ReverseID_WCh:
#     if s not in DataFrames:
#         print(f"Warning: {s} not in DataFrames — skipping.")
#         continue

#     # Convert Fake D⁰ BDT output to float32 if it exists
#     if "Ds_FakeD0BDT" in DataFrames[s].columns:
#         DataFrames[s]["Ds_FakeD0BDT"] = DataFrames[s]["Ds_FakeD0BDT"].astype(np.float32)

#     # Set original file name
#     if s == "Signal_ReverseID_WCh":
#         original_name = "Ds2D0enu-Signal_ReverseID_WCh.root"
#     else:
#         tag = s.replace("_ReverseID_WCh", "")
#         original_name = f"Ds2D0e-Generic_Ds_{Date_ReverseID_WCh}25_{Attempt_ReverseID_WCh}_{tag}.root"

#     # Build output path with _withBDT suffix
#     output_name = original_name.replace(".root", "_withBDT.root")
#     out_path = os.path.join(output_dir_ReverseID_WCh, output_name)

#     # Save DataFrame to ROOT
#     with uproot.recreate(out_path) as f:
#         f["Dstree"] = DataFrames[s]

#     print(f"Saved: {out_path}")