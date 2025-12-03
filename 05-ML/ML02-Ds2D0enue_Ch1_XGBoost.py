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
import uproot
import pandas as pd

# === Load only selected branches ===
with open("/home/belle2/amubarak/Ds2D0enue_Analysis/03-Grid/Save_var.txt") as f:
    variables_to_load = [
        line.strip().strip(",").strip('"').strip("'")
        for line in f
        if line.strip() and not line.strip().startswith("#")
    ]

# Make sure to include BDT output variable
if "Ds_FakeD0BDT" not in variables_to_load:
    variables_to_load.append("Ds_FakeD0BDT")

# === Sample list ===
samples = ["Signal", "BB", "ccbar", "ddbar", "ssbar", "taupair", "uubar"]
GenEvents = samples.copy()

# === Input configuration ===
Date = "0530"
Attempt = "0"
input_dir = "/group/belle/users/amubarak/03-ML/FakeD0/"

# === Load ROOT files into DataFrames ===
DataFrames = {}

for s in samples:
    if s == "Signal":
        file_path = os.path.join(input_dir, "Ds2D0enu-Signal_withBDT.root")
    else:
        file_path = os.path.join(
            input_dir, f"Ds2D0e-Generic_Ds_{Date}25_{Attempt}_{s}_withBDT.root"
        )

    print(f"Loading: {file_path}")
    DataFrames[s] = uproot.concatenate(
        f"{file_path}:Dstree",
        filter_name=variables_to_load,
        library="pd"
    )

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
import pandas as pd

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

# === Wrong-charge samples ===
samples_WCh = ["Signal_WCh", "BB_WCh", "ccbar_WCh", "ddbar_WCh", "ssbar_WCh", "taupair_WCh", "uubar_WCh", "Data_WCh"]
background_WCh = ["BB_WCh", "ccbar_WCh", "ddbar_WCh", "ssbar_WCh", "taupair_WCh", "uubar_WCh"]

Date_WCh = "0630"
Attempt_WCh = "0"
input_dir_WCh = "/group/belle/users/amubarak/03-ML/FakeD0_WCh/"

# === Load wrong-charge ROOT files into DataFrames ===
DataFrames = {} if "DataFrames" not in globals() else DataFrames

for s in samples_WCh:
    if s == "Signal_WCh":
        file_path = os.path.join(input_dir_WCh, "Ds2D0enu-Signal_WCh_withBDT.root")
    else:
        tag = s.replace("_WCh", "")
        file_path = os.path.join(
            input_dir_WCh,
            f"Ds2D0e-Generic_Ds_{Date_WCh}25_{Attempt_WCh}_{tag}_withBDT.root"
        )

    print(f"Loading: {file_path}")
    DataFrames[s] = uproot.concatenate(
        f"{file_path}:Dstree",
        filter_name=variables_to_load,
        library="pd"
    )

# === Combine wrong-charge backgrounds ===
DataFrames["All_WCh"] = pd.concat([DataFrames[s] for s in background_WCh], ignore_index=True)
DataFrames["uds_WCh"] = pd.concat(
    [DataFrames[s] for s in ["uubar_WCh", "ddbar_WCh", "ssbar_WCh"]],
    ignore_index=True
)

import os
import uproot
import pandas as pd

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

samples_ReverseID = ["Signal_ReverseID", "BB_ReverseID", "ccbar_ReverseID", "ddbar_ReverseID", "ssbar_ReverseID", "taupair_ReverseID", "uubar_ReverseID", "Data_ReverseID"]
background_ReverseID = ["BB_ReverseID", "ccbar_ReverseID", "ddbar_ReverseID", "ssbar_ReverseID", "taupair_ReverseID", "uubar_ReverseID"]

Date_ReverseID = "0626"
Attempt_ReverseID = "0"
input_dir_ReverseID = "/group/belle/users/amubarak/03-ML/FakeD0_ReverseID/"

DataFrames = {} if "DataFrames" not in globals() else DataFrames
step_size = 100_000  # Load in steps of 100k entries

for s in samples_ReverseID:
    if s == "Signal_ReverseID":
        file_path = os.path.join(input_dir_ReverseID, "Ds2D0enu-Signal_ReverseID_withBDT.root")
    else:
        tag = s.replace("_ReverseID", "")
        file_path = os.path.join(
            input_dir_ReverseID,
            f"Ds2D0e-Generic_Ds_{Date_ReverseID}25_{Attempt_ReverseID}_{tag}_withBDT.root"
        )

    print(f"Loading in chunks: {file_path}")
    
    chunks = []
    for chunk in uproot.iterate(
        f"{file_path}:Dstree",
        filter_name=variables_to_load,
        step_size=step_size,
        library="pd"
    ):
        chunks.append(chunk)

    DataFrames[s] = pd.concat(chunks, ignore_index=True)

# === Combine wrong-charge backgrounds ===
DataFrames["All_ReverseID"] = pd.concat([DataFrames[s] for s in background_ReverseID], ignore_index=True)
DataFrames["uds_ReverseID"] = pd.concat(
    [DataFrames[s] for s in ["uubar_ReverseID", "ddbar_ReverseID", "ssbar_ReverseID"]],
    ignore_index=True
)

# import os
# import uproot
# import pandas as pd

# # === Load only selected branches ===
# with open("/home/belle2/amubarak/Ds2D0enue_Analysis/03-Grid/Save_var.txt") as f:
#     variables_to_load = [
#         line.strip().strip(",").strip('"').strip("'")
#         for line in f
#         if line.strip() and not line.strip().startswith("#")
#     ]

# # Ensure BDT variable is included
# if "Ds_FakeD0BDT" not in variables_to_load:
#     variables_to_load.append("Ds_FakeD0BDT")

# # === _ReverseID_WCh samples ===
# samples_ReverseID_WCh = ["BB_ReverseID_WCh", "ccbar_ReverseID_WCh", "ddbar_ReverseID_WCh", "ssbar_ReverseID_WCh", "taupair_ReverseID_WCh", "uubar_ReverseID_WCh", "Data_ReverseID_WCh"]
# background_ReverseID_WCh = ["BB_ReverseID_WCh", "ccbar_ReverseID_WCh", "ddbar_ReverseID_WCh", "ssbar_ReverseID_WCh", "taupair_ReverseID_WCh", "uubar_ReverseID_WCh"]

# Date_ReverseID_WCh = "0708"
# Attempt_ReverseID_WCh = "0"
# input_dir_ReverseID_WCh = "/group/belle/users/amubarak/03-ML/FakeD0_ReverseID_WCh/"

# # === Load _ReverseID_WCh ROOT files into DataFrames ===
# DataFrames = {} if "DataFrames" not in globals() else DataFrames

# for s in samples_ReverseID_WCh:
#     if s == "Signal_ReverseID_WCh":
#         file_path = os.path.join(input_dir_ReverseID_WCh, "Ds2D0enu-Signal_ReverseID_WCh_withBDT.root")
#     else:
#         tag = s.replace("_ReverseID_WCh", "")
#         file_path = os.path.join(
#             input_dir_ReverseID_WCh,
#             f"Ds2D0e-Generic_Ds_{Date_ReverseID_WCh}25_{Attempt_ReverseID_WCh}_{tag}_withBDT.root"
#         )

#     print(f"Loading: {file_path}")
#     DataFrames[s] = uproot.concatenate(
#         f"{file_path}:Dstree",
#         filter_name=variables_to_load,
#         library="pd"
#     )

# # === Combine wrong-charge backgrounds ===
# DataFrames["All_ReverseID_WCh"] = pd.concat([DataFrames[s] for s in background_ReverseID_WCh], ignore_index=True)
# DataFrames["uds_ReverseID_WCh"] = pd.concat(
#     [DataFrames[s] for s in ["uubar_ReverseID_WCh", "ddbar_ReverseID_WCh", "ssbar_ReverseID_WCh"]],
#     ignore_index=True
# )

# === Inject derived variables ===
for key in DataFrames:
    df = DataFrames[key]

    DataFrames[key]["log_e_p"] = np.log(df["e_p"].replace(0, np.nan))
    DataFrames[key]["pe_ratio"] = df["D0_p"] / df["e_p"].replace(0, np.nan)

DataFrames["Signal"].columns.tolist()

print(DataFrames.keys())

# # Setup

# # Photon Conversion
# #------------------------
# for key in DataFrames.keys():
#     df = DataFrames[key]
#     if "Ds_gammaveto_M_Correction" in df.columns:
#         DataFrames[key] = df[
#             (df["Ds_gammaveto_M_Correction"] >= 0.1)
#         ]

# # Vertex Fit
# #----------------
# for key in DataFrames.keys():
#     df = DataFrames[key]
#     if "Ds_chiProb_Ds_rank" in df.columns:
#         DataFrames[key] = df[
#             (df['Ds_chiProb_Ds_rank'] == 1)
#         ]

DataFrames["All"]["D0_isSignal"] = DataFrames["All"]["D0_isSignal"].replace(np.nan, 0)

for s in GenEvents[0:]: # loop over samples
    DataFrames[s]["D0_isSignal"] = DataFrames[s]["D0_isSignal"].replace(np.nan, 0)

DataFrames["All"]["Ds_isSignal"] = DataFrames["All"]["Ds_isSignal"].replace(np.nan, 0)

for s in GenEvents[0:]: # loop over samples
    DataFrames[s]["Ds_isSignal"] = DataFrames[s]["Ds_isSignal"].replace(np.nan, 0)

Variables = [
    # Electron
    "e_cos_theta", 
    # "e_electronID",
    "e_dr", "e_dz",
    # Kaon
    "K_dr",
    # Pion
    "pi_dr",
    # D0
    # "D0_significanceOfDistance",
    "Ds_FakeD0BDT",
    "D0_useCMSFrame_p",
    # Ds
    # "Ds_Ds_starminusDs_M_Correction",
    # "log_e_p",
    # "pe_ratio",
    # "Ds_cos_theta",
    # "Ds_psi",
    # "Ds_pointangle","Ds_daughterDiffOf_0_1_cos_theta",
    # "Ds_daughterDiffOfPhi_0_1"
    # "Ds_flightDistance",
    # 'Ds_vertexDistanceOfDaughter_0',
    # 'Ds_vertexDistanceOfDaughterErr_0',
    # "log_e_p",
    # "pe_ratio"
]

# === Organize data for ML: features (X) and labels (y) ===

all_MC = []  # list of all MC feature arrays
all_y = []   # list of all MC label arrays

for s in GenEvents:
    if s == "data":
        continue  # skip data

    df = DataFrames[s]

    if "Signal" in s:
        # Use only true signal (Ds_isSignal == 1)
        true_signal = df[df["Ds_isSignal"] == 1]
        all_MC.append(true_signal[Variables])
        all_y.append(np.ones(true_signal.shape[0], dtype=np.int32))
    else:
        # All background MC
        all_MC.append(df[Variables])
        all_y.append(np.zeros(df.shape[0], dtype=np.int32))

# Concatenate into final training arrays
X = np.concatenate(all_MC)
y = np.concatenate(all_y)

#splitting with  Holdout method for eval_set
X_train, X_test, y_train, y_test = train_test_split(X, y,
                                                    test_size=0.20,
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

# xgbm_final = bdt

# Get feature importance scores
print(xgbm_final.feature_importances_)

# Make predictions on the test set
y_pred_proba = xgbm_final.predict_proba(X_test)[:, 1]

# Calculate the ROC AUC score
roc_auc = roc_auc_score(y_test, y_pred_proba)

print(f"ROC AUC Score: {roc_auc:.2f}")

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

# Apply BDT to all DataFrames that contain the required Variables
for key in DataFrames.keys():
    df = DataFrames[key]
    
    # Check: make sure all input BDT variables exist in this DataFrame
    if all(var in df.columns for var in Variables):
        # Apply BDT and store the result
        DataFrames[key]["Ds_BkgBDT"] = xgbm_final.predict_proba(df[Variables])[:, 1].astype(np.float32)

from Functions import optimize_cut, plot_save

cut = optimize_cut(
    df_sig=DataFrames["Signal"],                  # used for plotting signal vs background
    df_bkg=DataFrames["All"],
    Signal=DataFrames["Signal"],                  # used for FoM numerator (truth-matched signal)
    Background=DataFrames["All"],                 # used for FoM denominator (everything else)
    var="Ds_BkgBDT",                              # new classifier variable
    FoM="Ds_BkgBDT",                              # same as var here
    xlabel="Background Classifier Output",
    Bins=50,
    Range=[0, 1],
    varmin=0,
    varmax=0.98,
    select="right",                               # keep events with higher classifier output
    Width=False,
    query_signal="Ds_isSignal == 1"
)

print(f"Best cut is: {cut:.3f}")

import os
import uproot
import numpy as np  # Make sure this is imported if you're working in a standalone script

# === Samples to process ===
samples = ["Signal", "BB", "ccbar", "ddbar", "ssbar", "taupair", "uubar"]

# === Output directory for new Bkg BDT files ===
output_dir = "/group/belle/users/amubarak/03-ML/BkgBDT/"
os.makedirs(output_dir, exist_ok=True)

# === Base input info used to construct filenames ===
base_input_dir = "/group/belle/users/amubarak/02-Grid/Sample_Grid"
Date = "0530"
Attempt = "0"

# === Loop over samples and write output ROOT files ===
for s in samples:
    if s not in DataFrames:
        print(f"Warning: {s} not in DataFrames — skipping.")
        continue

    # Convert Ds_BkgBDT to float32 if present
    if "Ds_BkgBDT" in DataFrames[s].columns:
        DataFrames[s]["Ds_BkgBDT"] = DataFrames[s]["Ds_BkgBDT"].astype(np.float32)

    # Construct the original file name
    if s == "Signal":
        original_name = "Ds2D0enu-Signal.root"
    else:
        original_name = f"Ds2D0e-Generic_Ds_{Date}25_{Attempt}_{s}.root"

    # Build output file name with BkgBDT tag
    output_name = original_name.replace(".root", "_withBkgBDT.root")
    out_path = os.path.join(output_dir, output_name)

    # Save DataFrame to ROOT
    with uproot.recreate(out_path) as f:
        f["Dstree"] = DataFrames[s]

    print(f"Saved: {out_path}")

import os
import uproot
import numpy as np  # Required for dtype conversion

# === Define wrong-charge samples ===
samples_WCh = ["Signal_WCh", "BB_WCh", "ccbar_WCh", "ddbar_WCh", "ssbar_WCh", "taupair_WCh", "uubar_WCh", "Data_WCh"]

# === Output directory for BkgBDT files (wrong charge) ===
output_dir_WCh = "/group/belle/users/amubarak/03-ML/BkgBDT_WCh/"
os.makedirs(output_dir_WCh, exist_ok=True)

# === Base input path for wrong-charge files ===
base_input_dir_WCh = "/group/belle/users/amubarak/02-Grid/Sample_Grid_WCh"
Date_WCh = "0630"
Attempt_WCh = "0"

# === Save each wrong-charge DataFrame with updated Ds_BkgBDT ===
for s in samples_WCh:
    if s not in DataFrames:
        print(f"Warning: {s} not in DataFrames — skipping.")
        continue

    # Convert Ds_BkgBDT to float32 if present
    if "Ds_BkgBDT" in DataFrames[s].columns:
        DataFrames[s]["Ds_BkgBDT"] = DataFrames[s]["Ds_BkgBDT"].astype(np.float32)

    # Set original file name
    if s == "Signal_WCh":
        original_name = "Ds2D0enu-Signal_WCh.root"
    else:
        tag = s.replace("_WCh", "")
        original_name = f"Ds2D0e-Generic_Ds_{Date_WCh}25_{Attempt_WCh}_{tag}.root"

    # Build output path with _withBkgBDT suffix
    output_name = original_name.replace(".root", "_withBkgBDT.root")
    out_path = os.path.join(output_dir_WCh, output_name)

    # Save DataFrame to ROOT
    with uproot.recreate(out_path) as f:
        f["Dstree"] = DataFrames[s]

    print(f"Saved: {out_path}")

import os
import uproot
import numpy as np

# === Helper: Chunked ROOT writer ===
def save_df_in_chunks(df, output_path, chunk_size=500_000):
    # Define tree structure based on column data types
    tree_dict = {}
    for col in df.columns:
        dtype = df[col].dtype
        if np.issubdtype(dtype, np.integer):
            tree_dict[col] = "int32"
        elif np.issubdtype(dtype, np.floating):
            tree_dict[col] = "float64"
        else:
            continue  # skip strings/objects/etc.

    with uproot.recreate(output_path) as f:
        tree = f.mktree("Dstree", tree_dict)

        for start in range(0, len(df), chunk_size):
            end = min(start + chunk_size, len(df))
            chunk = df.iloc[start:end]
            tree.extend({col: chunk[col].values for col in tree_dict})


# === Define ReverseID samples ===
samples_ReverseID = [
    "Signal_ReverseID", "BB_ReverseID", "ccbar_ReverseID",
    "ddbar_ReverseID", "ssbar_ReverseID", "taupair_ReverseID",
    "uubar_ReverseID", "Data_ReverseID"
]

# === Output directory for BkgBDT ReverseID files ===
output_dir_ReverseID = "/group/belle/users/amubarak/03-ML/BkgBDT_ReverseID/"
os.makedirs(output_dir_ReverseID, exist_ok=True)

# === Grid metadata ===
Date_ReverseID = "0626"
Attempt_ReverseID = "0"

# === Save each ReverseID DataFrame with Ds_BkgBDT ===
for s in samples_ReverseID:
    if s not in DataFrames:
        print(f"Warning: {s} not in DataFrames — skipping.")
        continue

    # Convert Ds_BkgBDT to float32 if present
    if "Ds_BkgBDT" in DataFrames[s].columns:
        DataFrames[s]["Ds_BkgBDT"] = DataFrames[s]["Ds_BkgBDT"].astype(np.float32)

    # Determine original file name
    if s == "Signal_ReverseID":
        original_name = "Ds2D0enu-Signal_ReverseID.root"
    else:
        tag = s.replace("_ReverseID", "")
        original_name = f"Ds2D0e-Generic_Ds_{Date_ReverseID}25_{Attempt_ReverseID}_{tag}.root"

    # Output file with BDT tag
    output_name = original_name.replace(".root", "_withBkgBDT.root")
    out_path = os.path.join(output_dir_ReverseID, output_name)

    print(f"Saving in chunks: {out_path}")
    save_df_in_chunks(DataFrames[s], out_path)

# import os
# import uproot
# import numpy as np  # Ensure this is available for dtype conversion

# # === Define ReverseID samples ===
# samples_ReverseID_WCh = ["All_ReverseID_WCh", "BB_ReverseID_WCh", "ccbar_ReverseID_WCh", "ddbar_ReverseID_WCh", "ssbar_ReverseID_WCh", "taupair_ReverseID_WCh", "uubar_ReverseID_WCh", "uds_ReverseID_WCh", "Data_ReverseID_WCh"]

# # === Output directory for BkgBDT ReverseID files ===
# output_dir_ReverseID_WCh = "/group/belle/users/amubarak/03-ML/BkgBDT_ReverseID_WCh/"
# os.makedirs(output_dir_ReverseID_WCh, exist_ok=True)

# # === Base input path for ReverseID ===
# base_input_dir_ReverseID_WCh = "/group/belle/users/amubarak/02-Grid/Sample_Grid_ReverseID_WCh"
# Date_ReverseID_WCh = "0708"
# Attempt_ReverseID_WCh = "0"

# # === Save each ReverseID DataFrame with Ds_BkgBDT ===
# for s in samples_ReverseID_WCh:
#     if s not in DataFrames:
#         print(f"Warning: {s} not in DataFrames — skipping.")
#         continue

#     # Convert Ds_BkgBDT to float32 if present
#     if "Ds_BkgBDT" in DataFrames[s].columns:
#         DataFrames[s]["Ds_BkgBDT"] = DataFrames[s]["Ds_BkgBDT"].astype(np.float32)

#     # Set original file name
#     if s == "Signal_ReverseID_WCh":
#         original_name = "Ds2D0enu-Signal_ReverseID_WCh.root"
#     else:
#         tag = s.replace("_ReverseID_WCh", "")
#         original_name = f"Ds2D0e-Generic_Ds_{Date_ReverseID_WCh}25_{Attempt_ReverseID_WCh}_{tag}.root"

#     # Build output file name with _withBkgBDT suffix
#     output_name = original_name.replace(".root", "_withBkgBDT.root")
#     out_path = os.path.join(output_dir_ReverseID_WCh, output_name)

#     # Save DataFrame to ROOT
#     with uproot.recreate(out_path) as f:
#         f["Dstree"] = DataFrames[s]

#     print(f"Saved: {out_path}")