# %% [markdown]
# # Fake D⁰ BDT Training - Multi-Mode Implementation
# 
# ## Overview
# 
# This notebook trains an XGBoost BDT to suppress fake D⁰ backgrounds using the **corrected labeling procedure**:
# 
# - **Real D⁰**: `abs(D0_mcPDG) == 421` from BOTH signal MC AND generic MC
# - **Fake D⁰**: `abs(D0_mcPDG) != 421` or NaN from generic MC
# 
# **Key Features**:
# - Multi-mode support: kmpip, km3pi, kmpippi0_eff20_May2020
# - Uses ALL variables from `final_variables.py` (278, 478, 788 vars per mode)
# - Sequential processing to avoid memory bloat
# - Punzi FoM for cut optimization (90% CL, 95% CL, 3σ)
# - Configurable hyperparameter optimization toggle
# - Control sample support with toggle
# - Comprehensive plotting and diagnostics
# - Automatic saving to `/group/belle/users/amubarak/04-ML/`

# %% [markdown]
# ## System Information
# 
# Check available CPU cores and RAM before processing.

# %%
import os

os.cpu_count()

# %%
import psutil

mem = psutil.virtual_memory()
print(f"Total RAM: {mem.total / 1e9:.2f} GB")
print(f"Available: {mem.available / 1e9:.2f} GB")

# %% [markdown]
# ## Imports

# %%
import sys
import gc
import seaborn as sns

import numpy as np
import pandas as pd
import uproot
from matplotlib import pyplot as plt
from tqdm import tqdm

from sklearn.ensemble import GradientBoostingClassifier
from sklearn.metrics import auc, roc_curve, confusion_matrix, classification_report
from sklearn.metrics import precision_recall_curve, mean_squared_error, accuracy_score, roc_auc_score
from sklearn.model_selection import GridSearchCV, cross_validate, validation_curve
from sklearn.model_selection import train_test_split, KFold, learning_curve, cross_val_score
from sklearn.utils import compute_sample_weight
from scipy.stats import ks_2samp, randint, uniform

import xgboost
from xgboost import XGBClassifier
from sklearn.model_selection import RandomizedSearchCV

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
sys.path.append("/home/belle2/amubarak/Ds2D0enue_Analysis/08-Python_Functions/")
sys.path.append("/home/belle2/amubarak/Ds2D0enue_Analysis/05-ML/Variables/")

# Import custom functions
from Functions import optimize_cut, plot_save
from Ds2D0e_config import DECAY_CONFIG, BACKGROUND_SAMPLES, get_signal_file, get_generic_file
from final_variables import VARIABLES

print("Imported custom functions successfully")

# %% [markdown]
# ## Configuration
# 
# Set the parameters for BDT training, data loading, and output saving.
# 
# **Important toggles**:
# - `LOAD_CONTROL_SAMPLES`: Enable/disable control sample processing (WCh, ReverseID, etc.)
# - `RUN_OPTIMIZATION`: Enable/disable hyperparameter optimization via RandomizedSearchCV
# - `SAVE_OUTPUT`: Enable/disable ROOT file output saving
# - `SAVE_IMAGES`: Enable/disable plot saving

# %%
# ========================================================================
# CONFIGURATION SECTION
# ========================================================================

# === Mode Selection ===
# Set to specific mode: "kmpip", "km3pi", "kmpippi0_eff20_May2020"
# Set to "all" to train on all modes sequentially
TRAIN_MODE = "all"  # Change this to train different modes

# === Control Samples ===
# Toggle for loading control samples (noEID, wrongCharge, reverseID)
LOAD_CONTROL_SAMPLES = False  # Set to True to load control samples
CONTROL_SAMPLES = ["WCh", "ReverseID", "ReverseID_WCh"]  # Which control samples to load

# === Image Saving ===
SAVE_IMAGES = True  # Set to True to save all plots
OUTPUT_DIR = "/home/belle2/amubarak/Ds2D0enue_Analysis/05-ML/Figures/FakeD0_BDT"

# === ROOT File Output Saving ===
SAVE_OUTPUT = False  # Set to True to save DataFrames with BDT to ROOT files
OUTPUT_BASE_DIR = "/group/belle/users/amubarak/04-ML/"

# === Normalization Parameters for Punzi FoM ===
N_SIGNAL_EVENTS = 100_000  # Number of signal MC events generated
LUMINOSITY_FB = 200  # Integrated luminosity in fb^-1
PUNZI_A_VALUES = [1.64, 1.96, 3.0]  # 90% CL, 95% CL, 3σ discovery

# === BDT Training Parameters ===
RANDOM_STATE = 42
TEST_SIZE = 0.30
N_ESTIMATORS = 100

# === Hyperparameter Optimization Toggle ===
RUN_OPTIMIZATION = False  # Set to False to skip hyperparameter optimization
RANDOM_SEARCH_ITERS = 50  # Number of random search iterations (only if RUN_OPTIMIZATION=True)

print(f"Configuration loaded:")
print(f"  Train mode: {TRAIN_MODE}")
print(f"  Load control samples: {LOAD_CONTROL_SAMPLES}")
print(f"  Save images: {SAVE_IMAGES}")
print(f"  Save output: {SAVE_OUTPUT}")
print(f"  Run hyperparameter optimization: {RUN_OPTIMIZATION}")
print(f"  N signal events: {N_SIGNAL_EVENTS:,}")
print(f"  Luminosity: {LUMINOSITY_FB} fb^-1")

# %% [markdown]
# ## Import Data Configuration and Variables
# 
# Load decay mode configurations from `Ds2D0e_config.py` and the final variable lists from `final_variables.py`.
# 
# **Variable counts per mode**:
# - kmpip: 278 variables
# - km3pi: 478 variables  
# - kmpippi0_eff20_May2020: 788 variables
# 
# We use ALL variables from `final_variables.py` without reduction at this stage.

# %%
from Ds2D0e_config import DECAY_CONFIG, BACKGROUND_SAMPLES, get_signal_file, get_generic_file
from final_variables import VARIABLES

print("Available decay modes:")
for mode in DECAY_CONFIG.keys():
    print(f"  - {mode}")

print("\nVariable counts per mode:")
for mode in VARIABLES.keys():
    print(f"  {mode}: {len(VARIABLES[mode]['all_vars'])} variables")

# %% [markdown]
# ## Helper Functions
# 
# ### Punzi Figure of Merit
# 
# The Punzi FoM is used for cut optimization in the presence of systematic uncertainties:
# 
# $$\text{Punzi FoM} = \frac{\epsilon_S}{\frac{a}{2} + \sqrt{B}}$$
# 
# where:
# - $\epsilon_S$ is the signal efficiency
# - $B$ is the number of background events
# - $a$ is the confidence level parameter:
#   - $a = 1.64$ for 90% CL
#   - $a = 1.96$ for 95% CL
#   - $a = 3.0$ for 3σ discovery
# 
# This is preferred over the simple $S/\sqrt{S+B}$ FoM when systematic uncertainties are significant.

# %%
def punzi_fom(epsS, B, a=1.96):
    """
    Compute Punzi figure of merit.
    
    Parameters:
        epsS (float): Signal efficiency
        B (float): Number of background events
        a (float): Confidence level parameter (default 1.96 for 95% CL)
    
    Returns:
        float: Punzi FoM = ε_S / (a/2 + √B)
    """
    denom = (a / 2.0) + np.sqrt(max(B, 0.0))
    return 0.0 if denom <= 0.0 else epsS / denom


def compute_punzi_scan(bdt_scores_sig, bdt_scores_bkg, n_sig_total, n_thresholds=200, a=1.96):
    """
    Scan BDT thresholds and compute Punzi FoM.
    
    Parameters:
        bdt_scores_sig (array): BDT scores for signal events
        bdt_scores_bkg (array): BDT scores for background events
        n_sig_total (int): Total number of signal events (for efficiency calculation)
        n_thresholds (int): Number of thresholds to scan
        a (float): Confidence level parameter
    
    Returns:
        thresholds, foms, best_threshold, best_fom
    """
    thresholds = np.linspace(0, 1, n_thresholds)
    foms = []
    
    for t in thresholds:
        n_sig_pass = np.sum(bdt_scores_sig > t)
        n_bkg_pass = np.sum(bdt_scores_bkg > t)
        
        epsS = n_sig_pass / n_sig_total if n_sig_total > 0 else 0.0
        fom = punzi_fom(epsS, n_bkg_pass, a)
        foms.append(fom)
    
    foms = np.array(foms)
    best_idx = np.argmax(foms)
    return thresholds, foms, thresholds[best_idx], foms[best_idx]


def get_pulls(counts, errors, pdf):
    """
    Compute pull values for comparing test vs train distributions.
    
    Parameters:
        counts (array): Test counts
        errors (array): Test errors (uncertainties)
        pdf (array): Train probability density (normalized)
    
    Returns:
        array: Pull values = (counts - pdf) / errors
    """
    pull = (counts - pdf) / errors
    return pull

# %% [markdown]
# ### Data Loading Functions
# 
# Functions to load signal and generic MC for each decay mode, applying D⁰ mass window cuts from the configuration.
# 
# **Data sources**:
# - Signal MC: Ds → D⁰ e ν samples
# - Generic MC: Combined generic BB̄ backgrounds (charged, mixed, uubar, ddbar, ssbar, ccbar)

# %%
# Mode titles with proper LaTeX formatting for decay chains
MODE_TITLES = {
    "kmpip": r"$D_s^{+} \rightarrow [D^{0} \rightarrow K^{-} \pi^{+}] e^{+} \nu_{e}$",
    "kmpippi0_eff20_May2020": r"$D_s^{+} \rightarrow [D^{0} \rightarrow K^{-} \pi^{+} \pi^{0}] e^{+} \nu_{e}$",
    "km3pi": r"$D_s^{+} \rightarrow [D^{0} \rightarrow K^{-} 3\pi] e^{+} \nu_{e}$",
}

print("Mode titles loaded:")
for mode, title in MODE_TITLES.items():
    print(f"  {mode}: {title}")

# %%
def plot_feature_importance(model, features, mode, num=None, save_dir=None):
    """
    Plot feature importance from trained XGBoost model and save all importances to CSV.
    
    Parameters:
        model: Trained XGBoost classifier
        features: List of feature names
        mode: Decay mode (for title)
        num: Number of top features to show in plot (None = show all)
        save_dir: Directory to save plot and CSV (None = show plot)
    """
    feature_imp = pd.DataFrame({
        'Value': model.feature_importances_,
        'Feature': features
    })
    
    # Sort by importance
    feature_imp_sorted = feature_imp.sort_values(by="Value", ascending=False)
    
    # Save complete table to CSV
    if save_dir:
        csv_path = os.path.join(save_dir, "feature_importance_all.csv")
        feature_imp_sorted.to_csv(csv_path, index=False)
        print(f"\n✓ Saved complete feature importance table to: {csv_path}")
    
    # If num is specified, take top N for plotting, otherwise take all
    if num is not None:
        feature_imp_plot = feature_imp_sorted.head(num)
        title_suffix = f"Top {num} Feature Importances"
    else:
        feature_imp_plot = feature_imp_sorted
        title_suffix = f"Feature Importances (All {len(features)} variables)"
    
    # Adjust figure size based on number of features
    n_features = len(feature_imp_plot)
    fig_height = max(8, min(100, n_features * 0.3))  # At least 8, at most 100
    
    plt.figure(figsize=(16, fig_height))
    sns.barplot(
        x="Value", y="Feature",
        data=feature_imp_plot
    )
    
    title = MODE_TITLES.get(mode, mode)
    plt.title(f'{title}\n{title_suffix}', loc='left')
    plt.tight_layout()
    
    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
        if num is not None:
            plt.savefig(os.path.join(save_dir, f"feature_importance_top{num}.png"), dpi=150, bbox_inches='tight')
        else:
            plt.savefig(os.path.join(save_dir, "feature_importance_all.png"), dpi=150, bbox_inches='tight')
        plt.close()
    else:
        plt.show()
    
    return feature_imp_sorted

# %%
def compare_train_test(clf, X_train, y_train, X_test, y_test, mode):
    """
    Plot BDT output comparing train vs test with pull plots.
    
    Parameters:
        clf: Trained classifier
        X_train, y_train: Training data and labels
        X_test, y_test: Test data and labels
        mode: Decay mode (for title)
    
    Returns:
        decisions: List of [train_bkg, train_sig, test_bkg, test_sig] BDT outputs
    """
    decisions = []  # list to hold decisions of classifier
    for X, y in ((X_train, y_train), (X_test, y_test)):  # train and test
        if hasattr(clf, "predict_proba"):  # if predict_proba function exists
            d1 = clf.predict_proba(X[y<0.5])[:, 1]  # background
            d2 = clf.predict_proba(X[y>0.5])[:, 1]  # signal
        else:  # predict_proba function doesn't exist
            X_tensor = torch.as_tensor(X, dtype=torch.float)
            y_tensor = torch.as_tensor(y, dtype=torch.long)
            X_var, y_var = Variable(X_tensor), Variable(y_tensor)
            d1 = clf(X_var[y_var<0.5])[1][:, 1].cpu().detach().numpy()  # background
            d2 = clf(X_var[y_var>0.5])[1][:, 1].cpu().detach().numpy()  # signal
        decisions += [d1, d2]  # add to list of classifier decision

    lw = 3
    fig, axs = plt.subplots(3, 1, figsize=(10, 10), gridspec_kw={'height_ratios':[1, 0.2, 0.2]})

    bins = 50
    bin_edges = np.linspace(0, 1, bins)
    
    test_bkg_count_weight = bins / len(decisions[2])
    test_sig_count_weight = bins / len(decisions[3])
    test_bkg_counts, test_bkg_bins = np.histogram(decisions[2], bins=bins, range=(0, 1))
    test_sig_counts, test_sig_bins = np.histogram(decisions[3], bins=bins, range=(0, 1))

    train_bkg_counts, train_bkg_bins, _etc = axs[0].hist(
        decisions[0], color='tab:blue', histtype='step', bins=bins,
        density=True, range=(0, 1), linewidth=lw, label='Train Background'
    )
    train_sig_counts, train_sig_bins, _etc = axs[0].hist(
        decisions[1], color='tab:red', histtype='step', bins=bins,
        density=True, range=(0, 1), linewidth=lw, label=r'Train Signal'
    )
    axs[0].hist(decisions[0], color='tab:blue', histtype='stepfilled',
                alpha=0.4, bins=bins, density=True, range=(0, 1))
    axs[0].hist(decisions[1], color='tab:red', histtype='stepfilled',
                alpha=0.4, bins=bins, density=True, range=(0, 1))
    
    bin_width = test_bkg_bins[1] - test_bkg_bins[0]
    bin_centers = [el + (bin_width/2) for el in test_bkg_bins[:-1]]

    axs[0].errorbar(
        bin_centers, test_bkg_count_weight * test_bkg_counts,
        yerr=test_bkg_count_weight * np.sqrt(test_bkg_counts),
        label='Test Background', color='tab:blue', marker='o', linewidth=lw, ls=''
    )
    axs[0].errorbar(
        bin_centers, test_sig_count_weight * test_sig_counts,
        yerr=test_sig_count_weight * np.sqrt(test_sig_counts),
        label='Test Signal', color='tab:red', marker='o', linewidth=lw, ls=''
    )
    
    # Title with proper LaTeX decay chain
    title = MODE_TITLES.get(mode, mode)
    axs[0].set_title(title, loc='left')
    axs[0].set_xlim(0, 1)
    axs[0].set_ylim(0)
    axs[0].set_ylabel('Event Density')

    # K-S test scores
    ks_p_value_sig = ks_2samp(decisions[1], decisions[3])[1]
    ks_p_value_bkg = ks_2samp(decisions[0], decisions[2])[1]

    leg = axs[0].legend(
        loc='upper center',
        title=f"Sig K-S test score: {ks_p_value_sig:0.3f}\nBkg K-S test score: {ks_p_value_bkg:0.3f}"
    )
    leg._legend_box.align = "left"

    # Background pulls
    pulls = get_pulls(
        test_bkg_count_weight * test_bkg_counts,
        test_bkg_count_weight * np.sqrt(test_bkg_counts),
        np.array(train_bkg_counts)
    )
    axs[1].bar(bin_centers, pulls, width=bin_width)
    axs[1].set_xlim(0, 1)
    axs[1].set_ylabel('Pulls')
    axs[1].set_ylim(-5, 5)

    # Signal pulls
    pulls = get_pulls(
        test_sig_count_weight * test_sig_counts,
        test_sig_count_weight * np.sqrt(test_sig_counts),
        np.array(train_sig_counts)
    )
    axs[2].bar(bin_centers, pulls, width=bin_width, color='tab:red')
    axs[2].set_xlim(0, 1)
    axs[2].set_ylabel('Pulls')
    axs[2].set_ylim(-5, 5)
    axs[2].set_xlabel(r'BDT output')

    return decisions

# %%
def filter_valid_variables(df, variable_list, nan_threshold=0.5):
    """
    Filter variable list to only include valid columns from dataframe.
    
    Removes variables that:
    - Don't exist in the dataframe
    - Are not numeric
    - Have >50% NaN values (by default)
    - Have any infinity values
    
    Parameters:
        df: DataFrame to check
        variable_list: List of variable names to filter
        nan_threshold: Maximum fraction of NaN values allowed (0.5 = 50%)
    
    Returns:
        List of valid variable names
    """
    valid_vars = []
    
    for var in variable_list:
        # Check if variable exists
        if var not in df.columns:
            continue
        
        # Check if numeric
        if not pd.api.types.is_numeric_dtype(df[var]):
            continue
        
        col_data = df[var]
        
        # Check NaN fraction
        nan_fraction = col_data.isna().sum() / len(col_data)
        if nan_fraction > nan_threshold:
            print(f"  ⚠ Skipping {var}: {nan_fraction*100:.1f}% NaN values")
            continue
        
        # Check for infinity
        if np.isinf(col_data.dropna()).any():
            print(f"  ⚠ Skipping {var}: contains infinity values")
            continue
        
        valid_vars.append(var)
    
    return valid_vars


def load_mode_data(mode, control_sample=None):
    """
    Load signal and generic MC for a given decay mode.
    Applies D⁰ mass window cuts from config.
    
    Parameters:
        mode (str): Decay mode (kmpip, km3pi, kmpippi0_eff20_May2020)
        control_sample (str): Control sample tag ("WCh", "ReverseID", etc.) or None for nominal
    
    Returns:
        df_signal, df_generic (DataFrames)
    """
    config = DECAY_CONFIG[mode]
    tree_name = config["ds_tree"]
    mass_cut = config["cut"]
    
    print(f"\n{'='*80}")
    print(f"Loading data for mode: {mode}")
    if control_sample:
        print(f"Control sample: {control_sample}")
    print(f"Tree: {tree_name}")
    print(f"D⁰ mass cut: {mass_cut}")
    print(f"{'='*80}\n")
    
    # Load signal MC
    use_control = control_sample is not None
    control_tag = control_sample if control_sample else None
    
    signal_path = get_signal_file(mode, use_control_sample=use_control, control_sample_tag=control_tag)
    if isinstance(signal_path, list):
        signal_path = signal_path[0]  # Take first if multiple returned
    
    print(f"Loading Signal MC from: {signal_path}")
    try:
        df_signal = uproot.concatenate(f"{signal_path}:{tree_name}", library='pd')
        df_signal = df_signal.query(mass_cut)  # Apply D⁰ mass window
        print(f"  → Signal MC: {len(df_signal):,} events after D⁰ mass cut")
    except Exception as e:
        print(f"  ✗ Failed to load signal: {e}")
        df_signal = pd.DataFrame()
    
    # Load generic MC (all background samples combined)
    df_generic_list = []
    for sample in tqdm(BACKGROUND_SAMPLES, desc="Loading generic MC"):
        generic_path = get_generic_file(sample, mode, use_control_sample=use_control, 
                                       control_sample_tag=control_tag)
        if isinstance(generic_path, list):
            generic_path = generic_path[0]
        
        try:
            df = uproot.concatenate(f"{generic_path}:{tree_name}", library='pd')
            df = df.query(mass_cut)  # Apply D⁰ mass window
            df_generic_list.append(df)
            print(f"  → {sample}: {len(df):,} events after D⁰ mass cut")
        except Exception as e:
            print(f"  ✗ Failed to load {sample}: {e}")
    
    df_generic = pd.concat(df_generic_list, ignore_index=True) if df_generic_list else pd.DataFrame()
    print(f"\n  → Total Generic MC: {len(df_generic):,} events\n")
    
    return df_signal, df_generic

# %%
def get_pulls(test_data, test_errors, train_data):
    """
    Calculate pulls for overtraining check.
    
    Parameters:
        test_data: Test histogram values
        test_errors: Test histogram errors
        train_data: Train histogram values
    
    Returns:
        pulls: (test_data - train_data) / test_errors array
    """
    # Avoid division by zero
    pulls = np.zeros(len(test_data))
    mask = test_errors > 0
    pulls[mask] = (test_data[mask] - train_data[mask]) / test_errors[mask]
    return pulls


def compute_punzi_scan(signal_scores, background_scores, n_sig_events, n_thresholds=100, a=1.96):
    """
    Scan BDT cuts and compute Punzi Figure of Merit (FoM).
    
    Punzi FoM = n_sig / (a + sqrt(n_bkg))
    
    Parameters:
        signal_scores: Array of BDT scores for signal events
        background_scores: Array of BDT scores for background events
        n_sig_events: Number of signal events (for FoM calculation)
        n_thresholds: Number of threshold points to scan
        a: Confidence level parameter (default 1.96 for 95% CL)
    
    Returns:
        thresholds: Array of BDT cut values
        foms: Array of FoM values
        best_threshold: Optimal BDT cut
        best_fom: Maximum FoM value
    """
    thresholds = np.linspace(0, 1, n_thresholds)
    foms = []
    
    for cut in thresholds:
        # Count events passing cut
        n_sig = np.sum(signal_scores > cut)
        n_bkg = np.sum(background_scores > cut)
        
        # Compute FoM
        if n_bkg < 0:
            n_bkg = 0
        
        # Scale FoM to full dataset
        if n_sig > 0:
            fom = n_sig / np.sqrt(a**2 + n_bkg) if (a**2 + n_bkg) > 0 else 0
        else:
            fom = 0
        
        foms.append(fom)
    
    foms = np.array(foms)
    best_idx = np.argmax(foms)
    best_threshold = thresholds[best_idx]
    best_fom = foms[best_idx]
    
    return thresholds, foms, best_threshold, best_fom


def print_bdt_statistics(model, X_train, y_train, X_test, y_test):
    """
    Print BDT training and test statistics.
    
    Parameters:
        model: Trained XGBoost classifier
        X_train, y_train: Training features and labels
        X_test, y_test: Test features and labels
    """
    print("\n" + "="*80)
    print("BDT STATISTICS")
    print("="*80)
    
    # Get predictions
    y_pred_train = model.predict(X_train)
    y_pred_test = model.predict(X_test)
    
    y_proba_train = model.predict_proba(X_train)[:, 1]
    y_proba_test = model.predict_proba(X_test)[:, 1]
    
    # Training metrics
    train_acc = accuracy_score(y_train, y_pred_train)
    train_auc = roc_auc_score(y_train, y_proba_train)
    
    # Test metrics
    test_acc = accuracy_score(y_test, y_pred_test)
    test_auc = roc_auc_score(y_test, y_proba_test)
    
    # Compute precision and recall
    precision_train = confusion_matrix(y_train, y_pred_train)[1, 1] / (confusion_matrix(y_train, y_pred_train)[1, 1] + confusion_matrix(y_train, y_pred_train)[0, 1]) if (confusion_matrix(y_train, y_pred_train)[1, 1] + confusion_matrix(y_train, y_pred_train)[0, 1]) > 0 else 0
    precision_test = confusion_matrix(y_test, y_pred_test)[1, 1] / (confusion_matrix(y_test, y_pred_test)[1, 1] + confusion_matrix(y_test, y_pred_test)[0, 1]) if (confusion_matrix(y_test, y_pred_test)[1, 1] + confusion_matrix(y_test, y_pred_test)[0, 1]) > 0 else 0
    
    print(f"\nTraining Set Performance:")
    print(f"  Accuracy: {train_acc:.4f}")
    print(f"  ROC AUC:  {train_auc:.4f}")
    print(f"\nTest Set Performance:")
    print(f"  Accuracy: {test_acc:.4f}")
    print(f"  ROC AUC:  {test_auc:.4f}")
    print(f"\nOvertraining Check (Test - Train):")
    print(f"  ΔAccuracy: {test_acc - train_acc:+.4f}")
    print(f"  ΔAUC:      {test_auc - train_auc:+.4f}")
    
    # Check for overtraining
    if abs(test_auc - train_auc) < 0.01 and abs(test_acc - train_acc) < 0.01:
        print(f"  ✓ Good agreement (no overtraining detected)")
    elif test_auc < train_auc or test_acc < train_acc:
        print(f"  ⚠ Slight overtraining detected")
    else:
        print(f"  ✓ Test performance better than train (normal for some cases)")
    
    print("="*80 + "\n")

# %% [markdown]
# ### Training Set Creation with Corrected Labeling
# 
# **CORRECTED PROCEDURE**:
# 
# The key insight is that we want the BDT to learn "is this a real D⁰?" rather than "is this from signal MC?".
# 
# - **Real D⁰ (label = 1)**: Events with `abs(D0_mcPDG) == 421` from BOTH:
#   - Signal MC (Ds → D⁰ e ν samples)
#   - Generic MC (combinatoric backgrounds that happen to form real D⁰)
# 
# - **Fake D⁰ (label = 0)**: Events from generic MC with:
#   - `abs(D0_mcPDG) != 421` (combinatoric background)
#   - `D0_mcPDG` is NaN (no truth match)
# 
# This ensures the BDT learns real vs fake D⁰ topology, not signal vs background event characteristics.

# %%
def create_training_set(df_signal, df_generic, mode_variables):
    """
    Create training set with CORRECTED LABELING PROCEDURE.
    
    Real D⁰: abs(D0_mcPDG) == 421 from BOTH signal MC AND generic MC
    Fake D⁰: abs(D0_mcPDG) != 421 or NaN from generic MC
    
    Parameters:
        df_signal: Signal MC DataFrame
        df_generic: Generic MC DataFrame
        mode_variables: List of variable names to use
    
    Returns:
        X, y, df_train, df_real, df_fake, mode_variables_filtered
    """
    print("\n" + "="*80)
    print("Creating training set with CORRECTED labeling")
    print("="*80)
    
    # Real D⁰ from signal MC
    real_signal_mask = (abs(df_signal["D0_mcPDG"]) == 421)
    df_real_signal = df_signal[real_signal_mask]
    print(f"Real D⁰ from Signal MC: {len(df_real_signal):,} events")
    
    # Real D⁰ from generic MC
    real_generic_mask = (abs(df_generic["D0_mcPDG"]) == 421)
    df_real_generic = df_generic[real_generic_mask]
    print(f"Real D⁰ from Generic MC: {len(df_real_generic):,} events")
    
    # Combine real D⁰
    df_real = pd.concat([df_real_signal, df_real_generic], ignore_index=True)
    print(f"Total Real D⁰: {len(df_real):,} events")
    
    # Fake D⁰ from generic MC
    fake_mask = (abs(df_generic["D0_mcPDG"]) != 421) | df_generic["D0_mcPDG"].isna()
    df_fake = df_generic[fake_mask]
    print(f"Fake D⁰ from Generic MC: {len(df_fake):,} events")
    
    # Check variable availability
    missing_vars = [v for v in mode_variables if v not in df_real.columns]
    if missing_vars:
        print(f"\n⚠ Warning: {len(missing_vars)} variables not available in data:")
        for v in missing_vars[:10]:  # Show first 10
            print(f"  - {v}")
        if len(missing_vars) > 10:
            print(f"  ... and {len(missing_vars)-10} more")
        
        # Filter to available variables
        mode_variables_filtered = [v for v in mode_variables if v in df_real.columns]
        print(f"\nProceeding with {len(mode_variables_filtered)} available variables\n")
    else:
        mode_variables_filtered = mode_variables
    
    # Combine
    df_train = pd.concat([df_real, df_fake], axis=0, ignore_index=True)
    
    # Labels: 1 for real D⁰, 0 for fake D⁰
    labels = np.concatenate([
        np.ones(len(df_real), dtype=np.int64),
        np.zeros(len(df_fake), dtype=np.int64)
    ])
    
    # Features
    X = df_train[mode_variables_filtered].to_numpy(dtype=np.float32)
    y = labels
    
    # Clean data: replace inf with NaN
    X = np.where(np.isinf(X), np.nan, X)
    
    # Check for NaN/inf values
    n_nan = np.sum(np.isnan(X))
    n_inf = np.sum(np.isinf(X))
    
    if n_nan > 0 or n_inf > 0:
        print(f"\n⚠ Data cleaning:")
        print(f"  - NaN values found: {n_nan:,}")
        print(f"  - Inf values found: {n_inf:,}")
        print(f"  - NaN/Inf values will be handled as missing by XGBoost")
    
    print(f"\nTraining set: {len(df_train):,} events")
    print(f"  - Real D⁰ (label=1): {np.sum(y==1):,}")
    print(f"  - Fake D⁰ (label=0): {np.sum(y==0):,}")
    print(f"  - Class ratio (Real/Fake): {np.sum(y==1)/np.sum(y==0):.3f}")
    print(f"  - Number of features: {len(mode_variables_filtered)}")
    print("="*80 + "\n")
    
    return X, y, df_train, df_real, df_fake, mode_variables_filtered

# %%
def train_bdt(X, y, run_optimization=True):
    """
    Train XGBoost BDT with optional hyperparameter optimization.
    
    Automatically handles class imbalance using scale_pos_weight parameter.
    
    Parameters:
        X: Feature matrix
        y: Labels
        run_optimization: Whether to run RandomizedSearchCV
    
    Returns:
        bdt_final, X_train, X_test, y_train, y_test
    """
    print("\n" + "="*80)
    print("TRAINING BDT")
    print("="*80)
    
    # Split data
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=TEST_SIZE, random_state=RANDOM_STATE, stratify=y
    )
    
    print(f"Training set: {len(X_train):,} events")
    print(f"Test set: {len(X_test):,} events")
    
    # Calculate class weights to handle imbalance
    n_real = np.sum(y_train == 1)
    n_fake = np.sum(y_train == 0)
    scale_pos_weight = n_fake / n_real if n_real > 0 else 1.0
    
    print(f"\nClass balance in training set:")
    print(f"  Real D⁰ (positive, label=1): {n_real:,} ({100*n_real/len(y_train):.1f}%)")
    print(f"  Fake D⁰ (negative, label=0): {n_fake:,} ({100*n_fake/len(y_train):.1f}%)")
    print(f"  Class ratio (Fake/Real): {n_fake/n_real:.2f}")
    print(f"  XGBoost scale_pos_weight: {scale_pos_weight:.2f} (for balancing)")
    
    if run_optimization:
        print(f"\nRunning hyperparameter optimization with {RANDOM_SEARCH_ITERS} iterations...")
        
        # Define base estimator with scale_pos_weight for class imbalance
        bdt_base = XGBClassifier(
            objective="binary:logistic",
            eval_metric="logloss",
            max_delta_step=1,
            scale_pos_weight=scale_pos_weight,  # Handle class imbalance
            random_state=RANDOM_STATE,
            missing=np.nan  # Handle NaN values as missing
        )
        
        # Hyperparameter search space
        param_dist = {
            'learning_rate': uniform(0.01, 0.19),
            'max_depth': randint(1, 6),
            'n_estimators': randint(100, 201),
            'reg_lambda': uniform(1, 4),
            'gamma': uniform(0, 4),
            'subsample': uniform(0.5, 0.5),
            'min_child_weight': randint(1, 7),
            'colsample_bytree': uniform(0.3, 0.7)
        }
        
        # Random search with cross-validation
        random_search = RandomizedSearchCV(
            bdt_base,
            param_distributions=param_dist,
            n_iter=RANDOM_SEARCH_ITERS,
            cv=5,
            scoring='roc_auc',
            n_jobs=-1,
            random_state=RANDOM_STATE,
            verbose=1
        )
        
        random_search.fit(X_train, y_train)
        
        print("\nBest hyperparameters:")
        for param, value in random_search.best_params_.items():
            print(f"  {param}: {value}")
        
        print(f"\nBest CV score (ROC AUC): {random_search.best_score_:.4f}")
        
        bdt_final = random_search.best_estimator_
    else:
        print("\nTraining with default parameters (no optimization)...")
        print(f"Using scale_pos_weight={scale_pos_weight:.2f} to handle class imbalance")
        
        bdt_final = XGBClassifier(
            objective="binary:logistic",
            eval_metric="logloss",
            max_delta_step=1,
            scale_pos_weight=scale_pos_weight,  # Handle class imbalance
            random_state=RANDOM_STATE,
            n_estimators=N_ESTIMATORS,
            missing=np.nan  # Handle NaN values as missing
        )
        
        bdt_final.fit(X_train, y_train)
    
    print("\n✓ BDT training complete")
    print("="*80 + "\n")
    
    return bdt_final, X_train, X_test, y_train, y_test

# %%
def apply_bdt_and_save(bdt, mode_variables, mode, df_signal, df_generic):
    """
    Apply trained BDT to signal and generic MC, and save to ROOT files.
    
    Parameters:
        bdt: Trained XGBoost classifier
        mode_variables: List of variable names used in training
        mode: Decay mode
        df_signal: Signal MC DataFrame
        df_generic: Generic MC DataFrame
    """
    import uproot
    
    if not SAVE_OUTPUT:
        print("\nSAVE_OUTPUT is False. Skipping file saving.")
        return
    
    print("\n" + "="*80)
    print("SAVING OUTPUT FILES")
    print("="*80)
    
    # Create output directory
    output_dir = os.path.join(OUTPUT_BASE_DIR, "FakeD0_BDT", mode)
    os.makedirs(output_dir, exist_ok=True)
    
    # Save signal with BDT
    signal_path = os.path.join(output_dir, f"Ds2D0enu-Signal_{mode}_withBDT.root")
    print(f"\nSaving signal MC to: {signal_path}")
    
    config = DECAY_CONFIG[mode]
    tree_name = config["ds_tree"]
    
    with uproot.recreate(signal_path) as f:
        f[tree_name] = df_signal
    
    print(f"  ✓ Saved {len(df_signal):,} signal events")
    
    # Save generic with BDT
    generic_path = os.path.join(output_dir, f"Ds2D0e-Generic_{mode}_withBDT.root")
    print(f"\nSaving generic MC to: {generic_path}")
    
    with uproot.recreate(generic_path) as f:
        f[tree_name] = df_generic
    
    print(f"  ✓ Saved {len(df_generic):,} generic events")
    print("="*80 + "\n")

# %% [markdown]
# ## Model Training Function
# 
# Train XGBoost BDT with optional hyperparameter optimization.
# 
# **Hyperparameter search space** (when `RUN_OPTIMIZATION = True`):
# - `learning_rate`: [0.01, 0.2]
# - `max_depth`: [1, 5]
# - `n_estimators`: [100, 200]
# - `reg_lambda`: [1, 5] (L2 regularization)
# - `gamma`: [0, 4] (min loss reduction for split)
# - `subsample`: [0.5, 1.0]
# - `min_child_weight`: [1, 6]
# - `colsample_bytree`: [0.3, 1.0]
# 
# Uses `RandomizedSearchCV` with 5-fold cross-validation.

# %%
def process_mode(mode):
    """
    Complete BDT training pipeline for a single mode.
    
    Parameters:
        mode: Decay mode to process
    """
    print("\n" + "█"*80)
    print(f"PROCESSING MODE: {mode}")
    print("█"*80 + "\n")
    
    # Get proper title
    mode_title = MODE_TITLES.get(mode, mode)
    print(f"Decay chain: {mode_title}")
    
    # 1. Get variable list for this mode
    if mode not in VARIABLES:
        print(f"ERROR: No variables defined for mode {mode}")
        return
    
    mode_variables = VARIABLES[mode]["all_vars"]
    print(f"Using {len(mode_variables)} variables from final_variables.py")
    
    # 2. Load data
    df_signal, df_generic = load_mode_data(mode)
    
    if df_signal.empty or df_generic.empty:
        print(f"ERROR: Failed to load data for {mode}")
        return
    
    # 3. Create training set
    X, y, df_train, df_real, df_fake, mode_variables_filtered = create_training_set(
        df_signal, df_generic, mode_variables
    )
    
    # 4. Train BDT
    bdt_final, X_train, X_test, y_train, y_test = train_bdt(X, y, run_optimization=RUN_OPTIMIZATION)
    
    # 5. BDT Statistics and overtraining metrics
    print_bdt_statistics(bdt_final, X_train, y_train, X_test, y_test)
    
    # 6. Feature Importance
    if SAVE_IMAGES:
        save_dir = os.path.join(OUTPUT_DIR, mode)
        os.makedirs(save_dir, exist_ok=True)
    else:
        save_dir = None
    
    print("\n" + "="*80)
    print("FEATURE IMPORTANCE")
    print(f"\nTotal variables being used for training: {len(mode_variables_filtered)}")
    print("\nVerifying all variables from final_variables.py are included...")
    
    # Check if all expected variables are present
    expected_vars = set(VARIABLES[mode]["all_vars"])
    actual_vars = set(mode_variables_filtered)
    
    missing_vars = expected_vars - actual_vars
    if missing_vars:
        print(f"  ⚠ {len(missing_vars)} variables from final_variables.py not found in data:")
        for v in list(missing_vars)[:10]:
            print(f"    - {v}")
        if len(missing_vars) > 10:
            print(f"    ... and {len(missing_vars)-10} more")
    else:
        print(f"  ✓ All {len(expected_vars)} variables from final_variables.py are being used")
    print("="*80)
    feature_imp_df = plot_feature_importance(
        bdt_final, mode_variables_filtered, mode, num=20, save_dir=save_dir
    )
    print("\nTop 10 Most Important Features:")
    print(feature_imp_df.sort_values(by="Value", ascending=False).head(10))
    
    # 7. BDT Output Comparison with Pull Plots
    print("\n" + "="*80)
    print("BDT OUTPUT COMPARISON (Train vs Test with Pull Plots)")
    print("="*80)
    decisions = compare_train_test(bdt_final, X_train, y_train, X_test, y_test, mode)
    
    if save_dir:
        plt.savefig(os.path.join(save_dir, "bdt_output_comparison.png"), dpi=150, bbox_inches='tight')
        plt.close()
    else:
        plt.show()
    
    # 8. ROC Curve (keep existing Belle II style)
    print("\nPlotting ROC curve...")
    y_score_test = bdt_final.predict_proba(X_test)[:, 1]
    fpr_test, tpr_test, _ = roc_curve(y_test, y_score_test)
    area_test = auc(fpr_test, tpr_test)
    
    y_score_train = bdt_final.predict_proba(X_train)[:, 1]
    fpr_train, tpr_train, _ = roc_curve(y_train, y_score_train)
    area_train = auc(fpr_train, tpr_train)
    
    # Compute background rejection vs signal efficiency
    bdt_cuts = np.linspace(0, 1, 100)
    sig_train = y_score_train[y_train == 1]
    bkg_train = y_score_train[y_train == 0]
    sig_test = y_score_test[y_test == 1]
    bkg_test = y_score_test[y_test == 0]
    
    sig_eff_train = []
    bkg_rej_train = []
    sig_eff_test = []
    bkg_rej_test = []
    
    for cut in bdt_cuts:
        sig_eff_train.append(np.sum(sig_train > cut) / len(sig_train))
        bkg_rej_train.append(1 - (np.sum(bkg_train > cut) / len(bkg_train)))
        sig_eff_test.append(np.sum(sig_test > cut) / len(sig_test))
        bkg_rej_test.append(1 - (np.sum(bkg_test > cut) / len(bkg_test)))
    
    fig, axs = plt.subplots(1, 1, figsize=(7, 6))
    
    axs.plot(bkg_rej_train, sig_eff_train, color='tab:blue', linewidth=2,
            label=f'Train (AUC = {area_train:.2f})')
    axs.plot(bkg_rej_test, sig_eff_test, color='tab:red', linestyle='--', linewidth=2,
            label=f'Test (AUC = {area_test:.2f})')
    
    axs.fill_between(bkg_rej_test, sig_eff_train, sig_eff_test,
                     where=(np.array(sig_eff_train) > np.array(sig_eff_test)),
                     color='gray', alpha=0.2, label='Overfit Gap')
    
    axs.set_title(mode_title, loc='left')
    axs.set_ylim(0, 1.05)
    axs.set_xlim(0, 1.05)
    axs.set_xlabel('Background rejection')
    axs.set_ylabel('Signal efficiency')
    axs.legend(loc='lower left')
    axs.grid(True)
    
    if save_dir:
        plt.savefig(os.path.join(save_dir, "roc_curve_belle2.png"), dpi=150, bbox_inches='tight')
        plt.close()
    else:
        plt.show()
    
    # 9. Optimize cut using Punzi FoM
    print("\n" + "="*80)
    print("OPTIMIZING BDT CUT USING PUNZI FoM")
    print("="*80)
    
    # Apply BDT to signal and generic
    df_signal["Ds_FakeD0BDT"] = bdt_final.predict_proba(
        df_signal[mode_variables_filtered]
    )[:, 1].astype(np.float32)
    
    df_generic["Ds_FakeD0BDT"] = bdt_final.predict_proba(
        df_generic[mode_variables_filtered]
    )[:, 1].astype(np.float32)
    
    # Get real D0 from signal for Punzi FoM calculation
    sig_real_mask = (abs(df_signal["D0_mcPDG"]) == 421)
    bdt_scores_sig = df_signal.loc[sig_real_mask, "Ds_FakeD0BDT"].values
    bdt_scores_bkg = df_generic["Ds_FakeD0BDT"].values
    
    # Compute Punzi FoM for different confidence levels
    print("\nOptimizing BDT cut for different confidence levels:\n")
    
    best_cuts = {}
    best_foms = {}
    
    for a_val in PUNZI_A_VALUES:
        thresholds, foms, best_thresh, best_fom = compute_punzi_scan(
            bdt_scores_sig, bdt_scores_bkg, N_SIGNAL_EVENTS, n_thresholds=200, a=a_val
        )
        
        cl_name = {1.64: "90% CL", 1.96: "95% CL", 3.0: "3σ"}.get(a_val, f"a={a_val}")
        best_cuts[cl_name] = best_thresh
        best_foms[cl_name] = best_fom
        print(f"  {cl_name:8s}: Optimal cut = {best_thresh:.3f}, Punzi FoM = {best_fom:.6g}")
    
    # Plot Punzi FoM curves for all confidence levels
    if save_dir:
        fig, ax = plt.subplots(figsize=(10, 7))
        
        # Colors: purple for 95% CL (main), black for others
        colors = {"90% CL": "black", "95% CL": "purple", "3σ": "black"}
        linestyles = {"90% CL": "--", "95% CL": "-", "3σ": "-."}
        linewidths = {"90% CL": 2, "95% CL": 3, "3σ": 2}
        
        for a_val in PUNZI_A_VALUES:
            thresholds, foms, best_thresh, best_fom = compute_punzi_scan(
                bdt_scores_sig, bdt_scores_bkg, N_SIGNAL_EVENTS, n_thresholds=200, a=a_val
            )
            
            cl_name = {1.64: "90% CL", 1.96: "95% CL", 3.0: "3σ"}.get(a_val, f"a={a_val}")
            
            ax.plot(
                thresholds, foms,
                color=colors[cl_name],
                linestyle=linestyles[cl_name],
                linewidth=linewidths[cl_name],
                label=f"{cl_name} (cut={best_thresh:.3f})"
            )
            
            # Mark optimal point
            ax.plot(best_thresh, best_fom, "o", color=colors[cl_name], markersize=8)
        
        ax.set_xlabel("BDT Cut", color="black")
        ax.set_ylabel("Punzi FoM", color="black")
        ax.set_title(mode_title, loc="left")
        ax.legend(loc="best")
        ax.grid(True, alpha=0.3)
        
        # Make sure tick labels are black
        ax.tick_params(axis="x", colors="black")
        ax.tick_params(axis="y", colors="black")
        
        plt.tight_layout()
        plt.savefig(os.path.join(save_dir, "punzi_fom_optimization.png"), dpi=150, bbox_inches="tight")
        plt.close()  # Close to avoid display
        
        print(f"\n✓ Saved Punzi FoM plot to: {os.path.join(save_dir, 'punzi_fom_optimization.png')}")
    
    print("="*80 + "\n")

# %% [markdown]
# ## Run Training
# 
# Execute the training pipeline for the selected mode(s).
# 
# Set `TRAIN_MODE` in the configuration section to:
# - `"kmpip"` - Train only on K⁻ π⁺ mode
# - `"km3pi"` - Train only on K⁻ 3π mode
# - `"kmpippi0_eff20_May2020"` - Train only on K⁻ π⁺ π⁰ mode
# - `"all"` - Train on all modes sequentially

# %%
if TRAIN_MODE == "all":
    # Train on all modes sequentially
    for mode in DECAY_CONFIG.keys():
        process_mode(mode)
else:
    # Train on single mode
    if TRAIN_MODE not in DECAY_CONFIG:
        raise ValueError(f"Unknown mode: {TRAIN_MODE}. Choose from {list(DECAY_CONFIG.keys())} or 'all'")
    process_mode(TRAIN_MODE)

print("\n" + "="*80)
print("✓ ALL PROCESSING COMPLETE")
print("="*80)

# %% [markdown]
# ## Control Sample Processing (Optional)
# 
# Load and apply BDT to control samples for systematic studies.
# 
# **Available control samples**:
# - `WCh`: Wrong charge (opposite sign combinations)
# - `ReverseID`: Reversed particle ID (e.g., π identified as K)
# - `ReverseID_WCh`: Combination of both
# 
# **Note**: This section is currently disabled. Set `LOAD_CONTROL_SAMPLES = True` in the configuration to enable.
# 
# Control samples are useful for:
# - Validating BDT performance on data-driven backgrounds
# - Estimating systematic uncertainties
# - Cross-checking signal/background discrimination

# %%
if LOAD_CONTROL_SAMPLES:
    print("\n" + "="*80)
    print("PROCESSING CONTROL SAMPLES")
    print("="*80)
    
    # Get the trained BDT (need to reload or keep from previous cell)
    # For now, we'll process control samples separately per mode
    
    for control_sample in CONTROL_SAMPLES:
        print(f"\nProcessing control sample: {control_sample}")
        
        # Load data
        df_signal_ctrl, df_generic_ctrl = load_mode_data(TRAIN_MODE, control_sample=control_sample)
        
        if df_signal_ctrl.empty or df_generic_ctrl.empty:
            print(f"  ✗ No data for control sample {control_sample}")
            continue
        
        # Get variables for this mode
        mode_variables = VARIABLES[TRAIN_MODE]["all_vars"]
        mode_variables_filtered = [v for v in mode_variables if v in df_signal_ctrl.columns]
        
        # Apply BDT (would need to reload trained model from saved file)
        # For demonstration, assuming we have bdt_final from main training
        # In practice, you'd save the model and reload it here
        
        print(f"  Loaded {len(df_signal_ctrl):,} signal events")
        print(f"  Loaded {len(df_generic_ctrl):,} generic events")
        
        # Save outputs
        if SAVE_OUTPUT:
            output_dir_ctrl = os.path.join(OUTPUT_BASE_DIR, f"FakeD0_{control_sample}", TRAIN_MODE)
            os.makedirs(output_dir_ctrl, exist_ok=True)
            
            # Save signal
            signal_path = os.path.join(output_dir_ctrl, f"Ds2D0enu-Signal_{TRAIN_MODE}_withBDT.root")
            with uproot.recreate(signal_path) as f:
                f["Dstree"] = df_signal_ctrl
            print(f"  Saved: {signal_path}")
            
            # Save generic
            generic_path = os.path.join(output_dir_ctrl, f"Ds2D0e-Generic_{TRAIN_MODE}_withBDT.root")
            with uproot.recreate(generic_path) as f:
                f["Dstree"] = df_generic_ctrl
            print(f"  Saved: {generic_path}")
        
        # Clear memory
        del df_signal_ctrl, df_generic_ctrl
        gc.collect()
        
    print("\n✓ Control sample processing complete")
else:
    print("\nControl sample processing is disabled (LOAD_CONTROL_SAMPLES = False)")

# %% [markdown]
# ## Control Sample / Data-MC Comparison (Optional)
# 
# This section is for validating the BDT on control samples or sideband data.
# 
# **Use cases**:
# 1. **Control Samples**: Load control samples (WCh, ReverseID, etc.) and apply the trained BDT
# 2. **Sideband Data**: Load signal region sideband data for data-MC comparison
# 3. **BDT Output Comparison**: Compare BDT distributions between MC and data/control samples
# 4. **Variable Comparison**: Compare input variables between MC and control samples
# 
# **Available control samples**:
# - `WCh`: Wrong charge (opposite sign combinations)
# - `ReverseID`: Reversed particle ID (e.g., π identified as K)
# - `ReverseID_WCh`: Combination of both
# 
# **Note**: Set `LOAD_CONTROL_SAMPLES = True` in the configuration to enable.

# %%
# This cell provides a template for data-MC comparison using control samples
# Modify as needed for your specific analysis

if LOAD_CONTROL_SAMPLES:
    print("\n" + "="*80)
    print("DATA-MC COMPARISON WITH CONTROL SAMPLES")
    print("="*80)
    
    # TODO: Load the trained BDT model for the mode you want to validate
    # You may need to save/load the model from the training step
    # For now, assuming you have bdt_final and mode_variables_filtered from training
    
    for control_sample in CONTROL_SAMPLES:
        print(f"\n{'-'*80}")
        print(f"Processing control sample: {control_sample}")
        print(f"{'-'*80}")
        
        # Load control sample data
        df_signal_ctrl, df_generic_ctrl = load_mode_data(TRAIN_MODE, control_sample=control_sample)
        
        if df_signal_ctrl.empty or df_generic_ctrl.empty:
            print(f"  ✗ No data for control sample {control_sample}")
            continue
        
        # Filter variables (same as training)
        mode_variables_ctrl = filter_valid_variables(
            pd.concat([df_signal_ctrl, df_generic_ctrl], ignore_index=True),
            VARIABLES[TRAIN_MODE]["all_vars"],
            nan_threshold=0.5
        )
        
        print(f"  Valid variables: {len(mode_variables_ctrl)}")
        
        # Apply BDT to control samples (requires trained model)
        # df_signal_ctrl["Ds_FakeD0BDT"] = bdt_final.predict_proba(
        #     df_signal_ctrl[mode_variables_ctrl]
        # )[:, 1].astype(np.float32)
        
        # df_generic_ctrl["Ds_FakeD0BDT"] = bdt_final.predict_proba(
        #     df_generic_ctrl[mode_variables_ctrl]
        # )[:, 1].astype(np.float32)
        
        # Example: Compare variable distributions (MC vs control sample)
        # Select a few important variables to compare
        # vars_to_compare = mode_variables_ctrl[:5]  # First 5 variables
        
        # for var in vars_to_compare:
        #     plt.figure(figsize=(10, 6))
        #     
        #     # Plot MC
        #     plt.hist(df_signal[var].dropna(), bins=50, alpha=0.5, 
        #              label='MC', density=True, histtype='step', linewidth=2)
        #     
        #     # Plot control sample
        #     plt.hist(df_signal_ctrl[var].dropna(), bins=50, alpha=0.5,
        #              label=f'Control ({control_sample})', density=True, 
        #              histtype='step', linewidth=2)
        #     
        #     plt.xlabel(var)
        #     plt.ylabel('Event Density')
        #     plt.legend()
        #     plt.title(f'{MODE_TITLES[TRAIN_MODE]} - {var}', loc='left')
        #     plt.show()
        
        # Example: Compare BDT output distributions
        # plt.figure(figsize=(10, 6))
        # plt.hist(df_signal["Ds_FakeD0BDT"], bins=50, alpha=0.5,
        #          label='MC', density=True, range=(0, 1), histtype='step', linewidth=2)
        # plt.hist(df_signal_ctrl["Ds_FakeD0BDT"], bins=50, alpha=0.5,
        #          label=f'Control ({control_sample})', density=True, range=(0, 1),
        #          histtype='step', linewidth=2)
        # plt.xlabel('Fake D⁰ BDT Output')
        # plt.ylabel('Event Density')
        # plt.legend()
        # plt.title(f'{MODE_TITLES[TRAIN_MODE]} - BDT Output Comparison', loc='left')
        # plt.show()
        
        print(f"  Loaded {len(df_signal_ctrl):,} signal events")
        print(f"  Loaded {len(df_generic_ctrl):,} generic events")
        
        # Save outputs if requested
        if SAVE_OUTPUT:
            output_dir_ctrl = os.path.join(OUTPUT_BASE_DIR, f"FakeD0_{control_sample}", TRAIN_MODE)
            os.makedirs(output_dir_ctrl, exist_ok=True)
            
            # Save signal
            signal_path = os.path.join(output_dir_ctrl, f"Ds2D0enu-Signal_{TRAIN_MODE}_withBDT.root")
            with uproot.recreate(signal_path) as f:
                f["Dstree"] = df_signal_ctrl
            print(f"  Saved: {signal_path}")
            
            # Save generic
            generic_path = os.path.join(output_dir_ctrl, f"Ds2D0e-Generic_{TRAIN_MODE}_withBDT.root")
            with uproot.recreate(generic_path) as f:
                f["Dstree"] = df_generic_ctrl
            print(f"  Saved: {generic_path}")
        
        # Clear memory
        del df_signal_ctrl, df_generic_ctrl
        gc.collect()
        
    print("\n✓ Control sample processing complete")
else:
    print("\nControl sample processing is disabled (LOAD_CONTROL_SAMPLES = False)")
    print("Set LOAD_CONTROL_SAMPLES = True in the configuration to enable data-MC comparison.")


