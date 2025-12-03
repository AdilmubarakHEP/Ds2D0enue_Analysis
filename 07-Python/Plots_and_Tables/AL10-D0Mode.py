# %% [markdown]
#  # Unified D⁰ Mass Fits + π⁰ Tag Optimization with Punzi FoM
# 
#  **Workflow:**
#  1. Fit D⁰ invariant mass for all modes (kmpip, km3pi) and all π⁰ tags (kmpippi0)
#  2. Extract fitted σ values to define 3σ mass windows
#  3. Load Ds⁺ reconstruction data (DstreeCh2), apply 3σ D⁰ cuts
#  4. Compute Punzi FoM across π⁰ efficiency tags for multiple confidence levels

# %% [markdown]
#  ## Imports and Setup

# %%
# Core
import numpy as np
import pandas as pd
import uproot
import matplotlib.pyplot as plt
from IPython.display import display
import os
import time
import gc

# RooFit / ROOT
import ROOT

# plothist
try:
    import plothist as ph
    from plothist import make_hist, plot_hist, plot_data_model_comparison
except Exception as e:
    raise RuntimeError(
        "plothist is required. Install: pip install --user plothist"
    ) from e

# SciPy utilities
from scipy.interpolate import interp1d
from scipy.integrate import trapezoid

# Inline plots and ROOT batch mode
%matplotlib inline
ROOT.gROOT.SetBatch(True)

print(f"ROOT version: {ROOT.gROOT.GetVersion()}")
print(f"plothist version: {ph.__version__}")

# %% [markdown]
#  ## Configuration

# %%
# === D⁰ Signal MC Configuration ===
SIGNAL_BASE = "/home/belle2/amubarak/C01-Simulated_Events/Signal"

RUN_MODES = ["kmpip", "km3pi"]  # Fit these once
EFF_TAGS = ["eff10_May2020", "eff20_May2020", "eff30_May2020", 
            "eff40_May2020", "eff50_May2020"]  # kmpippi0 tags

mode_cfg = {
    "kmpip": {
        "tree": "D02kmpiptree",
        "file": "output_test_kmpip.root",
        "var": "D0_kmpip_dM",
        "isig": "D0_kmpip_isSignal",
        "bins": 50,
        "range": (-0.03, 0.03),
        "title": r"$D^{0}\!\to K^{-}\pi^{+}$",
        "sig_model": "gauss",
        "bkg_model": "flat",
    },
    "kmpippi0": {
        "tree": "D02kmpippi0tree",
        "file_fmt": "output_test_kmpippi0_{tag}.root",  # {tag} replaced per scan
        "var": "D0_kmpippi0_dM",
        "isig": "D0_kmpippi0_isSignal",
        "bins": 50,
        "range": (-0.15, 0.10),
        "title": r"$D^{0}\!\to K^{-}\pi^{+}\pi^{0}$",
        "sig_model": "bifurgauss",
        "bkg_model": "cheb2",
    },
    "km3pi": {
        "tree": "D02km3pitree",
        "file": "output_test_km3pi.root",
        "var": "D0_km3pi_dM",
        "isig": "D0_km3pi_isSignal",
        "bins": 50,
        "range": (-0.04, 0.04),
        "title": r"$D^{0}\!\to K^{-}\pi^{+}\pi^{-}\pi^{+}$",
        "sig_model": "gauss",
        "bkg_model": "flat",
    },
}

# Fit seeds per mode
fit_seeds = {
    "kmpip": {
        "mu0": 0.0, "mu_min": -0.01, "mu_max": 0.01,
        "sigma0": 0.005, "sigma_min": 1e-4, "sigma_max": 0.02,
        "Nsig0": 10000, "Nbkg0": 3000,
    },
    "kmpippi0": {
        "mu0": -0.010, "mu_min": -0.04, "mu_max": 0.04,
        "sigmaL0": 0.018, "sigmaL_min": 0.001, "sigmaL_max": 0.08,
        "sigmaR0": 0.023, "sigmaR_min": 0.001, "sigmaR_max": 0.08,
        "c1_0": -0.2, "c2_0": 0.0,
        "Nsig0": 15000, "Nbkg0": 35000,
    },
    "km3pi": {
        "mu0": 0.0, "mu_min": -0.01, "mu_max": 0.01,
        "sigma0": 0.004, "sigma_min": 1e-4, "sigma_max": 0.02,
        "Nsig0": 9000, "Nbkg0": 4000,
    },
}

NGEN_SIGNAL = 100000  # Generated signal events per sample

# === Ds⁺ Reconstruction Configuration (for Punzi FoM) ===
DS_TREE = "DstreeCh2"
VAR_D0_DM = "D0_dM"  # D⁰ mass variable in Ds tree
VAR_LABEL = "Ds_ifNANgiveX_isSignal_5"  # Truth label for Ds⁺

# Generic MC paths (for background)
BKG_BASE = "/group/belle2/users2022/amubarak/Sample_KEKCC"
BKG_EVENT_TYPES = ["ccbar", "ddbar", "uubar", "ssbar", "mixed", "charged"]

# Punzi FoM settings
PUNZI_A_VALUES = [1.64, 1.96, 3.0]  # 90% CL, 95% CL, 3σ discovery
TARGET_LUMI_FB = 200  # Target luminosity for reference (not used in Punzi calculation)

# Plotting style
LINEWIDTH = 2.0
COLORS_STACK = ["purple", "#fd7f6f", "#4C6EB1"]  # NaN, bkg, signal

# %% [markdown]
#  ## RooFit Utilities

# %%
def numpy_to_roodataset(x, var, name="data"):
    ds = ROOT.RooDataSet(name, name, ROOT.RooArgSet(var))
    aset = ROOT.RooArgSet(var)
    for v in x:
        var.setVal(float(v))
        ds.add(aset)
    return ds

def tmatrix_to_df(m, names):
    n = m.GetNrows()
    arr = np.zeros((n, n), float)
    for i in range(n):
        for j in range(n):
            try:
                arr[i, j] = m[i][j]
            except Exception:
                arr[i, j] = m(i, j)
    return pd.DataFrame(arr, index=names, columns=names)

def save_pdf(var, pdf, n_points=10000):
    xmin, xmax = var.getMin(), var.getMax()
    xs = np.linspace(xmin, xmax, int(n_points))
    ys = np.zeros_like(xs, float)
    for i, xv in enumerate(xs):
        var.setVal(float(xv))
        ys[i] = float(pdf.getVal(var))
    return interp1d(xs, ys, bounds_error=False, fill_value=0.0)

def renormalize_to_yield(pdf_func, x_range, n_bins, yield_events, n_int=5000):
    xmin, xmax = x_range
    bw = (xmax - xmin) / n_bins
    xs = np.linspace(xmin, xmax, int(n_int))
    ys = pdf_func(xs)
    area = trapezoid(ys, xs)
    if area <= 0:
        return lambda x: 0.0 * np.asarray(x)
    scale = (yield_events * bw) / area
    return lambda x: pdf_func(x) * scale

def _build_signal_pdf(dm, model_name, seeds):
    """Build signal PDF: 'gauss', 'bifurgauss', or 'doublegauss'."""
    mu = ROOT.RooRealVar(
        "mu", "mu",
        float(seeds.get("mu0", 0.0)),
        float(seeds.get("mu_min", dm.getMin())),
        float(seeds.get("mu_max", dm.getMax())),
    )

    if model_name == "gauss":
        sigma = ROOT.RooRealVar(
            "sigma", "sigma",
            float(seeds.get("sigma0", 0.005)),
            float(seeds.get("sigma_min", 1e-5)),
            float(seeds.get("sigma_max", 0.1)),
        )
        sig_pdf = ROOT.RooGaussian("sig_gauss", "sig_gauss", dm, mu, sigma)
        return sig_pdf, 2, {"mu": mu, "sigma": sigma}

    elif model_name == "bifurgauss":
        sigmaL = ROOT.RooRealVar(
            "sigmaL", "sigmaL",
            float(seeds.get("sigmaL0", 0.018)),
            float(seeds.get("sigmaL_min", 1e-5)),
            float(seeds.get("sigmaL_max", 0.2)),
        )
        sigmaR = ROOT.RooRealVar(
            "sigmaR", "sigmaR",
            float(seeds.get("sigmaR0", 0.023)),
            float(seeds.get("sigmaR_min", 1e-5)),
            float(seeds.get("sigmaR_max", 0.2)),
        )
        sig_pdf = ROOT.RooBifurGauss("sig_bifur", "sig_bifur", dm, mu, sigmaL, sigmaR)
        return sig_pdf, 3, {"mu": mu, "sigmaL": sigmaL, "sigmaR": sigmaR}
    
    raise ValueError(f"Unknown signal model '{model_name}'")

def _build_background_pdf(dm, model_name, seeds):
    """Build background PDF: 'flat', 'cheb1', 'cheb2', 'expo'."""
    if model_name == "flat":
        bkg = ROOT.RooUniform("bkg_flat", "bkg_flat", ROOT.RooArgSet(dm))
        return bkg, 0, {}

    elif model_name == "cheb1":
        c1 = ROOT.RooRealVar("c1", "c1", float(seeds.get("c1_0", 0.0)), -2.0, 2.0)
        bkg = ROOT.RooChebychev("bkg_cheb1", "bkg_cheb1", dm, ROOT.RooArgList(c1))
        return bkg, 1, {"c1": c1}

    elif model_name == "cheb2":
        c1 = ROOT.RooRealVar("c1", "c1", float(seeds.get("c1_0", 0.0)), -2.0, 2.0)
        c2 = ROOT.RooRealVar("c2", "c2", float(seeds.get("c2_0", 0.0)), -2.0, 2.0)
        bkg = ROOT.RooChebychev("bkg_cheb2", "bkg_cheb2", dm, ROOT.RooArgList(c1, c2))
        return bkg, 2, {"c1": c1, "c2": c2}
    
    raise ValueError(f"Unknown background model '{model_name}'")

def fit_signal_plus_bkg(x, x_range, bins, signal_model, background_model, seeds):
    """Unbinned extended maximum likelihood fit."""
    xmin, xmax = x_range
    dm = ROOT.RooRealVar("dm", "m(D^{0})-m_{PDG}(D^{0}) [GeV/c^{2}]", xmin, xmax)
    dm.setRange("fitRange", xmin, xmax)

    ds = numpy_to_roodataset(x, dm)
    n_tot = len(x)

    sig_pdf, npar_sig, sig_pars = _build_signal_pdf(dm, signal_model, seeds)
    bkg_pdf, npar_bkg, bkg_pars = _build_background_pdf(dm, background_model, seeds)

    Nsig = ROOT.RooRealVar("Nsig", "Nsig", float(seeds.get("Nsig0", 0.6*n_tot)), 0.0, 2.0*n_tot)
    Nbkg = ROOT.RooRealVar("Nbkg", "Nbkg", float(seeds.get("Nbkg0", 0.4*n_tot)), 0.0, 2.0*n_tot)

    model = ROOT.RooAddPdf("model", "sig+bkg",
                           ROOT.RooArgList(sig_pdf, bkg_pdf),
                           ROOT.RooArgList(Nsig, Nbkg))

    fitRes = model.fitTo(ds, ROOT.RooFit.Save(True), ROOT.RooFit.Extended(True), 
                         ROOT.RooFit.PrintLevel(-1))

    out = {
        "Nsig": (Nsig.getVal(), Nsig.getError()),
        "Nbkg": (Nbkg.getVal(), Nbkg.getError()),
        "Ndata": int(n_tot),
        "status": int(fitRes.status()),
        "edm": float(fitRes.edm()),
        "covQual": int(fitRes.covQual()),
        "signal_model": signal_model,
        "npar_signal": npar_sig,
        "npar_bkg": npar_bkg,
    }

    # Extract shape parameters
    if signal_model == "gauss":
        mu_v, sig_v = sig_pars["mu"], sig_pars["sigma"]
        out["mu"] = (mu_v.getVal(), mu_v.getError())
        out["sigma"] = (sig_v.getVal(), sig_v.getError())
        out["sigmaL"] = out["sigmaR"] = out["sigma"]  # For 3σ window calc
    elif signal_model == "bifurgauss":
        mu_v, sL_v, sR_v = sig_pars["mu"], sig_pars["sigmaL"], sig_pars["sigmaR"]
        out["mu"] = (mu_v.getVal(), mu_v.getError())
        out["sigmaL"] = (sL_v.getVal(), sL_v.getError())
        out["sigmaR"] = (sR_v.getVal(), sR_v.getError())

    # Covariance & correlation
    pars = fitRes.floatParsFinal()
    names = [pars[i].GetName() for i in range(pars.getSize())]
    out["covariance_df"] = tmatrix_to_df(fitRes.covarianceMatrix(), names)
    out["correlation_df"] = tmatrix_to_df(fitRes.correlationMatrix(), names)

    # Callables for plotting
    f_sig = save_pdf(dm, sig_pdf)
    f_bkg = save_pdf(dm, bkg_pdf)
    out["sig_callable"] = renormalize_to_yield(f_sig, x_range, bins, out["Nsig"][0])
    out["bkg_callable"] = renormalize_to_yield(f_bkg, x_range, bins, out["Nbkg"][0])
    out["pdf_callable"] = lambda xs: out["sig_callable"](xs) + out["bkg_callable"](xs)

    return out

def entries_per_bin_mev(x_range, bins):
    return (x_range[1] - x_range[0]) / bins * 1000.0

def chisq_ndf_from_hist_and_func(h_data, f_model, npar):
    centers = h_data.axes[0].centers
    data = h_data.values()
    var = h_data.variances()
    if var is None:
        var = np.maximum(data, 1.0)
    var = np.asarray(var, dtype=float)
    model = np.asarray(f_model(centers), dtype=float)
    mask = var > 0
    chi2 = np.sum((data[mask] - model[mask])**2 / var[mask])
    ndf = int(mask.sum() - npar)
    return chi2, ndf, chi2 / max(ndf, 1)

# %% [markdown]
#  ## Phase 1: Fit D⁰ Mass Peaks and Extract σ

# %% [markdown]
#  ### Fit kmpip and km3pi (once each)

# %%
fit_results = {}

for mode in RUN_MODES:
    cfg = mode_cfg[mode]
    path = f"{SIGNAL_BASE}/{cfg['file']}:{cfg['tree']}"
    
    print(f"\n{'='*80}\n[{mode}] Loading {path}")
    df = uproot.concatenate(path, library="pd")
    print(f"Loaded {len(df):,} rows")
    
    x_all = df[cfg['var']].dropna().to_numpy(dtype=float)
    
    # Fit
    res = fit_signal_plus_bkg(
        x_all, 
        x_range=cfg['range'], 
        bins=cfg['bins'],
        signal_model=cfg['sig_model'],
        background_model=cfg['bkg_model'],
        seeds=fit_seeds[mode]
    )
    
    fit_results[mode] = res
    
    # Print summary
    Nsig_val, Nsig_err = res["Nsig"]
    eff = Nsig_val / NGEN_SIGNAL
    eff_err = Nsig_err / NGEN_SIGNAL
    
    print(f"Mode: {mode} | Ndata: {res['Ndata']} | status: {res['status']} | EDM: {res['edm']:.2e}")
    print(f"Nsig: {Nsig_val:.0f} ± {Nsig_err:.0f}")
    print(f"Reco eff: {eff:.5f} ± {eff_err:.5f} ({eff*100:.2f} ± {eff_err*100:.2f}%)")
    print(f"μ: {res['mu'][0]*1e3:.2f} ± {res['mu'][1]*1e3:.2f} MeV/c²")
    print(f"σ: {res['sigma'][0]*1e3:.2f} ± {res['sigma'][1]*1e3:.2f} MeV/c²")
    
    # Calculate 3σ window (symmetric for Gaussian)
    mu = res['mu'][0]
    sigma = res['sigma'][0]
    lower_3sig = mu - 3*sigma
    upper_3sig = mu + 3*sigma
    print(f"3σ window: [{lower_3sig:.5f}, {upper_3sig:.5f}] GeV/c²")
    
    # Plot
    x_range = cfg['range']
    xmin, xmax = x_range
    x_plot = x_all[(x_all >= xmin) & (x_all <= xmax)]
    h_data = make_hist(x_plot, bins=cfg['bins'], range=x_range, weights=1)
    
    npar = 2 + res["npar_signal"] + res["npar_bkg"]
    chi2, ndf, chi2_ndf = chisq_ndf_from_hist_and_func(h_data, res["pdf_callable"], npar)
    
    per_bin_mev = entries_per_bin_mev(x_range, cfg['bins'])
    xlabel = r"$m(D^{0}) - m_{\mathrm{PDG}}(D^{0})\;[\mathrm{GeV}/c^{2}]$"
    ylabel = rf"Entries/({per_bin_mev:.2f} MeV/$c^2$)"
    
    fig, ax_main, ax_comp = plot_data_model_comparison(
        data_hist=h_data,
        stacked_components=[res["bkg_callable"]],
        stacked_labels=["Background"],
        unstacked_components=[res["sig_callable"]],
        unstacked_labels=["Signal"],
        xlabel=xlabel,
        ylabel=ylabel,
        model_sum_kwargs={"show": True, "label": "Model"},
        comparison="pull",
    )
    
    # === ADD 3σ BOUNDARY VISUALIZATION ===
    # Shade excluded regions (gray)
    ax_main.axvspan(xmin, lower_3sig, color='gray', alpha=0.2, zorder=0)
    ax_main.axvspan(upper_3sig, xmax, color='gray', alpha=0.2, zorder=0)
    # Vertical lines at 3σ boundaries
    ax_main.axvline(lower_3sig, ls='--', color='gray', lw=1.5, zorder=3)
    ax_main.axvline(upper_3sig, ls='--', color='gray', lw=1.5, zorder=3)
    
    # Fix legend order: move Data to top
    handles, labels = ax_main.get_legend_handles_labels()
    if "Data" in labels:
        idx_data = labels.index("Data")
        handles = [handles[idx_data]] + [h for i, h in enumerate(handles) if i != idx_data]
        labels = [labels[idx_data]] + [l for i, l in enumerate(labels) if i != idx_data]
    
    # Annotation (fit parameters at top-right)
    s, se = res["sigma"]
    txt = (
        rf"$\sigma = {s*1e3:.2f} \pm {se*1e3:.2f}\ \mathrm{{MeV}}/c^2$" + "\n" +
        rf"$\chi^2/\mathrm{{ndf}} = {chi2_ndf:.3f}$"
    )
    ax_main.text(0.98, 0.98, txt, transform=ax_main.transAxes, va="top", ha="right", fontsize=10)
    
    ax_main.set_title(cfg['title'], loc="left")
    ax_main.set_title("Signal MC", loc="right")
    ax_main.legend(handles, labels, loc="upper left", frameon=False)
    
    plt.tight_layout()
    plt.show()
    
    # Matrices
    print("\nCovariance Matrix:")
    display(res["covariance_df"].round(6))
    print("\nCorrelation Matrix:")
    display(res["correlation_df"].round(3))

# %% [markdown]
#  ### Fit all kmpippi0 π⁰ efficiency tags

# %%
cfg_kmpippi0 = mode_cfg["kmpippi0"]

for tag in EFF_TAGS:
    path = f"{SIGNAL_BASE}/{cfg_kmpippi0['file_fmt'].format(tag=tag)}:{cfg_kmpippi0['tree']}"
    
    print(f"\n{'='*80}\n[kmpippi0-{tag}] Loading {path}")
    try:
        df = uproot.concatenate(path, library="pd")
        print(f"Loaded {len(df):,} rows")
    except Exception as e:
        print(f"[ERROR] Failed to load: {e}")
        continue
    
    x_all = df[cfg_kmpippi0['var']].dropna().to_numpy(dtype=float)
    
    # Fit
    res = fit_signal_plus_bkg(
        x_all,
        x_range=cfg_kmpippi0['range'],
        bins=cfg_kmpippi0['bins'],
        signal_model=cfg_kmpippi0['sig_model'],
        background_model=cfg_kmpippi0['bkg_model'],
        seeds=fit_seeds["kmpippi0"]
    )
    
    fit_results[f"kmpippi0_{tag}"] = res
    
    # Print summary
    Nsig_val, Nsig_err = res["Nsig"]
    eff = Nsig_val / NGEN_SIGNAL
    eff_err = Nsig_err / NGEN_SIGNAL
    
    print(f"Tag: {tag} | Ndata: {res['Ndata']} | status: {res['status']} | EDM: {res['edm']:.2e}")
    print(f"Nsig: {Nsig_val:.0f} ± {Nsig_err:.0f}")
    print(f"Reco eff: {eff:.5f} ± {eff_err:.5f} ({eff*100:.2f} ± {eff_err*100:.2f}%)")
    
    mu, muErr = res['mu']
    sL, sLErr = res['sigmaL']
    sR, sRErr = res['sigmaR']
    
    print(f"μ: {mu*1e3:.2f} ± {muErr*1e3:.2f} MeV/c²")
    print(f"σ_L: {sL*1e3:.2f} ± {sLErr*1e3:.2f} MeV/c²")
    print(f"σ_R: {sR*1e3:.2f} ± {sRErr*1e3:.2f} MeV/c²")
    
    # Calculate 3σ window (asymmetric for BifurGauss)
    lower_3sig = mu - 3*sL
    upper_3sig = mu + 3*sR
    print(f"3σ window: [{lower_3sig:.5f}, {upper_3sig:.5f}] GeV/c²")
    
    # Plot
    x_range = cfg_kmpippi0['range']
    xmin, xmax = x_range
    x_plot = x_all[(x_all >= xmin) & (x_all <= xmax)]
    h_data = make_hist(x_plot, bins=cfg_kmpippi0['bins'], range=x_range, weights=1)
    
    npar = 2 + res["npar_signal"] + res["npar_bkg"]
    chi2, ndf, chi2_ndf = chisq_ndf_from_hist_and_func(h_data, res["pdf_callable"], npar)
    
    per_bin_mev = entries_per_bin_mev(x_range, cfg_kmpippi0['bins'])
    xlabel = r"$m(D^{0}) - m_{\mathrm{PDG}}(D^{0})\;[\mathrm{GeV}/c^{2}]$"
    ylabel = rf"Entries/({per_bin_mev:.2f} MeV/$c^2$)"
    
    fig, ax_main, ax_comp = plot_data_model_comparison(
        data_hist=h_data,
        stacked_components=[res["bkg_callable"]],
        stacked_labels=["Background"],
        unstacked_components=[res["sig_callable"]],
        unstacked_labels=["Signal (BifurGauss)"],
        xlabel=xlabel,
        ylabel=ylabel,
        model_sum_kwargs={"show": True, "label": "Model"},
        comparison="pull",
    )
    
    # === ADD 3σ BOUNDARY VISUALIZATION (ASYMMETRIC) ===
    # Shade excluded regions (gray)
    ax_main.axvspan(xmin, lower_3sig, color='gray', alpha=0.2, zorder=0)
    ax_main.axvspan(upper_3sig, xmax, color='gray', alpha=0.2, zorder=0)
    # Vertical lines at 3σ boundaries
    ax_main.axvline(lower_3sig, ls='--', color='gray', lw=1.5, zorder=3)
    ax_main.axvline(upper_3sig, ls='--', color='gray', lw=1.5, zorder=3)
    
    # Fix legend order
    handles, labels = ax_main.get_legend_handles_labels()
    if "Data" in labels:
        idx_data = labels.index("Data")
        handles = [handles[idx_data]] + [h for i, h in enumerate(handles) if i != idx_data]
        labels = [labels[idx_data]] + [l for i, l in enumerate(labels) if i != idx_data]
    
    # Annotation (fit parameters at top-right)
    txt = (
        rf"$\sigma_L = {sL*1e3:.2f}\ \mathrm{{MeV}}/c^2$" + "\n" +
        rf"$\sigma_R = {sR*1e3:.2f}\ \mathrm{{MeV}}/c^2$" + "\n" +
        rf"$\chi^2/\mathrm{{ndf}} = {chi2_ndf:.3f}$"
    )
    ax_main.text(0.98, 0.98, txt, transform=ax_main.transAxes, va="top", ha="right", fontsize=10)
    
    ax_main.set_title(cfg_kmpippi0['title'], loc="left")
    ax_main.set_title(tag, loc="right")
    ax_main.legend(handles, labels, loc="upper left", frameon=False)
    
    plt.tight_layout()
    plt.show()
    
    # Matrices
    print("\nCovariance Matrix:")
    display(res["covariance_df"].round(6))
    print("\nCorrelation Matrix:")
    display(res["correlation_df"].round(3))

# %% [markdown]
#  ### Dalitz Plot for kmpippi0 (using eff50 data)

# %%
tag_dalitz = "eff50_May2020"
path_dalitz = f"{SIGNAL_BASE}/{cfg_kmpippi0['file_fmt'].format(tag=tag_dalitz)}:{cfg_kmpippi0['tree']}"

print(f"\n{'='*80}\n[Dalitz Plot] Loading {path_dalitz}")
df_dalitz = uproot.concatenate(path_dalitz, library="pd")

needed_cols = ['daughterInvM__bo0__cm__sp1__bc',  # m(K- pi+)
               'daughterInvM__bo1__cm__sp2__bc',  # m(pi+ pi0)
               cfg_kmpippi0['isig']]

missing = [c for c in needed_cols if c not in df_dalitz.columns]
if missing:
    print(f"[WARNING] Missing columns for Dalitz plot: {missing}")
else:
    df_sig_dalitz = df_dalitz.loc[df_dalitz[cfg_kmpippi0['isig']] == 1, needed_cols].dropna()
    
    m_kpi = df_sig_dalitz['daughterInvM__bo0__cm__sp1__bc'].to_numpy(dtype=float)
    m_pipi0 = df_sig_dalitz['daughterInvM__bo1__cm__sp2__bc'].to_numpy(dtype=float)
    
    x_dalitz = m_kpi**2
    y_dalitz = m_pipi0**2
    
    fig, ax = plt.subplots(figsize=(6, 5))
    hb = ax.hexbin(x_dalitz, y_dalitz, gridsize=60, bins="log", cmap="viridis")
    cb = fig.colorbar(hb, ax=ax)
    cb.set_label("Counts")
    
    ax.set_xlabel(r"$m^{2}(K^{-}\pi^{+})\;[\mathrm{GeV}^{2}/c^{4}]$")
    ax.set_ylabel(r"$m^{2}(\pi^{+}\pi^{0})\;[\mathrm{GeV}^{2}/c^{4}]$")
    ax.set_title(r"$D^{0}\!\to K^{-}\pi^{+}\pi^{0}$", loc="left")
    ax.set_title(f"Signal only — {tag_dalitz}", loc="right")
    
    plt.tight_layout()
    plt.show()

# %% [markdown]
#  ## Phase 2: Punzi FoM Calculation with 3σ D⁰ Mass Cut (CORRECTED)

# %% [markdown]
#  ### Helper Functions for Punzi FoM (OPTIMIZED)

# %%
def punzi_fom(epsS, B, a):
    """Punzi figure of merit: ε_S / (a/2 + sqrt(B))."""
    denom = (a / 2.0) + np.sqrt(max(B, 0.0))
    return 0.0 if denom <= 0.0 else epsS / denom

def load_ds_tree(tag, is_signal=True):
    """Load DstreeCh2 (Ds⁺ reconstruction) for given π⁰ tag - OPTIMIZED."""
    # ONLY load columns we actually need
    filter_cols = [VAR_D0_DM, VAR_LABEL]
    
    if is_signal:
        path = f"{SIGNAL_BASE}/output_test_kmpippi0_{tag}.root:{DS_TREE}"
        try:
            print(f"  Loading signal (columns: {filter_cols})... ", end="", flush=True)
            df = uproot.concatenate(path, filter_name=filter_cols, library="pd")
            print(f"{len(df):,} events")
            return df
        except Exception as e:
            print(f"\n[ERROR] Signal load failed for {tag}: {e}")
            return pd.DataFrame()
    else:
        # Load generic MC (all event types, merged)
        specs = []
        for ev_type in BKG_EVENT_TYPES:
            mono_path = f"{BKG_BASE}/Ds2D0e-Generic_M2_{tag}_Ds_101525_0_{ev_type}.root"
            if os.path.isfile(mono_path):
                specs.append(f"{mono_path}:{DS_TREE}")
            else:
                shard_dir = f"{BKG_BASE}/Ds2D0e-Generic_M2_{tag}_Ds_101525_0_{ev_type}_shards"
                if os.path.isdir(shard_dir):
                    import glob
                    parts = sorted(glob.glob(os.path.join(shard_dir, "part_*.root")))
                    specs.extend([f"{p}:{DS_TREE}" for p in parts])
        
        if not specs:
            print(f"[WARNING] No background files found for {tag}")
            return pd.DataFrame()
        
        try:
            print(f"  Loading background from {len(specs)} file(s) (columns: {filter_cols})... ", end="", flush=True)
            df = uproot.concatenate(specs, filter_name=filter_cols, library="pd")
            print(f"{len(df):,} events")
            return df
        except Exception as e:
            print(f"\n[ERROR] Background load failed for {tag}: {e}")
            return pd.DataFrame()

# %% [markdown]
#  ### Run Punzi FoM Scan (CORRECTED - Fixed Efficiency Denominator)

# %%
punzi_results = []
start_time = time.time()

for idx, tag in enumerate(EFF_TAGS, 1):
    print(f"\n{'='*80}")
    print(f"[Punzi FoM {idx}/{len(EFF_TAGS)}] Tag: {tag}")
    print(f"{'='*80}")
    
    # Get fitted D⁰ mass parameters (3σ window)
    fit_key = f"kmpippi0_{tag}"
    if fit_key not in fit_results:
        print(f"[SKIP] No fit results for {tag}")
        continue
    
    res_fit = fit_results[fit_key]
    mu = res_fit['mu'][0]
    sL = res_fit['sigmaL'][0]
    sR = res_fit['sigmaR'][0]
    d0_window = (mu - 3*sL, mu + 3*sR)
    
    print(f"3σ D⁰ mass window: [{d0_window[0]:.5f}, {d0_window[1]:.5f}] GeV/c²")
    
    # === Load FULL datasets (before cuts) - ONLY columns we need ===
    tag_start = time.time()
    
    sig_df_full = load_ds_tree(tag, is_signal=True)
    bkg_df_full = load_ds_tree(tag, is_signal=False)
    
    load_time = time.time() - tag_start
    print(f"  → Load time: {load_time:.1f}s")
    
    if sig_df_full.empty or bkg_df_full.empty:
        print(f"[SKIP] Empty signal or background for {tag}")
        continue
    
    print(f"\n[SANITY CHECK] Raw counts BEFORE 3σ cut:")
    print(f"  Signal MC reconstructed in Ds tree: {len(sig_df_full):,} events")
    print(f"  Background MC: {len(bkg_df_full):,} total events")
    
    # Check truth label exists
    if VAR_LABEL not in sig_df_full.columns:
        print(f"[SKIP] Missing truth label column: {VAR_LABEL}")
        continue
    
    # === SIGNAL EFFICIENCY (CORRECTED) ===
    # CRITICAL FIX: Use CONSTANT denominator (generated events)
    N_sig_total_generated = NGEN_SIGNAL  # 50,000 - CONSTANT for all tags!
    
    # Apply 3σ D⁰ mass cut
    sig_df_cut = sig_df_full[(sig_df_full[VAR_D0_DM] >= d0_window[0]) & 
                             (sig_df_full[VAR_D0_DM] <= d0_window[1])]
    
    # Count TRUE signal events passing cut
    sig_pure_pass = sig_df_cut[sig_df_cut[VAR_LABEL] == 1]
    N_sig_pass = len(sig_pure_pass)
    
    # TRUE total reconstruction efficiency (includes π⁰ reco losses + 3σ cut)
    ε_S = N_sig_pass / N_sig_total_generated  # Denominator NOW constant!
    
    # === BACKGROUND COUNT (no scaling) ===
    # All generic events after 3σ D⁰ cut (already at 200 fb⁻¹)
    bkg_df_cut = bkg_df_full[(bkg_df_full[VAR_D0_DM] >= d0_window[0]) & 
                             (bkg_df_full[VAR_D0_DM] <= d0_window[1])]
    B = len(bkg_df_cut)
    
    # Diagnostic breakdown
    N_sig_reconstructed = len(sig_df_full[sig_df_full[VAR_LABEL] == 1])
    ε_pi0_reco = N_sig_reconstructed / N_sig_total_generated
    ε_3sigma_acceptance = N_sig_pass / N_sig_reconstructed if N_sig_reconstructed > 0 else 0.0
    
    print(f"\n[EFFICIENCY BREAKDOWN]:")
    print(f"  Generated signal events: {N_sig_total_generated:,} (constant)")
    print(f"  Reconstructed in Ds tree: {N_sig_reconstructed:,}")
    print(f"  Passing 3σ D⁰ cut: {N_sig_pass:,}")
    print(f"  π⁰ reconstruction efficiency: {ε_pi0_reco:.4f} = {ε_pi0_reco*100:.2f}%")
    print(f"  3σ cut acceptance (of reconstructed): {ε_3sigma_acceptance:.4f} = {ε_3sigma_acceptance*100:.2f}%")
    print(f"  TOTAL efficiency ε_S: {ε_S:.6f} = {ε_S*100:.2f}%")
    print(f"\n[BACKGROUND]:")
    print(f"  Background events passing 3σ cut: {B:,}")
    
    # === PUNZI FoM ===
    for a_val in PUNZI_A_VALUES:
        P = punzi_fom(ε_S, B, a_val)
        punzi_results.append({
            'tag': tag,
            'a': a_val,
            'Punzi': P,
            'epsilonS': ε_S,
            'B': B,
            'N_sig_generated': N_sig_total_generated,
            'N_sig_reconstructed': N_sig_reconstructed,
            'N_sig_pass': N_sig_pass,
            'eps_pi0_reco': ε_pi0_reco,
            'eps_3sigma': ε_3sigma_acceptance,
            'd0_window': d0_window
        })
        print(f"  a={a_val:.2f}: Punzi={P:.6g}")
    
    # Clear memory
    del sig_df_full, bkg_df_full, sig_df_cut, bkg_df_cut
    gc.collect()

total_time = time.time() - start_time
print(f"\n{'='*80}")
print(f"Total Punzi scan time: {total_time:.1f}s ({total_time/60:.1f} min)")
print(f"{'='*80}")

# Convert to DataFrame for plotting
df_punzi = pd.DataFrame(punzi_results)

# %% [markdown]
#  ### Plot Punzi FoM Comparison

# %%
if not df_punzi.empty:
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    
    for i, a_val in enumerate(PUNZI_A_VALUES):
        df_a = df_punzi[df_punzi['a'] == a_val].sort_values('tag')
        
        ax = axes[i]
        bars = ax.bar(df_a['tag'], df_a['Punzi'], color='#4C6EB1', alpha=0.85)
        
        # Annotate bars with total efficiency
        for bar, row in zip(bars, df_a.itertuples()):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2, height,
                    f"ε={row.epsilonS:.4f}",
                    ha='center', va='bottom', fontsize=9)
        
        ax.set_ylabel(rf"Punzi FoM ($a={a_val:.2f}$)")
        ax.set_xlabel("π⁰ Efficiency Tag")
        ax.set_title(f"a = {a_val:.2f}", loc='center')
        ax.grid(axis='y', alpha=0.3)
        ax.set_ylim(0, 1.1*df_a['Punzi'].max() if len(df_a) > 0 else 1.0)
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=25, ha='right')
    
    plt.suptitle(r"Punzi FoM: $\epsilon_S / (a/2 + \sqrt{B})$ (Denominator = 100k Generated)", 
                 fontsize=14, y=1.02)
    plt.tight_layout()
    plt.show()
    
    # Print summary table
    print("\n" + "="*100)
    print("Punzi FoM Summary (CORRECTED: Denominator = 50,000 Generated Events)")
    print("="*100)
    for a_val in PUNZI_A_VALUES:
        print(f"\na = {a_val:.2f}:")
        df_a = df_punzi[df_punzi['a'] == a_val][['tag', 'Punzi', 'epsilonS', 'eps_pi0_reco', 'eps_3sigma', 'B', 'N_sig_pass']]
        df_a_display = df_a.copy()
        df_a_display.columns = ['tag', 'Punzi', 'ε_total', 'ε_π⁰_reco', 'ε_3σ_cut', 'B', 'N_pass']
        display(df_a_display.set_index('tag'))
else:
    print("[WARNING] No Punzi results to plot")

# %% [markdown]
#  ## 3σ Mass Window Summary for Steering Scripts

# %%
print("\n" + "="*100)
print("3σ D⁰ MASS WINDOWS FOR STEERING SCRIPTS")
print("="*100)
print("\nBased on fits to signal MC (50k events per sample)")
print("Apply these cuts to m(D⁰) - m_PDG(D⁰) variable\n")

# Collect all windows with simple column names
windows_data = []

# kmpip and km3pi (symmetric Gaussian)
for mode in RUN_MODES:
    if mode in fit_results:
        res = fit_results[mode]
        mu = res['mu'][0]
        sigma = res['sigma'][0]
        lower = mu - 3*sigma
        upper = mu + 3*sigma
        windows_data.append({
            'Mode': mode,
            'Tag': 'N/A',
            'μ': mu,
            'σ_L': sigma,
            'σ_R': sigma,
            'Lower': lower,
            'Upper': upper,
        })

# kmpippi0 (bifurcated Gaussian for each tag)
for tag in EFF_TAGS:
    fit_key = f"kmpippi0_{tag}"
    if fit_key in fit_results:
        res = fit_results[fit_key]
        mu = res['mu'][0]
        sL = res['sigmaL'][0]
        sR = res['sigmaR'][0]
        lower = mu - 3*sL
        upper = mu + 3*sR
        windows_data.append({
            'Mode': 'kmpippi0',
            'Tag': tag,
            'μ': mu,
            'σ_L': sL,
            'σ_R': sR,
            'Lower': lower,
            'Upper': upper,
        })

df_windows = pd.DataFrame(windows_data)

# Display formatted table
df_windows_display = df_windows.copy()
df_windows_display['μ (GeV/c²)'] = df_windows_display['μ'].apply(lambda x: f"{x:.6f}")
df_windows_display['σ_L (GeV/c²)'] = df_windows_display['σ_L'].apply(lambda x: f"{x:.6f}")
df_windows_display['σ_R (GeV/c²)'] = df_windows_display['σ_R'].apply(lambda x: f"{x:.6f}")
df_windows_display['Lower (GeV/c²)'] = df_windows_display['Lower'].apply(lambda x: f"{x:.6f}")
df_windows_display['Upper (GeV/c²)'] = df_windows_display['Upper'].apply(lambda x: f"{x:.6f}")

display(df_windows_display[['Mode', 'Tag', 'μ (GeV/c²)', 'σ_L (GeV/c²)', 'σ_R (GeV/c²)', 'Lower (GeV/c²)', 'Upper (GeV/c²)']])

print("\n" + "="*100)
print("STEERING SCRIPT CUTS (Python format - copy-paste ready)")
print("="*100)
print()

# Print in steering script format
for _, row in df_windows.iterrows():
    mode = row['Mode']
    tag = row['Tag']
    lower = row['Lower']
    upper = row['Upper']
    
    if tag == 'N/A':
        var_name = f"Ds_D0{mode}_Cut"
        print(f'{var_name} = "massDifference(0) <= 0.5 and {lower:.6f} <= daughter(0,dM) <= {upper:.6f}"')
    else:
        var_name = f"Ds_D0kmpippi0_{tag}_Cut"
        print(f'{var_name} = "massDifference(0) <= 0.5 and {lower:.6f} <= daughter(0,dM) <= {upper:.6f}"')

print("\n" + "="*100)

# %% [markdown]
#  ### Background Composition Analysis (D⁰ Level, After 3σ Cut)

# %%
# === USER CONFIGURATION: Choose π⁰ list for kmpippi0 analysis ===
SELECTED_PI0_TAG = "eff10_May2020"  # Change this to analyze different π⁰ lists

print(f"\n{'='*100}")
print(f"BACKGROUND COMPOSITION ANALYSIS (D⁰ Level, After 3σ Mass Cut)")
print(f"{'='*100}")
print(f"Selected π⁰ list for kmpippi0: {SELECTED_PI0_TAG}\n")

def load_generic_mc_d0_tree(mode, tag=None):
    """
    Load generic MC for D⁰ tree (background only, all event types merged).
    
    Parameters:
    -----------
    mode : str
        'kmpip', 'kmpippi0', or 'km3pi'
    tag : str, optional
        π⁰ efficiency tag (only for kmpippi0)
    
    Returns:
    --------
    pd.DataFrame : Merged background data from all event types
    """
    cfg = mode_cfg[mode]
    tree_name = cfg['tree']
    
    # For kmpippi0, we need the tag
    if mode == 'kmpippi0' and tag is None:
        raise ValueError("Must specify π⁰ tag for kmpippi0 mode")
    
    # Build file specifications (same as DstreeCh2 but different tree)
    specs = []
    for ev_type in BKG_EVENT_TYPES:
        if mode == 'kmpippi0':
            # Use tag-specific files for kmpippi0
            mono_path = f"{BKG_BASE}/Ds2D0e-Generic_M2_{tag}_Ds_101525_0_{ev_type}.root"
        else:
            # For kmpip and km3pi, use a default tag (e.g., eff10_May2020)
            # These modes don't depend on π⁰ list, so any file should have the tree
            mono_path = f"{BKG_BASE}/Ds2D0e-Generic_M2_eff10_May2020_Ds_101525_0_{ev_type}.root"
        
        if os.path.isfile(mono_path):
            specs.append(f"{mono_path}:{tree_name}")
        else:
            # Try sharded directory
            if mode == 'kmpippi0':
                shard_dir = f"{BKG_BASE}/Ds2D0e-Generic_M2_{tag}_Ds_101525_0_{ev_type}_shards"
            else:
                shard_dir = f"{BKG_BASE}/Ds2D0e-Generic_M2_eff10_May2020_Ds_101525_0_{ev_type}_shards"
            
            if os.path.isdir(shard_dir):
                import glob
                parts = sorted(glob.glob(os.path.join(shard_dir, "part_*.root")))
                specs.extend([f"{p}:{tree_name}" for p in parts])
    
    if not specs:
        print(f"[WARNING] No background files found for mode {mode}")
        return pd.DataFrame()
    
    # Load only the columns we need
    mass_var = cfg['var']
    mcpdg_var = f"{mode_cfg[mode]['var'].rsplit('_', 1)[0]}_mcPDG"  # e.g., D0_kmpip_mcPDG
    filter_cols = [mass_var, mcpdg_var]
    
    try:
        print(f"  Loading {len(specs)} file(s) from {len(BKG_EVENT_TYPES)} event types... ", end="", flush=True)
        df = uproot.concatenate(specs, filter_name=filter_cols, library="pd")
        print(f"{len(df):,} events")
        return df
    except Exception as e:
        print(f"\n[ERROR] Failed to load: {e}")
        return pd.DataFrame()

# === ANALYZE EACH MODE ===
for mode in ["kmpip", "kmpippi0", "km3pi"]:
    print(f"\n{'='*100}")
    print(f"MODE: {mode.upper()}")
    print(f"{'='*100}")
    
    cfg = mode_cfg[mode]
    
    # Get 3σ window from fit results
    if mode == "kmpippi0":
        fit_key = f"kmpippi0_{SELECTED_PI0_TAG}"
        if fit_key not in fit_results:
            print(f"[SKIP] No fit results for {fit_key}")
            continue
        res_fit = fit_results[fit_key]
        tag_to_use = SELECTED_PI0_TAG
    else:
        if mode not in fit_results:
            print(f"[SKIP] No fit results for {mode}")
            continue
        res_fit = fit_results[mode]
        tag_to_use = None
    
    # Extract 3σ window
    mu = res_fit['mu'][0]
    sL = res_fit['sigmaL'][0]
    sR = res_fit['sigmaR'][0]
    lower_3sig = mu - 3*sL
    upper_3sig = mu + 3*sR
    
    print(f"3σ D⁰ mass window: [{lower_3sig:.5f}, {upper_3sig:.5f}] GeV/c²")
    
    # Load background MC
    print(f"\nLoading background MC for {mode}:")
    df_bkg = load_generic_mc_d0_tree(mode, tag=tag_to_use)
    
    if df_bkg.empty:
        print(f"[SKIP] No background data for {mode}")
        continue
    
    mass_var = cfg['var']
    mcpdg_var = f"{cfg['var'].rsplit('_', 1)[0]}_mcPDG"
    
    print(f"\n[BEFORE 3σ CUT]")
    print(f"Total background events: {len(df_bkg):,}")
    
    # Apply 3σ cut
    df_bkg_cut = df_bkg[(df_bkg[mass_var] >= lower_3sig) & 
                        (df_bkg[mass_var] <= upper_3sig)]
    
    print(f"\n[AFTER 3σ CUT]")
    print(f"Background events passing cut: {len(df_bkg_cut):,}")
    print(f"Cut efficiency: {len(df_bkg_cut)/len(df_bkg)*100:.2f}%")
    
    # Show mcPDG composition
    print(f"\n[D⁰ mcPDG COMPOSITION] (Background only, after 3σ cut)")
    print(f"Variable: {mcpdg_var}\n")
    
    composition = df_bkg_cut[mcpdg_var].value_counts(normalize=True, dropna=False)
    composition_formatted = composition.apply(lambda x: f"{x:.6f}")
    
    # Display as DataFrame for better formatting
    df_comp = pd.DataFrame({
        'mcPDG': composition_formatted.index,
        'Fraction': composition_formatted.values,
        'Count': df_bkg_cut[mcpdg_var].value_counts(dropna=False)[composition_formatted.index].values
    })
    
    display(df_comp)
    
    # Also print in simple format
    print("\nFormatted output:")
    for pdg, frac in composition_formatted.items():
        count = df_bkg_cut[mcpdg_var].value_counts(dropna=False)[pdg]
        print(f"  mcPDG = {pdg:>10} : {frac} ({count:>8,} events)")
    
    # Cleanup
    del df_bkg, df_bkg_cut
    gc.collect()

print(f"\n{'='*100}")
print("BACKGROUND COMPOSITION ANALYSIS COMPLETE")
print(f"{'='*100}")

# %% [markdown]
#  ## Summary
# 
#  **Phase 1 Complete**: Fitted D⁰ mass for all modes and extracted σ values
#  **Phase 2 Complete**: Applied 3σ D⁰ cuts and computed Punzi FoM across π⁰ efficiency tags
# 
#  **CRITICAL FIX APPLIED:**
#  - Signal efficiency now computed with CONSTANT denominator: ε_S = N_sig_pass / 50,000
#  - This gives TRUE total reconstruction efficiency (π⁰ reco × 3σ cut acceptance)
#  - Background count uses raw events: B = len(bkg_df_cut) (no scaling, already at 200 fb⁻¹)
#  - Sanity checks print efficiency breakdown: π⁰ reco efficiency, 3σ acceptance, total efficiency
#  - Legend positioning fixed: Data at top-left, fit parameters at top-right
#  - Mass windows formatted as Python steering script cuts
#  - **NEW:** 3σ boundaries shown with vertical dashed lines + gray shading on excluded regions
# 
#  **Next Steps**:
#  - Validate 3σ windows on sideband data
#  - Optimize BDT for background suppression in Ds⁺ domain
#  - Systematic uncertainties (shape modeling, efficiency variations)