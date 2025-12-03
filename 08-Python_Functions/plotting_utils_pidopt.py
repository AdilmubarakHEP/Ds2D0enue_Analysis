# === plotting_utils_pidopt.py ===
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# === Color Settings (User Preferred) ===
color_signal = '#007C91'     # Teal
color_background = '#2E2E2E' # Black/Dark gray
color_highlight = '#4C6EB1'

plt.rcParams.update({
    "axes.labelsize": 16,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 16,
    "figure.titlesize": 20
})

# === Optional helper ===
def get_normalization_constant():
    return (1444e15) * (1.329e-9) * 2 * 0.10 * 3e-8 * 0.04 / 2e6

# === Plotting Comparison (used in FoM scan) ===
def plot_var_comparison(ax, df_sig, df_bkg, var, var_cut=None, Bins=30, Range=[0,10], select='right',
                        colors=[color_signal, color_background], labels=['Signal', 'Background']):
    ax.hist([df_sig[var], df_bkg[var]], 
            color=colors, 
            label=labels, 
            density=True, bins=Bins, linewidth=2.5, alpha=1,
            range=Range, stacked=False, histtype='step')
    var_min, var_max = Range
    if var_cut is not None:
        ax.axvline(var_cut, ls='--', color='gray')
        if select == 'left':
            ax.axvspan(var_cut, var_max, color='gray', alpha=0.2)
        elif select == 'right':
            ax.axvspan(var_min, var_cut, color='gray', alpha=0.2)

# === Unified optimization function ===
def optimize_cut(
    df_sig, df_bkg,
    Signal, Background,
    var, FoM,
    xlabel='', Bins=30, Range=[0, 1],
    varname=None, varmin=None, varmax=None,
    select='right', Width=True,
    query_signal='Ds_isSignal==1',
    plot_func=plot_var_comparison,
    plot_labels=['Signal', 'Background'],
    plot_colors=[color_signal, color_background],
    title=r'$D_s^{+} \rightarrow [D^{0} \rightarrow K^{-} \pi^{+}]\;e^+ \nu_e$'
):
    if varname is None:
        varname = var
    operator = '>=' if select == 'right' else '<=' if select == 'left' else None
    if operator is None:
        print('Invalid operator')
        return

    df_sig_plot = df_sig.query(query_signal) if query_signal else df_sig
    if varmin is None:
        varmin = df_sig_plot[var].min()
    if varmax is None:
        varmax = df_sig_plot[var].max()

    cuts = np.linspace(varmin, varmax, Bins)
    FoMs = np.zeros_like(cuts)

    C = get_normalization_constant()
    Scale = 2

    for i, cut in enumerate(cuts):
        nsig = Scale * C * len(Signal.query(f'Ds_isSignal==1 and {FoM} {operator} {cut}'))
        nbkg = Scale * len(Background.query(f'{FoM} {operator} {cut}'))
        FoMs[i] = nsig / np.sqrt(nsig + nbkg) if (nsig + nbkg) > 0 else 0

    cut_optimal = cuts[np.argmax(FoMs)]

    fig, axs = plt.subplots(2, 1, sharex=True, figsize=(8, 12))
    fig.subplots_adjust(hspace=0)

    plot_func(axs[0], df_sig_plot, df_bkg, var, cut_optimal, Bins, Range, select,
              colors=plot_colors, labels=plot_labels)
    axs[0].legend()
    axs[0].set_title(title, loc="left")
    perBin = ((Range[1] - Range[0]) / Bins) * 1000
    axs[0].set_ylabel(
        r'$Entries/(\; {width:.2f}\;MeV/c^2)$'.format(width=perBin) if Width else 'Entries'
    )

    axs[1].plot(cuts, FoMs, color=color_highlight, lw=2)
    axs[1].axvline(cut_optimal, ls='--', lw=2, color='gray')
    axs[1].set_ylabel(r'FoM=$S/\sqrt{S+B}$')
    axs[1].set_xlabel(xlabel)
    axs[1].grid(axis='x')
    plt.xlim(*Range)
    plt.show()

    print(f'FoM maximized at {cut_optimal:.3f}')
    return cut_optimal

# === ROC Curve with Threshold Marker ===
def plot_roc_curve(df_sig, df_bkg, var, query_signal='Ds_isSignal==1', Range=[0, 1], Bins=100, select='right'):
    df_sig_plot = df_sig.query(query_signal) if query_signal else df_sig
    operator = '>=' if select == 'right' else '<='

    cuts = np.linspace(Range[0], Range[1], Bins)
    sig_eff = []
    bkg_rej = []
    fom_vals = []

    for cut in cuts:
        nsig_total = len(df_sig_plot)
        nbkg_total = len(df_bkg)
        nsig_pass = len(df_sig_plot.query(f"{var} {operator} {cut}"))
        nbkg_pass = len(df_bkg.query(f"{var} {operator} {cut}"))

        sig_eff_val = nsig_pass / nsig_total if nsig_total > 0 else 0
        bkg_rej_val = 1 - (nbkg_pass / nbkg_total) if nbkg_total > 0 else 0
        fom = nsig_pass / np.sqrt(nsig_pass + nbkg_pass) if (nsig_pass + nbkg_pass) > 0 else 0

        sig_eff.append(sig_eff_val)
        bkg_rej.append(bkg_rej_val)
        fom_vals.append(fom)

    best_idx = np.argmax(fom_vals)
    best_cut = cuts[best_idx]

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(bkg_rej, sig_eff, lw=2, color=color_signal, label='ROC Curve')
    ax.plot([0, 1], [0, 1], 'k--', lw=1)

    ax.axhline(sig_eff[best_idx], color='black', ls='--', lw=1.5)
    ax.axvline(bkg_rej[best_idx], color='black', ls='--', lw=1.5)
    ax.scatter(bkg_rej[best_idx], sig_eff[best_idx], color='green', s=50,
               label=f'Best Cut = {best_cut:.3f}')

    ax.set_xlabel("Background Rejection")
    ax.set_ylabel("Signal Efficiency")
    ax.set_title("ROC Curve with Best Threshold", loc='left')
    ax.grid(True)
    ax.legend(loc='lower left')
    plt.tight_layout()
    plt.show()

    print(f"ROC Best Threshold = {best_cut:.3f} (FoM max)")
    return best_cut