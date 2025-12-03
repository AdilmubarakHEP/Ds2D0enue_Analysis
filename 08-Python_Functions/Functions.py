import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

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

# === Generalized plotting function ===
def plot_var_comparison(ax, df_sig, df_bkg, var, var_cut=None, Bins=30, Range=[0,10], select='right',
                        colors=['#fd7f6f', '#7eb0d5'], labels=['Signal', 'Background']):
    ax.hist([df_sig[var], df_bkg[var]], 
            color=colors, 
            label=labels, 
            density=True, bins=Bins, linewidth=2, alpha=1,
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
    plot_colors=['#fd7f6f', '#7eb0d5'],
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

    axs[1].plot(cuts, FoMs, color='C1', lw=2)
    axs[1].axvline(cut_optimal, ls='--', lw=2, color='gray')
    axs[1].set_ylabel(r'FoM=$S/\sqrt{S+B}$')
    axs[1].set_xlabel(xlabel)
    axs[1].grid(axis='x')
    plt.xlim(*Range)
    plt.show()

    print(f'FoM maximized at {cut_optimal:.3f}')
    return cut_optimal

# === Plot saving helper ===
def plot_save(df_sig, df_bkg, var, var_cut, select='right', xlabel='', title='', Bins=30, Range=[0,10], Width=True, filename="Plot.png"):
    fig, ax = plt.subplots(figsize=(8, 6))
    plot_var_comparison(ax, df_sig, df_bkg, var, var_cut, Bins, Range, select)
    perBin = ((Range[1] - Range[0]) / Bins) * 1000
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r'$Entries/(\; {width:.2f}\;MeV/c^2)$'.format(width=perBin) if Width else 'Entries')
    ax.set_title(title)
    ax.legend()
    fig.savefig(filename)