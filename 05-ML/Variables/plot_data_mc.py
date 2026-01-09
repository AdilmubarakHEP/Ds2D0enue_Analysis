"""
Publication-quality plotting module for Data vs MC comparison.

This module provides functions to create:
1. Individual plots for each variable with ratio panels (Data/MC)
2. Grid plots with configurable layout

Key features:
- Luminosity scaling of MC to match data
- Ratio panels showing Data/MC agreement
- No MC truth splitting (MC shown as single category)

Usage:
    from plot_data_mc import DataMCPlotManager

    # Initialize with data
    pm = DataMCPlotManager(df_data, df_mc, mode="kmpip",
                           lum_data=364.093, lum_MC=1443.999)

    # Generate individual plots (with ratio)
    pm.plot_all_individual(variables, output_dir="plots/individual")

    # Generate grid plots (with ratio)
    pm.plot_grid(variables, grid_shape=(3, 4), output_path="plots/grid.png")
"""

import os
import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from pathlib import Path


# Load configuration from the same directory as plot_fake_real_d0.py
CONFIG_PATH = Path(__file__).parent / "plot_config.json"


def load_config():
    """Load plot configuration from JSON file."""
    with open(CONFIG_PATH, 'r') as f:
        return json.load(f)


class DataMCPlotManager:
    """Manager for creating publication-quality Data vs MC comparison plots."""

    def __init__(self, df_data, df_mc, mode, lum_data=364.093, lum_MC=1443.999, config_path=None):
        """
        Initialize DataMCPlotManager.

        Parameters:
            df_data: DataFrame with real data
            df_mc: DataFrame with MC (all backgrounds combined, no truth split)
            mode: Decay mode (kmpip, km3pi, kmpippi0_eff20_May2020)
            lum_data: Data luminosity in fb^-1 (default: 364.093)
            lum_MC: MC luminosity in fb^-1 (default: 1443.999)
            config_path: Optional path to config JSON (default: plot_config.json)
        """
        self.df_data = df_data
        self.df_mc = df_mc
        self.mode = mode
        self.lum_data = lum_data
        self.lum_MC = lum_MC
        self.scale_factor = lum_data / lum_MC

        if config_path:
            with open(config_path, 'r') as f:
                self.config = json.load(f)
        else:
            self.config = load_config()

        self.global_settings = self.config["global_settings"]
        self.grid_settings = self.config["grid_settings"]
        self.mode_title = self.config["mode_titles"].get(mode, mode)

        # Load mode-specific variable configuration
        self.var_config = self.config.get(mode, {})

        # Colors for Data/MC plots
        self.color_data = "black"
        self.color_mc = "#4C6EB1"  # Blue for MC

        # Set up matplotlib defaults
        self._setup_matplotlib()

    def _setup_matplotlib(self):
        """Configure matplotlib for publication-quality plots."""
        plt.rcParams.update({
            "font.family": "serif",
            "font.size": self.global_settings["label_fontsize"],
            "axes.labelsize": self.global_settings["label_fontsize"],
            "xtick.labelsize": self.global_settings["tick_fontsize"],
            "ytick.labelsize": self.global_settings["tick_fontsize"],
            "legend.fontsize": self.global_settings["legend_fontsize"],
            "figure.titlesize": self.global_settings["title_fontsize"],
            "axes.linewidth": 1.2,
            "xtick.major.width": 1.0,
            "ytick.major.width": 1.0,
            "xtick.minor.width": 0.6,
            "ytick.minor.width": 0.6,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "axes.grid": False,
        })

    def _get_var_config(self, var_name):
        """Get configuration for a variable, with defaults."""
        defaults = {
            "nbins": 50,
            "logy": False,
            "logx": False,
            "xrange": None
        }

        var_cfg = self.var_config.get(var_name, {})

        return {
            "xlabel": var_cfg.get("xlabel", var_name),
            "xrange": var_cfg.get("xrange", defaults["xrange"]),
            "logy": var_cfg.get("logy", defaults["logy"]),
            "logx": var_cfg.get("logx", defaults["logx"]),
            "nbins": var_cfg.get("nbins", defaults["nbins"]),
        }

    def _prepare_data(self, var_name, cfg):
        """Prepare data arrays for plotting, applying range and cleaning."""
        x_data = self.df_data[var_name].to_numpy()
        x_mc = self.df_mc[var_name].to_numpy()

        # Remove non-finite values
        x_data = x_data[np.isfinite(x_data)]
        x_mc = x_mc[np.isfinite(x_mc)]

        # Determine range
        if cfg["xrange"] is not None:
            low, high = cfg["xrange"]
        else:
            # Use 1-99 percentile of combined data
            x_combined = np.concatenate([x_data, x_mc])
            low = np.nanpercentile(x_combined, 1)
            high = np.nanpercentile(x_combined, 99)

        # For log x scale, ensure positive range
        if cfg["logx"]:
            if low <= 0:
                all_positive = np.concatenate([x_data[x_data > 0], x_mc[x_mc > 0]])
                if len(all_positive) > 0:
                    low = np.nanpercentile(all_positive, 1)
                else:
                    low = 1e-6
            x_data = x_data[x_data > 0]
            x_mc = x_mc[x_mc > 0]

        return x_data, x_mc, low, high

    def _compute_ratio_and_errors(self, hist_data, hist_mc):
        """Compute Data/MC ratio and propagated errors."""
        err_data = np.sqrt(hist_data)

        # Scale MC
        hist_mc_scaled = self.scale_factor * hist_mc
        err_mc_scaled = self.scale_factor * np.sqrt(hist_mc)

        # Compute ratio
        with np.errstate(divide='ignore', invalid='ignore'):
            ratio = hist_data / hist_mc_scaled
            ratio[np.isnan(ratio)] = 0
            ratio[np.isinf(ratio)] = 0

            # Propagate errors
            err_ratio = ratio * np.sqrt(
                (err_data / np.maximum(hist_data, 1e-10))**2 +
                (err_mc_scaled / np.maximum(hist_mc_scaled, 1e-10))**2
            )
            err_ratio[np.isnan(err_ratio)] = 0
            err_ratio[np.isinf(err_ratio)] = 0

        return ratio, err_ratio, hist_mc_scaled, err_data, err_mc_scaled

    def plot_individual(self, var_name, output_path=None, show=True):
        """
        Create a single Data/MC comparison plot with ratio panel.

        Parameters:
            var_name: Variable name
            output_path: Path to save figure (optional)
            show: Whether to display the plot

        Returns:
            matplotlib figure
        """
        cfg = self._get_var_config(var_name)

        # Check if variable exists
        if var_name not in self.df_data.columns or var_name not in self.df_mc.columns:
            print(f"  Warning: Variable {var_name} not found in data")
            return None

        try:
            x_data, x_mc, low, high = self._prepare_data(var_name, cfg)
        except Exception as e:
            print(f"  Warning: Could not prepare data for {var_name}: {e}")
            return None

        if len(x_data) == 0 or len(x_mc) == 0 or low == high:
            print(f"  Warning: No valid data for {var_name}")
            return None

        # Create figure with ratio panel
        gs = self.global_settings
        fig, (ax1, ax2) = plt.subplots(
            nrows=2, sharex=True, figsize=(6, 6),
            gridspec_kw={"height_ratios": [3, 1]}
        )

        # Create bins
        nbins = cfg["nbins"]
        if cfg["logx"]:
            bins = np.logspace(np.log10(low), np.log10(high), nbins + 1)
        else:
            bins = np.linspace(low, high, nbins + 1)

        bin_centers = 0.5 * (bins[:-1] + bins[1:])
        bin_width = bins[1] - bins[0]

        # Compute histograms
        hist_data, _ = np.histogram(x_data, bins=bins)
        hist_mc, _ = np.histogram(x_mc, bins=bins)

        # Compute ratio and errors
        ratio, err_ratio, hist_mc_scaled, err_data, err_mc_scaled = \
            self._compute_ratio_and_errors(hist_data, hist_mc)

        # Upper panel: MC histogram + Data points
        ax1.hist(x_mc, bins=bins, histtype='step', linewidth=2,
                 color=self.color_mc, label="MC",
                 weights=np.full(len(x_mc), self.scale_factor))

        ax1.errorbar(bin_centers, hist_data, yerr=err_data, fmt='o',
                     color=self.color_data, markersize=3, label="Data",
                     capsize=1, elinewidth=1)

        # Format upper panel - with luminosity on right title
        ax1.set_ylabel(r'Entries / ({:.3g})'.format(bin_width), fontsize=12)
        ax1.set_xlim([low, high])
        ax1.set_title(r"$\int\mathcal{L}dt =\,$" + f"{self.lum_data:.1f}" + r" fb$^{-1}$",
                      loc="right", fontsize=10)
        ax1.legend(loc="upper right", fontsize=11)
        plt.setp(ax1.get_xticklabels(), visible=False)

        if cfg["logy"]:
            ax1.set_yscale("log")
        if cfg["logx"]:
            ax1.set_xscale("log")

        # Lower panel: Ratio
        ax2.axhline(1.0, color='black', lw=1)
        ax2.axhline(1.1, color='gray', lw=1, ls='dashed')
        ax2.axhline(0.9, color='gray', lw=1, ls='dashed')
        ax2.errorbar(bin_centers, ratio, yerr=err_ratio, fmt='o',
                     color=self.color_data, markersize=3)

        ax2.set_ylabel("Data / MC", fontsize=12)
        ax2.set_xlabel(cfg["xlabel"], fontsize=12)
        ax2.set_xlim([low, high])
        ax2.set_ylim(0.5, 1.5)

        plt.tight_layout()
        plt.subplots_adjust(hspace=0.05)

        if output_path:
            os.makedirs(os.path.dirname(output_path), exist_ok=True)
            fig.savefig(output_path, dpi=gs["dpi"], bbox_inches="tight",
                        facecolor="white", edgecolor="none")
            print(f"  Saved: {output_path}")

        if show:
            plt.show()
        else:
            plt.close(fig)

        return fig

    def plot_all_individual(self, variables, output_dir, show=False):
        """
        Create individual plots for all variables.

        Parameters:
            variables: List of variable names
            output_dir: Directory to save plots
            show: Whether to display each plot
        """
        os.makedirs(output_dir, exist_ok=True)
        print(f"\nGenerating {len(variables)} Data/MC plots for mode: {self.mode}")

        for var_name in variables:
            output_path = os.path.join(output_dir, f"{var_name}_DataMC.png")
            self.plot_individual(var_name, output_path=output_path, show=show)

        print(f"Data/MC plots saved to: {output_dir}")

    def _plot_single_with_ratio(self, axes, var_name, show_legend=True, compact=True):
        """
        Plot a single variable with ratio panel on given axes pair.

        Parameters:
            axes: Tuple of (ax_main, ax_ratio) from GridSpec
            var_name: Variable name
            show_legend: Whether to show legend
            compact: If True, reduce font sizes

        Returns:
            bool: Success status
        """
        ax1, ax2 = axes
        cfg = self._get_var_config(var_name)

        # Check if variable exists
        if var_name not in self.df_data.columns or var_name not in self.df_mc.columns:
            ax1.text(0.5, 0.5, f"Not found:\n{var_name}",
                     transform=ax1.transAxes, ha='center', va='center', fontsize=8)
            ax2.axis('off')
            return False

        try:
            x_data, x_mc, low, high = self._prepare_data(var_name, cfg)
        except Exception as e:
            ax1.text(0.5, 0.5, f"Error:\n{str(e)[:20]}",
                     transform=ax1.transAxes, ha='center', va='center', fontsize=8)
            ax2.axis('off')
            return False

        if len(x_data) == 0 or len(x_mc) == 0 or low == high:
            ax1.text(0.5, 0.5, "No valid data",
                     transform=ax1.transAxes, ha='center', va='center', fontsize=8)
            ax2.axis('off')
            return False

        # Create bins
        nbins = cfg["nbins"]
        if cfg["logx"]:
            bins = np.logspace(np.log10(low), np.log10(high), nbins + 1)
        else:
            bins = np.linspace(low, high, nbins + 1)

        bin_centers = 0.5 * (bins[:-1] + bins[1:])

        # Compute histograms
        hist_data, _ = np.histogram(x_data, bins=bins)
        hist_mc, _ = np.histogram(x_mc, bins=bins)

        # Compute ratio and errors
        ratio, err_ratio, hist_mc_scaled, err_data, err_mc_scaled = \
            self._compute_ratio_and_errors(hist_data, hist_mc)

        # Font sizes for compact mode (grid plots) - 1.25x increase for readability
        label_fs = 11 if compact else 12
        tick_fs = 10 if compact else 10
        legend_fs = 10 if compact else 10

        # Compute bin width for y-label
        bin_width = bins[1] - bins[0]

        # Upper panel: MC histogram + Data points
        ax1.hist(x_mc, bins=bins, histtype='step', linewidth=1.5,
                 color=self.color_mc, label="MC",
                 weights=np.full(len(x_mc), self.scale_factor))

        ax1.errorbar(bin_centers, hist_data, yerr=err_data, fmt='o',
                     color=self.color_data, markersize=2, label="Data",
                     capsize=0.5, elinewidth=0.5)

        ax1.set_xlim([low, high])
        ax1.set_ylabel(f"Entries / ({bin_width:.3g})", fontsize=label_fs)
        ax1.tick_params(labelsize=tick_fs)
        plt.setp(ax1.get_xticklabels(), visible=False)

        if show_legend:
            ax1.legend(loc="upper right", fontsize=legend_fs)

        if cfg["logy"]:
            ax1.set_yscale("log")
        if cfg["logx"]:
            ax1.set_xscale("log")

        # Lower panel: Ratio
        ax2.axhline(1.0, color='black', lw=0.8)
        ax2.axhline(1.1, color='gray', lw=0.5, ls='dashed')
        ax2.axhline(0.9, color='gray', lw=0.5, ls='dashed')
        ax2.errorbar(bin_centers, ratio, yerr=err_ratio, fmt='o',
                     color=self.color_data, markersize=2, capsize=0.5, elinewidth=0.5)

        ax2.set_xlabel(cfg["xlabel"], fontsize=label_fs)
        ax2.set_ylabel("Data/MC", fontsize=label_fs)
        ax2.set_xlim([low, high])
        ax2.set_ylim(0.5, 1.5)
        ax2.tick_params(labelsize=tick_fs)
        ax2.yaxis.set_major_locator(MaxNLocator(3))

        return True

    def plot_grid(self, variables, grid_shape=(3, 4), output_path=None,
                  show=True, title=None):
        """
        Create a grid of Data/MC comparison plots with ratio panels.

        Parameters:
            variables: List of variable names
            grid_shape: Tuple (nrows, ncols)
            output_path: Base path for output
            show: Whether to display the plots
            title: Optional super title

        Returns:
            List of matplotlib figures
        """
        nrows, ncols = grid_shape
        slots_per_page = nrows * ncols

        n_vars = len(variables)
        n_pages = (n_vars + slots_per_page - 1) // slots_per_page

        gs_cfg = self.grid_settings
        subplot_size = gs_cfg.get("subplot_size", 3.5)

        figures = []

        for page_idx in range(n_pages):
            start_idx = page_idx * slots_per_page
            end_idx = min(start_idx + slots_per_page, n_vars)
            page_vars = variables[start_idx:end_idx]

            # Create figure with GridSpec for ratio panels
            # Each plot row has 2 sub-rows: main (height 3) + ratio (height 1)
            figsize = (ncols * subplot_size * 1.2, nrows * subplot_size * 1.4)
            fig = plt.figure(figsize=figsize, dpi=self.global_settings["dpi"])

            # Use nested GridSpec: outer grid for variables, inner for main+ratio
            # This allows different spacing between main/ratio vs between variables
            from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec

            # Outer grid: one cell per variable, with space for x-axis labels between rows
            outer_gs = GridSpec(nrows, ncols, figure=fig, hspace=0.18, wspace=0.40)

            # Plot each variable
            for i, var_name in enumerate(page_vars):
                row = i // ncols
                col = i % ncols

                # Inner grid for this variable: main plot + ratio with minimal spacing
                inner_gs = GridSpecFromSubplotSpec(2, 1, subplot_spec=outer_gs[row, col],
                                                   height_ratios=[3, 1], hspace=0.05)

                ax_main = fig.add_subplot(inner_gs[0])
                ax_ratio = fig.add_subplot(inner_gs[1], sharex=ax_main)

                # Hide x-labels on main plot
                plt.setp(ax_main.get_xticklabels(), visible=False)

                # Only show legend on first plot of first page
                show_legend = (i == 0 and page_idx == 0)

                self._plot_single_with_ratio(
                    (ax_main, ax_ratio), var_name,
                    show_legend=show_legend, compact=True
                )

            # Handle empty slots
            for j in range(len(page_vars), slots_per_page):
                row = j // ncols
                col = j % ncols
                inner_gs = GridSpecFromSubplotSpec(2, 1, subplot_spec=outer_gs[row, col],
                                                   height_ratios=[3, 1], hspace=0.05)
                ax_main = fig.add_subplot(inner_gs[0])
                ax_ratio = fig.add_subplot(inner_gs[1])
                ax_main.axis('off')
                ax_ratio.axis('off')

            # Layout first, then add title close to plots
            plt.tight_layout()
            plt.subplots_adjust(top=0.93)  # Push plots up, leaving only 7% at top

            # Super title with luminosity - placed just above the plots
            if page_idx == 0:
                lum_str = r"$\int\mathcal{L}dt =\,$" + f"{self.lum_data:.1f}" + r" fb$^{-1}$"
                plot_title = title if title else f"{self.mode_title} - Data/MC Comparison"
                full_title = f"{plot_title}\n{lum_str}"
                fig.suptitle(full_title, fontsize=self.global_settings["title_fontsize"],
                             y=0.98)

            # Save
            if output_path:
                if page_idx == 0:
                    save_path = output_path
                else:
                    base, ext = os.path.splitext(output_path)
                    suffix = "_continued" if page_idx == 1 else f"_continued{page_idx}"
                    save_path = f"{base}{suffix}{ext}"

                os.makedirs(os.path.dirname(save_path), exist_ok=True)
                fig.savefig(save_path, dpi=self.global_settings["dpi"],
                            bbox_inches="tight", facecolor="white", edgecolor="none")
                print(f"  Saved: {save_path}")

            if show:
                plt.show()
            else:
                plt.close(fig)

            figures.append(fig)

        return figures

    def plot_bdt_output(self, output_path=None, show=True):
        """
        Create a dedicated BDT output comparison plot.

        Parameters:
            output_path: Path to save figure (optional)
            show: Whether to display the plot

        Returns:
            matplotlib figure
        """
        var_name = "Ds_FakeD0BDT"

        if var_name not in self.df_data.columns or var_name not in self.df_mc.columns:
            print(f"  Warning: BDT output {var_name} not found")
            return None

        # Use fixed range for BDT output
        cfg = {
            "xlabel": "Fake D$^0$ BDT Output",
            "xrange": (0, 1),
            "nbins": 50,
            "logy": False,
            "logx": False,
        }

        x_data = self.df_data[var_name].to_numpy()
        x_mc = self.df_mc[var_name].to_numpy()

        x_data = x_data[np.isfinite(x_data)]
        x_mc = x_mc[np.isfinite(x_mc)]

        if len(x_data) == 0 or len(x_mc) == 0:
            print(f"  Warning: No valid BDT data")
            return None

        # Create figure
        gs = self.global_settings
        fig, (ax1, ax2) = plt.subplots(
            nrows=2, sharex=True, figsize=(8, 7),
            gridspec_kw={"height_ratios": [3, 1]}
        )

        bins = np.linspace(0, 1, cfg["nbins"] + 1)
        bin_centers = 0.5 * (bins[:-1] + bins[1:])
        bin_width = bins[1] - bins[0]

        # Compute histograms
        hist_data, _ = np.histogram(x_data, bins=bins)
        hist_mc, _ = np.histogram(x_mc, bins=bins)

        # Compute ratio and errors
        ratio, err_ratio, hist_mc_scaled, err_data, err_mc_scaled = \
            self._compute_ratio_and_errors(hist_data, hist_mc)

        # Upper panel
        ax1.hist(x_mc, bins=bins, histtype='step', linewidth=2.5,
                 color=self.color_mc, label="MC",
                 weights=np.full(len(x_mc), self.scale_factor))

        ax1.errorbar(bin_centers, hist_data, yerr=err_data, fmt='o',
                     color=self.color_data, markersize=4, label="Data",
                     capsize=2, elinewidth=1)

        ax1.set_ylabel(r'Entries / ({:.3g})'.format(bin_width), fontsize=12)
        ax1.set_xlim([0, 1])
        ax1.set_title(r"$\int\mathcal{L}dt =\,$" + f"{self.lum_data:.1f}" + r" fb$^{-1}$",
                      loc="right", fontsize=10)
        ax1.legend(loc="upper right", fontsize=12)
        ax1.tick_params(labelsize=10)
        plt.setp(ax1.get_xticklabels(), visible=False)

        # Lower panel
        ax2.axhline(1.0, color='black', lw=1)
        ax2.axhline(1.1, color='gray', lw=1, ls='dashed')
        ax2.axhline(0.9, color='gray', lw=1, ls='dashed')
        ax2.errorbar(bin_centers, ratio, yerr=err_ratio, fmt='o',
                     color=self.color_data, markersize=4, capsize=2, elinewidth=1)

        ax2.set_ylabel("Data / MC", fontsize=12)
        ax2.set_xlabel(cfg["xlabel"], fontsize=12)
        ax2.set_xlim([0, 1])
        ax2.set_ylim(0.5, 1.5)
        ax2.tick_params(labelsize=10)

        plt.tight_layout()
        plt.subplots_adjust(hspace=0.05)

        if output_path:
            os.makedirs(os.path.dirname(output_path), exist_ok=True)
            fig.savefig(output_path, dpi=gs["dpi"], bbox_inches="tight",
                        facecolor="white", edgecolor="none")
            print(f"  Saved BDT plot: {output_path}")

        if show:
            plt.show()
        else:
            plt.close(fig)

        return fig
