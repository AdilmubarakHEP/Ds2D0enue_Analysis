"""
Publication-quality plotting module for Fake vs Real D0 distributions.

This module provides functions to create:
1. Individual plots for each variable (saved to individual_plots/ folder)
2. Grid plots with configurable layout (e.g., 3x3, 4x4)

Configuration is loaded from plot_config.json which controls:
- x-axis labels (LaTeX)
- x-axis ranges
- Log scale toggles for both x and y axes
- Number of bins

Usage:
    from plot_fake_real_d0 import PlotManager

    # Initialize with data
    pm = PlotManager(df_real, df_fake, mode="kmpip")

    # Generate individual plots
    pm.plot_all_individual(output_dir="plots/individual")

    # Generate grid plots
    pm.plot_grid(variables=var_list, grid_shape=(3, 3), output_path="plots/grid.png")
"""

import os
import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, LogLocator, LogFormatterSciNotation
import matplotlib.patches as mpatches
from pathlib import Path


# Load configuration
CONFIG_PATH = Path(__file__).parent / "plot_config.json"


def load_config():
    """Load plot configuration from JSON file."""
    with open(CONFIG_PATH, 'r') as f:
        return json.load(f)


class PlotManager:
    """Manager for creating publication-quality fake vs real D0 plots."""

    def __init__(self, df_real, df_fake, mode, config_path=None):
        """
        Initialize PlotManager.

        Parameters:
            df_real: DataFrame with real D0 candidates
            df_fake: DataFrame with fake D0 candidates
            mode: Decay mode (kmpip, km3pi, kmpippi0_eff20_May2020)
            config_path: Optional path to config JSON (default: plot_config.json)
        """
        self.df_real = df_real
        self.df_fake = df_fake
        self.mode = mode

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

        # Set up matplotlib defaults for publication quality
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
        # Default values if variable not found in config
        defaults = {
            "nbins": 50,
            "logy": False,
            "logx": False,
            "xrange": None  # Will use percentile-based auto range
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
        x_real = self.df_real[var_name].to_numpy()
        x_fake = self.df_fake[var_name].to_numpy()

        # Remove non-finite values
        x_real = x_real[np.isfinite(x_real)]
        x_fake = x_fake[np.isfinite(x_fake)]

        # Determine range
        if cfg["xrange"] is not None:
            low, high = cfg["xrange"]
        else:
            # Use 1-99 percentile of combined data
            x_combined = np.concatenate([x_real, x_fake])
            low = np.nanpercentile(x_combined, 1)
            high = np.nanpercentile(x_combined, 99)

        # For log x scale, ensure positive range
        if cfg["logx"]:
            if low <= 0:
                # Find minimum positive value
                all_positive = np.concatenate([x_real[x_real > 0], x_fake[x_fake > 0]])
                if len(all_positive) > 0:
                    low = np.nanpercentile(all_positive, 1)
                else:
                    low = 1e-6
            # Filter to positive values
            x_real = x_real[x_real > 0]
            x_fake = x_fake[x_fake > 0]

        return x_real, x_fake, low, high

    def _plot_single_histogram(self, ax, var_name, show_legend=True, show_xlabel=True,
                               show_ylabel=True, compact=False):
        """
        Plot a single variable histogram on the given axes.

        Parameters:
            ax: matplotlib axes
            var_name: Variable name
            show_legend: Whether to show legend
            show_xlabel: Whether to show x-axis label
            show_ylabel: Whether to show y-axis label
            compact: If True, reduce font sizes for grid plots
        """
        cfg = self._get_var_config(var_name)

        # Check if variable exists in both dataframes
        if var_name not in self.df_real.columns or var_name not in self.df_fake.columns:
            ax.text(0.5, 0.5, f"Variable not found:\n{var_name}",
                   transform=ax.transAxes, ha='center', va='center')
            ax.set_frame_on(True)
            return False

        try:
            x_real, x_fake, low, high = self._prepare_data(var_name, cfg)
        except Exception as e:
            ax.text(0.5, 0.5, f"Data error:\n{str(e)[:30]}",
                   transform=ax.transAxes, ha='center', va='center', fontsize=8)
            return False

        if len(x_real) == 0 or len(x_fake) == 0 or low == high:
            ax.text(0.5, 0.5, "No valid data",
                   transform=ax.transAxes, ha='center', va='center')
            return False

        # Create bins (linear or log-spaced)
        nbins = cfg["nbins"]
        if cfg["logx"]:
            bins = np.logspace(np.log10(low), np.log10(high), nbins + 1)
        else:
            bins = np.linspace(low, high, nbins + 1)

        # Calculate bin width for y-label
        if cfg["logx"]:
            bin_width = "variable"
        else:
            bin_width = (high - low) / nbins

        # Plot settings
        gs = self.global_settings
        histtype = gs.get("histtype", "step")
        use_density = gs.get("use_density", True)

        # Plot real D0 (step histogram - outline only, no fill)
        ax.hist(x_real, bins=bins, range=[low, high],
                histtype=histtype, color=gs["color_real"],
                linewidth=gs["linewidth"], density=use_density,
                label=r"Real $D^{0}$ ($|$mcPDG$|$ = 421)")

        # Plot fake D0 (step histogram - outline only, no fill)
        ax.hist(x_fake, bins=bins, range=[low, high],
                histtype=histtype, color=gs["color_fake"],
                linewidth=gs["linewidth"], density=use_density,
                label=r"Fake $D^{0}$")

        # Apply log scales
        if cfg["logy"]:
            ax.set_yscale("log")
        if cfg["logx"]:
            ax.set_xscale("log")

        # Set labels
        if show_xlabel:
            fontsize = 10 if compact else gs["label_fontsize"]
            ax.set_xlabel(cfg["xlabel"], fontsize=fontsize)

        if show_ylabel:
            fontsize = 10 if compact else gs["label_fontsize"]
            if isinstance(bin_width, str):
                # Log-spaced bins have variable width
                ylabel = "Normalized" if use_density else "Entries"
            else:
                # Include bin width in y-axis label
                ylabel = f"Normalized / ({bin_width:.3g})" if use_density else f"Entries / ({bin_width:.3g})"
            ax.set_ylabel(ylabel, fontsize=fontsize)

        # Set x limits
        ax.set_xlim(low, high)

        # Legend
        if show_legend:
            fontsize = 9 if compact else gs["legend_fontsize"]
            ax.legend(loc=self.grid_settings.get("legend_position", "upper right"),
                     fontsize=fontsize, framealpha=0.9)

        # Tick formatting
        tick_fontsize = 9 if compact else gs["tick_fontsize"]
        ax.tick_params(labelsize=tick_fontsize)

        return True

    def plot_individual(self, var_name, output_path=None, show=True):
        """
        Create a single publication-quality plot for one variable.

        Parameters:
            var_name: Variable name
            output_path: Path to save figure (optional)
            show: Whether to display the plot

        Returns:
            matplotlib figure
        """
        gs = self.global_settings
        figsize = gs.get("figsize_individual", [6, 6])

        fig, ax = plt.subplots(figsize=figsize, dpi=gs["dpi"])

        success = self._plot_single_histogram(ax, var_name, show_legend=True,
                                              show_xlabel=True, show_ylabel=True)

        if success:
            ax.set_title(self.mode_title, loc="left", fontsize=gs["title_fontsize"])

        plt.tight_layout()

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
        print(f"\nGenerating {len(variables)} individual plots for mode: {self.mode}")

        for var_name in variables:
            output_path = os.path.join(output_dir, f"{var_name}.png")
            self.plot_individual(var_name, output_path=output_path, show=show)

        print(f"Individual plots saved to: {output_dir}")

    def plot_grid(self, variables, grid_shape=(3, 3), output_path=None,
                  show=True, title=None):
        """
        Create a grid of plots for multiple variables.

        If there are more variables than grid slots, multiple figures are created.
        If there are fewer variables, empty subplots show the shared legend.

        Parameters:
            variables: List of variable names
            grid_shape: Tuple (nrows, ncols)
            output_path: Base path for output (will append _1.png, _2.png for multiple)
            show: Whether to display the plots
            title: Optional super title (defaults to mode title)

        Returns:
            List of matplotlib figures
        """
        nrows, ncols = grid_shape
        slots_per_page = nrows * ncols

        # Split variables into pages
        n_vars = len(variables)
        n_pages = (n_vars + slots_per_page - 1) // slots_per_page

        gs_cfg = self.grid_settings
        subplot_size = gs_cfg.get("subplot_size", 4.0)

        figures = []

        for page_idx in range(n_pages):
            start_idx = page_idx * slots_per_page
            end_idx = min(start_idx + slots_per_page, n_vars)
            page_vars = variables[start_idx:end_idx]

            # Create figure with square subplots
            figsize = (ncols * subplot_size, nrows * subplot_size)
            fig, axes = plt.subplots(nrows, ncols, figsize=figsize,
                                     dpi=self.global_settings["dpi"])

            # Flatten axes for easy iteration
            if nrows == 1 and ncols == 1:
                axes = np.array([[axes]])
            elif nrows == 1:
                axes = axes.reshape(1, -1)
            elif ncols == 1:
                axes = axes.reshape(-1, 1)

            axes_flat = axes.flatten()

            # Track if we need shared legend
            legend_handles = None
            legend_labels = None

            # Plot each variable
            for i, var_name in enumerate(page_vars):
                ax = axes_flat[i]

                # Only show legend on the FIRST plot of FIRST page
                show_legend = (i == 0 and page_idx == 0)

                success = self._plot_single_histogram(
                    ax, var_name,
                    show_legend=show_legend,
                    show_xlabel=True,
                    show_ylabel=True,
                    compact=True
                )

                # Capture legend handles from first successful plot
                if success and legend_handles is None:
                    legend_handles, legend_labels = ax.get_legend_handles_labels()

            # Handle empty subplots - just turn them off (legend is on first plot)
            n_empty = slots_per_page - len(page_vars)
            if n_empty > 0:
                for j in range(len(page_vars), slots_per_page):
                    ax = axes_flat[j]
                    ax.axis('off')

            # Super title - only on first page, no numbering
            if page_idx == 0:
                plot_title = title if title else self.mode_title
                fig.suptitle(plot_title, fontsize=self.global_settings["title_fontsize"] + 2,
                            y=1.01)

            # Adjust spacing
            plt.tight_layout()
            fig.subplots_adjust(
                wspace=gs_cfg.get("wspace", 0.30),
                hspace=gs_cfg.get("hspace", 0.40),
                top=0.93
            )

            # Save - first page uses original name, subsequent pages add "_continued", "_continued2", etc.
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


def make_all_plots(df_real, df_fake, mode, variables, output_base_dir,
                   grid_shape=(3, 4), show=False):
    """
    Convenience function to generate all plots for a mode.

    Creates:
    - Individual plots in {output_base_dir}/individual/{mode}/
    - Grid plots in {output_base_dir}/grid/{mode}/

    Parameters:
        df_real: DataFrame with real D0
        df_fake: DataFrame with fake D0
        mode: Decay mode string
        variables: List of variable names
        output_base_dir: Base output directory
        grid_shape: Grid dimensions (default 3x4, i.e., 3 rows x 4 cols)
        show: Whether to display plots
    """
    pm = PlotManager(df_real, df_fake, mode)

    # Individual plots
    individual_dir = os.path.join(output_base_dir, "individual", mode)
    pm.plot_all_individual(variables, individual_dir, show=show)

    # Grid plots
    grid_dir = os.path.join(output_base_dir, "grid")
    os.makedirs(grid_dir, exist_ok=True)
    grid_path = os.path.join(grid_dir, f"{mode}_grid.png")
    pm.plot_grid(variables, grid_shape=grid_shape, output_path=grid_path, show=show)

    print(f"\nAll plots generated for mode: {mode}")
    print(f"  Individual: {individual_dir}")
    print(f"  Grid: {grid_path}")
