# Ds2D0e_config.py

import os
from typing import Iterable, Optional, Sequence, Union, Dict, TYPE_CHECKING

# Lazy imports for uproot/pandas - only loaded when loading functions are called
if TYPE_CHECKING:
    import uproot
    import pandas as pd

# ======================================================================
# Base directories
# ======================================================================

BASE_SIGNAL_DIR = "/gpfs/home/belle2/amubarak/C01-Simulated_Events/Signal"
BASE_GENERIC_DIR = "/group/belle/users/amubarak/02-Grid"

# ======================================================================
# Control sample toggles
# ======================================================================

USE_CONTROL_SAMPLE = False  # whether to load control-sample ROOT files
CONTROL_SAMPLE_TAG: Union[str, Sequence[str]] = "noEID"  # "noEID", "wrongCharge", ["noEID", "wrongCharge"], or "all"
CONTROL_SAMPLE_CONFIG = {
    None: {"subdir": "Signal_Region", "suffix": ""},
    "nominal": {"subdir": "Signal_Region", "suffix": ""},
    "noEID": {"subdir": "Full_eID_Range", "suffix": "_CS_noEID"},
    "wrongCharge": {"subdir": "Wrong_Charge", "suffix": "_CS_wrongCharge"},
}
# Legacy alias for backwards compatibility
CONTROL_SAMPLE_SUFFIXES = {k: v["suffix"] for k, v in CONTROL_SAMPLE_CONFIG.items()}

# ======================================================================
# File name patterns
# ======================================================================

SIGNAL_PATTERN = "Ds2D0e-Signal_122225_0_{mode}{suffix}.root"
GENERIC_PATTERN = "Ds2D0e-Generic_{sample}_{mode}_122225_0{suffix}.root"
DATA_PATTERN = "Ds2D0e-Data_122225_0_{mode}{suffix}.root"

# ======================================================================
# Decay modes, trees, and D0 mass window cuts
# ======================================================================

DECAY_CONFIG = {
    "kmpip": {
        "ds_tree": "DstreeCh1",
        "cut": "(-0.014291 <= D0_dM) & (D0_dM <= 0.014152)",
    },
    "km3pi": {
        "ds_tree": "DstreeCh3",
        "cut": "(-0.013093 <= D0_dM) & (D0_dM <= 0.012520)",
    },
    "kmpippi0_eff20_May2020": {
        "ds_tree": "DstreeCh2",
        "cut": "(-0.052152 <= D0_dM) & (D0_dM <= 0.024237)",
    },
}

# ======================================================================
# D0 mass sideband configuration for Data/MC comparison
# full_range: from steering script (broader cut applied during reconstruction)
# signal_window: tight cut from DECAY_CONFIG (used for signal selection)
# sideband_cut: pandas query string to select events OUTSIDE the signal window
# ======================================================================

D0_SIDEBAND_CONFIG = {
    "kmpip": {
        "full_range": (-0.04, 0.04),  # GeV, from steering script
        "signal_window": (-0.014291, 0.014152),  # GeV, tight cut
        "sideband_cut": "(D0_dM < -0.014291) | (D0_dM > 0.014152)",
    },
    "km3pi": {
        "full_range": (-0.04, 0.04),  # GeV, from steering script
        "signal_window": (-0.013093, 0.012520),  # GeV, tight cut
        "sideband_cut": "(D0_dM < -0.013093) | (D0_dM > 0.012520)",
    },
    "kmpippi0_eff20_May2020": {
        "full_range": (-0.20, 0.15),  # GeV, from steering script (broader for pi0)
        "signal_window": (-0.052152, 0.024237),  # GeV, tight cut
        "sideband_cut": "(D0_dM < -0.052152) | (D0_dM > 0.024237)",
    },
}

# ======================================================================
# Background samples (BB = charged + mixed combined)
# ======================================================================

BACKGROUND_SAMPLES = ["BB", "ccbar", "ddbar", "ssbar", "taupair", "uubar"]

# ======================================================================
# Helper functions
# ======================================================================


def _normalize_control_tags(control_sample_tag: Optional[Union[str, Iterable[str]]]):
    """Return a list of control-sample tag strings from a user input."""
    tag = CONTROL_SAMPLE_TAG if control_sample_tag is None else control_sample_tag

    if isinstance(tag, str):
        if tag.lower() == "all":
            return [t for t in CONTROL_SAMPLE_SUFFIXES.keys() if t not in (None, "nominal")]
        return [tag]

    try:
        tags = list(tag)
    except TypeError:
        tags = [tag]

    if len(tags) == 0:
        raise ValueError("control_sample_tag cannot be empty")

    return tags


def _resolve_control_configs(
    use_control_sample: Optional[bool] = None,
    control_sample_tag: Optional[Union[str, Iterable[str]]] = None,
) -> list:
    """
    Decide which control-sample config(s) to use, returning a list of dicts with 'subdir' and 'suffix'.
    If use_control_sample is False, returns the nominal config only.
    If control_sample_tag is "all", both control samples are returned.
    """
    use_control = USE_CONTROL_SAMPLE if use_control_sample is None else use_control_sample
    if not use_control:
        return [CONTROL_SAMPLE_CONFIG[None]]

    tags = _normalize_control_tags(control_sample_tag)
    configs = []
    for tag in tags:
        config = CONTROL_SAMPLE_CONFIG.get(tag)
        if config is None:
            valid_tags = [k for k in CONTROL_SAMPLE_CONFIG.keys() if k]
            valid_tags.append("all")
            raise ValueError(
                f"Unknown control sample tag '{tag}'. Choose from {valid_tags}."
            )
        configs.append(config)
    return configs


def _resolve_control_suffixes(
    use_control_sample: Optional[bool] = None,
    control_sample_tag: Optional[Union[str, Iterable[str]]] = None,
) -> list:
    """Legacy function: returns just the suffixes for backwards compatibility."""
    configs = _resolve_control_configs(use_control_sample, control_sample_tag)
    return [c["suffix"] for c in configs]


def get_signal_file(
    mode: str,
    use_control_sample: Optional[bool] = None,
    control_sample_tag: Optional[Union[str, Iterable[str]]] = None,
) -> Union[str, list[str]]:
    """
    Absolute path(s) to the signal ROOT file(s) for a given mode.
    Returns a string for a single file, or a list of strings when multiple control tags are requested.
    """
    suffixes = _resolve_control_suffixes(use_control_sample, control_sample_tag)
    paths = [
        os.path.join(BASE_SIGNAL_DIR, SIGNAL_PATTERN.format(mode=mode, suffix=suffix))
        for suffix in suffixes
    ]
    return paths if len(paths) > 1 else paths[0]


def get_generic_file(
    sample: str,
    mode: str,
    use_control_sample: Optional[bool] = None,
    control_sample_tag: Optional[Union[str, Iterable[str]]] = None,
) -> Union[str, list[str]]:
    """
    Absolute path(s) to the generic ROOT file(s) for a given (sample, mode).
    Returns a string for a single file, or a list of strings when multiple control tags are requested.
    Files are organized in subdirectories: Signal_Region, Full_eID_Range, Wrong_Charge.
    """
    configs = _resolve_control_configs(use_control_sample, control_sample_tag)
    paths = [
        os.path.join(
            BASE_GENERIC_DIR,
            cfg["subdir"],
            GENERIC_PATTERN.format(sample=sample, mode=mode, suffix=cfg["suffix"]),
        )
        for cfg in configs
    ]
    return paths if len(paths) > 1 else paths[0]


def get_data_file(
    mode: str,
    use_control_sample: Optional[bool] = None,
    control_sample_tag: Optional[Union[str, Iterable[str]]] = None,
) -> Union[str, list[str]]:
    """
    Absolute path(s) to the real data ROOT file(s) for a given mode.
    Returns a string for a single file, or a list of strings when multiple control tags are requested.
    Files are organized in subdirectories: Signal_Region, Full_eID_Range, Wrong_Charge.
    """
    configs = _resolve_control_configs(use_control_sample, control_sample_tag)
    paths = [
        os.path.join(
            BASE_GENERIC_DIR,
            cfg["subdir"],
            DATA_PATTERN.format(mode=mode, suffix=cfg["suffix"]),
        )
        for cfg in configs
    ]
    return paths if len(paths) > 1 else paths[0]


# ======================================================================
# Uproot loading settings
# ======================================================================

UPROOT_NUM_WORKERS = 16

DEFAULT_BRANCH_FILTER = [
    "Ds_mcPDG",
    "D0_*",                  # includes D0_dM and D0 truth
    "Ds_massDifference_0",   # Δm_e
    "Ds_diff_D0pi",          # Δm_pi
    "Ds_isSignal",
    "D0_genMotherPDG",
    "e_mcPDG",
    "e_genMotherPDG",
    "K_*",
    "pi_*",
    "pi0_*",
    "pi1_*",                 # for km3pi mode
    "pi2_*",                 # for km3pi mode
    "pi3_*",                 # for km3pi mode
]


# ======================================================================
# Data loading functions
# ======================================================================


def _build_tree_paths(file_or_files: Union[str, list], tree_name: str) -> list:
    """Build uproot-compatible paths with tree names."""
    if isinstance(file_or_files, (list, tuple, set)):
        return [f"{f}:{tree_name}" for f in file_or_files]
    return [f"{file_or_files}:{tree_name}"]


def _load_root_file(
    tree_paths: list,
    cut: Optional[str] = None,
    branch_filter: Optional[list] = None,
    num_workers: Optional[int] = None,
) -> "pd.DataFrame":
    """
    Load ROOT file(s) into a pandas DataFrame using uproot.

    Parameters
    ----------
    tree_paths : list
        List of "filepath:treename" strings
    cut : str, optional
        Query string to apply after loading (e.g., D0 mass cut)
    branch_filter : list, optional
        Branch name patterns to load (default: DEFAULT_BRANCH_FILTER)
    num_workers : int, optional
        Number of parallel workers (default: UPROOT_NUM_WORKERS)

    Returns
    -------
    pd.DataFrame
        Loaded and optionally filtered data
    """
    import uproot
    import pandas as pd

    if branch_filter is None:
        branch_filter = DEFAULT_BRANCH_FILTER
    if num_workers is None:
        num_workers = UPROOT_NUM_WORKERS

    try:
        df = uproot.concatenate(
            tree_paths,
            library="pd",
            filter_name=branch_filter,
            num_workers=num_workers,
        )
    except TypeError:
        # fallback if uproot version does not support num_workers
        df = uproot.concatenate(
            tree_paths,
            library="pd",
            filter_name=branch_filter,
        )

    if cut is not None:
        df = df.query(cut)

    return df


def load_signal(
    mode: str,
    apply_d0_cut: bool = True,
    use_control_sample: Optional[bool] = None,
    control_sample_tag: Optional[Union[str, Iterable[str]]] = None,
    branch_filter: Optional[list] = None,
    num_workers: Optional[int] = None,
) -> "pd.DataFrame":
    """
    Load signal ROOT file(s) for a given decay mode.

    Parameters
    ----------
    mode : str
        Decay mode key (e.g., "kmpip", "km3pi", "kmpippi0_eff20_May2020")
    apply_d0_cut : bool
        Whether to apply the D0 mass window cut (default: True)
    use_control_sample : bool, optional
        Override global USE_CONTROL_SAMPLE setting
    control_sample_tag : str or list, optional
        Override global CONTROL_SAMPLE_TAG setting
    branch_filter : list, optional
        Branch patterns to load (default: DEFAULT_BRANCH_FILTER)
    num_workers : int, optional
        Number of parallel workers (default: UPROOT_NUM_WORKERS)

    Returns
    -------
    pd.DataFrame
        Loaded signal data
    """
    config = DECAY_CONFIG[mode]
    signal_file = get_signal_file(mode, use_control_sample, control_sample_tag)
    tree_paths = _build_tree_paths(signal_file, config["ds_tree"])
    cut = config["cut"] if apply_d0_cut else None

    return _load_root_file(tree_paths, cut, branch_filter, num_workers)


def load_generic(
    sample: str,
    mode: str,
    apply_d0_cut: bool = True,
    use_control_sample: Optional[bool] = None,
    control_sample_tag: Optional[Union[str, Iterable[str]]] = None,
    branch_filter: Optional[list] = None,
    num_workers: Optional[int] = None,
) -> "pd.DataFrame":
    """
    Load a single generic (background) ROOT file for a given sample and mode.

    Parameters
    ----------
    sample : str
        Background sample name (e.g., "BB", "ccbar", "taupair")
    mode : str
        Decay mode key (e.g., "kmpip", "km3pi", "kmpippi0_eff20_May2020")
    apply_d0_cut : bool
        Whether to apply the D0 mass window cut (default: True)
    use_control_sample : bool, optional
        Override global USE_CONTROL_SAMPLE setting
    control_sample_tag : str or list, optional
        Override global CONTROL_SAMPLE_TAG setting
    branch_filter : list, optional
        Branch patterns to load (default: DEFAULT_BRANCH_FILTER)
    num_workers : int, optional
        Number of parallel workers (default: UPROOT_NUM_WORKERS)

    Returns
    -------
    pd.DataFrame
        Loaded background data
    """
    config = DECAY_CONFIG[mode]
    generic_file = get_generic_file(sample, mode, use_control_sample, control_sample_tag)
    tree_paths = _build_tree_paths(generic_file, config["ds_tree"])
    cut = config["cut"] if apply_d0_cut else None

    return _load_root_file(tree_paths, cut, branch_filter, num_workers)


def load_all_generic(
    mode: str,
    samples: Optional[list] = None,
    apply_d0_cut: bool = True,
    use_control_sample: Optional[bool] = None,
    control_sample_tag: Optional[Union[str, Iterable[str]]] = None,
    branch_filter: Optional[list] = None,
    num_workers: Optional[int] = None,
) -> "Dict[str, pd.DataFrame]":
    """
    Load all background samples for a given mode.

    Parameters
    ----------
    mode : str
        Decay mode key (e.g., "kmpip", "km3pi", "kmpippi0_eff20_May2020")
    samples : list, optional
        List of sample names to load (default: BACKGROUND_SAMPLES)
    apply_d0_cut : bool
        Whether to apply the D0 mass window cut (default: True)
    use_control_sample : bool, optional
        Override global USE_CONTROL_SAMPLE setting
    control_sample_tag : str or list, optional
        Override global CONTROL_SAMPLE_TAG setting
    branch_filter : list, optional
        Branch patterns to load (default: DEFAULT_BRANCH_FILTER)
    num_workers : int, optional
        Number of parallel workers (default: UPROOT_NUM_WORKERS)

    Returns
    -------
    dict
        Dictionary mapping sample names to DataFrames
    """
    if samples is None:
        samples = BACKGROUND_SAMPLES

    result = {}
    for sample in samples:
        result[sample] = load_generic(
            sample, mode, apply_d0_cut,
            use_control_sample, control_sample_tag,
            branch_filter, num_workers
        )
    return result


def load_data(
    mode: str,
    apply_d0_cut: bool = True,
    use_control_sample: Optional[bool] = None,
    control_sample_tag: Optional[Union[str, Iterable[str]]] = None,
    branch_filter: Optional[list] = None,
    num_workers: Optional[int] = None,
) -> "pd.DataFrame":
    """
    Load real data ROOT file(s) for a given decay mode.

    Parameters
    ----------
    mode : str
        Decay mode key (e.g., "kmpip", "km3pi", "kmpippi0_eff20_May2020")
    apply_d0_cut : bool
        Whether to apply the D0 mass window cut (default: True)
    use_control_sample : bool, optional
        Override global USE_CONTROL_SAMPLE setting
    control_sample_tag : str or list, optional
        Override global CONTROL_SAMPLE_TAG setting
    branch_filter : list, optional
        Branch patterns to load (default: DEFAULT_BRANCH_FILTER)
    num_workers : int, optional
        Number of parallel workers (default: UPROOT_NUM_WORKERS)

    Returns
    -------
    pd.DataFrame
        Loaded real data
    """
    config = DECAY_CONFIG[mode]
    data_file = get_data_file(mode, use_control_sample, control_sample_tag)
    tree_paths = _build_tree_paths(data_file, config["ds_tree"])
    cut = config["cut"] if apply_d0_cut else None

    return _load_root_file(tree_paths, cut, branch_filter, num_workers)
