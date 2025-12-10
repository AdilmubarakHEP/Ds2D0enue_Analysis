# Ds2D0e_config.py

import os
from typing import Iterable, Optional, Sequence, Union

# ======================================================================
# Base directories
# ======================================================================

BASE_SIGNAL_DIR = "/home/belle2/amubarak/C01-Simulated_Events/Signal"
BASE_GENERIC_DIR = "/group/belle/users/amubarak/03-KEKCC"

# ======================================================================
# Control sample toggles
# ======================================================================

USE_CONTROL_SAMPLE = False  # whether to load control-sample ROOT files
CONTROL_SAMPLE_TAG: Union[str, Sequence[str]] = "noEID"  # "noEID", "wrongCharge", ["noEID", "wrongCharge"], or "all"
CONTROL_SAMPLE_SUFFIXES = {
    None: "",
    "nominal": "",
    "noEID": "_CS_noEID",
    "wrongCharge": "_CS_wrongCharge",
}

# ======================================================================
# File name patterns
# ======================================================================

SIGNAL_PATTERN = "Ds2D0e-Signal_1_{mode}{suffix}.root"
GENERIC_PATTERN = "Ds2D0e-Generic_Ds_120925_1_{sample}_{mode}{suffix}.root"

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
# Background samples (taupair removed)
# ======================================================================

BACKGROUND_SAMPLES = ["ccbar", "charged", "ddbar", "mixed", "ssbar", "uubar"]

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


def _resolve_control_suffixes(
    use_control_sample: Optional[bool] = None,
    control_sample_tag: Optional[Union[str, Iterable[str]]] = None,
) -> list:
    """
    Decide whether to append a control-sample suffix, and allow multiple control tags.
    If use_control_sample is False, returns the nominal suffix only.
    If control_sample_tag is "all", both control samples are returned.
    """
    use_control = USE_CONTROL_SAMPLE if use_control_sample is None else use_control_sample
    if not use_control:
        return [CONTROL_SAMPLE_SUFFIXES[None]]

    tags = _normalize_control_tags(control_sample_tag)
    suffixes = []
    for tag in tags:
        suffix = CONTROL_SAMPLE_SUFFIXES.get(tag)
        if suffix is None:
            valid_tags = [k for k in CONTROL_SAMPLE_SUFFIXES.keys() if k]
            valid_tags.append("all")
            raise ValueError(
                f"Unknown control sample tag '{tag}'. Choose from {valid_tags}."
            )
        suffixes.append(suffix)
    return suffixes


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
    """
    suffixes = _resolve_control_suffixes(use_control_sample, control_sample_tag)
    paths = [
        os.path.join(
            BASE_GENERIC_DIR,
            GENERIC_PATTERN.format(sample=sample, mode=mode, suffix=suffix),
        )
        for suffix in suffixes
    ]
    return paths if len(paths) > 1 else paths[0]
