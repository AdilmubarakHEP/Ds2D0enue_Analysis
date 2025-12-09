# Ds2D0e_config.py

import os

# ======================================================================
# Base directories
# ======================================================================

BASE_SIGNAL_DIR = "/home/belle2/amubarak/C01-Simulated_Events/Signal"
BASE_GENERIC_DIR = "/group/belle/users/amubarak/03-KEKCC"

# ======================================================================
# File name patterns
# ======================================================================

SIGNAL_PATTERN = "output_test_{mode}.root"
GENERIC_PATTERN = "Ds2D0e-Generic_Ds_120725_1_{sample}_{mode}.root"

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

def get_signal_file(mode: str) -> str:
    """
    Absolute path to the signal ROOT file for a given mode.
    """
    return os.path.join(BASE_SIGNAL_DIR, SIGNAL_PATTERN.format(mode=mode))


def get_generic_file(sample: str, mode: str) -> str:
    """
    Absolute path to the generic ROOT file for a given (sample, mode).
    """
    return os.path.join(
        BASE_GENERIC_DIR,
        GENERIC_PATTERN.format(sample=sample, mode=mode),
    )
