#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
import json
import shutil
import subprocess
import argparse
from typing import List, Optional, Tuple, Dict
import logging
from datetime import datetime

logging.basicConfig(level=logging.INFO)

import b2luigi
from b2luigi.basf2_helper import Basf2PathTask
from b2luigi.batch.processes.gbasf2 import Gbasf2GridProjectTarget, get_unique_project_name

# Disable automatic dataset download - we'll use explicit DsDownloadTask instead
# This gives us proper Luigi parallelism for downloads
b2luigi.set_setting("gbasf2_download_dataset", False)

# ============================================================
# Paths and configuration
# ============================================================

# Script directory (where this file and settings.json live - in home directory)
SCRIPT_DIR = "/home/belle2/amubarak/Ds2D0enue_Analysis/04-b2luigi/gbasf2"

# Working directory for downloaded ROOT files (use GROUP SPACE with large quota!)
# This is where temporary downloaded chunk files will be saved
WORK_DIR = "/group/belle2/users2022/amubarak/gbasf2_workspace"

SETTINGS_PATH = os.path.join(SCRIPT_DIR, "settings.json")

# Logs in home directory for easy access
LOG_ROOT = os.path.join(SCRIPT_DIR, "logs")

# Final merged outputs
MERGE_ROOT = "/group/belle/users/amubarak/02-Grid"

# Signal MC outputs and inputs
SIGNAL_OUTPUT_ROOT = "/home/belle2/amubarak/C01-Simulated_Events/Signal"
SIGNAL_INPUT_DIR = "/group/belle/users/amubarak/00-Generation/Signal"
SIGNAL_INPUTS = {
    "kmpip": "ccbarDs+EventGeneration_Mode1.root",
    "kmpippi0": "ccbarDs+EventGeneration_Mode2.root",
    "km3pi": "ccbarDs+EventGeneration_Mode3.root",
}

# Modes
M_KMPIP = "kmpip"
M_KMPIPPI0 = "kmpippi0"
M_KM3PI = "km3pi"

# Global, filled after reading settings.json
SELECTED_PI0_LISTS: List[str] = []
ENABLED_CONTROL_SAMPLES: dict = {}
RUN_ON_DATA: bool = False


# ============================================================
# Lazy import steering: create_reconstruction_path
# ============================================================

# Steering script path
STEERING_FILE = "/home/belle2/amubarak/Ds2D0enue_Analysis/02-Reconstruction_Scripts/Ds2D0e-Reconstruction.py"

# Cache for lazy loading
_steering_module = None


def get_create_reconstruction_path():
    """
    Lazily import the steering script. This allows the b2luigi script to run
    without basf2 being sourced (e.g., for --dry-run), while still working
    when basf2 is available (on the grid or after sourcing locally).
    """
    global _steering_module
    if _steering_module is None:
        import importlib.util
        _spec = importlib.util.spec_from_file_location("ds_reco", STEERING_FILE)
        if _spec is None or _spec.loader is None:
            raise RuntimeError(f"Cannot import steering from {STEERING_FILE}")
        _steering_module = importlib.util.module_from_spec(_spec)
        _spec.loader.exec_module(_steering_module)
    return _steering_module.create_reconstruction_path

# ============================================================
# Config helpers
# ============================================================

def get_auto_date() -> str:
    """Generate date string in MMDDYY format from current date."""
    from datetime import datetime
    return datetime.now().strftime("%m%d%y")


def load_settings(path: str) -> dict:
    with open(path, "r") as f:
        cfg = json.load(f)

    # Auto-generate date if set to "auto" or not present
    if cfg.get("date", "auto").lower() == "auto":
        cfg["date"] = get_auto_date()
        logging.info(f"Auto-generated date: {cfg['date']}")

    if "attempt" not in cfg:
        raise ValueError("settings.json must contain 'attempt'")
    if "pi0_lists" not in cfg:
        raise ValueError("settings.json must contain 'pi0_lists'")
    # Default control samples: signal only
    if "control_samples" not in cfg:
        cfg["control_samples"] = {"signal": True, "no_electron_id": False, "wrong_charge": False, "both": False}
    # Default: run on MC
    if "run_on_data" not in cfg:
        cfg["run_on_data"] = False
    return cfg


def control_sample_flags(sample_type: str) -> Tuple[bool, bool]:
    """
    Convert sample_type string to (no_electron_id, wrong_charge) flags.

    Args:
        sample_type: one of "signal", "no_electron_id", "wrong_charge", "both"

    Returns:
        Tuple of (no_electron_id, wrong_charge) booleans
    """
    if sample_type == "signal":
        return (False, False)
    elif sample_type == "no_electron_id":
        return (True, False)
    elif sample_type == "wrong_charge":
        return (False, True)
    elif sample_type == "both":
        return (True, True)
    else:
        raise ValueError(f"Unknown sample_type: {sample_type}")


def control_sample_suffix(sample_type: str) -> str:
    """
    Get output filename suffix for a given sample type.

    Returns:
        String suffix like "" (signal), "_CS_noEID", "_CS_wrongCharge", "_CS_noEID_wrongCharge"
    """
    if sample_type == "signal":
        return ""
    elif sample_type == "no_electron_id":
        return "_CS_noEID"
    elif sample_type == "wrong_charge":
        return "_CS_wrongCharge"
    elif sample_type == "both":
        return "_CS_noEID_wrongCharge"
    else:
        raise ValueError(f"Unknown sample_type: {sample_type}")


def control_sample_folder(sample_type: str) -> str:
    """
    Get output folder name for a given sample type.

    Returns:
        Folder name like "Signal_Region", "Full_eID_Range", "Wrong_Charge", "No_eID_Wrong_Charge"
    """
    if sample_type == "signal":
        return "Signal_Region"
    elif sample_type == "no_electron_id":
        return "Full_eID_Range"
    elif sample_type == "wrong_charge":
        return "Wrong_Charge"
    elif sample_type == "both":
        return "No_eID_Wrong_Charge"
    else:
        raise ValueError(f"Unknown sample_type: {sample_type}")

# ============================================================
# Helpers
# ============================================================

def get_project_name_prefix(date: str, attempt: str, mode: str,
                            sample_type: str, is_data: bool) -> str:
    """
    Create gbasf2 project name prefix (< 22 chars, b2luigi adds hash).
    Examples:
      Ds_120725_1_kmpip_sig    (data=False, sample_type=signal)
      Ds_120725_1_kmpip_noEID  (data=False, sample_type=no_electron_id)
      DsD_120725_1_kmpip       (data=True, sample_type=signal)
    """
    prefix = "DsD" if is_data else "Ds"
    cs_code = ""
    if sample_type == "signal":
        cs_code = "sig"
    elif sample_type == "no_electron_id":
        cs_code = "noEID"
    elif sample_type == "wrong_charge":
        cs_code = "wch"
    elif sample_type == "both":
        cs_code = "both"

    # Format: Ds_DDMMYY_A_mode_cs
    return f"{prefix}_{date}_{attempt}_{mode}_{cs_code}"[:21]


def merged_output_path(date: str, attempt: str, mode: str,
                      pi0: Optional[str], sample_type: str, is_data: bool,
                      campaign: Optional[str] = None, background: Optional[str] = None) -> str:
    """
    /group/belle/users/amubarak/02-Grid/<sample_folder>/
        Ds2D0e-Data_120725_1_kmpip[_CS_suffix].root  (if is_data=True)
        Ds2D0e-Generic_ccbar_kmpip_120725_1[_CS_suffix].root    (if is_data=False)

    Sample folders:
        Signal_Region, Full_eID_Range, Wrong_Charge, No_eID_Wrong_Charge
    """
    cs_suffix = control_sample_suffix(sample_type)
    cs_folder = control_sample_folder(sample_type)

    if is_data:
        # Data naming
        if mode == M_KMPIPPI0 and pi0:
            name = f"Ds2D0e-Data_{date}_{attempt}_{mode}_{pi0}{cs_suffix}.root"
        else:
            name = f"Ds2D0e-Data_{date}_{attempt}_{mode}{cs_suffix}.root"
    else:
        # MC naming: include background
        bg_str = background if background else "MC"
        if mode == M_KMPIPPI0 and pi0:
            name = f"Ds2D0e-Generic_{bg_str}_{mode}_{pi0}_{date}_{attempt}{cs_suffix}.root"
        else:
            name = f"Ds2D0e-Generic_{bg_str}_{mode}_{date}_{attempt}{cs_suffix}.root"

    return os.path.join(MERGE_ROOT, cs_folder, name)


def signal_output_path(date: str, attempt: str, mode: str,
                       pi0: Optional[str], sample_type: str) -> str:
    """
    Signal MC output path:
    /group/belle/users/amubarak/02-Grid/Ds2D0e-Signal_<date>_<attempt>_<mode>[_<pi0>][_CS_suffix].root
    """
    cs_suffix = control_sample_suffix(sample_type)
    if mode == M_KMPIPPI0 and pi0:
        name = f"Ds2D0e-Signal_{date}_{attempt}_{mode}_{pi0}{cs_suffix}.root"
    else:
        name = f"Ds2D0e-Signal_{date}_{attempt}_{mode}{cs_suffix}.root"
    return os.path.join(SIGNAL_OUTPUT_ROOT, name)


def get_signal_input(mode: str) -> str:
    """Get the full path to the signal input file for a given mode."""
    filename = SIGNAL_INPUTS.get(mode)
    if not filename:
        raise ValueError(f"Unknown mode for signal: {mode}")
    return os.path.join(SIGNAL_INPUT_DIR, filename)


def _which_hadd() -> str:
    exe = shutil.which("hadd")
    if not exe:
        raise RuntimeError("Cannot find 'hadd' in PATH. Source basf2/ROOT first.")
    return exe


def _run_hadd(out_file: str, inputs: List[str]) -> str:
    """
    Run hadd to merge ROOT files. Uses batch merging if there are many files
    to avoid ROOT TBuffer overflow (limit ~1GB).

    Returns: Path to batch directory if batches were created, None otherwise.
    """
    # If few files, merge directly
    if len(inputs) <= 100:
        cmd = [_which_hadd(), "-f", out_file] + inputs
        subprocess.check_call(cmd)
        return None

    # Batch merging for large file counts
    import tempfile
    batch_size = 100
    # Create official folder for batch files (not hidden)
    batch_dir_name = os.path.basename(out_file).replace(".root", "_batches")
    batch_dir = os.path.join(os.path.dirname(out_file), batch_dir_name)
    os.makedirs(batch_dir, exist_ok=True)

    # Step 1: Merge in batches
    intermediate_files = []
    for i in range(0, len(inputs), batch_size):
        batch = inputs[i:i+batch_size]
        intermediate = os.path.join(batch_dir, f"batch_{i//batch_size:04d}.root")
        cmd = [_which_hadd(), "-f", intermediate] + batch
        subprocess.check_call(cmd)
        intermediate_files.append(intermediate)
        logging.info(f"Merged batch {i//batch_size + 1}/{(len(inputs)-1)//batch_size + 1} ({len(batch)} files)")

    # Step 2: Merge intermediate files (if <= 100, otherwise keep them as-is)
    if len(intermediate_files) <= 100:
        logging.info(f"Merging {len(intermediate_files)} intermediate files into final output")
        cmd = [_which_hadd(), "-f", out_file] + intermediate_files
        subprocess.check_call(cmd)
        # Don't delete batches yet - caller should delete after validation
        logging.info(f"Batch merge complete. Batches saved in {batch_dir} (will be deleted after validation)")
        return batch_dir
    else:
        # Too many intermediate files - keep them as-is
        logging.warning(f"Too many intermediate files ({len(intermediate_files)}). Keeping them in {batch_dir}")
        logging.warning(f"You can manually merge later or re-run with smaller batch_size")
        raise RuntimeError(f"Cannot merge {len(intermediate_files)} intermediate files (limit 100). Check {batch_dir}")


def _expected_trees(mode: str) -> Tuple[str, str]:
    """
    Tree names as produced by your steering.
    """
    if mode == M_KMPIP:
        return "DstreeCh1", "D02kmpiptree"
    if mode == M_KM3PI:
        return "DstreeCh3", "D02km3pitree"
    return "DstreeCh2", "D02kmpippi0tree"


def _validate_root_has_two_trees(path: str, mode: str) -> bool:
    try:
        import uproot
        if (not os.path.isfile(path)) or os.path.getsize(path) == 0:
            return False
        t_ds, t_d0 = _expected_trees(mode)
        with uproot.open(path) as f:
            if t_ds not in f or t_d0 not in f:
                return False
            t1, t2 = f[t_ds], f[t_d0]
            n1 = int(getattr(t1, "num_entries", 0) or len(t1))
            n2 = int(getattr(t2, "num_entries", 0) or len(t2))
            if n1 > 0:
                _ = t1.arrays(entry_start=0, entry_stop=min(5, n1), library="np")
            if n2 > 0:
                _ = t2.arrays(entry_start=0, entry_stop=min(5, n2), library="np")
        return True
    except Exception:
        return False


def _keep_only_stderr(log_dir: str) -> None:
    """
    For a given log directory, keep only files that look like stderr,
    strip stdout and other junk. Remove empty dirs.
    """
    if not os.path.isdir(log_dir):
        return
    try:
        for root, dirs, files in os.walk(log_dir, topdown=False):
            for name in files:
                lower = name.lower()
                if "stderr" in lower:
                    continue
                try:
                    os.remove(os.path.join(root, name))
                except Exception:
                    pass
            for d in dirs:
                p = os.path.join(root, d)
                try:
                    os.rmdir(p)
                except OSError:
                    pass
        try:
            if not os.listdir(log_dir):
                os.rmdir(log_dir)
        except Exception:
            pass
    except Exception:
        pass

# ============================================================
# Grid Helper Functions
# ============================================================

# Directory for download logs
DOWNLOAD_LOG_DIR = os.path.join(SCRIPT_DIR, "worker_logs")

# Cache for gbasf2 environment
_gbasf2_env_cache = None

def get_gbasf2_env_cached() -> dict:
    """
    Get gbasf2 environment, caching it for reuse.
    Sources gbasf2 bashrc and captures only proper environment variables.
    """
    global _gbasf2_env_cache
    if _gbasf2_env_cache is not None:
        return _gbasf2_env_cache

    try:
        # Source gbasf2 and capture environment (using env -0 for null-separated output)
        import shlex
        cmd = shlex.split("bash -c 'source /cvmfs/belle.kek.jp/grid/gbasf2/pro/bashrc && env -0'")
        result = subprocess.run(cmd, capture_output=True, timeout=60)
        if result.returncode == 0:
            # Parse null-separated environment, filtering out bash functions
            env = {}
            for item in result.stdout.split(b'\x00'):
                if b'=' in item:
                    try:
                        line = item.decode('utf-8', errors='replace')
                        key, value = line.split('=', 1)
                        # Skip bash function definitions (they start with BASH_FUNC_)
                        # and any keys with problematic characters
                        if not key.startswith('BASH_FUNC_') and key.replace('_', '').isalnum():
                            env[key] = value
                    except (ValueError, UnicodeDecodeError):
                        pass
            _gbasf2_env_cache = env
            logging.info(f"Loaded gbasf2 environment with {len(env)} variables")
            return _gbasf2_env_cache
    except Exception as e:
        logging.warning(f"Could not source gbasf2: {e}")

    # Fallback: return current environment
    return os.environ.copy()

def get_gbasf2_projects(user: str = "adilmub") -> List[Dict]:
    """Query grid for all projects belonging to user."""
    try:
        env = get_gbasf2_env_cached()
        result = subprocess.run(
            ["gb2_project_summary", "-u", user],
            capture_output=True, text=True, timeout=120,
            env=env
        )
        projects = []
        for line in result.stdout.splitlines():
            parts = line.split()
            if len(parts) >= 6 and parts[2] == "Good":
                projects.append({
                    "name": parts[0],
                    "owner": parts[1],
                    "status": parts[2],
                    "done": int(parts[3]),
                    "fail": int(parts[4]),
                })
        return projects
    except Exception as e:
        logging.error(f"Failed to query grid projects: {e}")
        return []


# ============================================================
# Tasks
# ============================================================

class DsGbasf2Task(Basf2PathTask):
    date = b2luigi.Parameter()
    attempt = b2luigi.Parameter()
    mode = b2luigi.Parameter()      # "kmpip", "kmpippi0", "km3pi"
    pi0 = b2luigi.Parameter(default="", significant=True)
    sample_type = b2luigi.Parameter(default="signal")  # "signal", "no_electron_id", "wrong_charge", "both"
    is_data = b2luigi.BoolParameter(default=False)
    background = b2luigi.Parameter(default="")  # MC background type (e.g., "ccbar", "taupair"), empty for data

    batch_system = "gbasf2"

    @property
    def result_dir(self) -> str:
        """
        Directory for downloaded gbasf2 ROOT files.
        CRITICAL: Use WORK_DIR (in group space) to avoid filling home directory!
        """
        return WORK_DIR

    @property
    def log_dir(self) -> str:
        """
        Base log directory:
        logs/<date>_<attempt>/<data_or_mc>/<sample_type>/<mode>[/<pi0_list>]
        b2luigi will append parameter subfolders and TaskName under this.
        """
        data_tag = "Data" if self.is_data else "MC"
        base = os.path.join(LOG_ROOT, f"{self.date}_{self.attempt}", data_tag, self.sample_type, self.mode)
        if self.mode == M_KMPIPPI0 and self.pi0:
            base = os.path.join(base, self.pi0)
        return base

    def get_log_file_dir(self) -> str:
        # Hook used by b2luigi for stdout/stderr paths (per docs)
        return self.log_dir

    @property
    def gbasf2_project_name_prefix(self) -> str:
        """Short project name (b2luigi adds unique hash)."""
        return get_project_name_prefix(self.date, self.attempt, self.mode,
                                       self.sample_type, self.is_data)

    @property
    def gbasf2_input_dataset(self) -> str:
        """Input dataset path from settings."""
        cfg = load_settings(SETTINGS_PATH)
        if self.is_data:
            return cfg.get("data_collection", "")
        else:
            # MC: use campaign and background from task parameters
            mc_campaign = cfg.get("mc_campaign", "MC15rd")
            background = self.background

            # Get campaign data
            if mc_campaign not in cfg:
                raise ValueError(f"Campaign '{mc_campaign}' not found in settings.json")

            campaign_data = cfg[mc_campaign]
            datasets = campaign_data.get("datasets", {})

            # Validate background exists for this campaign
            if background not in datasets:
                available = ", ".join(datasets.keys())
                raise ValueError(
                    f"Background '{background}' not available for {mc_campaign}. "
                    f"Available backgrounds: {available}"
                )

            return datasets[background]

    @property
    def gbasf2_release(self) -> str:
        """basf2 release version."""
        cfg = load_settings(SETTINGS_PATH)
        return cfg.get("gbasf2_release", "light-2409-toyger")

    @property
    def gbasf2_max_retries(self) -> int:
        cfg = load_settings(SETTINGS_PATH)
        return int(cfg.get("gbasf2_max_retries", 3))

    @property
    def gbasf2_download_logs(self) -> bool:
        cfg = load_settings(SETTINGS_PATH)
        return cfg.get("gbasf2_download_logs", False)

    @property
    def gbasf2_cputime(self) -> int:
        """
        CPU time in normalized minutes for gbasf2 jobs.
        Rule of thumb: cputime = 20 * <expected time on KEKCC in minutes>
        If 0, gbasf2 uses its default (which may be overestimated).
        """
        cfg = load_settings(SETTINGS_PATH)
        return int(cfg.get("gbasf2_cputime", 0))

    def _out_filename(self) -> str:
        """Output filename for this gbasf2 task."""
        if self.mode == M_KMPIPPI0:
            return f"Ds2D0e_{self.mode}_{self.pi0}.root"
        else:
            return f"Ds2D0e_{self.mode}.root"

    def output(self):
        """
        Return a Gbasf2GridProjectTarget - indicates job is complete on grid.
        Downloads are handled by DsDownloadTask, not here.
        """
        # Get the unique project name (with hash) for this task
        project_name = get_unique_project_name(self)
        return Gbasf2GridProjectTarget(project_name)

    def create_path(self):
        """Create basf2 analysis path."""
        if self.mode == M_KMPIPPI0:
            pi0_list = self.pi0
        else:
            # Steering default for non pi0 modes
            pi0_list = "eff50_May2020"

        # Get control sample flags
        no_electron_id, wrong_charge = control_sample_flags(self.sample_type)

        # Note: gbasf2 will provide the input files automatically
        # For gbasf2, use just the filename - grid jobs write to current directory
        # and gbasf2 handles transferring outputs to grid storage
        out_full = self._out_filename()

        return get_create_reconstruction_path()(
            mode=self.mode,
            infile="",  # gbasf2 fills this automatically
            outfile=out_full,
            pi0_list=pi0_list,
            truth=False,
            data=self.is_data,
            no_electron_id=no_electron_id,
            wrong_charge=wrong_charge,
        )


class DsDownloadTask(b2luigi.Task):
    """
    Download task that runs gb2_ds_get to download files from a completed grid project.

    This is a proper Luigi task, so b2luigi's scheduler manages parallelism correctly.
    Can optionally run on LSF worker nodes instead of login node.
    """
    date = b2luigi.Parameter()
    attempt = b2luigi.Parameter()
    mode = b2luigi.Parameter()
    pi0 = b2luigi.Parameter(default="", significant=True)
    sample_type = b2luigi.Parameter(default="signal")
    is_data = b2luigi.BoolParameter(default=False)
    background = b2luigi.Parameter(default="")

    # Set to "lsf" to run downloads on worker nodes instead of login node
    # batch_system = "lsf"

    @property
    def download_dir(self) -> str:
        """Directory where downloaded files will be stored."""
        # Build b2luigi-style parameter path
        params = f"date={self.date}/attempt={self.attempt}/mode={self.mode}/pi0={self.pi0}"
        params += f"/sample_type={self.sample_type}/is_data={self.is_data}"
        if not self.is_data and self.background:
            params += f"/background={self.background}"
        else:
            params += "/background="
        return os.path.join(WORK_DIR, params)

    @property
    def log_dir(self) -> str:
        data_tag = "Data" if self.is_data else "MC"
        base = os.path.join(LOG_ROOT, f"{self.date}_{self.attempt}", "download", data_tag, self.sample_type, self.mode)
        if self.mode == M_KMPIPPI0 and self.pi0:
            base = os.path.join(base, self.pi0)
        if not self.is_data and self.background:
            base = os.path.join(base, self.background)
        return base

    def get_log_file_dir(self) -> str:
        return self.log_dir

    def requires(self):
        """Require the grid job to be complete."""
        return DsGbasf2Task(
            date=self.date,
            attempt=self.attempt,
            mode=self.mode,
            pi0=self.pi0,
            sample_type=self.sample_type,
            is_data=self.is_data,
            background=self.background,
        )

    def output(self):
        """Output is a marker file indicating download is complete."""
        return b2luigi.LocalTarget(os.path.join(self.download_dir, ".download_complete"))

    def run(self):
        """Download files from grid using gb2_ds_get."""
        # Get the grid project name from the upstream task
        gbasf2_task = self.requires()
        project_name = get_unique_project_name(gbasf2_task)

        # Create directories
        os.makedirs(self.download_dir, exist_ok=True)
        os.makedirs(DOWNLOAD_LOG_DIR, exist_ok=True)

        # Check if download was already completed (marker exists)
        # If marker exists, skip. Otherwise always run gb2_ds_get to validate/complete download.
        marker_path = self.output().path
        if os.path.exists(marker_path):
            logging.info(f"Download already complete: {marker_path}")
            return

        # Create log file for this download
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        log_file = os.path.join(DOWNLOAD_LOG_DIR, f"{project_name}_{timestamp}.log")

        # Get output filename pattern
        out_filename = gbasf2_task._out_filename()
        base_name = out_filename.replace(".root", "")

        # Build gb2_ds_get command
        # Download pattern: /belle/user/<user>/<project>/sub*/<basename>_*.root
        cfg = load_settings(SETTINGS_PATH)
        grid_user = cfg.get("grid_user", "adilmub")

        download_pattern = f"/belle/user/{grid_user}/{project_name}/sub*/{base_name}_*.root"
        failed_lfns_file = os.path.join(self.download_dir, "failed_files.txt")

        cmd = [
            "gb2_ds_get",
            "--new",
            "--force",
            download_pattern,
            "--failed_lfns", failed_lfns_file,
        ]

        logging.info(f"Downloading project {project_name} to {self.download_dir}")
        logging.info(f"Log: {log_file}")

        # Write initial log entry
        with open(log_file, "w") as f:
            f.write(f"Download started: {datetime.now().isoformat()}\n")
            f.write(f"Project: {project_name}\n")
            f.write(f"Target: {self.download_dir}\n")
            f.write(f"Command: {' '.join(cmd)}\n")
            f.write(f"\n{'='*60}\n\n")

        # Change to download directory and run with gbasf2 environment
        original_dir = os.getcwd()
        downloaded_files = []
        env = get_gbasf2_env_cached()  # Get gbasf2 environment (gb2_ds_get needs this)
        try:
            os.chdir(self.download_dir)

            # Run command and capture output to log file
            with open(log_file, "a") as f:
                result = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, env=env)

            if result.returncode != 0:
                with open(log_file, "a") as f:
                    f.write(f"\n\nDownload FAILED with code {result.returncode}\n")
                logging.error(f"gb2_ds_get failed with code {result.returncode}")
                raise RuntimeError(f"Download failed for {project_name} (code {result.returncode})")

            # Count downloaded files
            downloaded_files = glob.glob(os.path.join(self.download_dir, "**", "*.root"), recursive=True)

            with open(log_file, "a") as f:
                f.write(f"\n\n{'='*60}\n")
                f.write(f"Download COMPLETED: {datetime.now().isoformat()}\n")
                f.write(f"Files downloaded: {len(downloaded_files)}\n")

            logging.info(f"Download complete for {project_name}: {len(downloaded_files)} files")

        finally:
            os.chdir(original_dir)

        # Create marker file to indicate success
        with open(self.output().path, "w") as f:
            f.write(f"Download complete: {project_name}\n")
            f.write(f"Files: {len(downloaded_files)}\n")
            f.write(f"Log: {log_file}\n")
            f.write(f"Timestamp: {datetime.now().isoformat()}\n")


class DsMergeTask(b2luigi.Task):
    date = b2luigi.Parameter()
    attempt = b2luigi.Parameter()
    mode = b2luigi.Parameter()
    pi0 = b2luigi.Parameter(default="", significant=True)
    sample_type = b2luigi.Parameter(default="signal")
    is_data = b2luigi.BoolParameter(default=False)
    background = b2luigi.Parameter(default="")  # MC background type, empty for data

    @property
    def log_dir(self) -> str:
        # logs/<date>_<attempt>/merge/<data_or_mc>/<sample_type>/<mode>
        data_tag = "Data" if self.is_data else "MC"
        return os.path.join(LOG_ROOT, f"{self.date}_{self.attempt}", "merge", data_tag, self.sample_type, self.mode)

    def get_log_file_dir(self) -> str:
        return self.log_dir

    @property
    def out_root(self) -> str:
        cfg = load_settings(SETTINGS_PATH)
        campaign = cfg.get("mc_campaign", None) if not self.is_data else None
        background = self.background if self.background else None
        return merged_output_path(self.date, self.attempt, self.mode,
                                 (self.pi0 or None), self.sample_type, self.is_data,
                                 campaign, background)

    def requires(self):
        """Require the download task to be complete."""
        return DsDownloadTask(
            date=self.date,
            attempt=self.attempt,
            mode=self.mode,
            pi0=self.pi0,
            sample_type=self.sample_type,
            is_data=self.is_data,
            background=self.background,
        )

    def output(self):
        """
        Output is either:
        1. The merged ROOT file (if merge succeeded), OR
        2. The batches folder with README.txt (if merge failed but batches are valid)
        """
        class MergeOrBatchesTarget(b2luigi.Target):
            def __init__(self, merge_path):
                self.merge_path = merge_path
                self.batch_dir = merge_path.replace(".root", "_batches")

            def exists(self):
                # Check if merged file exists and is valid
                if os.path.isfile(self.merge_path):
                    return True

                # Check if batches folder exists with README (indicates completed batch merge)
                readme_path = os.path.join(self.batch_dir, "README.txt")
                if os.path.isdir(self.batch_dir) and os.path.isfile(readme_path):
                    return True

                return False

            def makedirs(self):
                """Create parent directory for the merge output."""
                os.makedirs(os.path.dirname(self.merge_path), exist_ok=True)

            @property
            def path(self):
                # Return merged file path if it exists, otherwise batch dir
                if os.path.isfile(self.merge_path):
                    return self.merge_path
                return self.batch_dir

        return MergeOrBatchesTarget(self.out_root)

    def run(self):
        # Get the download directory from the download task
        download_task = self.requires()
        download_dir = download_task.download_dir

        logging.info(f"Looking for chunk files in: {download_dir}")

        # Find all ROOT files in the download directory (recursively, as gb2_ds_get creates subdirs)
        chunks = sorted(glob.glob(os.path.join(download_dir, "**", "*.root"), recursive=True))

        if not chunks:
            raise RuntimeError(f"No chunk files found for merge in: {download_dir}")

        os.makedirs(os.path.dirname(self.out_root), exist_ok=True)

        # If previous invalid file exists, wipe it
        if os.path.isfile(self.out_root) and not _validate_root_has_two_trees(self.out_root, self.mode):
            try:
                os.remove(self.out_root)
            except OSError:
                pass

        batch_dir = _run_hadd(self.out_root, chunks)

        if not _validate_root_has_two_trees(self.out_root, self.mode):
            # Validation failed
            if batch_dir and os.path.isdir(batch_dir):
                # Delete corrupted merge file, keep batches
                try:
                    os.remove(self.out_root)
                    logging.warning(f"Deleted corrupted merge file: {self.out_root}")
                except OSError as e:
                    logging.warning(f"Could not delete corrupted file: {e}")

                # Create a README in batches folder
                readme_path = os.path.join(batch_dir, "README.txt")
                with open(readme_path, "w") as f:
                    f.write(f"Batch files for: {os.path.basename(self.out_root)}\n")
                    f.write(f"Created: {datetime.now().isoformat()}\n")
                    f.write(f"Reason: Final merge validation failed (likely TBuffer overflow)\n")
                    f.write(f"Use these {len(os.listdir(batch_dir))-1} batch files directly for analysis.\n")

                logging.info(f"Merge validation failed for single file.")
                logging.info(f"Using batch files in: {batch_dir}")
                logging.info(f"Task marked as SUCCESS - batches are valid output.")
                return  # SUCCESS - batches are the valid output

            # No batches available - this is a real failure
            raise RuntimeError(f"Merged ROOT file failed validation and no batches available.")

        # Validation passed - clean up batch files if they exist
        if batch_dir and os.path.isdir(batch_dir):
            try:
                shutil.rmtree(batch_dir)
                logging.info(f"Cleaned up batch files: {batch_dir}")
            except Exception as e:
                logging.warning(f"Failed to clean up batch files {batch_dir}: {e}")

        # Optionally remove downloaded chunks to save space
        cfg = load_settings(SETTINGS_PATH)
        if cfg.get("remove_downloads_after_merge", False):
            try:
                if os.path.isdir(download_dir):
                    shutil.rmtree(download_dir)
                    logging.info(f"Cleaned up downloaded chunks: {download_dir}")
            except Exception as e:
                logging.warning(f"Failed to clean up {download_dir}: {e}")


class SignalTask(Basf2PathTask):
    """
    Process signal MC locally (single input file per mode).
    Runs with batch_system="local" - b2luigi workers handle parallelism.
    """
    date = b2luigi.Parameter()
    attempt = b2luigi.Parameter()
    mode = b2luigi.Parameter()      # "kmpip", "kmpippi0", "km3pi"
    pi0 = b2luigi.Parameter(default="", significant=True)
    sample_type = b2luigi.Parameter(default="signal")

    batch_system = "local"  # Run locally, not on grid or LSF
    max_event = 0  # Process all events

    @property
    def log_dir(self) -> str:
        """Log directory for signal tasks."""
        base = os.path.join(LOG_ROOT, f"{self.date}_{self.attempt}", "Signal", self.sample_type, self.mode)
        if self.mode == M_KMPIPPI0 and self.pi0:
            base = os.path.join(base, self.pi0)
        return base

    def get_log_file_dir(self) -> str:
        return self.log_dir

    def _out_fullpath(self) -> str:
        return signal_output_path(self.date, self.attempt, self.mode,
                                  (self.pi0 or None), self.sample_type)

    def output(self):
        return b2luigi.LocalTarget(self._out_fullpath())

    def create_path(self):
        out_full = self._out_fullpath()
        os.makedirs(os.path.dirname(out_full), exist_ok=True)

        # Get signal input file for this mode
        infile = get_signal_input(self.mode)

        if self.mode == M_KMPIPPI0:
            pi0_list = self.pi0
        else:
            pi0_list = "eff50_May2020"

        # Get control sample flags
        no_electron_id, wrong_charge = control_sample_flags(self.sample_type)

        return get_create_reconstruction_path()(
            mode=self.mode,
            infile=infile,
            outfile=out_full,
            pi0_list=pi0_list,
            truth=False,
            data=False,
            no_electron_id=no_electron_id,
            wrong_charge=wrong_charge,
        )


class SignalCampaign(b2luigi.WrapperTask):
    """Generate all signal MC samples based on settings.json configuration."""
    date = b2luigi.Parameter()
    attempt = b2luigi.Parameter()

    def requires(self):
        enabled_types = [st for st, enabled in ENABLED_CONTROL_SAMPLES.items() if enabled]
        if not enabled_types:
            enabled_types = ["signal"]

        for sample_type in enabled_types:
            # kmpip and km3pi modes
            yield SignalTask(date=self.date, attempt=self.attempt, mode=M_KMPIP, sample_type=sample_type)
            yield SignalTask(date=self.date, attempt=self.attempt, mode=M_KM3PI, sample_type=sample_type)
            # kmpippi0 with each pi0 list
            for pi0 in SELECTED_PI0_LISTS:
                yield SignalTask(
                    date=self.date,
                    attempt=self.attempt,
                    mode=M_KMPIPPI0,
                    pi0=pi0,
                    sample_type=sample_type,
                )


class DsCampaign(b2luigi.WrapperTask):
    date = b2luigi.Parameter()
    attempt = b2luigi.Parameter()

    def requires(self):
        # Get enabled control sample types
        enabled_types = [st for st, enabled in ENABLED_CONTROL_SAMPLES.items() if enabled]
        if not enabled_types:
            logging.warning("No control samples enabled! Defaulting to signal only.")
            enabled_types = ["signal"]

        cfg = load_settings(SETTINGS_PATH)
        backgrounds = cfg.get("backgrounds_to_run", [])

        for sample_type in enabled_types:
            # Generic MC backgrounds: ALWAYS run
            for background in backgrounds:
                yield DsMergeTask(date=self.date, attempt=self.attempt, mode=M_KMPIP,
                                sample_type=sample_type, is_data=False, background=background)
                yield DsMergeTask(date=self.date, attempt=self.attempt, mode=M_KM3PI,
                                sample_type=sample_type, is_data=False, background=background)
                for pi0 in SELECTED_PI0_LISTS:
                    yield DsMergeTask(
                        date=self.date,
                        attempt=self.attempt,
                        mode=M_KMPIPPI0,
                        pi0=pi0,
                        sample_type=sample_type,
                        is_data=False,
                        background=background,
                    )

            # Data: only run if enabled
            if RUN_ON_DATA:
                yield DsMergeTask(date=self.date, attempt=self.attempt, mode=M_KMPIP,
                                sample_type=sample_type, is_data=True, background="")
                yield DsMergeTask(date=self.date, attempt=self.attempt, mode=M_KM3PI,
                                sample_type=sample_type, is_data=True, background="")
                for pi0 in SELECTED_PI0_LISTS:
                    yield DsMergeTask(
                        date=self.date,
                        attempt=self.attempt,
                        mode=M_KMPIPPI0,
                        pi0=pi0,
                        sample_type=sample_type,
                        is_data=True,
                        background="",
                    )


class CleanupLogsTask(b2luigi.Task):
    date = b2luigi.Parameter()
    attempt = b2luigi.Parameter()

    def requires(self):
        # Always run both:
        # 1. Generic MC/Data via gbasf2
        yield DsCampaign(date=self.date, attempt=self.attempt)
        # 2. Signal MC locally (runs in parallel with gbasf2 jobs)
        yield SignalCampaign(date=self.date, attempt=self.attempt)

    def output(self):
        # Keep small marker files in home directory
        marker = os.path.join(SCRIPT_DIR, f".cleanup_{self.date}_{self.attempt}.done")
        return b2luigi.LocalTarget(marker)

    def run(self):
        # Final sweep through LOG_ROOT: keep only stderr files, remove stdout and junk
        if os.path.isdir(LOG_ROOT):
            for root, dirs, files in os.walk(LOG_ROOT, topdown=False):
                for name in files:
                    lower = name.lower()
                    if "stderr" in lower:
                        continue
                    try:
                        os.remove(os.path.join(root, name))
                    except Exception:
                        pass
                try:
                    if not os.listdir(root):
                        os.rmdir(root)
                except Exception:
                    pass

        os.makedirs(os.path.dirname(self.output().path), exist_ok=True)
        with self.output().open("w") as f:
            f.write("ok\n")

# ============================================================
# Event handlers: keep only stderr on failure
# ============================================================

@DsGbasf2Task.event_handler(b2luigi.Event.FAILURE)
def _gbasf2_on_failure(task, *args, **kwargs):
    try:
        _keep_only_stderr(task.get_log_file_dir())
    except Exception:
        pass


@DsMergeTask.event_handler(b2luigi.Event.FAILURE)
def _merge_on_failure(task, *args, **kwargs):
    try:
        _keep_only_stderr(task.get_log_file_dir())
    except Exception:
        pass


@SignalTask.event_handler(b2luigi.Event.FAILURE)
def _signal_on_failure(task, *args, **kwargs):
    try:
        _keep_only_stderr(task.get_log_file_dir())
    except Exception:
        pass


BROKEN_EVENT = getattr(b2luigi.Event, "BROKEN_TASK", None)
if BROKEN_EVENT:
    @DsGbasf2Task.event_handler(BROKEN_EVENT)
    def _gbasf2_on_broken(task, *a, **k):
        _gbasf2_on_failure(task, *a, **k)

    @DsMergeTask.event_handler(BROKEN_EVENT)
    def _merge_on_broken(task, *a, **k):
        _merge_on_failure(task, *a, **k)

    @SignalTask.event_handler(BROKEN_EVENT)
    def _signal_on_broken(task, *a, **k):
        _signal_on_failure(task, *a, **k)

# ============================================================
# Main
# ============================================================

def _ensure_dirs():
    # Ensure script directory exists (should already exist)
    os.makedirs(SCRIPT_DIR, exist_ok=True)
    # Ensure working directory exists (in group space - for downloads)
    os.makedirs(WORK_DIR, exist_ok=True)
    os.makedirs(LOG_ROOT, exist_ok=True)
    # Ensure output directories exist
    os.makedirs(MERGE_ROOT, exist_ok=True)
    os.makedirs(SIGNAL_OUTPUT_ROOT, exist_ok=True)


def _has_flag(args, name: str) -> bool:
    """
    Check whether a luigi style flag (like --scheduler-port) is present
    in args, either as '--flag' or '--flag=VALUE'.
    """
    for a in args:
        if a == name:
            return True
        if a.startswith(name + "="):
            return True
    return False


if __name__ == "__main__":
    _ensure_dirs()

    parser = argparse.ArgumentParser(description="Ds -> D0 e nu gbasf2 campaign runner")
    parser.add_argument("--config", default=SETTINGS_PATH,
                        help="Path to settings.json (default: %(default)s)")
    parser.add_argument("--workers", type=int, default=None,
                        help="Override workers (otherwise read from settings.json)")

    args, unknown = parser.parse_known_args()

    cfg = load_settings(args.config)
    date = cfg["date"]
    attempt = cfg["attempt"]
    pi0_lists = cfg.get("pi0_lists", [])
    control_samples = cfg.get("control_samples", {"signal": True, "no_electron_id": False, "wrong_charge": False, "both": False})
    run_on_data = cfg.get("run_on_data", False)
    workers_cfg = int(cfg.get("workers", 10))
    scheduler_port_cfg = int(cfg.get("scheduler_port", 8082))
    scheduler_host_cfg = cfg.get("scheduler_host", "localhost")

    SELECTED_PI0_LISTS[:] = list(pi0_lists)
    ENABLED_CONTROL_SAMPLES.clear()
    ENABLED_CONTROL_SAMPLES.update(control_samples)
    RUN_ON_DATA = run_on_data

    # Log which control samples are enabled
    enabled = [k for k, v in ENABLED_CONTROL_SAMPLES.items() if v]
    logging.info(f"Enabled control samples: {enabled}")

    if RUN_ON_DATA:
        logging.info(f"Running on: Data")
        logging.info(f"Data collection: {cfg.get('data_collection', 'Not specified')}")
    else:
        mc_campaign = cfg.get("mc_campaign", "MC15rd")
        backgrounds = cfg.get("backgrounds_to_run", [])
        int_lumi = cfg.get(mc_campaign, {}).get("int_luminosity_fb", "unknown")

        # Validate backgrounds for the selected campaign
        if mc_campaign not in cfg:
            raise ValueError(f"Campaign '{mc_campaign}' not found in settings.json")

        available_backgrounds = list(cfg[mc_campaign].get("datasets", {}).keys())

        # Validate each background
        invalid_backgrounds = [bg for bg in backgrounds if bg not in available_backgrounds]
        if invalid_backgrounds:
            raise ValueError(
                f"Invalid backgrounds for {mc_campaign}: {invalid_backgrounds}. "
                f"Available: {available_backgrounds}"
            )

        logging.info(f"Running on: MC")
        logging.info(f"Campaign: {mc_campaign} ({int_lumi} fb^-1)")
        logging.info(f"Backgrounds: {', '.join(backgrounds)} ({len(backgrounds)} total)")

    # Always run both signal (local) and generic (gbasf2)
    logging.info(f"Signal MC: local processing (runs in parallel)")

    workers = args.workers if args.workers is not None else workers_cfg

    # Build luigi argv: keep unknown args, inject scheduler flags from settings.json
    luigi_args = list(unknown)

    # Only add scheduler flags if not doing dry-run or test mode
    # (dry-run/test require local scheduler)
    is_test_mode = _has_flag(luigi_args, "--dry-run") or _has_flag(luigi_args, "--test")

    if not is_test_mode:
        if not _has_flag(luigi_args, "--scheduler-port"):
            luigi_args.append(f"--scheduler-port={scheduler_port_cfg}")

        if not _has_flag(luigi_args, "--scheduler-host"):
            luigi_args.append(f"--scheduler-host={scheduler_host_cfg}")

    # Let luigi see only its own args
    sys.argv = [sys.argv[0]] + luigi_args

    # Define the root task (depends on all other tasks)
    root_task = CleanupLogsTask(date=date, attempt=attempt)

    # Use central scheduler configured via CLI flags above
    b2luigi.process(
        root_task,
        workers=workers,
    )
