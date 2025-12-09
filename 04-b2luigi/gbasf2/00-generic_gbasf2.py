#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
import json
import shutil
import subprocess
import argparse
from typing import List, Optional, Tuple
import logging

logging.basicConfig(level=logging.INFO)

import b2luigi
from b2luigi.basf2_helper import Basf2PathTask

# ============================================================
# Paths and configuration
# ============================================================

# Workflow directory (this script lives here)
WORK_DIR = "/home/belle2/amubarak/Ds2D0enue_Analysis/04-b2luigi/gbasf2"

SETTINGS_PATH = os.path.join(WORK_DIR, "settings.json")

LOG_ROOT = os.path.join(WORK_DIR, "logs")

# Final merged outputs
MERGE_ROOT = "/group/belle/users/amubarak/03-Grid"

# Modes
M_KMPIP = "kmpip"
M_KMPIPPI0 = "kmpippi0"
M_KM3PI = "km3pi"

# Global, filled after reading settings.json
SELECTED_PI0_LISTS: List[str] = []
ENABLED_CONTROL_SAMPLES: dict = {}
RUN_ON_DATA: bool = False

# ============================================================
# Import steering: create_reconstruction_path
# ============================================================

# Steering script path
STEERING_FILE = "/home/belle2/amubarak/Ds2D0enue_Analysis/02-Reconstruction_Scripts/Ds2D0e-Reconstruction.py"

import importlib.util

_spec = importlib.util.spec_from_file_location("ds_reco", STEERING_FILE)
if _spec is None or _spec.loader is None:
    raise RuntimeError(f"Cannot import steering from {STEERING_FILE}")
_ds_reco = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_ds_reco)

create_reconstruction_path = _ds_reco.create_reconstruction_path

# ============================================================
# Config helpers
# ============================================================

def load_settings(path: str) -> dict:
    with open(path, "r") as f:
        cfg = json.load(f)
    if "date" not in cfg or "attempt" not in cfg:
        raise ValueError("settings.json must contain 'date' and 'attempt'")
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
    /group/belle/users/amubarak/03-Grid/
        Ds2D0e-Data_120725_1_kmpip[_CS_suffix].root  (if is_data=True)
        Ds2D0e-Generic_ccbar_kmpip_120725_1[_CS_suffix].root    (if is_data=False)
    """
    cs_suffix = control_sample_suffix(sample_type)

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

    return os.path.join(MERGE_ROOT, name)


def _which_hadd() -> str:
    exe = shutil.which("hadd")
    if not exe:
        raise RuntimeError("Cannot find 'hadd' in PATH. Source basf2/ROOT first.")
    return exe


def _run_hadd(out_file: str, inputs: List[str]) -> None:
    cmd = [_which_hadd(), "-f", out_file] + inputs
    subprocess.check_call(cmd)


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

    def _out_filename(self) -> str:
        """Output filename for this gbasf2 task."""
        if self.mode == M_KMPIPPI0:
            return f"Ds2D0e_{self.mode}_{self.pi0}.root"
        else:
            return f"Ds2D0e_{self.mode}.root"

    def output(self):
        return self.add_to_output(self._out_filename())

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
        # We create a dummy output filename here; gbasf2 handles actual output naming
        out_full = self.get_output_file_name(self._out_filename())

        return create_reconstruction_path(
            mode=self.mode,
            infile="",  # gbasf2 fills this automatically
            outfile=out_full,
            pi0_list=pi0_list,
            truth=False,
            data=self.is_data,
            no_electron_id=no_electron_id,
            wrong_charge=wrong_charge,
        )


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
        return b2luigi.LocalTarget(self.out_root)

    def run(self):
        # Get the gbasf2 task output
        gbasf2_task = self.requires()
        gbasf2_output_target = gbasf2_task.output()

        # The output from gbasf2 is a directory containing multiple ROOT files
        # Get the path to this directory
        if hasattr(gbasf2_output_target, 'path'):
            output_dir = gbasf2_output_target.path
        else:
            # Handle if it's wrapped differently
            output_dir = str(gbasf2_output_target)

        # Find all ROOT files in the output directory
        if os.path.isdir(output_dir):
            chunk_glob = os.path.join(output_dir, "*.root")
            chunks = sorted(glob.glob(chunk_glob))
        else:
            # If it's a single file (shouldn't happen with gbasf2 but handle anyway)
            chunks = [output_dir] if os.path.isfile(output_dir) else []

        if not chunks:
            raise RuntimeError(f"No chunk files found for merge in: {output_dir}")

        os.makedirs(os.path.dirname(self.out_root), exist_ok=True)

        # If previous invalid file exists, wipe it
        if os.path.isfile(self.out_root) and not _validate_root_has_two_trees(self.out_root, self.mode):
            try:
                os.remove(self.out_root)
            except OSError:
                pass

        _run_hadd(self.out_root, chunks)

        if not _validate_root_has_two_trees(self.out_root, self.mode):
            raise RuntimeError(f"Merged ROOT file failed validation: {self.out_root}")

        # Optionally remove downloaded chunks to save space
        cfg = load_settings(SETTINGS_PATH)
        if cfg.get("remove_downloads_after_merge", False):
            try:
                if os.path.isdir(output_dir):
                    shutil.rmtree(output_dir)
                    logging.info(f"Cleaned up downloaded chunks: {output_dir}")
            except Exception as e:
                logging.warning(f"Failed to clean up {output_dir}: {e}")


class DsCampaign(b2luigi.WrapperTask):
    date = b2luigi.Parameter()
    attempt = b2luigi.Parameter()

    def requires(self):
        # Get enabled control sample types
        enabled_types = [st for st, enabled in ENABLED_CONTROL_SAMPLES.items() if enabled]
        if not enabled_types:
            logging.warning("No control samples enabled! Defaulting to signal only.")
            enabled_types = ["signal"]

        # Determine if we're running on data or MC
        for sample_type in enabled_types:
            if RUN_ON_DATA:
                # Data: run all modes for this control sample type
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
            else:
                # MC: run all modes for each background in the list
                cfg = load_settings(SETTINGS_PATH)
                backgrounds = cfg.get("backgrounds_to_run", [])

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


class CleanupLogsTask(b2luigi.Task):
    date = b2luigi.Parameter()
    attempt = b2luigi.Parameter()

    def requires(self):
        return DsCampaign(date=self.date, attempt=self.attempt)

    def output(self):
        marker = os.path.join(WORK_DIR, f".cleanup_{self.date}_{self.attempt}.done")
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


BROKEN_EVENT = getattr(b2luigi.Event, "BROKEN_TASK", None)
if BROKEN_EVENT:
    @DsGbasf2Task.event_handler(BROKEN_EVENT)
    def _gbasf2_on_broken(task, *a, **k):
        _gbasf2_on_failure(task, *a, **k)

    @DsMergeTask.event_handler(BROKEN_EVENT)
    def _merge_on_broken(task, *a, **k):
        _merge_on_failure(task, *a, **k)

# ============================================================
# Main
# ============================================================

def _ensure_dirs():
    os.makedirs(WORK_DIR, exist_ok=True)
    os.makedirs(LOG_ROOT, exist_ok=True)
    os.makedirs(MERGE_ROOT, exist_ok=True)


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

    workers = args.workers if args.workers is not None else workers_cfg

    root_task = CleanupLogsTask(date=date, attempt=attempt)

    # Build luigi argv: keep unknown args, inject scheduler flags from settings.json
    luigi_args = list(unknown)

    if not _has_flag(luigi_args, "--scheduler-port"):
        luigi_args.append(f"--scheduler-port={scheduler_port_cfg}")

    if not _has_flag(luigi_args, "--scheduler-host"):
        luigi_args.append(f"--scheduler-host={scheduler_host_cfg}")

    # Let luigi see only its own args
    sys.argv = [sys.argv[0]] + luigi_args

    # Use central scheduler configured via CLI flags above
    b2luigi.process(
        root_task,
        workers=workers,
    )
