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
import luigi
from luigi import LocalTarget, Event
from b2luigi.basf2_helper import Basf2PathTask

# ============================================================
# Paths and configuration
# ============================================================

# Workflow directory (this script lives here)
WORK_DIR = "/home/belle2/amubarak/Ds2D0enue_Analysis/04-b2luigi/KEKCC"

SETTINGS_PATH = os.path.join(WORK_DIR, "settings.json")

LOG_ROOT = os.path.join(WORK_DIR, "logs")

# Temporary chunk outputs
TEMP_BASE = "/group/belle2/users2022/amubarak"

# Final merged outputs
MERGE_ROOT = "/group/belle/users/amubarak/03-KEKCC"

# Input mDSTs
MDST_ROOT = "/group/belle2/dataprod/MC/MC15ri"
SAMPLES = ["ccbar", "charged", "ddbar", "mixed", "ssbar", "uubar"]

# Modes
M_KMPIP = "kmpip"
M_KMPIPPI0 = "kmpippi0"
M_KM3PI = "km3pi"

# Global, filled after reading settings.json
SELECTED_PI0_LISTS: List[str] = []
ENABLED_CONTROL_SAMPLES: dict = {}

DEFAULT_BATCH_SYSTEM = os.environ.get("B2L_BATCH", "lsf").strip().lower()

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

def discover_mdsts(sample: str) -> List[str]:
    pattern = os.path.join(MDST_ROOT, sample, "sub*", "*.root")
    return sorted(glob.glob(pattern))


def sample_subdir(infile: str) -> str:
    return os.path.basename(os.path.dirname(infile))


def chunk_output_path(date: str, attempt: str, mode: str,
                      sample: str, subdir: str, index: int,
                      pi0: Optional[str], sample_type: str = "signal") -> str:
    """
    /group/belle2/users2022/amubarak/temp_<date>_<attempt>/
        <sample>_<mode>[_<pi0_list>][_CS_suffix]/<subdir>/ntuple_<index>.root
    """
    temp_root = os.path.join(TEMP_BASE, f"temp_{date}_{attempt}")
    base_name = f"{sample}_{mode}"
    if mode == M_KMPIPPI0 and pi0:
        base_name = f"{sample}_{mode}_{pi0}"
    # Add control sample suffix
    cs_suffix = control_sample_suffix(sample_type)
    if cs_suffix:
        base_name = base_name + cs_suffix
    out_dir = os.path.join(temp_root, base_name, subdir)
    return os.path.join(out_dir, f"ntuple_{index}.root")


def merged_output_path(date: str, attempt: str, mode: str,
                       sample: str, pi0: Optional[str], sample_type: str = "signal") -> str:
    """
    /group/belle/users/amubarak/03-KEKCC/
        Ds2D0e-Generic_Ds_<date>_<attempt>_<sample>_<mode>[_<pi0>][_CS_suffix].root
    """
    cs_suffix = control_sample_suffix(sample_type)
    if mode == M_KMPIPPI0 and pi0:
        name = f"Ds2D0e-Generic_Ds_{date}_{attempt}_{sample}_{mode}_{pi0}{cs_suffix}.root"
    else:
        name = f"Ds2D0e-Generic_Ds_{date}_{attempt}_{sample}_{mode}{cs_suffix}.root"
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

class DsChunkTask(Basf2PathTask):
    sample = b2luigi.Parameter()
    infile = b2luigi.Parameter(hashed=True)
    index = b2luigi.IntParameter()
    subdir = b2luigi.Parameter()
    date = b2luigi.Parameter()
    attempt = b2luigi.Parameter()
    mode = b2luigi.Parameter()      # "kmpip", "kmpippi0", "km3pi"
    pi0 = b2luigi.Parameter(default="", significant=True)
    sample_type = b2luigi.Parameter(default="signal")  # "signal", "no_electron_id", "wrong_charge", "both"

    batch_system = DEFAULT_BATCH_SYSTEM
    queue = "l"
    max_event = 0  # full file

    @property
    def log_dir(self) -> str:
        """
        Base log directory:
        logs/<date>_<attempt>/<sample_type>/<mode>/<sample>[/<pi0_list>]/<subdir>
        b2luigi will append parameter subfolders and TaskName under this.
        """
        base = os.path.join(LOG_ROOT, f"{self.date}_{self.attempt}", self.sample_type, self.mode, self.sample)
        if self.mode == M_KMPIPPI0 and self.pi0:
            base = os.path.join(base, self.pi0)
        return os.path.join(base, self.subdir)

    def get_log_file_dir(self) -> str:
        # Hook used by b2luigi for stdout/stderr paths (per docs)
        return self.log_dir

    def _out_fullpath(self) -> str:
        return chunk_output_path(
            self.date, self.attempt, self.mode,
            self.sample, self.subdir, self.index,
            (self.pi0 or None),
            sample_type=self.sample_type,
        )

    def output(self):
        return [LocalTarget(self._out_fullpath())]

    def create_path(self):
        out_full = self._out_fullpath()
        os.makedirs(os.path.dirname(out_full), exist_ok=True)

        if self.mode == M_KMPIPPI0:
            pi0_list = self.pi0
        else:
            # Steering default for non pi0 modes
            pi0_list = "eff50_May2020"

        # Get control sample flags
        no_electron_id, wrong_charge = control_sample_flags(self.sample_type)

        return create_reconstruction_path(
            mode=self.mode,
            infile=self.infile,
            outfile=out_full,
            pi0_list=pi0_list,
            truth=False,
            data=False,
            no_electron_id=no_electron_id,
            wrong_charge=wrong_charge,
        )


class DsMergeTask(b2luigi.Task):
    sample = b2luigi.Parameter()
    date = b2luigi.Parameter()
    attempt = b2luigi.Parameter()
    mode = b2luigi.Parameter()
    pi0 = b2luigi.Parameter(default="", significant=True)
    sample_type = b2luigi.Parameter(default="signal")

    @property
    def log_dir(self) -> str:
        # logs/<date>_<attempt>/merge/<sample_type>/<mode>/<sample>
        return os.path.join(LOG_ROOT, f"{self.date}_{self.attempt}", "merge", self.sample_type, self.mode, self.sample)

    def get_log_file_dir(self) -> str:
        return self.log_dir

    @property
    def out_root(self) -> str:
        return merged_output_path(self.date, self.attempt, self.mode, self.sample, (self.pi0 or None), sample_type=self.sample_type)

    def requires(self):
        mdsts = discover_mdsts(self.sample)
        for i, mdst in enumerate(mdsts):
            yield DsChunkTask(
                sample=self.sample,
                infile=mdst,
                index=i,
                subdir=sample_subdir(mdst),
                date=self.date,
                attempt=self.attempt,
                mode=self.mode,
                pi0=self.pi0,
                sample_type=self.sample_type,
            )

    def output(self):
        return [LocalTarget(self.out_root)]

    def run(self):
        # Find chunks:
        temp_root = os.path.join(TEMP_BASE, f"temp_{self.date}_{self.attempt}")
        base_name = f"{self.sample}_{self.mode}"
        if self.mode == M_KMPIPPI0 and self.pi0:
            base_name = f"{self.sample}_{self.mode}_{self.pi0}"
        # Add control sample suffix
        cs_suffix = control_sample_suffix(self.sample_type)
        if cs_suffix:
            base_name = base_name + cs_suffix
        base_dir = os.path.join(temp_root, base_name)

        # any subdir (sub00, sub01, ...)
        chunk_glob = os.path.join(base_dir, "*", "ntuple_*.root")
        chunks = sorted(glob.glob(chunk_glob))

        if not chunks:
            raise RuntimeError(f"No chunk files found for merge: {chunk_glob}")

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

        # Remove chunk files
        for f in chunks:
            try:
                os.remove(f)
            except OSError:
                pass

        # Clean empty dirs under temp_<date>_<attempt>/sample_mode...
        if os.path.isdir(base_dir):
            for root, dirs, files in os.walk(base_dir, topdown=False):
                if not dirs and not files:
                    try:
                        os.rmdir(root)
                    except OSError:
                        pass


class DsCampaign(b2luigi.WrapperTask):
    date = b2luigi.Parameter()
    attempt = b2luigi.Parameter()

    def requires(self):
        # Get enabled control sample types
        enabled_types = [st for st, enabled in ENABLED_CONTROL_SAMPLES.items() if enabled]
        if not enabled_types:
            logging.warning("No control samples enabled! Defaulting to signal only.")
            enabled_types = ["signal"]

        for s in SAMPLES:
            for sample_type in enabled_types:
                # No pi0 modes
                yield DsMergeTask(sample=s, date=self.date, attempt=self.attempt, mode=M_KMPIP, sample_type=sample_type)
                yield DsMergeTask(sample=s, date=self.date, attempt=self.attempt, mode=M_KM3PI, sample_type=sample_type)
                # With pi0 lists from settings.json
                for pi0 in SELECTED_PI0_LISTS:
                    yield DsMergeTask(
                        sample=s,
                        date=self.date,
                        attempt=self.attempt,
                        mode=M_KMPIPPI0,
                        pi0=pi0,
                        sample_type=sample_type,
                    )


class CleanupLogsTask(b2luigi.Task):
    date = b2luigi.Parameter()
    attempt = b2luigi.Parameter()

    def requires(self):
        return DsCampaign(date=self.date, attempt=self.attempt)

    def output(self):
        marker = os.path.join(WORK_DIR, f".cleanup_{self.date}_{self.attempt}.done")
        return [LocalTarget(marker)]

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

        # Remove full temp directory for this campaign
        temp_root = os.path.join(TEMP_BASE, f"temp_{self.date}_{self.attempt}")
        try:
            if os.path.isdir(temp_root):
                shutil.rmtree(temp_root, ignore_errors=True)
        except Exception:
            pass

        os.makedirs(os.path.dirname(self.output()[0].path), exist_ok=True)
        with self.output()[0].open("w") as f:
            f.write("ok\n")

# ============================================================
# Event handlers: keep only stderr on failure
# ============================================================

@DsChunkTask.event_handler(Event.FAILURE)
def _chunk_on_failure(task, *args, **kwargs):
    try:
        _keep_only_stderr(task.get_log_file_dir())
    except Exception:
        pass


@DsMergeTask.event_handler(Event.FAILURE)
def _merge_on_failure(task, *args, **kwargs):
    try:
        _keep_only_stderr(task.get_log_file_dir())
    except Exception:
        pass


BROKEN_EVENT = getattr(Event, "BROKEN_TASK", None)
if BROKEN_EVENT:
    @DsChunkTask.event_handler(BROKEN_EVENT)
    def _chunk_on_broken(task, *a, **k):
        _chunk_on_failure(task, *a, **k)

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
    os.makedirs(TEMP_BASE, exist_ok=True)


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

    parser = argparse.ArgumentParser(description="Ds -> D0 e nu b2luigi campaign runner")
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
    workers_cfg = int(cfg.get("workers", 50))
    scheduler_port_cfg = int(cfg.get("scheduler_port", 8082))
    scheduler_host_cfg = cfg.get("scheduler_host", "localhost")

    SELECTED_PI0_LISTS[:] = list(pi0_lists)
    ENABLED_CONTROL_SAMPLES.clear()
    ENABLED_CONTROL_SAMPLES.update(control_samples)

    # Log which control samples are enabled
    enabled = [k for k, v in ENABLED_CONTROL_SAMPLES.items() if v]
    logging.info(f"Enabled control samples: {enabled}")

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
