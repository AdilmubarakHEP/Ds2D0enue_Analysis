#!/usr/bin/env python3
"""
b2luigi workflow for Ds+ -> D0 e+ nu reconstruction with gbasf2.

Usage:
    python gbasf2_analysis_task.py --workers 4
    python gbasf2_analysis_task.py --dry-run
    
Author: Adil Mubarak
Date: December 2025
"""

import b2luigi
from b2luigi.basf2_helper.tasks import Basf2PathTask
import json
from pathlib import Path
from datetime import datetime


#==============================================================================
# LOAD CONFIGURATION
#==============================================================================

SCRIPT_DIR = Path("/home/belle2/amubarak/Ds2D0enue_Analysis/04-b2luigi/gbasf2")
SETTINGS_FILE = SCRIPT_DIR / "settings.json"

if not SETTINGS_FILE.exists():
    raise FileNotFoundError(f"Settings file not found: {SETTINGS_FILE}")

with open(SETTINGS_FILE, "r") as f:
    settings = json.load(f)

# Get active MC campaign
MC_CAMPAIGN = settings["mc_campaign"]
CAMPAIGN_DATA = settings[MC_CAMPAIGN]
DATASETS = CAMPAIGN_DATA["datasets"]
INT_LUMINOSITY = CAMPAIGN_DATA["int_luminosity_fb"]

# Project naming
DATE = settings["date"]
ATTEMPT = settings["attempt"]
YEAR = datetime.now().strftime("%y")

# Settings
MODES = settings["modes_to_run"]
BACKGROUNDS = settings["backgrounds_to_run"]
BASF2_RELEASE = settings["basf2_release"]
PI0_LIST = settings["pi0_list"]

# Background code mapping (shortened for gbasf2 project name)
BACKGROUND_CODES = {
    "taupair": "tau",
    "ccbar": "cc",
    "ddbar": "dd",
    "ssbar": "ss",
    "uubar": "uu",
    "BB": "BB",
    "charged": "chg",
    "mixed": "mix"
}

# Validate backgrounds
INVALID_BACKGROUNDS = [bg for bg in BACKGROUNDS if bg not in DATASETS]
if INVALID_BACKGROUNDS:
    raise ValueError(f"Invalid backgrounds for {MC_CAMPAIGN}: {INVALID_BACKGROUNDS}")

# Directories
MERGE_OUTPUT_DIR = Path(settings["merge_output_dir"])
TERMINAL_LOGS_DIR = Path(settings["terminal_logs_dir"])
REMOVE_DOWNLOADS = settings.get("remove_downloads_after_merge", False)

# Central scheduler
USE_CENTRAL_SCHEDULER = settings.get("use_central_scheduler", False)
CENTRAL_SCHEDULER_HOST = settings.get("central_scheduler_host", "localhost")
CENTRAL_SCHEDULER_PORT = settings.get("central_scheduler_port", 8082)

# Path to steering script
STEERING_SCRIPT = SCRIPT_DIR / "Ds2D0enue_b2luigi_steering.py"


#==============================================================================
# HELPER FUNCTIONS
#==============================================================================

def get_project_name(background, mode):
    """
    Generate short gbasf2 project name (b2luigi adds unique hash automatically).
    Example: Ds_tau_kmpippi0 (18 chars) + b2luigi hash (10 chars) = 28 chars
    """
    bg_code = BACKGROUND_CODES.get(background, background[:3])
    template = settings.get("gbasf2_project_name_template", "Ds_{background}_{mode}")
    return template.format(background=bg_code, mode=mode)


def get_output_filename(background, mode):
    """
    Generate merged output filename with full date and attempt for tracking.
    Example: Ds2D0e-Generic_taupair_kmpippi0_251201_1.root
    """
    return f"Ds2D0e-Generic_{background}_{mode}_{DATE}_{ATTEMPT}.root"


#==============================================================================
# B2LUIGI TASKS
#==============================================================================

class AnalysisTask(Basf2PathTask):
    """Grid analysis task using b2luigi's Basf2PathTask."""
    
    batch_system = "gbasf2"
    
    background = b2luigi.Parameter()
    mode = b2luigi.Parameter()
    
    @property
    def gbasf2_project_name_prefix(self):
        """Short project name (b2luigi adds unique hash)."""
        return get_project_name(self.background, self.mode)
    
    @property
    def gbasf2_input_dataset(self):
        """Input dataset path."""
        return DATASETS[self.background]
    
    @property
    def gbasf2_release(self):
        """basf2 release version."""
        return BASF2_RELEASE
    
    @property
    def gbasf2_max_retries(self):
        return settings.get("gbasf2_max_retries", 3)
    
    @property
    def gbasf2_download_logs(self):
        return settings.get("gbasf2_download_logs", False)
    
    def output(self):
        """Output ROOT file."""
        yield self.add_to_output(f"output_{self.mode}.root")
    
    def create_path(self):
        """Create basf2 analysis path."""
        # Import inside method (when basf2 environment is available)
        import sys
        sys.path.insert(0, str(SCRIPT_DIR))
        from Ds2D0enue_b2luigi_steering import create_analysis_path
        
        return create_analysis_path(
            mode=self.mode,
            pi0_list=PI0_LIST,
            output_filename=self.get_output_file_name(f"output_{self.mode}.root")
        )


class HaddMergeTask(b2luigi.Task):
    """Merge grid outputs with hadd."""
    
    background = b2luigi.Parameter()
    mode = b2luigi.Parameter()
    
    def requires(self):
        return AnalysisTask(background=self.background, mode=self.mode)
    
    def output(self):
        output_file = MERGE_OUTPUT_DIR / get_output_filename(self.background, self.mode)
        return b2luigi.LocalTarget(str(output_file))
    
    def run(self):
        import subprocess
        
        project_name = get_project_name(self.background, self.mode)
        output_filename = get_output_filename(self.background, self.mode)
        
        # Get input directory from AnalysisTask
        input_task = AnalysisTask(background=self.background, mode=self.mode)
        input_file_path = Path(input_task.get_output_file_name(f"output_{self.mode}.root"))
        input_dir = input_file_path.parent
        
        # Find all ROOT files
        root_files = sorted(input_dir.rglob("*.root"))
        
        if not root_files:
            raise FileNotFoundError(f"No ROOT files found in {input_dir}")
        
        MERGE_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        output_path = self.output().path
        
        print(f"\n{'='*70}")
        print(f"[{project_name}] MERGING OUTPUT")
        print(f"{'='*70}")
        print(f"  Campaign:    {MC_CAMPAIGN}")
        print(f"  Background:  {self.background}")
        print(f"  Mode:        {self.mode}")
        print(f"  Date:        {DATE}")
        print(f"  Attempt:     {ATTEMPT}")
        print(f"  Input files: {len(root_files)}")
        print(f"  Output:      {output_filename}")
        print(f"{'='*70}\n")
        
        cmd = ["hadd", "-f", str(output_path)] + [str(f) for f in root_files]
        
        TERMINAL_LOGS_DIR.mkdir(parents=True, exist_ok=True)
        log_file = TERMINAL_LOGS_DIR / f"{project_name}_merge.log"
        
        with open(log_file, "w") as f:
            result = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, text=True)
        
        if result.returncode != 0:
            raise RuntimeError(f"hadd failed. Check: {log_file}")
        
        file_size_mb = Path(output_path).stat().st_size / (1024 * 1024)
        print(f"[{project_name}] Merge complete ({file_size_mb:.1f} MB)")
        print(f"[{project_name}] Saved: {output_path}\n")
        
        if REMOVE_DOWNLOADS:
            import shutil
            try:
                shutil.rmtree(input_dir)
                print(f"[{project_name}] Cleaned up downloads")
            except Exception as e:
                print(f"[{project_name}] WARNING: Cleanup failed: {e}")


class WrapperTask(b2luigi.WrapperTask):
    """Top-level workflow orchestrator."""
    
    def requires(self):
        tasks = []
        for bg in BACKGROUNDS:
            for mode in MODES:
                tasks.append(HaddMergeTask(background=bg, mode=mode))
        
        print(f"\n{'='*70}")
        print(f"Ds+ -> D0 e+ nu GRID ANALYSIS WORKFLOW")
        print(f"{'='*70}")
        print(f"MC Campaign:        {MC_CAMPAIGN}")
        print(f"Integrated Lumi:    {INT_LUMINOSITY:.2f} fb^-1")
        print(f"Date:               {DATE}")
        print(f"Attempt:            {ATTEMPT}")
        print(f"Backgrounds:        {', '.join(BACKGROUNDS)}")
        print(f"D0 Modes:           {', '.join(MODES)}")
        print(f"Total tasks:        {len(tasks)} grid jobs")
        print(f"")
        print(f"gbasf2 project names (b2luigi adds unique hash):")
        for bg in BACKGROUNDS:
            for mode in MODES:
                pname = get_project_name(bg, mode)
                print(f"  {pname}")
        print(f"")
        print(f"Output files (with full date/attempt):")
        for bg in BACKGROUNDS:
            for mode in MODES:
                fname = get_output_filename(bg, mode)
                print(f"  {fname}")
        print(f"{'='*70}\n")
        
        return tasks


#==============================================================================
# MAIN
#==============================================================================

if __name__ == "__main__":
    import sys
    
    # Parse --workers manually
    workers = 4
    filtered_args = [sys.argv[0]]
    i = 1
    while i < len(sys.argv):
        if sys.argv[i] == '--workers' and i + 1 < len(sys.argv):
            workers = int(sys.argv[i + 1])
            i += 2
        else:
            filtered_args.append(sys.argv[i])
            i += 1
    sys.argv = filtered_args
    
    # Configure b2luigi
    b2luigi.set_setting("result_dir", "/home/belle2/amubarak/C02-Grid/b2luigi")
    b2luigi.set_setting("log_file_dir", "/home/belle2/amubarak/C02-Grid/b2luigi_logs")
    b2luigi.set_setting("workers", workers)
    
    # Force gbasf2 submission without prompts (for automation)
    b2luigi.set_setting("gbasf2_additional_params", "--force")
    
    if USE_CENTRAL_SCHEDULER:
        b2luigi.set_setting("central_scheduler_host", CENTRAL_SCHEDULER_HOST)
        b2luigi.set_setting("central_scheduler_port", CENTRAL_SCHEDULER_PORT)
    
    print(f"Using {workers} workers\n")
    
    # Run workflow
    try:
        b2luigi.process(WrapperTask())
        print(f"\n✅ WORKFLOW COMPLETE\nOutput: {MERGE_OUTPUT_DIR}\n")
    except Exception as e:
        print(f"\n❌ WORKFLOW FAILED: {e}\n")
        raise