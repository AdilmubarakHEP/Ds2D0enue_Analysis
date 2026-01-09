# b2luigi Workflow for Ds+ -> D0 e+ nu Analysis

This directory contains b2luigi workflows for running the Ds+ -> D0 e+ nu reconstruction analysis on the Belle II computing infrastructure.

## Directory Structure

```
04-b2luigi/
├── gbasf2/                    # Grid (gbasf2) workflow
│   ├── 00-generic_gbasf2.py   # Main workflow script
│   ├── settings.json          # Configuration file
│   └── logs/                  # Log files from runs
├── KEKCC/                     # KEKCC batch (LSF) workflow
│   ├── 00-generic_KEKCC.py    # Main workflow script
│   ├── settings.json          # Configuration file
│   └── logs/                  # Log files from runs
└── README.md                  # This file
```

## Critical: Environment Setup

**The most important requirement** is using the correct Python environment. b2luigi and basf2 must use the **same Python installation**.

### The Problem

Belle II's basf2 framework comes with its own Python environment (currently Python 3.11). If you have:
- A conda environment
- A Python virtual environment (venv)
- Any other custom Python setup

...these will **conflict** with basf2 and cause errors like:
```
ModuleNotFoundError: No module named 'libcppyy3_9'
ImportError: Failed to import libcppyy3_9. Please check that ROOT has been built for Python 3.9
```

### The Solution

**Always use basf2's Python environment.** Never activate conda, venv, or any other Python environment when running these workflows.

## Setup Instructions

### One-Time Setup

1. **Remove/rename any conflicting environments:**
   ```bash
   # If you have a venv that might interfere
   mv ~/env ~/env_backup

   # If conda auto-activates in your .bashrc, disable it:
   conda config --set auto_activate_base false
   ```

2. **Install b2luigi in basf2's Python (one-time):**
   ```bash
   # Start clean shell
   exec bash --norc --noprofile

   # Source basf2
   source /cvmfs/belle.cern.ch/tools/b2setup release-09-00-03

   # Install b2luigi for this Python version
   pip install --user b2luigi
   ```

### Running the Workflow

**Every time** you want to run the workflow:

```bash
# 1. Start a CLEAN shell (no conda, no venv, no custom environments)
exec bash --norc --noprofile

# 2. Source basf2 (this sets up Python, ROOT, and all dependencies)
source /cvmfs/belle.cern.ch/tools/b2setup release-09-00-03

# 3. Verify you're using basf2's Python
which python
# Should show: /cvmfs/belle.cern.ch/el9/externals/.../python

# 4. Verify basf2 works
python -c "import basf2; print('basf2 OK')"

# 5. Get a grid proxy (for gbasf2 workflow)
gb2_proxy_init

# 6. Run the workflow
cd ~/Ds2D0enue_Analysis/04-b2luigi/gbasf2
python 00-generic_gbasf2.py
```

## Configuration (settings.json)

### gbasf2/settings.json

```json
{
  "date": "auto",                    // Auto-generates DDMMYY format
  "attempt": "1",                    // Attempt number for this run

  "run_on_data": false,              // true for data, false for MC
  "mc_campaign": "MC15rd",           // "MC15rd" (1444 fb^-1) or "MC16rd" (215 fb^-1)
  "backgrounds_to_run": ["taupair"], // Which MC backgrounds to process

  "pi0_lists": ["eff20_May2020"],    // Pi0 reconstruction lists for kmpippi0 mode

  "control_samples": {
    "signal": true,                  // Standard signal selection
    "no_electron_id": false,         // Fake rate control sample
    "wrong_charge": false,           // D*0 background control sample
    "both": false                    // Combined control sample
  },

  "gbasf2_release": "light-2409-toyger",  // basf2 release for grid jobs
  "gbasf2_cputime": 600,                   // CPU time (20 * KEKCC minutes)
  "gbasf2_max_retries": 3,

  "scheduler_host": "localhost",
  "scheduler_port": 9882,
  "workers": 10
}
```

### Key Settings Explained

| Setting | Description |
|---------|-------------|
| `date` | Set to `"auto"` for today's date (DDMMYY), or specify manually like `"181225"` |
| `attempt` | Increment this for re-runs with the same date |
| `mc_campaign` | MC15rd has 1444 fb^-1, MC16rd has 215 fb^-1 |
| `backgrounds_to_run` | List of backgrounds: `["ccbar", "BB", "ddbar", "ssbar", "taupair", "uubar"]` for MC15rd |
| `gbasf2_cputime` | Rule of thumb: `cputime = 20 * <expected KEKCC time in minutes>` |

## Workflow Components

### Tasks

1. **DsGbasf2Task**: Submits reconstruction jobs to the grid
2. **DsMergeTask**: Merges downloaded ROOT files using `hadd`
3. **SignalTask**: Processes signal MC locally (runs in parallel with grid jobs)
4. **DsCampaign**: Wrapper that runs all gbasf2 tasks for selected backgrounds
5. **SignalCampaign**: Wrapper that runs all signal MC tasks
6. **CleanupLogsTask**: Final task that cleans up log files after completion

### Output Locations

- **Grid MC outputs**: `/group/belle/users/amubarak/02-Grid/`
- **Signal MC outputs**: `/home/belle2/amubarak/C01-Simulated_Events/Signal/`
- **Signal MC inputs**: `/group/belle/users/amubarak/00-Generation/Signal/`

### Output File Naming

Grid MC:
```
Ds2D0e-Generic_{background}_{mode}_{date}_{attempt}.root
# Example: Ds2D0e-Generic_taupair_kmpip_181225_1.root
```

Signal MC:
```
Ds2D0e-Signal_{date}_{attempt}_{mode}.root
# Example: Ds2D0e-Signal_181225_1_kmpip.root
```

## Command Line Options

```bash
# Normal run (connects to central scheduler)
python 00-generic_gbasf2.py

# Dry run (shows what would be executed, uses local scheduler)
python 00-generic_gbasf2.py --dry-run

# Override number of workers
python 00-generic_gbasf2.py --workers 5

# Use custom config file
python 00-generic_gbasf2.py --config /path/to/custom_settings.json
```

## Monitoring

### Central Scheduler Web Interface

When running without `--dry-run`, the workflow connects to a central scheduler that provides a web interface:

```bash
# Start the scheduler (in a separate terminal, also with basf2 sourced)
luigid --port 9882

# Access the web interface
# http://localhost:9882
```

### Checking Grid Job Status

```bash
# Check all your gbasf2 projects
gb2_job_status

# Check specific project
gb2_job_status -p Ds_181225_1_kmpip_xxx
```

## Troubleshooting

### "No module named 'libcppyy3_9'" or similar cppyy errors

**Cause**: Python environment mismatch. A conda/venv environment is interfering with basf2.

**Solution**:
```bash
# Remove/rename conflicting environments
mv ~/env ~/env_backup
conda deactivate  # if conda is active

# Start fresh
exec bash --norc --noprofile
source /cvmfs/belle.cern.ch/tools/b2setup release-09-00-03
```

### "Can not find basf2. Can not use the basf2 task."

**Cause**: basf2 was not sourced before running the workflow.

**Solution**: Always source basf2 before running:
```bash
source /cvmfs/belle.cern.ch/tools/b2setup release-09-00-03
```

### "Can not test while using a central scheduler!"

**Cause**: Using `--dry-run` with central scheduler flags.

**Solution**: The script automatically handles this. If you still see this error, ensure you're using the latest version of `00-generic_gbasf2.py`.

### Worker processes fail but main process works

**Cause**: b2luigi spawns subprocess workers that might pick up a different Python environment.

**Solution**: Completely remove any venv/conda that might be in your PATH:
```bash
mv ~/env ~/env_backup
exec bash --norc --noprofile
source /cvmfs/belle.cern.ch/tools/b2setup release-09-00-03
```

### Grid jobs submitted but no output

**Cause**: Jobs may still be running, or there was an error.

**Solution**: Check job status:
```bash
gb2_job_status -p <project_name>
```

## Quick Reference

### Complete Run Sequence

```bash
# 1. Clean environment
exec bash --norc --noprofile

# 2. Setup basf2
source /cvmfs/belle.cern.ch/tools/b2setup release-09-00-03

# 3. Get grid proxy
gb2_proxy_init

# 4. Navigate to workflow directory
cd ~/Ds2D0enue_Analysis/04-b2luigi/gbasf2

# 5. Edit settings if needed
vim settings.json

# 6. Test with dry-run first
python 00-generic_gbasf2.py --dry-run

# 7. Run for real
python 00-generic_gbasf2.py
```

### Useful Commands

```bash
# Check basf2 version
basf2 --version

# Check grid proxy validity
gb2_proxy_info

# List grid projects
gb2_job_status

# Kill a grid project
gb2_job_kill -p <project_name>

# Download logs from grid job
gb2_job_output -p <project_name> --logs
```

## Support

- Belle II questions: https://questions.belle2.org
- b2luigi documentation: https://b2luigi.readthedocs.io
- basf2 documentation: https://software.belle2.org
