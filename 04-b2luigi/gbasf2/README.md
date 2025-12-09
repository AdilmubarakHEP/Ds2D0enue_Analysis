# Ds+ → D0 e+ nu gbasf2 Analysis Workflow

## Overview

This directory contains a b2luigi workflow for grid-based (gbasf2) analysis of Ds+ → D0 e+ nu decays with control sample support.

**Main script:** `00-generic_gbasf2.py`
**Configuration:** `settings.json`

## Features

- **Control Sample Support**: Signal, no electron ID, wrong charge, or combined control samples
- **Data/MC Toggle**: Process either collision data or Monte Carlo via `run_on_data` flag
- **Multiple D0 Modes**: kmpip, kmpippi0 (with pi0 list selection), km3pi
- **Clean Directory Structure**: Organized logs by date/attempt/data-or-MC/sample-type/mode
- **Central Scheduler**: Web-based monitoring via luigi scheduler
- **Automatic Merging**: Uses hadd to merge grid outputs
- **Validation**: ROOT file validation with tree existence checks

## Quick Start

### 1. Setup Environment

```bash
cd /home/belle2/amubarak/Ds2D0enue_Analysis/04-b2luigi/gbasf2
source /cvmfs/belle.cern.ch/tools/b2setup light-2409-toyger
```

### 2. Configure Your Run

Edit `settings.json`:

```json
{
  "date": "120725",           // DDMMYY format
  "attempt": "1",             // Attempt number
  "run_on_data": false,       // true = Data, false = MC

  // MC Campaign Selection (only used when run_on_data=false)
  "mc_campaign": "MC15rd",    // "MC15rd" (1444 fb^-1) or "MC16rd" (215 fb^-1)
  "backgrounds_to_run": ["ccbar", "BB", "ddbar", "ssbar", "taupair", "uubar"]  // All backgrounds to process

  "pi0_lists": ["eff20_May2020"],
  "control_samples": {
    "signal": true,           // Standard selection
    "no_electron_id": false,  // Remove electron ID (fake rate study)
    "wrong_charge": false,    // Flip e+ → e- (D*0 background)
    "both": false             // Both modifications combined
  }
}
```

### 3. Run the Workflow

**Basic run:**
```bash
python 00-generic_gbasf2.py
```

**With custom workers:**
```bash
python 00-generic_gbasf2.py --workers 20
```

**With central scheduler (web monitoring):**
```bash
# Terminal 1: Start scheduler
luigid --port 8082

# Terminal 2: Run workflow
python 00-generic_gbasf2.py --scheduler-host localhost --scheduler-port 8082
```

**Dry run (check tasks without submitting):**
```bash
python 00-generic_gbasf2.py --dry-run
```

## Configuration Details

### Control Samples

| Sample Type | `no_electron_id` | `wrong_charge` | Physics Purpose |
|------------|------------------|----------------|-----------------|
| `signal` | false | false | Standard signal region |
| `no_electron_id` | true | false | D*+ peaking background & fake rate study |
| `wrong_charge` | false | true | D*0 background control (photon conversion) |
| `both` | true | true | Combined modifications |

### Data vs MC

**For Data** (`run_on_data: true`):
- Uses `data_collection` path from settings
- Output files: `Ds2D0e-Data_<date>_<attempt>_<mode>[_CS_suffix].root`
- Project names: `DsD_<date>_<attempt>_<mode>_<cs_code>`

**For MC** (`run_on_data: false`):
- Uses campaign and background from settings (`mc_campaign` and `background_to_run`)
- Output files: `Ds2D0e-Generic_<background>_<mode>_<date>_<attempt>[_CS_suffix].root`
- Project names: `Ds_<date>_<attempt>_<mode>_<cs_code>`
- Example: `Ds2D0e-Generic_ccbar_kmpip_120725_1.root` (MC15rd ccbar)

### MC Campaigns and Backgrounds

The workflow processes **all backgrounds in the list** for the selected campaign. Edit `settings.json`:

```json
{
  "mc_campaign": "MC15rd",      // "MC15rd" or "MC16rd"
  "backgrounds_to_run": ["ccbar", "BB", "ddbar", "ssbar", "taupair", "uubar"]  // All backgrounds
}
```

**Available backgrounds per campaign:**

**MC15rd** (1444 fb^-1):
- `ccbar`, `BB`, `ddbar`, `ssbar`, `taupair`, `uubar`

**MC16rd** (215 fb^-1):
- `ccbar`, `charged`, `ddbar`, `mixed`, `ssbar`, `taupair`, `uubar`

**Examples:**

```json
// Run ALL backgrounds for MC15rd (typical usage)
"mc_campaign": "MC15rd",
"backgrounds_to_run": ["ccbar", "BB", "ddbar", "ssbar", "taupair", "uubar"]

// Run only specific backgrounds for MC15rd
"mc_campaign": "MC15rd",
"backgrounds_to_run": ["ccbar", "taupair"]

// Run ALL backgrounds for MC16rd
"mc_campaign": "MC16rd",
"backgrounds_to_run": ["ccbar", "charged", "ddbar", "mixed", "ssbar", "taupair", "uubar"]
```

The workflow will create separate output files for each background, control sample, and mode combination.

## Output Structure

### Merged Outputs
```
/group/belle/users/amubarak/03-Grid/
# MC outputs (include background name):
├── Ds2D0e-Generic_ccbar_kmpip_120725_1.root                     # Signal region
├── Ds2D0e-Generic_ccbar_kmpip_120725_1_CS_noEID.root           # No electron ID
├── Ds2D0e-Generic_ccbar_kmpip_120725_1_CS_wrongCharge.root     # Wrong charge
├── Ds2D0e-Generic_taupair_kmpippi0_eff20_May2020_120725_1.root # With pi0 list
│
# Data outputs:
├── Ds2D0e-Data_120725_1_kmpip.root                              # Signal region
├── Ds2D0e-Data_120725_1_kmpip_CS_noEID.root                     # No electron ID
└── ...
```

### Log Structure
```
logs/
└── 120725_1/                 # date_attempt
    ├── MC/                   # or Data/
    │   ├── signal/
    │   │   ├── kmpip/
    │   │   ├── kmpippi0/
    │   │   │   └── eff20_May2020/
    │   │   └── km3pi/
    │   ├── no_electron_id/
    │   ├── wrong_charge/
    │   └── both/
    └── merge/                # Merge task logs
```

## Workflow Tasks

The workflow consists of three main task types:

1. **DsGbasf2Task**: Submits jobs to grid, downloads results
2. **DsMergeTask**: Merges downloaded ROOT files with hadd
3. **DsCampaign**: Wrapper that orchestrates all enabled control samples
4. **CleanupLogsTask**: Removes stdout, keeps only stderr logs

## Advanced Usage

### Multiple Control Samples

Run signal + both control samples:
```json
"control_samples": {
  "signal": true,
  "no_electron_id": true,
  "wrong_charge": true,
  "both": false
}
```

This will create 3 output files per mode (signal, CS_noEID, CS_wrongCharge).

### Custom Configuration File

```bash
python 00-generic_gbasf2.py --config /path/to/custom_settings.json
```

### Monitoring Progress

If using central scheduler, open in browser:
```
http://localhost:8082/
```

### Checking gbasf2 Project Status

```bash
gb2_job_status
```

## Troubleshooting

**Issue:** `ModuleNotFoundError: No module named 'b2luigi'`
**Solution:** Source basf2 environment first:
```bash
source /cvmfs/belle.cern.ch/tools/b2setup light-2409-toyger
```

**Issue:** No ROOT files found for merge
**Solution:** Check gbasf2 job status with `gb2_job_status`, ensure jobs completed successfully

**Issue:** Validation failed after merge
**Solution:** Check that expected trees exist in output ROOT files. For kmpip mode expects: `DstreeCh1` and `D02kmpiptree`

**Issue:** Grid submission fails
**Solution:** Verify dataset path in settings.json is correct and accessible

## Comparison to Previous Workflow

### Old Workflow (gbasf2_analysis_task.py)
- Separate tasks per background type
- Fixed MC campaign selection
- No control sample support
- Complex nested output structure

### New Workflow (00-generic_gbasf2.py)
- Single generic task for data or MC
- Flexible collection specification
- Built-in control sample support
- Clean, flat output structure
- Matches KEKCC workflow design

## Files

- `00-generic_gbasf2.py` - Main workflow script (this replaces manual .sh scripts)
- `settings.json` - Configuration file
- `Ds2D0enue_b2luigi_steering.py` - Steering file (created by gbasf2 tasks)
- `README.md` - This file

## Contact

For issues or questions about this workflow, contact the developer or refer to b2luigi documentation:
https://b2luigi.readthedocs.io/

---

**Author:** Adil Mubarak
**Last Updated:** December 2025
