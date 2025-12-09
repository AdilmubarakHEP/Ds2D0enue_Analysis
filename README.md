# Ds → D0 e ν Analysis

**Belle II Rare Semileptonic Decay Analysis**

## Physics Process

Study of the rare charm semileptonic decay:

**D<sub>s</sub><sup>+</sup> → D<sup>0</sup> e<sup>+</sup> ν<sub>e</sub>**

with subsequent D<sup>0</sup> decays to:
- **Mode 1 (kmpip)**: D<sup>0</sup> → K<sup>-</sup> π<sup>+</sup>
- **Mode 2 (kmpippi0)**: D<sup>0</sup> → K<sup>-</sup> π<sup>+</sup> π<sup>0</sup>
- **Mode 3 (km3pi)**: D<sup>0</sup> → K<sup>-</sup> π<sup>+</sup> π<sup>+</sup> π<sup>-</sup>

## Analysis Workflow

1. **Signal Generation** → Generate signal MC events using EvtGen
2. **Reconstruction** → Process signal and generic MC samples through basf2 reconstruction
3. **BDT Training** → Train machine learning classifiers for background suppression
4. **Fitting** → Extract signal yields and perform systematic studies
5. **Results** → Publication plots and tables

## Directory Structure

### Core Analysis Directories

#### `00-Setup/`
Environment setup scripts and configuration files for the Belle II software framework.

#### `01-Event_Generation/`
Signal MC generation using EvtGen:
- `.dec` files: Decay descriptors for signal and peaking backgrounds
- Steering files for event generation
- Output: Signal MC samples stored in `C00-Generation/`

#### `02-Reconstruction_Scripts/`
basf2 steering files for event reconstruction:
- **`Ds2D0e-Reconstruction.py`**: Main reconstruction script supporting all three D<sup>0</sup> decay modes
- Background simulation scripts for peaking backgrounds
- Control sample generation (no electron ID, wrong charge)
- Output variables: Kinematic, vertex, PID, and MC truth information

#### `03-Grid/`
Scripts for large-scale production on Belle II distributed computing grid.

#### `04-b2luigi/`
**Batch processing workflow using b2luigi (Luigi-based workflow manager)**:

**`KEKCC/`**: Production scripts for KEK Computing Center
- **`00-generic_KEKCC.py`**: Main b2luigi workflow script
  - Processes generic MC samples (ccbar, charged, ddbar, mixed, ssbar, uubar)
  - Optionally generates signal samples
  - Supports multiple π<sup>0</sup> lists (eff20_May2020, eff30_May2020, eff40_May2020, eff50_May2020, eff60_May2020)
  - Control sample variations (signal, no_electron_id, wrong_charge, both)
  - Automatic chunking, parallel processing, and merging
  - Configurable via `settings.json`
- **`settings.json`**: Campaign configuration
  - `date`: Campaign identifier (DDMMYY)
  - `attempt`: Version number
  - `pi0_lists`: Which π<sup>0</sup> reconstruction lists to use
  - `control_samples`: Which control samples to generate
  - `generate_signal`: Toggle signal MC processing
  - `workers`, `scheduler_port`, `scheduler_host`: Batch system configuration
- **`logs/`**: Batch job logs organized by date/attempt/sample/mode
- Output locations:
  - Generic MC: `/group/belle/users/amubarak/03-KEKCC/`
  - Signal: `/home/belle2/amubarak/C01-Simulated_Events/Signal/`

**Usage**:
```bash
# Dry-run to preview tasks
python3 00-generic_KEKCC.py --dry-run

# Run full campaign
python3 00-generic_KEKCC.py

# Override number of workers
python3 00-generic_KEKCC.py --workers 100
```

#### `05-ML/`
Machine learning for background suppression:

**basf2/**: BDT training within basf2 framework
- Fake D<sup>0</sup> suppression BDT
- General background suppression BDT

**Local/**: External ML training (XGBoost, scikit-learn)
- Feature engineering and selection
- Model optimization
- Performance evaluation

#### `06-ROOT/`
ROOT-based analysis:
- Invariant mass fitting
- Signal extraction
- Systematic uncertainty evaluation
- Efficiency calculations

#### `07-Python/`
Python analysis scripts:

**Machine_Learning/**:
- Variable correlation studies
- Feature discrimination analysis
- ROC curves and performance metrics

**Plots_and_Tables/**:
- Jupyter notebooks for publication plots
- Tables for paper and presentations
- Kinematic distributions
- Efficiency studies

#### `08-Python_Functions/`
Custom Python utilities and helper functions to avoid code duplication across analysis scripts.

#### `09-Images/`
Generated plots, figures, and diagrams for presentations and publications.

#### `10-Sources/`
References, papers, and documentation relevant to the analysis.

#### `11-Save/`
Training exercises and tutorials for Belle II software framework.

#### `Unused/`
Archived scripts that may be useful for future reference.

---

## Key Data Locations

### Input Data
- **Generic MC**: `/group/belle2/dataprod/MC/MC15ri/{ccbar,charged,ddbar,mixed,ssbar,uubar}/`
- **Generated Signal**: `/home/belle2/amubarak/C00-Generation/`

### Output Data
- **Reconstructed Generic MC**: `/group/belle/users/amubarak/03-KEKCC/`
- **Reconstructed Signal**: `/home/belle2/amubarak/C01-Simulated_Events/Signal/`
- **Temporary Processing**: `/group/belle2/users2022/amubarak/temp_*/`

---

## Control Samples

The analysis includes control samples to validate electron identification and charge requirements:

1. **signal**: Standard signal selection
2. **no_electron_id**: Signal selection without electron ID requirement (for PID studies)
3. **wrong_charge**: Wrong-sign D<sub>s</sub> candidates (D<sup>0</sup> e<sup>-</sup> combinations for D*<sup>0</sup> background)
4. **both**: Combined no electron ID + wrong charge

---

## π<sup>0</sup> Lists

Multiple π<sup>0</sup> reconstruction configurations tested for Mode 2 (kmpippi0):
- `eff20_May2020`, `eff30_May2020`, `eff40_May2020`, `eff50_May2020`, `eff60_May2020`

These represent different efficiency/purity trade-offs in π<sup>0</sup> reconstruction.

---

## Software Requirements

- **basf2**: Belle II Analysis Software Framework (release-09-00-03 or later)
- **b2luigi**: Workflow management (for batch processing)
- **ROOT**: For data analysis and fitting
- **Python 3.x** with:
  - numpy, pandas, matplotlib, seaborn
  - uproot (ROOT file I/O)
  - scikit-learn, XGBoost (ML)
  - Jupyter (notebooks)

---

*Last updated: December 2025*
