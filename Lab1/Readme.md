# GNSS Single Point Positioning (SPP) Analysis Pipeline

A comprehensive toolkit for processing GNSS UBX data files and performing Single Point Positioning (SPP) analysis with weighted least squares estimation.

---

## Table of Contents

1. [Overview](#overview)  
2. [Lab Materials (2026/27)](#lab-materials-202627)  
3. [Datasets](#datasets)  
4. [Lab Requirements](#lab-requirements)  
5. [Prerequisites](#prerequisites)  
6. [Installation](#installation)  
7. [Workflow Overview](#workflow-overview)
8. [Step 0: GNSS Data Collection](#step-0-gnss-data-collection)
9. [Step 1: UBX to RINEX Conversion](#step-1-ubx-to-rinex-conversion)
10. [Step 2: RINEX Data Processing](#step-2-rinex-data-processing)  
11. [Step 3: UBX to CSV Conversion](#step-3-ubx-to-csv-conversion)  
12. [Step 4: SPP Analysis](#step-4-spp-analysis)  
13. [Step 5: Accuracy Challenge and Leaderboard](#step-5-accuracy-challenge-and-leaderboard)  
14. [Results and Visualization](#results-and-visualization)  
15. [Output Files](#output-files)  
16. [Troubleshooting](#troubleshooting)  
17. [Technical Notes](#technical-notes)  
18. [Contributing](#contributing)  
19. [License](#license)  
20. [References](#references)  
21. [Contact](#contact)  

---

## Overview

This pipeline provides a complete workflow for GNSS data processing, from raw UBX file collection to advanced positioning analysis with error evaluation. The workflow includes data collection, data conversion, visualization, and statistical analysis of positioning accuracy.

**New in 2026/27:** (i) an urban benchmark dataset from UrbanNav, (ii) a CSV submission ranked on a class leaderboard by horizontal RMSE, and (iii) explicit permission to use generative AI tools for coding, subject to the rules in the tutorial (Section 3.3).

---

## Lab Materials (2026/27)

| File | Description |
|------|-------------|
| [AAE4203_26_27_S1_Lab_Tutorial.pdf](AAE4203_26_27_S1_Lab_Tutorial.pdf) | Full lab tutorial: theory, worked examples, data collection, processing, challenge rules (the group list is on Blackboard) |
| [AAE4203_Week9_GNSS_Lab_Slides.pdf](AAE4203_Week9_GNSS_Lab_Slides.pdf) | Week 9 lab session slides |
| [PolyU_AAE4203_Lab_report_template.docx](PolyU_AAE4203_Lab_report_template.docx) | Lab report template (Word) |
| [tools/check_submission.py](tools/check_submission.py) | Checks the format of your leaderboard CSV before upload |
| [rinex2csv.m](rinex2csv.m), [matRTKLIB/](matRTKLIB/) | RINEX -> CSV conversion (MATLAB) |
| [pyubx2_csv_converter_gui.py](pyubx2_csv_converter_gui.py) | UBX -> CSV conversion (Python GUI) |
| [COM6___9600_250916_020025_blockz.ubx](COM6___9600_250916_020025_blockz.ubx) | Example raw dataset (Block X, 2025) |

---

## Datasets

| Dataset | Content | Ground truth | Purpose |
|---------|---------|--------------|---------|
| Own field data | UBX log recorded by your group (15–20 min, static) near Block X | none | Data collection, own-data validation and discussion in the report |
| Example (Block X, 2025) | UBX log, RINEX .obs/.nav, CSV files from `rinex2csv.m`, receiver solution (NAV-HPPOSECEF), example Python code | reference results provided | Warm-up and regression test: check that your code reproduces the example results |
| UrbanNav benchmark | `UrbanNav-HK-Deep-Urban-1.ublox.f9p.obs` and `UrbanNav_20210521_multiGNSS_clean.nav` | kept by the teaching team | Accuracy challenge and leaderboard |

The example dataset, the example SPP code (`WLSE_spp.py`) and the benchmark files are linked on Blackboard.

---

## Lab Requirements

Each group submits the following deliverables:

| Deliverable | Content | Where |
|-------------|---------|-------|
| Leaderboard entry | `GXX_benchmark.csv`: your trajectory for the benchmark dataset | TA's Google Drive folder (link on Blackboard) |
| Own-data solution | `GXX_owndata.csv`: your trajectory for your own data (same format, not ranked) | Same folder |
| Lab report | PDF, at most 10 pages excluding the cover; use the [report template](PolyU_AAE4203_Lab_report_template.docx) | Blackboard |
| Code and data | GitHub repository with raw data, scripts and README.md; it must reproduce both CSV files | Link given in the report |

The report must contain the following parts (code listings in the report: at most 2 pages):

| Part | Required content |
|------|------------------|
| **1. Data collection** | Hardware and software set-up, site photos, time, environment and sky view near Block X. |
| **2. Raw data analysis** | Satellites tracked per constellation, C/N₀, elevation and sky plot, DOP, other quality indicators. |
| **3. SPP algorithm** | Measurement model, linearisation and weighted least squares; your Python implementation. |
| **4. Own-data validation** | Your SPP compared with RTKPOST and with the receiver solution (mean and standard deviation of the E/N/U differences), scatter about the mean, input checks. |
| **5. Benchmark results** | Each technique you tried, why, and its measured effect (ablation table); final 2D/3D RMSE, 95th percentile, availability and rank; discussion of urban effects. |
| **6. AI-use statement** | Tools used, what they were used for, key prompts (appendix), and how the group verified the output. |

📅 **Submission Deadline:** announced on Blackboard.

---


## Prerequisites

### Software Requirements
- **[u-center](https://www.u-blox.com/zh/product/u-center)** for data collection
- **[MATLAB](https://www.mathworks.com/products/matlab.html)** (optional, the CSV files of the example dataset are provided) with the following toolboxes:
  - Mapping Toolbox
- **[Python 3.9+](https://www.python.org/downloads/)** with the following packages:
  - numpy
  - matplotlib
  - pandas
  - pyubx2
  - cartopy (optional, for enhanced map visualization)
- **[RTKLIB 2.5](https://github.com/rtklibexplorer/RTKLIB/releases/download/v2.5.0/RTKLIB_EX_2.5.0.zip)** (RTKconv and RTKPOST utilities)
- **[Google Earth](https://earth.google.com/web/)** for .kml trajectory visualization

### Hardware Requirements

- Minimum 4GB RAM
- 1GB free disk space for data processing

## Installation

### 1. Clone the Repository

```bash
git clone https://github.com/weisongwen/AAE4203-S12627.git
cd AAE4203-S12627/Lab1
```

### 2. Install Python Dependencies

```bash
pip install numpy matplotlib pandas pyubx2 cartopy
```

**Note**: If you encounter issues with `cartopy`, you can install it using conda:
```bash
conda install -c conda-forge cartopy
```

For systems without cartopy, the code will automatically use fallback plotting methods.

### 3. Install MATLAB Dependencies

Ensure you have [MATLAB](https://www.mathworks.com/products/matlab.html) installed for rinex2csv convertion. Add the matRTKLIB directory to your MATLAB path.

### 4. Install RTKLIB

Download and install [RTKLIB_EX_2.5.0.zip](https://github.com/rtklibexplorer/RTKLIB/releases/download/v2.5.0/RTKLIB_EX_2.5.0.zip) from the official repository. Ensure `RTKCONV` and `RTKPOST` utilities are accessible from your system PATH.

## Workflow Overview

The complete pipeline consists of four main steps:

![Lab1 Workflow](examples/flow.png)

## Step 0: GNSS Data Collection
For detailed GNSS data collection procedures, please refer to Section 5 of the [AAE4203_26_27_S1_Lab_Tutorial.pdf](AAE4203_26_27_S1_Lab_Tutorial.pdf). Log 15–20 minutes of static data near Block X and name the file e.g. `G05_owndata_20261030.ubx`.

## Step 1: UBX to RINEX Conversion

### 1.1 Using rtkconv

1. Launch RTKconv from RTKLIB
2. Load your UBX file as input
3. Configure output settings:
   - Format: RINEX
   - Generate both observation (.obs) and navigation (.nav) files
4. Run conversion to generate RINEX files

### 1.2 Generate Reference Solution

1. Launch RTKPOST from RTKLIB
2. Load the generated RINEX observation and navigation files
3. Configure positioning mode (typically SPP for reference)
4. Generate reference positioning solution for accuracy comparison

## Step 2: RINEX Data Processing

### 2.1 Configure MATLAB Script

Open `rinex2csv.m` in MATLAB and update the file paths:

```matlab
% Change the obs and nav files to yours
obs = gt.Gobs(fullfile('COM6___9600_250916_020025_blockz.obs'));
nav = gt.Gnav(fullfile('COM6___9600_250916_020025_blockz.nav'));
```
```matlab
% Update these paths to match your data location
obs = gt.Gobs(fullfile('path/to/your/***.obs'));
nav = gt.Gnav(fullfile('path/to/your/***.nav'));
```

### 2.2 Run MATLAB Processing

Execute the script in MATLAB:

```matlab
run('rinex2csv.m')
```

### 2.3 Expected Outputs and Visualizations

The script will generate several visualizations and CSV files. Below are the examples:

#### Signal Quality Analysis
<p  align="center">
<img src="examples/rinex2csv1.png" width="50%">
</p>

*SNR threshold analysis showing signal quality distribution*

#### Satellite Observation Count
<p  align="center">
<img src="examples/rinex2csv2.png" width="50%">
</p>

*L1 observation count over time*

#### Sky Plot
<p  align="center">
<img src="examples/rinex2csv3.png" width="50%">
</p>

*Satellite visibility and elevation angles*

#### SPP Results
<p  align="center">
<img src="examples/rinex2csv4.png" width="50%">
</p>

*Single Point Positioning results from MATLAB processing*

### 2.4 Generated CSV Files

The following CSV files will be created for subsequent processing:
- `pseudoranges_meas.csv` - Raw pseudorange measurements
- `satellite_clock_bias.csv` - Satellite clock corrections
- `ionospheric_delay.csv` - Ionospheric delay corrections
- `tropospheric_delay.csv` - Tropospheric delay corrections
- `satellite_positions.csv` - Computed satellite positions

## Step 3: UBX to CSV Conversion

### 3.1 Install Dependencies

Ensure pyubx2 is installed:

```bash
pip install pyubx2
```

### 3.2 Run GUI Converter

Execute the Python GUI converter:

```bash
python pyubx2_csv_converter_gui.py
```

### 3.3 Using the GUI

1. **Select Input UBX File**: Browse and select your raw UBX data file
2. **Choose Output Directory**: Select where to save the converted CSV files
3. **Start Conversion**: Click "Start" to begin the conversion process

The GUI provides:
- Real-time progress indication
- Background processing to keep UI responsive
- Error handling and status updates

### 3.4 Key Output: Receiver Reference Solution

The conversion generates `NAV-HPPOSECEF.csv`, which contains the receiver's own high-precision ECEF solution. Without RTK corrections this is a standalone real-time solution that can be wrong by several metres in an urban canyon: use it as a reference for comparison, not as ground truth. Accuracy against a true trajectory is measured only on the benchmark dataset.

## Step 4: SPP Analysis

### 4.1 Configure Python Script

The example script `WLSE_spp.py` is provided with the example dataset (link on Blackboard). Open it and verify the CSV file paths match your generated files:

```python
# Verify these file paths
satellite_positions = np.loadtxt('satellite_positions.csv', delimiter=',')
pseudoranges_meas = np.loadtxt('pseudoranges_meas.csv', delimiter=',')
satellite_clock_bias = np.loadtxt('satellite_clock_bias.csv', delimiter=',')
ionospheric_delay = np.loadtxt('ionospheric_delay.csv', delimiter=',')
tropospheric_delay = np.loadtxt('tropospheric_delay.csv', delimiter=',')
ground_truth = np.genfromtxt('NAV-HPPOSECEF.csv', delimiter=',', skip_header=1, usecols=[2,3,4,5])
```

### 4.2 Run SPP Analysis

Execute the weighted least squares SPP analysis:

```bash
python WLSE_spp.py
```

### 4.3 Algorithm Features

The SPP implementation includes:
- **Weighted Least Squares**: Elevation-angle dependent weighting
- **Multi-iteration Convergence**: Iterative position refinement
- **Comprehensive Error Correction**: Ionospheric, tropospheric, and clock corrections
- **Robust Filtering**: Low elevation and poor signal quality masking

## Step 5: Accuracy Challenge and Leaderboard

Process the UrbanNav benchmark dataset with your best algorithm and submit one position per epoch. Verify every new version on the example dataset first.

### 5.1 CSV File Format

One header line, then one row per benchmark epoch in time order. Positions are WGS-84 ECEF coordinates in metres; no latitude/longitude and no extra columns.

| Column | Meaning | Unit / format |
|--------|---------|---------------|
| `gps_week` | GPS week number | integer |
| `gps_tow_s` | GPS time of week of the epoch | s, at least 1 decimal |
| `x_ecef_m` | Estimated receiver position, ECEF X | m, 3 decimals |
| `y_ecef_m` | Estimated receiver position, ECEF Y | m, 3 decimals |
| `z_ecef_m` | Estimated receiver position, ECEF Z | m, 3 decimals |

```text
gps_week,gps_tow_s,x_ecef_m,y_ecef_m,z_ecef_m
2158,455346.0,-2418202.586,5385847.524,2405648.799
2158,455347.0,-2418202.174,5385847.238,2405648.952
```

File names: `GXX_benchmark.csv` and `GXX_owndata.csv`, where `XX` is your two-digit group number (G01–G20).

### 5.2 Check and Upload

```bash
python tools/check_submission.py G05_benchmark.csv
```

Fix all errors, then upload both CSV files to your group's sub-folder in the TA's Google Drive folder (link on Blackboard). Replace the file whenever you have a better solution; only the latest `GXX_benchmark.csv` at the deadline is used for the final ranking. The leaderboard platform source code is available at [PolyU-TASLAB/AAE4203_ASN_2026](https://github.com/PolyU-TASLAB/AAE4203_ASN_2026).

### 5.3 Ranking Rules

- Errors are computed in the local East-North-Up frame at the ground-truth position, for epochs matched by GPS week and time of week (to 0.1 s).
- Only submissions with availability ≥ 90 % are ranked.
- Ranked submissions are ordered by 2D RMSE (smaller is better); ties are broken by 3D RMSE. The leaderboard also reports the 95th-percentile 2D error and availability.
- Your repository must reproduce the submitted file; the teaching team may re-run the code of any group.

## Results and Visualization

### Geographic and 3D Trajectory Analysis
![Trajectory Analysis](examples/Figure_1.png)
*Left: Geographic trajectory with map overlay; Right: 3D ENU trajectory visualization*

### Positioning Error Analysis
![Positioning Errors](examples/Figure_2.png)
*Time series analysis of East, North, Up, 2D, and 3D positioning errors*

### Statistical Error Distribution
![Error Distribution](examples/Figure_3.png)
*Histogram analysis showing 2D and 3D error distributions with statistical metrics*

### Error Correlation Analysis
![Error Scatter Analysis](examples/Figure_4.png)
*Left: Error evolution over time with smoothed trends; Right: East vs North error correlation*

## Output Files

### Primary Results
- `estimated_positions.csv`: Complete positioning solution with ECEF, LLA, and ENU coordinates
- `positioning_errors.csv`: Detailed error analysis for each epoch

### Generated Visualizations
- Geographic trajectory plots with cartographic features
- 3D trajectory visualization in ENU frame
- Comprehensive error analysis plots
- Statistical distribution analysis

### CSV File Structure

#### estimated_positions.csv
| Column | Description |
|--------|-------------|
| Epoch | Time epoch number |
| Latitude_deg | Latitude in decimal degrees |
| Longitude_deg | Longitude in decimal degrees |
| Altitude_m | Altitude in meters |
| Clock_Bias_m | Receiver clock bias in meters |
| East_m, North_m, Up_m | ENU coordinates |
| ECEF_X_m, ECEF_Y_m, ECEF_Z_m | ECEF coordinates |

#### positioning_errors.csv
| Column | Description |
|--------|-------------|
| Epoch | Time epoch number |
| Error_East_m | East positioning error |
| Error_North_m | North positioning error |
| Error_Up_m | Up positioning error |
| Error_2D_m | 2D horizontal error |
| Error_3D_m | 3D total error |

## Troubleshooting

### Common Issues

#### MATLAB Issues
- **matRTKLIB not found**: Ensure the matRTKLIB directory is added to MATLAB path
- **File path errors**: Use absolute paths or ensure files are in MATLAB working directory
- **Memory issues**: For large datasets, consider processing in smaller segments

#### Python Issues
- **pyubx2 import error**: Install using `pip install pyubx2`
- **cartopy installation issues**: Use conda installation or rely on fallback plotting
- **CSV file not found**: Verify file paths and ensure Step 2 completed successfully

#### Data Quality Issues
- **Insufficient satellites**: Check SNR thresholds and elevation angle masks
- **Poor positioning accuracy**: Verify correction models and satellite geometry
- **Conversion failures**: Ensure UBX files are valid and not corrupted

### Performance Optimization

#### For Large Datasets
- Process data in segments if memory becomes an issue
- Adjust visualization parameters for better performance
- Consider using more restrictive satellite selection criteria

#### Accuracy Improvements
- Adjust SNR and elevation angle thresholds
- Fine-tune weighting functions
- Consider additional error correction models

## Technical Notes

### Coordinate Systems
- **ECEF**: Earth-Centered, Earth-Fixed coordinates
- **LLA**: Latitude, Longitude, Altitude
- **ENU**: East, North, Up local coordinate system

### Error Correction Models
- Ionospheric delay correction using broadcast parameters
- Tropospheric delay using standard atmospheric models
- Satellite clock corrections from navigation message
- Elevation-angle dependent weighting for measurement quality

### Algorithm Convergence
- Maximum 100 iterations per epoch
- Convergence threshold: 1e-4 meters
- Automatic fallback to unweighted least squares if needed

## Contributing

Contributions are welcome! Please ensure:
- Code follows existing style conventions
- New features include appropriate documentation
- Test with various datasets before submission

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## References

- [RTKLIB-2.5](https://github.com/rtklibexplorer/RTKLIB/releases): An Open Source Program Package for GNSS Positioning 
- pyubx2: Python library for UBX protocol parsing
- [matRTKLIB](https://github.com/taroz/MatRTKLIB.git): MATLAB interface for RTKLIB functions

## Contact

For questions or support, please open an issue on the GitHub repository or contact TA dirrectly.
