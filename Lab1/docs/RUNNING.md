# Lab1 Running Instructions

## 1 Open the right folder

Open a terminal (for example PowerShell) in the `Lab1` folder of the course repository.
All relative paths below assume that `Lab1` is the current folder. Paths may contain spaces.

```powershell
python --version
python -m pip install -r requirements.txt
```

`numpy` is the only required third-party package for the numerical programs. Matplotlib and Cartopy are not required.
If `python` is not found, try `py -3` or the full path of the Python interpreter in your environment.

## 2 Run the example dataset

```powershell
python WLSE_spp.py --reference dataset/NAV-HPPOSECEF.csv
```

Outputs:

| File | Content |
| --- | --- |
| `outputs/example/G01_owndata.csv` | Standard five-column ECEF submission; times on a 0.1 s grid |
| `outputs/example/solution_details.csv` | Coordinates, receiver clock and number of measurements, with the full OBS time precision |
| `outputs/example/omitted_epochs.csv` | Epochs without a solution and the reason; only the header for the normal example |
| `outputs/example/summary.json` | Number of solved epochs and the optional comparison with the receiver reference |

Expected: 186/186 epochs, 2D RMSE 37.9107584162 m, 3D RMSE 51.0428443069 m.
The reference coordinates are stored in centimetres and `iTOW` in milliseconds. Epochs are matched by time, not by row number.
These numbers are differences from the receiver solution; the receiver solution is not a ground truth.

Check the submission format (the original one-argument call still works; the options below also check the actual observation times):

```powershell
python tools/check_submission.py outputs/example/G01_owndata.csv --epoch-file dataset/epochs.csv
```

This file is a demonstration of the example dataset, not the field data of any group.

Additional check: if the raw example OBS/NAV files are converted again with the new `rinex2csv.m` and solved in `gps` mode, 113/186 epochs are solved and 73 epochs are omitted under the current GPS masks. The 2D/3D RMS differences from the receiver reference are about 15.393/82.106 m. This GPS route uses different inputs and models from the archived CSV files, so it is not expected to reproduce the numbers in the slides.

## 3 Run the UrbanNav benchmark

Quick route without MATLAB:

```powershell
python WLSE_spp.py --mode gps --input-dir benchmark/preprocessed_gps --output-dir outputs/benchmark --submission G01_benchmark.csv
python tools/check_submission.py outputs/benchmark/G01_benchmark.csv --week 2158 --epoch-file benchmark/preprocessed_gps/epochs.csv
```

The program reads only the measurement CSV files, satellite positions, corrections and OBS times. No ground truth is needed.
`coverage_percent` in `summary.json` is the share of OBS epochs that were solved; it is not a measure of accuracy.

To repeat the RINEX conversion, open MATLAB with `Lab1` as the current folder and run:

```matlab
rinex2csv('benchmark/UrbanNav-HK-Deep-Urban-1.ublox.f9p.obs', ...
          'benchmark/UrbanNav_20210521_multiGNSS_clean.nav', ...
          'outputs/benchmark_inputs');
```

The conversion opens no figures. It uses the MatRTKLIB MEX files and does not need the Mapping Toolbox.
After the conversion:

```powershell
python WLSE_spp.py --mode gps --input-dir outputs/benchmark_inputs --output-dir outputs/benchmark_regenerated --submission G01_benchmark.csv
python tools/check_submission.py outputs/benchmark_regenerated/G01_benchmark.csv --week 2158 --epoch-file outputs/benchmark_inputs/epochs.csv
```

GPS baseline settings: GPS L1, 15 degree elevation mask, 15 dB-Hz minimum signal strength, broadcast ephemeris, model-based atmospheric corrections, TGD and Earth-rotation correction. Each failed epoch is omitted individually; the previous position is never copied.
The converter is written for the RINEX 3 GPS baseline of this lab; check it again before using other receivers or data.

## 4 Process your own data

1. Record a UBX log with u-center and back up the raw file.
2. Convert the UBX file to OBS/NAV with RTKCONV; note the software version and settings.
3. In MATLAB, call `rinex2csv` with your OBS file, NAV file and an output folder:

```matlab
rinex2csv('mydata/my.obs','mydata/my.nav','outputs/own_inputs');
```

4. Solve with Python, replacing G01 with your group number:

```powershell
python WLSE_spp.py --mode gps --input-dir outputs/own_inputs --output-dir outputs/own --submission G01_owndata.csv
python tools/check_submission.py outputs/own/G01_owndata.csv --epoch-file outputs/own_inputs/epochs.csv
```

5. Compare your solution with RTKPOST and with the receiver solution after aligning them in time. Do not call any of the three a ground truth.
6. Prepare the figures, report, code description and AI-use statement as required by the course.

## 5 Extract the receiver reference from UBX

```powershell
python -m pip install -r requirements-optional.txt
python pyubx2_csv_converter_gui.py
```

Select `COM6___9600_250916_020025_blockz.ubx` in the `Lab1` folder and write the output to `outputs/receiver`.
A conversion of this file without the GUI produced a 186-row `NAV-HPPOSECEF.csv`; the GUI buttons have not been tested separately.

## 6 Troubleshooting

| Message or situation | What to do |
| --- | --- |
| CSV file not found | Check the current folder and `--input-dir` |
| numpy missing | Run `-m pip install -r requirements.txt` with the same Python that runs the program |
| MATLAB cannot find `gt` | Keep `matRTKLIB/` next to `rinex2csv.m`; the function adds it to the path itself |
| MEX file cannot be loaded | Check that MATLAB and the operating system match the precompiled files; see the toolbox README |
| Fewer than 4 usable satellites | The epoch may be omitted; see the conversion summary and `omitted_epochs.csv`; never invent a position |
| Benchmark file only passes the format check | The error can only be computed by the TA server; this package contains no ground truth and cannot predict the rank |
| Output differs from the slides | Make sure you used `historical` mode with the 186-epoch inputs in `dataset/`; the new GPS baseline does not reproduce the historical numbers |

`legacy/` keeps the old scripts for comparison. They contain fixed file paths, plotting dependencies and a risk of mixing datasets; use the entry points in the `Lab1` folder for normal work.

## 7 Notes on the course versions

The marking weights, page limits, AI-use rules, submission locations and ranking rules of the course are kept in `Readme.md`; this file only adds the running steps.
The new MATLAB entry point is `rinex2csv(obs_file, nav_file, out_dir)`; do not call it with `run('rinex2csv.m')`.
The checker uses only the Python standard library. It supports the original one-file call as well as `--week`, `--epochs` and `--epoch-file`.
Its basic CSV format checks follow the online platform, with additional checks of the group file name and of the position near the Earth's surface; passing the checker does not guarantee acceptance by the server.
Verification results: [VERIFICATION.md](VERIFICATION.md).
