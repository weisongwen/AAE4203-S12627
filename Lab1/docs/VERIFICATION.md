# Lab1 Verification

Verification of the runnable Lab1 package (October 2026).

## Course files and rules

- The lab report contents, page limits, AI-use rules, submission locations, CSV column definitions, ranking threshold and rules, marking weights, course policies and schedule are unchanged.
- The PDFs, report template, licence, MatRTKLIB, workflow figure and the materials of other weeks are unchanged.
- The original course versions of the converter and the checker are kept in `legacy/`.

## Results

| Check | Result |
| --- | --- |
| Example dataset | 186/186 epochs; ECEF positions identical to the archived results (maximum difference 0 m) |
| Difference from the receiver reference | 2D RMS 37.9107584162 m, 3D RMS 51.0428443069 m; not an accuracy against an independent ground truth |
| Included GPS cache to Python | 1437/1534 epochs, availability 93.6766623207 %; 97 epochs explicitly omitted |
| MATLAB conversion repeated | The 10 files generated from the included OBS/NAV are byte-identical to the verified cache |
| Python solution after reconversion | 1437/1534 epochs; CSV format and OBS time matching checks passed |
| Original one-file checker call | `python tools/check_submission.py file.csv` works |
| Checker with time tables | `--week`, `--epochs` and `--epoch-file` work as described in the README |
| Numerical regression | Synthetic GPS position and clock recovered within 1 mm; historical numbers unchanged |
| Missing epochs | Omitted with a reason; never filled with the previous position |
| Malformed CSV | Duplicate times, NaN, wrong number of columns, wrong coordinate units, empty data and wrong order are rejected |
| Online format comparison | 22 representative valid and invalid inputs give the same decision as the current server parser |

The checker additionally keeps the course checks of the group file name and of the position near the Earth's surface.
Agreement on representative inputs does not prove equivalence for every input; acceptance, time matching and scoring are decided by the online system.
With `--epoch-file`, the reported coverage is the share of OBS epochs matched; a format pass does not replace an accuracy assessment.

## How to run the checks

In the `Lab1` folder:

```powershell
python -m pip install -r requirements.txt
python WLSE_spp.py --reference example/NAV-HPPOSECEF.csv
python WLSE_spp.py --mode gps --input-dir benchmark/gps_csv --output-dir outputs/benchmark --submission G05_benchmark.csv
python tools/check_submission.py outputs/benchmark/G05_benchmark.csv --week 2158 --epoch-file benchmark/gps_csv/epochs.csv
python tools/verify_package.py
```

In MATLAB, with `Lab1` as the current folder:

```matlab
rinex2csv('benchmark/urbannav_deep1.obs', ...
          'benchmark/urbannav_deep1.nav', ...
          'outputs/benchmark_regenerated_inputs');
```

Then in the terminal:

```powershell
python WLSE_spp.py --mode gps --input-dir outputs/benchmark_regenerated_inputs --output-dir outputs/benchmark_regenerated --submission G05_benchmark.csv
python tools/check_submission.py outputs/benchmark_regenerated/G05_benchmark.csv --week 2158 --epoch-file outputs/benchmark_regenerated_inputs/epochs.csv
```

## Scope

Environment: Windows x64, Python 3.10, NumPy 2.2, MATLAB R2023a. The numerical command-line tools need only NumPy; the receiver converter is optional.

The GPS baseline does not implement multi-constellation positioning. This package contains no UrbanNav scoring ground truth, so benchmark ranks cannot be computed offline.
The archived example CSV files and the new GPS conversion use different inputs and masks; the 113/186-epoch result of reconverting the raw example must not be confused with the 186-epoch historical regression.
Field logging with a receiver, the RTKCONV GUI, campus network access, Google Drive synchronisation and online scoring were not tested.
Generated files are written to `outputs/`, which is excluded by `.gitignore`.
