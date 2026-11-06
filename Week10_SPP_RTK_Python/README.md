# Week 10: Pure-Python SPP and Carrier-Phase RTK Demo

This example processes one hour of real rover and HKSC reference-station data.
The positioning algorithms are implemented in Python and do not invoke or read
RTKLIB.

## What is included

```text
Week10_SPP_RTK_Python/
├── README.md
├── requirements.txt
├── run_demo.ps1
├── python_gnss.py
├── python_batch_rtk.py
├── evaluate_python.py
└── data/
    ├── rover_gps_20210610_0300_0400.obs
    ├── hksc_gps_20210610_0300_0400.obs
    ├── hksc_gps_20210610.nav
    └── ground_truth_20210610_0300_0400.txt
```

Only the GPS observations and the 03:00–04:00 GPST interval used by this
example are included. The larger multi-GNSS source logs are intentionally not
duplicated here.

## Algorithms implemented

`python_gnss.py` contains:

- RINEX 3 observation and GPS navigation parsing
- GPS broadcast ephemeris and satellite-clock computation
- relativistic and Earth-rotation corrections
- Klobuchar ionosphere and Saastamoinen troposphere models
- iterative weighted least-squares SPP
- rover/base carrier-phase measurement construction
- integer least-squares candidate search

`python_batch_rtk.py` jointly estimates:

- rover position and receiver-clock difference at each epoch
- persistent single-difference L1 carrier ambiguities
- robust code and carrier residuals

The carrier result is labelled **RTK float**. The code reports a fixed solution
only if integer validation passes. It does not pass for this GPS-only subset,
so no unreliable centimetre-level result is claimed.

## Installation

Python 3.10 or newer is recommended.

```bash
python -m pip install -r requirements.txt
```

## Run the complete example

On Windows PowerShell:

```powershell
.\run_demo.ps1
```

Or run the two steps directly:

```bash
python python_batch_rtk.py \
  --rover data/rover_gps_20210610_0300_0400.obs \
  --base data/hksc_gps_20210610_0300_0400.obs \
  --nav data/hksc_gps_20210610.nav \
  --out results

python evaluate_python.py
```

The run creates:

- `results/spp_python.csv`
- `results/rtk_float_python.csv`
- `results/ambiguities_python.csv`
- `results/run_summary.json`
- `results/python_stats.json`
- three PNG figures

## Expected result

For the supplied one-hour data window:

- SPP solves approximately 3,593 epochs.
- Median SPP horizontal error is about 4.3 m.
- GPS-only RTK float is available for approximately 29% of the hour.
- The integer best/second-best ratio is close to 1, below the acceptance
  threshold of 3; therefore FIX is rejected.

The ground-truth point is used only after positioning to evaluate error. It is
not used to initialise or constrain either solution. The truth reference point
is not demonstrably the rover antenna phase centre, so a constant lever-arm
offset remains in the raw error.

## Data notes

- Rover: u-blox RINEX 3.03 observation data.
- Base: HKSC reference-station RINEX 3.02 observation data.
- Navigation: GPS broadcast navigation data associated with HKSC.
- Time system: GPST, 10 June 2021, 03:00–04:00.
- The rover uses GPS L2L while HKSC uses L2S. This signal mismatch and sparse
  common GPS geometry make integer ambiguity validation difficult.
