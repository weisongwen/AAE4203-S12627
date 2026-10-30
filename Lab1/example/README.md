# Example Dataset (186 Epochs, Block X, 2025)

These inputs belong to the example shown in the lab slides and the tutorial. All example files, including the raw UBX log, are in this folder.
The five measurement CSV files use one column per epoch (satellite positions use three adjacent columns per epoch); they are not one row per epoch.

| Content | Files |
| --- | --- |
| Raw UBX log | `blockx_2025.ubx` |
| Raw observations and navigation | `blockx_2025.obs` and the matching `.nav` (RTKCONV output) |
| HK SatRef base station for DGNSS | `hksc_2025.obs` (station HKSC, 1 Hz, 02:00:30 to 02:04:30 GPST) |
| Measurement inputs | `pseudoranges_meas.csv`, `satellite_positions.csv`, `satellite_clock_bias.csv`, `ionospheric_delay.csv`, `tropospheric_delay.csv` |
| Actual observation times | `epochs.csv`, extracted from the matching OBS file |
| Receiver reference | `NAV-HPPOSECEF.csv` (receiver solution, not a ground truth) |
| Archived results | `reference_results/estimated_positions.csv`, `positioning_errors.csv` |
| Other original example results | `blockx_2025.pos` (RTKPOST SPP) and `blockx_2025.kml` |

The default `historical` mode reads only the inputs in this folder and writes to a separate `outputs/example/` folder.
The original CSV files carry no satellite identifiers, so results of the historical mode and of the new GPS-only baseline should be interpreted separately.

For a DGNSS solution in RTKPOST, use `blockx_2025.obs` as the rover, `hksc_2025.obs` as the base station and `blockx_2025.nav` as the navigation file (see `docs/RUNNING.md`).
