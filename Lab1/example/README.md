# Example Dataset (186 Epochs, Block X, 2025)

These inputs belong to the example shown in the lab slides. The raw UBX log is in the parent folder.
The five measurement CSV files use one column per epoch (satellite positions use three adjacent columns per epoch); they are not one row per epoch.

| Content | Files |
| --- | --- |
| Raw observations and navigation | `blockx_2025.obs` and the matching `.nav` |
| Measurement inputs | `pseudoranges_meas.csv`, `satellite_positions.csv`, `satellite_clock_bias.csv`, `ionospheric_delay.csv`, `tropospheric_delay.csv` |
| Actual observation times | `epochs.csv`, extracted from the matching OBS file |
| Receiver reference | `NAV-HPPOSECEF.csv` (receiver solution, not a ground truth) |
| Archived results | `reference_results/estimated_positions.csv`, `positioning_errors.csv` |
| Other original example results | `.pos` and `.kml` files with the same base name |

The default `historical` mode reads only the inputs in this folder and writes to a separate `outputs/example/` folder.
The original CSV files carry no satellite identifiers, so results of the historical mode and of the new GPS-only baseline should be interpreted separately.
