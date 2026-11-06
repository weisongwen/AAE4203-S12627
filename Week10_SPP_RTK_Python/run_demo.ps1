$ErrorActionPreference = "Stop"
$Root = Split-Path -Parent $MyInvocation.MyCommand.Path

Push-Location $Root
try {
    python python_batch_rtk.py `
        --rover data\rover_gps_20210610_0300_0400.obs `
        --base data\hksc_gps_20210610_0300_0400.obs `
        --nav data\hksc_gps_20210610.nav `
        --out results
    python evaluate_python.py
} finally {
    Pop-Location
}
