"""AAE4203 Week 9 lab (2026/27): check a leaderboard CSV before uploading it.

Usage:  python check_submission.py G05_benchmark.csv
"""
import csv
import math
import re
import sys

HEADER = ["gps_week", "gps_tow_s", "x_ecef_m", "y_ecef_m", "z_ecef_m"]


def check(path):
    errors, warnings = [], []
    name = path.replace("\\", "/").split("/")[-1]
    if not re.fullmatch(r"G(0[1-9]|1[0-9]|20)_(benchmark|owndata)\.csv", name):
        errors.append(f"file name '{name}' must be GXX_benchmark.csv or GXX_owndata.csv (XX = 01..20)")
    with open(path, newline="", encoding="utf-8-sig") as f:
        rows = list(csv.reader(f))
    if not rows or [c.strip() for c in rows[0]] != HEADER:
        errors.append("first line must be exactly: " + ",".join(HEADER))
        return errors, warnings, 0
    seen, last = set(), -math.inf
    for i, r in enumerate(rows[1:], start=2):
        if len(r) != 5:
            errors.append(f"line {i}: expected 5 values, found {len(r)}"); continue
        try:
            wk, tow, x, y, z = int(r[0]), float(r[1]), float(r[2]), float(r[3]), float(r[4])
        except ValueError:
            errors.append(f"line {i}: non-numeric value"); continue
        if not 0 <= tow < 604800:
            errors.append(f"line {i}: gps_tow_s must be in [0, 604800)")
        t = round(wk * 604800 + tow, 1)
        if t in seen:
            errors.append(f"line {i}: duplicate epoch")
        if t < last:
            warnings.append(f"line {i}: epochs are not in time order")
        seen.add(t); last = t
        rad = math.sqrt(x * x + y * y + z * z)
        if not 6.30e6 < rad < 6.42e6:
            errors.append(f"line {i}: |r| = {rad:.0f} m is not on the Earth's surface (ECEF in metres?)")
    return errors, warnings, len(rows) - 1


if __name__ == "__main__":
    errs, warns, n = check(sys.argv[1])
    for w in warns[:10]: print("WARNING:", w)
    for e in errs[:20]: print("ERROR:", e)
    print(f"{n} epochs read; " + ("format OK - ready to upload." if not errs else f"{len(errs)} error(s) - please fix before uploading."))
    sys.exit(1 if errs else 0)
