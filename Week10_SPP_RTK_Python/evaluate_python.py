"""Evaluate and plot the from-scratch Python SPP and RTK-float solutions."""
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import python_gnss as g

ROOT = Path(__file__).resolve().parent
TRUTH = ROOT / "data" / "ground_truth_20210610_0300_0400.txt"
if not TRUTH.exists():
    TRUTH = ROOT.parent / "20210610_1hz_converted.txt"
if not TRUTH.exists():
    TRUTH = Path(r"C:\Users\user\Dropbox\03_PolyU Teaching & Service\Teaching\PolyU AAE 4203 Guidance and Navigation"
                 r"\2026-2027 Semester 1\Subject lecture Notes\06 Week 10 GNSS RTK\ExampleData"
                 r"\20210610_1hz_converted.txt")
IN = ROOT / "results"
if not IN.exists():
    IN = ROOT / "python_full"
T0 = 356400


def read_solution(path):
    t, p = [], []
    for r in csv.DictReader(path.open()):
        t.append(int(r["gps_tow_s"]))
        p.append([float(r[k]) for k in ("x_ecef_m", "y_ecef_m", "z_ecef_m")])
    return np.array(t), np.array(p)


raw = np.loadtxt(TRUTH)
truth = {int(round(r[2])): g.llh2xyz(*r[3:6]) for r in raw}
ref = truth[T0]
lat, lon, _ = g.xyz2llh(ref)
R = np.array([[-np.sin(lon), np.cos(lon), 0],
              [-np.sin(lat)*np.cos(lon), -np.sin(lat)*np.sin(lon), np.cos(lat)],
              [np.cos(lat)*np.cos(lon), np.cos(lat)*np.sin(lon), np.sin(lat)]])


def evaluate(path):
    t, p = read_solution(path)
    tr = np.array([truth[k] for k in t])
    e = (p - tr) @ R.T
    return t, p, e


ts, ps, es = evaluate(IN / "spp_python.csv")
tr, pr, er = evaluate(IN / "rtk_float_python.csv")
truth_s = np.array([truth[k] for k in ts])
truth_r = np.array([truth[k] for k in tr])
origin = truth_s[0]
track = (truth_s - origin) @ R.T
spp_track = (ps - origin) @ R.T
rtk_track = (pr - origin) @ R.T


def metrics(e, total):
    h = np.linalg.norm(e[:, :2], axis=1)
    v = abs(e[:, 2])
    centred = e - np.median(e, axis=0)
    hc = np.linalg.norm(centred[:, :2], axis=1)
    return {"epochs": len(e), "availability_percent": 100 * len(e) / total,
            "median_horizontal_m": float(np.median(h)),
            "p95_horizontal_m": float(np.percentile(h, 95)),
            "median_vertical_m": float(np.median(v)),
            "p95_vertical_m": float(np.percentile(v, 95)),
            "centred_median_horizontal_m": float(np.median(hc)),
            "centred_p95_horizontal_m": float(np.percentile(hc, 95))}


stats = {"SPP": metrics(es, 3600), "RTK float": metrics(er, 3600),
         "notes": ["Truth reference point is not the GNSS antenna phase centre.",
                   "Centred statistics remove only the median 3-D offset for precision comparison.",
                   "No integer-fixed solution is claimed."]}
(IN / "python_stats.json").write_text(json.dumps(stats, indent=2) + "\n")
print(json.dumps(stats, indent=2))

plt.rcParams.update({"font.family": "Arial", "font.size": 10, "axes.grid": True, "grid.alpha": .3,
                     "savefig.dpi": 200, "savefig.bbox": "tight"})

fig, (a, b) = plt.subplots(1, 2, figsize=(10.5, 5.2))
a.plot(spp_track[:, 0], spp_track[:, 1], ".", ms=2, color="#C0392B", label="Python SPP")
a.plot(rtk_track[:, 0], rtk_track[:, 1], ".", ms=2, color="#6C3483", label="Python RTK float")
a.plot(track[:, 0], track[:, 1], "k-", lw=1, label="Ground truth")
a.set_aspect("equal"); a.set_xlabel("East (m)"); a.set_ylabel("North (m)"); a.set_title("One-hour trajectory")
a.legend(fontsize=9)
m = (ts >= T0 + 60) & (ts < T0 + 600)
c = track[m].mean(0)
b.plot(spp_track[m, 0] - c[0], spp_track[m, 1] - c[1], ".", ms=3, color="#C0392B", label="SPP")
mr = (tr >= T0 + 60) & (tr < T0 + 600)
b.plot(rtk_track[mr, 0] - c[0], rtk_track[mr, 1] - c[1], ".", ms=3, color="#6C3483", label="RTK float")
b.plot(0, 0, "k+", ms=12, mew=2, label="Truth reference")
b.set_aspect("equal"); b.set_xlabel("East (m)"); b.set_ylabel("North (m)")
b.set_title("Parked segment, 03:01-03:10"); b.legend(fontsize=9)
fig.tight_layout(); fig.savefig(IN / "python_trajectory.png"); plt.close(fig)

fig, a = plt.subplots(figsize=(11.5, 4.3))
a.semilogy((ts-T0)/60, np.maximum(np.linalg.norm(es[:, :2], axis=1), .01), ".", ms=2.2,
           color="#C0392B", label="Python SPP")
a.semilogy((tr-T0)/60, np.maximum(np.linalg.norm(er[:, :2], axis=1), .01), ".", ms=2.5,
           color="#6C3483", label="Python RTK float")
a.set_xlabel("Minutes after 03:00 GPST"); a.set_ylabel("Horizontal error to truth reference (m)")
a.set_xlim(0, 60); a.legend(); fig.tight_layout(); fig.savefig(IN / "python_error_time.png"); plt.close(fig)

fig, a = plt.subplots(figsize=(6.2, 4.3))
for label, e, col in (("SPP", es, "#C0392B"), ("RTK float", er, "#6C3483")):
    h = np.sort(np.linalg.norm(e[:, :2], axis=1))
    a.semilogx(h, 100*np.arange(1, len(h)+1)/len(h), lw=2, label=label, color=col)
a.set_xlabel("Horizontal error to truth reference (m)"); a.set_ylabel("Epochs below x (%)")
a.set_ylim(0, 100); a.legend(); fig.tight_layout(); fig.savefig(IN / "python_cdf.png"); plt.close(fig)
