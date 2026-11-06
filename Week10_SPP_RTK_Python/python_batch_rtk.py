"""Pure-Python batch GPS carrier-phase differential positioning (RTK float).

The estimator jointly solves each epoch's rover position/receiver-clock
difference and persistent single-difference L1 carrier ambiguities.  It uses
SciPy only for sparse least squares; no RTKLIB executable, library, output, or
position is read.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import lsqr

import python_gnss as g


def build_records(rover, base, nav, alpha, beta, spp, base_xyz):
    records = []
    arc_no = {}
    last_seen = {}
    for tow in sorted(set(rover) & set(base) & set(spp)):
        rows = g.common_measurements(rover[tow], base[tow], tow, nav, spp[tow][0],
                                     base_xyz, alpha, beta)
        if len(rows) < 4:
            continue
        for r in rows:
            sat = r["sat"]
            if sat not in arc_no:
                arc_no[sat] = 0
            if (r["lli"] & 1) or (sat in last_seen and tow - last_seen[sat] > 120):
                arc_no[sat] += 1
            r["arc"] = (sat, arc_no[sat])
            last_seen[sat] = tow
        records.append((tow, rows))
    return records


def solve(records, spp, base_xyz, iterations=5):
    n_epoch = len(records)
    arcs = sorted({r["arc"] for _, rows in records for r in rows})
    arc_col = {a: 4 * n_epoch + i for i, a in enumerate(arcs)}
    n_state = 4 * n_epoch + len(arcs)
    x = np.zeros((n_epoch, 4))
    for i, (tow, rows) in enumerate(records):
        x[i, :3] = spp[tow][0]
        clocks = [r["Pr"] - (r["rr"] - r["rb"]) - r["code_corr"] for r in rows]
        x[i, 3] = np.median(clocks)
    amb = np.zeros(len(arcs))
    vals = {a: [] for a in arcs}
    for i, (_, rows) in enumerate(records):
        for r in rows:
            vals[r["arc"]].append(r["Lr"] - (r["rr"] - r["rb"]) - x[i, 3] - r["phase_corr"])
    for a, v in vals.items():
        amb[arc_col[a] - 4 * n_epoch] = np.median(v)

    code_res = phase_res = np.array([])
    for it in range(iterations):
        rr_idx, cc_idx, data, rhs = [], [], [], []
        code_residuals, phase_residuals = [], []
        row = 0
        for i, (tow, obs) in enumerate(records):
            pos, clk = x[i, :3], x[i, 3]
            for r in obs:
                rho_r = np.linalg.norm(r["sr"] - pos)
                rho_b = np.linalg.norm(r["sb"] - base_xyz)
                los = (pos - r["sr"]) / rho_r
                common = rho_r - rho_b + clk
                residuals = (
                    r["Pr"] - (common + r["code_corr"]),
                    r["Lr"] - (common + r["phase_corr"] + amb[arc_col[r["arc"]] - 4 * n_epoch]),
                )
                sigmas = (max(1.2 / max(np.sin(r["elev"]), .25), 1.2),
                          max(.012 / max(np.sin(r["elev"]), .25), .012))
                for is_phase, (res, sigma) in enumerate(zip(residuals, sigmas)):
                    # Huber downweighting prevents a code outlier or slip from
                    # pulling the whole ambiguity arc.
                    norm = abs(res) / sigma
                    robust = 1.0 if norm <= 3 else np.sqrt(3 / norm)
                    w = robust / sigma
                    cols = [4 * i, 4 * i + 1, 4 * i + 2, 4 * i + 3]
                    values = [*los, 1.0]
                    if is_phase:
                        cols.append(arc_col[r["arc"]])
                        values.append(1.0)
                        phase_residuals.append(res)
                    else:
                        code_residuals.append(res)
                    for c, v in zip(cols, values):
                        rr_idx.append(row); cc_idx.append(c); data.append(w * v)
                    rhs.append(w * res)
                    row += 1
            # Weak SPP prior only stabilises epochs with marginal geometry.
            prior = spp[tow][0]
            for j in range(3):
                rr_idx.append(row); cc_idx.append(4 * i + j); data.append(1 / 30.0)
                rhs.append((prior[j] - pos[j]) / 30.0)
                row += 1
        A = coo_matrix((data, (rr_idx, cc_idx)), shape=(row, n_state)).tocsr()
        d = lsqr(A, np.asarray(rhs), atol=1e-7, btol=1e-7, iter_lim=1000)[0]
        x += d[:4 * n_epoch].reshape(-1, 4)
        amb += d[4 * n_epoch:]
        code_res = np.asarray(code_residuals)
        phase_res = np.asarray(phase_residuals)
        print(f"iteration {it + 1}: step position RMS "
              f"{np.sqrt(np.mean(d[:4*n_epoch].reshape(-1,4)[:,:3]**2)):.4f} m; "
              f"code median {np.median(abs(code_res)):.3f} m; "
              f"phase median {np.median(abs(phase_res)):.4f} m")
        if np.sqrt(np.mean(d[:4*n_epoch].reshape(-1, 4)[:, :3] ** 2)) < 1e-3:
            break
    # Float ambiguity covariance after eliminating each epoch's independent
    # position/clock block (Schur complement).
    N = (A.T @ A).tocsc()
    Naa = N[4 * n_epoch:, 4 * n_epoch:].toarray()
    schur = Naa.copy()
    for i in range(n_epoch):
        sl = slice(4 * i, 4 * i + 4)
        B = N[sl, sl].toarray()
        C = N[sl, 4 * n_epoch:].toarray()
        try:
            schur -= C.T @ np.linalg.solve(B, C)
        except np.linalg.LinAlgError:
            pass
    Qamb = np.linalg.pinv((schur + schur.T) / 2, rcond=1e-12)
    return x, arcs, amb, Qamb, code_res, phase_res


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--rover", required=True, type=Path)
    ap.add_argument("--base", required=True, type=Path)
    ap.add_argument("--nav", required=True, type=Path)
    ap.add_argument("--out", type=Path, default=Path("python_batch_results"))
    ap.add_argument("--start", type=int, default=356400)
    ap.add_argument("--end", type=int, default=359999)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    nav, alpha, beta = g.read_nav(args.nav)
    rover, _, _ = g.read_obs(args.rover, args.start, args.end)
    base, base_xyz, _ = g.read_obs(args.base, args.start, args.end)
    spp = g.solve_spp(rover, nav, alpha, beta)
    records = build_records(rover, base, nav, alpha, beta, spp, base_xyz)
    x, arcs, amb, Qamb, cr, pr = solve(records, spp, base_xyz)
    # Integer least-squares diagnostic.  A common single-difference phase
    # offset is removed by differencing all arcs against the best-observed arc.
    counts = {a: sum(r["arc"] == a for _, rows in records for r in rows) for a in arcs}
    ref = max(arcs, key=counts.get)
    ir = arcs.index(ref)
    keep = [i for i in range(len(arcs)) if i != ir and counts[arcs[i]] >= 20]
    T = np.zeros((len(keep), len(arcs)))
    for k, i in enumerate(keep):
        T[k, i] = 1 / g.L1
        T[k, ir] = -1 / g.L1
    af = T @ amb
    Qdd = T @ Qamb @ T.T
    ils = g.integer_least_squares(af, Qdd, max_nodes=1000000) if len(keep) >= 3 else None
    integer_diag = None
    if ils:
        n1, s1, n2, s2, nodes = ils
        integer_diag = {"reference_arc": list(ref), "ambiguities": len(n1), "best_score": s1,
                        "second_score": s2, "ratio": s2 / max(s1, 1e-12), "search_nodes": nodes,
                        "accepted": bool(s2 / max(s1, 1e-12) >= 3)}
        print("integer diagnostic", integer_diag)
    with (args.out / "spp_python.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f); w.writerow(["gps_tow_s", "x_ecef_m", "y_ecef_m", "z_ecef_m"])
        for t, v in spp.items():
            w.writerow([t, *v[0]])
    with (args.out / "rtk_float_python.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["gps_tow_s", "x_ecef_m", "y_ecef_m", "z_ecef_m", "receiver_clock_sd_m", "satellites"])
        for i, (t, rows) in enumerate(records):
            w.writerow([t, *x[i, :3], x[i, 3], len(rows)])
    with (args.out / "ambiguities_python.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f); w.writerow(["satellite", "arc", "float_ambiguity_m", "float_ambiguity_L1_cycles"])
        for a, b in zip(arcs, amb):
            w.writerow([a[0], a[1], b, b / g.L1])
    summary = {
        "implementation": "pure Python/NumPy/SciPy; no RTKLIB input",
        "method": "batch single-difference GPS L1 code/carrier RTK float",
        "spp_epochs": len(spp), "rtk_float_epochs": len(records), "ambiguity_arcs": len(arcs),
        "median_absolute_code_residual_m": float(np.median(abs(cr))),
        "median_absolute_phase_residual_m": float(np.median(abs(pr))),
        "integer_fix_claimed": False, "integer_diagnostic": integer_diag,
    }
    (args.out / "run_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
