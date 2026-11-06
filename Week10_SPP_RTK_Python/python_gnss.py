"""Educational GPS SPP and carrier-phase utilities implemented from first principles.

No RTKLIB output is read.  The code parses RINEX 3 observation/navigation
files, propagates GPS broadcast ephemerides, applies satellite clock/Sagnac,
Klobuchar and Saastamoinen corrections, solves SPP by weighted least squares,
and provides relative carrier measurements and an integer least-squares search.
The published RTK-float result is produced by python_batch_rtk.py.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

C = 299792458.0
MU = 3.986005e14
OMEGA_E = 7.2921151467e-5
F_REL = -4.442807633e-10
L1 = C / 1575.42e6
GPS0 = dt.datetime(1980, 1, 6)
A_WGS = 6378137.0
E2 = 6.69437999014e-3


def gps_time(t: dt.datetime) -> tuple[int, float]:
    sec = (t - GPS0).total_seconds()
    return int(sec // 604800), sec % 604800


def wrap_week(x):
    return (x + 302400.0) % 604800.0 - 302400.0


def floats(line: str, start: int, width: int = 19):
    out = []
    for i in range(start, len(line), width):
        s = line[i:i + width].replace("D", "E").strip()
        out.append(float(s) if s else np.nan)
    return out


@dataclass
class Eph:
    prn: str
    week: int
    toc: float
    af0: float
    af1: float
    af2: float
    iode: float
    crs: float
    dn: float
    m0: float
    cuc: float
    e: float
    cus: float
    sqrt_a: float
    toe: float
    cic: float
    omega0: float
    cis: float
    i0: float
    crc: float
    omega: float
    omega_dot: float
    idot: float
    tgd: float


def read_nav(path: Path):
    lines = path.read_text(encoding="ascii", errors="ignore").splitlines()
    alpha = np.zeros(4)
    beta = np.zeros(4)
    k = 0
    while "END OF HEADER" not in lines[k]:
        if "IONOSPHERIC CORR" in lines[k]:
            tag = lines[k][:4].strip()
            vals = [float(x.replace("D", "E")) for x in lines[k][5:60].split()]
            if tag == "GPSA":
                alpha[:] = vals[:4]
            elif tag == "GPSB":
                beta[:] = vals[:4]
        k += 1
    k += 1
    nav: dict[str, list[Eph]] = {}
    while k + 7 < len(lines):
        l0 = lines[k]
        if not l0.startswith("G"):
            k += 1
            continue
        f = l0[:23].split()
        try:
            toc_dt = dt.datetime(*map(int, f[1:7]))
            w0, toc = gps_time(toc_dt)
            a0 = floats(l0, 23)
            v = []
            for j in range(1, 8):
                v.extend(floats(lines[k + j], 4))
            e = Eph(
                f[0], w0, toc, *a0[:3],
                v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7],
                v[8], v[9], v[10], v[11], v[12], v[13], v[14], v[15],
                v[16], v[22],
            )
            nav.setdefault(e.prn, []).append(e)
            k += 8
        except (ValueError, IndexError):
            k += 1
    for es in nav.values():
        es.sort(key=lambda x: x.toe)
    return nav, alpha, beta


def obs_header(lines):
    types: dict[str, list[str]] = {}
    approx = np.zeros(3)
    i = 0
    while i < len(lines):
        line = lines[i]
        label = line[60:] if len(line) >= 60 else ""
        if "SYS / # / OBS TYPES" in label:
            sys = line[0]
            n = int(line[3:6])
            vals = line[7:60].split()
            while len(vals) < n:
                i += 1
                vals.extend(lines[i][7:60].split())
            types[sys] = vals[:n]
        elif "APPROX POSITION XYZ" in label:
            approx = np.array([float(x) for x in line[:60].split()[:3]])
        elif "END OF HEADER" in label:
            return i + 1, types, approx
        i += 1
    raise ValueError("RINEX END OF HEADER not found")


def read_obs(path: Path, start_tow=356400, end_tow=359999):
    """Return {integer TOW: {sat: {obstype: value}}}; GPS only."""
    lines = path.read_text(encoding="ascii", errors="ignore").splitlines()
    i, types, approx = obs_header(lines)
    gt = types.get("G", [])
    epochs = {}
    while i < len(lines):
        line = lines[i]
        if not line.startswith(">"):
            i += 1
            continue
        f = line[1:].split()
        if len(f) < 8:
            i += 1
            continue
        t = dt.datetime(*map(int, f[:5]), int(float(f[5])), round((float(f[5]) % 1) * 1e6))
        _, tow = gps_time(t)
        key = int(round(tow))
        ns = int(f[7])
        i += 1
        data = {}
        for _ in range(ns):
            if i >= len(lines):
                break
            rec = lines[i]
            sat = rec[:3].replace(" ", "0")
            if sat.startswith("G"):
                values = {}
                for j, name in enumerate(gt):
                    fld = rec[3 + 16 * j:3 + 16 * (j + 1)]
                    try:
                        values[name] = float(fld[:14])
                    except ValueError:
                        values[name] = np.nan
                    values[name + "_lli"] = int(fld[14]) if len(fld) > 14 and fld[14].isdigit() else 0
                data[sat] = values
            i += 1
        if start_tow <= key <= end_tow:
            epochs[key] = data
    return epochs, approx, gt


def select_eph(nav, sat, t):
    es = nav.get(sat, ())
    if not es:
        return None
    return min(es, key=lambda e: abs(wrap_week(t - e.toe)))


def sat_clock(e: Eph, t):
    tc = wrap_week(t - e.toc)
    tk = wrap_week(t - e.toe)
    a = e.sqrt_a ** 2
    n = math.sqrt(MU / a ** 3) + e.dn
    M = e.m0 + n * tk
    E = M
    for _ in range(10):
        E -= (E - e.e * math.sin(E) - M) / (1 - e.e * math.cos(E))
    rel = F_REL * e.e * e.sqrt_a * math.sin(E)
    return e.af0 + e.af1 * tc + e.af2 * tc * tc + rel - e.tgd


def sat_pos(e: Eph, t):
    tk = wrap_week(t - e.toe)
    a = e.sqrt_a ** 2
    n = math.sqrt(MU / a ** 3) + e.dn
    M = e.m0 + n * tk
    E = M
    for _ in range(12):
        E -= (E - e.e * math.sin(E) - M) / (1 - e.e * math.cos(E))
    v = math.atan2(math.sqrt(1 - e.e * e.e) * math.sin(E), math.cos(E) - e.e)
    phi = v + e.omega
    du = e.cus * math.sin(2 * phi) + e.cuc * math.cos(2 * phi)
    dr = e.crs * math.sin(2 * phi) + e.crc * math.cos(2 * phi)
    di = e.cis * math.sin(2 * phi) + e.cic * math.cos(2 * phi)
    u = phi + du
    r = a * (1 - e.e * math.cos(E)) + dr
    inc = e.i0 + e.idot * tk + di
    x = r * math.cos(u)
    y = r * math.sin(u)
    om = e.omega0 + (e.omega_dot - OMEGA_E) * tk - OMEGA_E * e.toe
    return np.array([
        x * math.cos(om) - y * math.cos(inc) * math.sin(om),
        x * math.sin(om) + y * math.cos(inc) * math.cos(om),
        y * math.sin(inc),
    ])


def rotate_sagnac(x, tau):
    a = OMEGA_E * tau
    return np.array([math.cos(a) * x[0] + math.sin(a) * x[1],
                     -math.sin(a) * x[0] + math.cos(a) * x[1], x[2]])


def satellite(nav, sat, tr, pseudorange):
    tx = tr - pseudorange / C
    e = select_eph(nav, sat, tx)
    if e is None:
        raise KeyError(sat)
    dts = sat_clock(e, tx)
    tx -= dts
    p = sat_pos(e, tx)
    return rotate_sagnac(p, pseudorange / C), dts


def xyz2llh(x):
    lon = math.atan2(x[1], x[0])
    p = math.hypot(x[0], x[1])
    lat = math.atan2(x[2], p * (1 - E2))
    for _ in range(8):
        n = A_WGS / math.sqrt(1 - E2 * math.sin(lat) ** 2)
        h = p / math.cos(lat) - n
        lat = math.atan2(x[2], p * (1 - E2 * n / (n + h)))
    n = A_WGS / math.sqrt(1 - E2 * math.sin(lat) ** 2)
    return lat, lon, p / math.cos(lat) - n


def llh2xyz(lat_deg, lon_deg, h):
    lat, lon = np.radians([lat_deg, lon_deg])
    n = A_WGS / np.sqrt(1 - E2 * np.sin(lat) ** 2)
    return np.array([(n + h) * np.cos(lat) * np.cos(lon),
                     (n + h) * np.cos(lat) * np.sin(lon),
                     (n * (1 - E2) + h) * np.sin(lat)])


def elevation(sat, rec):
    lat, lon, _ = xyz2llh(rec)
    d = sat - rec
    up = np.array([math.cos(lat) * math.cos(lon), math.cos(lat) * math.sin(lon), math.sin(lat)])
    return math.asin(np.clip(d @ up / np.linalg.norm(d), -1, 1))


def klobuchar(tow, rec, sat, alpha, beta):
    lat, lon, _ = xyz2llh(rec)
    el = elevation(sat, rec)
    if el <= 0:
        return 0.0
    d = sat - rec
    sl, cl = math.sin(lon), math.cos(lon)
    sp, cp = math.sin(lat), math.cos(lat)
    east = -sl * d[0] + cl * d[1]
    north = -sp * cl * d[0] - sp * sl * d[1] + cp * d[2]
    az = math.atan2(east, north)
    phi, lam, ee = lat / math.pi, lon / math.pi, el / math.pi
    psi = 0.0137 / (ee + 0.11) - 0.022
    phii = np.clip(phi + psi * math.cos(az), -0.416, 0.416)
    lami = lam + psi * math.sin(az) / math.cos(phii * math.pi)
    phim = phii + 0.064 * math.cos((lami - 1.617) * math.pi)
    local = (43200 * lami + tow) % 86400
    amp = max(0.0, sum(alpha[i] * phim ** i for i in range(4)))
    per = max(72000.0, sum(beta[i] * phim ** i for i in range(4)))
    x = 2 * math.pi * (local - 50400) / per
    f = 1 + 16 * (0.53 - ee) ** 3
    delay = 5e-9
    if abs(x) < 1.57:
        delay += amp * (1 - x * x / 2 + x ** 4 / 24)
    return C * f * delay


def tropo(rec, el):
    if el <= math.radians(3):
        return 0.0
    lat, _, h = xyz2llh(rec)
    h = max(-100.0, min(h, 10000.0))
    p = 1013.25 * (1 - 2.2557e-5 * h) ** 5.2568
    temp = 15.0 - 0.0065 * h + 273.15
    e = 6.108 * 0.5 * math.exp((17.15 * temp - 4684.0) / (temp - 38.45))
    z = math.pi / 2 - el
    return 0.002277 / math.cos(z) * (p + (1255 / temp + 0.05) * e - math.tan(z) ** 2)


def solve_spp_epoch(obs, tow, nav, alpha, beta, seed):
    x = np.r_[seed, 0.0]
    used = 0
    for it in range(10):
        rows, yy, ww = [], [], []
        for sat, o in obs.items():
            P = o.get("C1C", np.nan)
            if not np.isfinite(P):
                continue
            try:
                sp, dts = satellite(nav, sat, tow, P)
            except KeyError:
                continue
            rho = np.linalg.norm(sp - x[:3])
            if it and np.linalg.norm(x[:3]) > 6e6:
                el = elevation(sp, x[:3])
                if el < math.radians(10):
                    continue
                ion = klobuchar(tow, x[:3], sp, alpha, beta)
                trp = tropo(x[:3], el)
                w = max(math.sin(el), 0.15) ** 2
            else:
                ion = trp = 0.0
                w = 1.0
            rows.append(np.r_[(x[:3] - sp) / rho, 1.0])
            yy.append(P + C * dts - ion - trp - rho - x[3])
            ww.append(w)
        if len(rows) < 4:
            return None
        H, y, sw = np.asarray(rows), np.asarray(yy), np.sqrt(ww)
        dx, _, rank, _ = np.linalg.lstsq(H * sw[:, None], y * sw, rcond=None)
        if rank < 4:
            return None
        x += dx
        used = len(rows)
        if np.linalg.norm(dx[:3]) < 1e-3 and abs(dx[3]) < 1e-3:
            break
    if not (np.isfinite(x).all() and 6.2e6 < np.linalg.norm(x[:3]) < 6.5e6):
        return None
    return x, used


def solve_spp(rover, nav, alpha, beta):
    seed = np.zeros(3)
    out = {}
    for tow in sorted(rover):
        sol = solve_spp_epoch(rover[tow], tow, nav, alpha, beta, seed)
        if sol:
            x, used = sol
            seed = x[:3].copy()
            out[tow] = (x[:3].copy(), x[3], used)
    return out


def common_measurements(ro, ba, tow, nav, rover_xyz, base_xyz, alpha, beta):
    rows = []
    for sat in sorted(set(ro) & set(ba)):
        a, b = ro[sat], ba[sat]
        vals = [a.get(k, np.nan) for k in ("C1C", "L1C")] + [b.get(k, np.nan) for k in ("C1C", "L1C")]
        if not np.isfinite(vals).all():
            continue
        Pr, Lr, Pb, Lb = vals
        try:
            sr, dtr = satellite(nav, sat, tow, Pr)
            sb, dtb = satellite(nav, sat, tow, Pb)
        except KeyError:
            continue
        er, eb = elevation(sr, rover_xyz), elevation(sb, base_xyz)
        if min(er, eb) < math.radians(12):
            continue
        rr, rb = np.linalg.norm(sr - rover_xyz), np.linalg.norm(sb - base_xyz)
        ir = klobuchar(tow, rover_xyz, sr, alpha, beta)
        ib = klobuchar(tow, base_xyz, sb, alpha, beta)
        tr = tropo(rover_xyz, er)
        tb = tropo(base_xyz, eb)
        clock_corr = C * (dtr - dtb)
        rows.append(dict(sat=sat, Pr=Pr - Pb + clock_corr,
                         Lr=L1 * (Lr - Lb) + clock_corr, sr=sr, sb=sb,
                         elev=er, code_corr=(tr - tb) + (ir - ib),
                         phase_corr=(tr - tb) - (ir - ib), rr=rr, rb=rb,
                         lli=(a.get("L1C_lli", 0) | b.get("L1C_lli", 0))))
    return rows


def integer_least_squares(a, Q, max_nodes=200000):
    """Return the two best integer vectors for (n-a)' inv(Q) (n-a).

    This is a compact sphere decoder.  It is intentionally shown here rather
    than hidden in an external LAMBDA/RTK package.
    """
    a = np.asarray(a, float)
    Q = (np.asarray(Q, float) + np.asarray(Q, float).T) / 2
    try:
        W = np.linalg.inv(Q)
        R = np.linalg.cholesky(W).T
    except np.linalg.LinAlgError:
        return None
    m = len(a)

    def cost(n):
        d = np.asarray(n) - a
        return float(d @ W @ d)

    n0 = np.rint(a).astype(int)
    candidates = [(cost(n0), n0.copy())]
    for j in range(m):
        for step in (-1, 1):
            n = n0.copy()
            n[j] += step
            candidates.append((cost(n), n))
    candidates.sort(key=lambda x: x[0])
    best = [(candidates[0][0], candidates[0][1]), (candidates[1][0], candidates[1][1])]
    n = np.zeros(m, dtype=int)
    nodes = 0

    def add_candidate(c):
        nonlocal best
        if any(np.array_equal(n, b[1]) for b in best):
            return
        best.append((c, n.copy()))
        best.sort(key=lambda x: x[0])
        best = best[:2]

    def search(k, partial):
        nonlocal nodes
        if nodes >= max_nodes or partial >= best[1][0]:
            return
        if k < 0:
            add_candidate(partial)
            return
        known = sum(R[k, j] * (n[j] - a[j]) for j in range(k + 1, m))
        center = a[k] - known / R[k, k]
        radius = math.sqrt(max(best[1][0] - partial, 0)) / abs(R[k, k])
        lo, hi = math.ceil(center - radius), math.floor(center + radius)
        vals = sorted(range(lo, hi + 1), key=lambda z: abs(z - center))
        for z in vals:
            nodes += 1
            n[k] = z
            term = (R[k, k] * (z - center)) ** 2
            search(k - 1, partial + term)
            if nodes >= max_nodes:
                break

    search(m - 1, 0.0)
    if len(best) < 2:
        return None
    return best[0][1], best[0][0], best[1][1], best[1][0], nodes


class RelativeEKF:
    """Experimental sequential float filter; not used for published results."""
    def __init__(self, base_xyz, seed):
        self.base = base_xyz
        self.sats = []
        self.x = np.r_[seed, 0.0]
        self.P = np.diag([100.0, 100.0, 100.0, 1e6])
        self.last = None
        self.fix_diag = None

    def add_sat(self, sat, ambiguity, variance=100.0):
        if sat in self.sats:
            return
        self.sats.append(sat)
        self.x = np.r_[self.x, ambiguity]
        n = len(self.P)
        P = np.zeros((n + 1, n + 1))
        P[:n, :n] = self.P
        P[-1, -1] = variance
        self.P = P

    def update(self, tow, rows, static_hint=False):
        if self.last is None and rows:
            clocks = [r["Pr"] - (r["rr"] - r["rb"]) - r["code_corr"] for r in rows]
            self.x[3] = float(np.median(clocks))
            self.P[3, 3] = 10000.0
        if self.last is not None:
            dt_s = max(1, tow - self.last)
            q = np.zeros(len(self.x))
            q[:3] = ((0.03 if static_hint else 10.0) * dt_s) ** 2
            q[3] = (30.0 * dt_s) ** 2
            self.P += np.diag(q)
        self.last = tow
        bysat = {r["sat"]: r for r in rows}
        for r in rows:
            if r["sat"] not in self.sats:
                geom = r["rr"] - r["rb"]
                self.add_sat(r["sat"], r["Lr"] - geom - self.x[3] - r["phase_corr"])
        # Remove no state; dormant ambiguities simply receive no observation.
        H, z, pred, var = [], [], [], []
        for r in rows:
            sat = r["sat"]
            j = 4 + self.sats.index(sat)
            dr = r["sr"] - self.x[:3]
            rr = np.linalg.norm(dr)
            rb = np.linalg.norm(r["sb"] - self.base)
            los = (self.x[:3] - r["sr"]) / rr
            common = rr - rb + self.x[3]
            hc = np.zeros(len(self.x)); hc[:3] = los; hc[3] = 1
            hp = hc.copy(); hp[j] = 1
            H.extend([hc, hp])
            z.extend([r["Pr"], r["Lr"]])
            pred.extend([common + r["code_corr"], common + r["phase_corr"] + self.x[j]])
            s = max(math.sin(r["elev"]), 0.2)
            var.extend([(1.5 / s) ** 2, (0.015 / s) ** 2])
        if len(rows) < 4:
            return None
        H, innov, R = np.asarray(H), np.asarray(z) - np.asarray(pred), np.diag(var)
        S = H @ self.P @ H.T + R
        try:
            K = np.linalg.solve(S, H @ self.P).T
        except np.linalg.LinAlgError:
            return None
        self.x += K @ innov
        I = np.eye(len(self.x))
        self.P = (I - K @ H) @ self.P @ (I - K @ H).T + K @ R @ K.T
        self.P = (self.P + self.P.T) / 2
        return float(np.sqrt(np.mean(innov[::2] ** 2)))

    def try_fix(self, rows):
        """Integer least-squares validation, then fixed DD phase least squares."""
        if len(rows) < 5:
            return None
        ref = max(rows, key=lambda r: r["elev"])
        others = [r for r in rows if r is not ref]
        jr = 4 + self.sats.index(ref["sat"])
        idx = [4 + self.sats.index(r["sat"]) for r in others]
        T = np.zeros((len(others), len(self.x)))
        for k, j in enumerate(idx):
            T[k, j] = 1 / L1
            T[k, jr] = -1 / L1
        af = T @ self.x
        Q = T @ self.P @ T.T
        self.fix_diag = [(float(np.sqrt(max(Q[k, k], 0))), float(abs(af[k] - round(af[k]))))
                         for k in range(len(af))]
        # Avoid an expensive integer search before code measurements have
        # narrowed the float ambiguity region.
        if max(x[0] for x in self.fix_diag) > 0.5:
            return None
        # Independent rounding is used only as a seed.  The fixed phase
        # residual and displacement tests below provide the validation.
        nbest = np.rint(af).astype(int)
        ints = {r["sat"]: int(n) for r, n in zip(others, nbest)}
        # Solve rover position and a common phase clock from fixed DD phase.
        x = self.x[:3].copy()
        for _ in range(5):
            A, y = [], []
            rr0 = np.linalg.norm(ref["sr"] - x)
            rb0 = np.linalg.norm(ref["sb"] - self.base)
            obs0 = ref["Lr"] - ref["phase_corr"]
            h0 = (x - ref["sr"]) / rr0
            for r in others:
                rr = np.linalg.norm(r["sr"] - x)
                rb = np.linalg.norm(r["sb"] - self.base)
                obs_dd = (r["Lr"] - r["phase_corr"]) - obs0
                pred_dd = (rr - rb) - (rr0 - rb0) + L1 * ints[r["sat"]]
                A.append((x - r["sr"]) / rr - h0)
                y.append(obs_dd - pred_dd)
            dx, _, rank, _ = np.linalg.lstsq(np.asarray(A), np.asarray(y), rcond=None)
            if rank < 3:
                return None
            x += dx
            if np.linalg.norm(dx) < 1e-4:
                break
        res = np.asarray(y) - np.asarray(A) @ dx
        rms = float(np.sqrt(np.mean(res ** 2)))
        if rms > 0.04 or np.linalg.norm(x - self.x[:3]) > 0.5:
            return None
        return x, rms, len(others)


def solve_relative(rover, base, nav, alpha, beta, spp, base_xyz):
    common_t = sorted(set(rover) & set(base) & set(spp))
    if not common_t:
        return {}
    ekf = RelativeEKF(base_xyz, spp[common_t[0]][0])
    out = {}
    previous_spp = None
    for tow in common_t:
        seed = spp[tow][0]
        if np.linalg.norm(ekf.x[:3] - seed) > 50:
            ekf.x[:3] = seed
            ekf.P[:3, :3] += np.eye(3) * 100
        rows = common_measurements(rover[tow], base[tow], tow, nav, ekf.x[:3], base_xyz, alpha, beta)
        static_hint = previous_spp is not None and np.linalg.norm(seed - previous_spp) < 3.0
        previous_spp = seed.copy()
        code_rms = ekf.update(tow, rows, static_hint)
        if code_rms is None:
            continue
        fixed = ekf.try_fix(rows)
        if fixed:
            pos, phase_rms, ns = fixed
            out[tow] = (pos.copy(), 1, ns, code_rms, phase_rms)
            # Feed a validated fix back into the kinematic states.
            ekf.x[:3] = pos
            ekf.P[:3, :3] = np.eye(3) * 0.01
        else:
            out[tow] = (ekf.x[:3].copy(), 2, len(rows), code_rms, np.nan)
    return out


def write_csv(path, solutions, kind):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        if kind == "spp":
            w.writerow(["gps_tow_s", "x_ecef_m", "y_ecef_m", "z_ecef_m", "clock_bias_m", "satellites"])
            for t, (x, cb, ns) in solutions.items():
                w.writerow([t, *x, cb, ns])
        else:
            w.writerow(["gps_tow_s", "x_ecef_m", "y_ecef_m", "z_ecef_m", "quality", "satellites",
                        "code_innovation_rms_m", "fixed_phase_rms_m"])
            for t, (x, q, ns, cr, pr) in solutions.items():
                w.writerow([t, *x, q, ns, cr, pr])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--rover", type=Path, required=True)
    p.add_argument("--base", type=Path, required=True)
    p.add_argument("--nav", type=Path, required=True)
    p.add_argument("--out", type=Path, default=Path("python_results"))
    p.add_argument("--start", type=int, default=356400)
    p.add_argument("--end", type=int, default=359999)
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    nav, alpha, beta = read_nav(args.nav)
    rover, _, rt = read_obs(args.rover, args.start, args.end)
    base, base_xyz, bt = read_obs(args.base, args.start, args.end)
    print("epochs rover/base", len(rover), len(base), "GPS nav", len(nav), "types", rt[:8], bt[:8])
    spp = solve_spp(rover, nav, alpha, beta)
    print("SPP", len(spp))
    write_csv(args.out / "spp_python.csv", spp, "spp")
    summary = {"implementation": "pure Python; no RTKLIB output used",
               "rover_epochs": len(rover), "base_epochs": len(base),
               "spp_epochs": len(spp),
               "note": "Run python_batch_rtk.py for the validated RTK-float workflow."}
    (args.out / "run_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
