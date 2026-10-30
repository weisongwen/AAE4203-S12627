#!/usr/bin/env python3
"""Offline course CSV checker. It does not calculate positioning accuracy.

Supports the original one-file CLI and optional OBS epoch checking. Basic
numeric formatting/order checks follow the current leaderboard contract;
course-specific filename and Earth-surface checks are retained.
"""
import argparse
import csv
import math
from pathlib import Path
import re

HEADER = ['gps_week', 'gps_tow_s', 'x_ecef_m', 'y_ecef_m', 'z_ecef_m']


def check(path, week=0, epochs=0, epoch_file=None):
    path = Path(path)
    errors = []
    if not re.fullmatch(r'G(0[1-9]|1[0-9]|20)_(benchmark|owndata)\.csv', path.name):
        errors.append('Use filename G01..G20_benchmark.csv or G01..G20_owndata.csv')
    if path.stat().st_size > 20*1024*1024:
        return ['File exceeds the 20 MiB submission limit'], 0, None
    with path.open(encoding='utf-8-sig', newline='') as stream:
        data = list(csv.reader(stream))
    if not data or data[0] != HEADER:
        return ['Header must be exactly '+','.join(HEADER)], 0, None
    previous, seen, millisecond_keys = None, set(), set()
    count = 0
    for line, row in enumerate(data[1:], 2):
        if not row or all(not c.strip() for c in row):
            continue
        count += 1
        if count > 200000:
            errors.append('More than 200000 data rows'); break
        row = [c.strip() for c in row]
        if len(row) != 5:
            errors.append(f'Line {line}: expected exactly five columns'); continue
        if not re.fullmatch(r'[0-9]+', row[0]):
            errors.append(f'Line {line}: GPS week must be an integer'); continue
        if not re.fullmatch(r'[0-9]+\.[0-9]+', row[1]):
            errors.append(f'Line {line}: time must be decimal seconds (no scientific notation)'); continue
        if any(not re.fullmatch(r'-?[0-9]+\.[0-9]{3}', x) for x in row[2:]):
            errors.append(f'Line {line}: ECEF coordinates must have exactly three decimals'); continue
        w, tow = int(row[0]), float(row[1])
        xyz = list(map(float, row[2:]))
        if not all(math.isfinite(v) for v in [tow, *xyz]) or not 0 <= tow < 604800:
            errors.append(f'Line {line}: invalid/non-finite values'); continue
        if not 6.30e6 < math.sqrt(sum(v*v for v in xyz)) < 6.42e6:
            errors.append(f'Line {line}: ECEF radius is outside the surface range')
        if week and w != week:
            errors.append(f'Line {line}: expected GPS week {week}, got {w}')
        # Server parsing rejects duplicates at millisecond precision. The
        # separate rounded 0.1 s key is only for optional OBS membership checks.
        exact_key = (w,tow)
        ms_key = (w,round(tow,3))
        if ms_key in millisecond_keys:
            errors.append(f'Line {line}: duplicate epoch')
        if previous is not None and exact_key <= previous:
            errors.append(f'Line {line}: epochs are not strictly increasing')
        millisecond_keys.add(ms_key)
        seen.add((w,round(tow,1)))
        previous = exact_key
    if not seen:
        errors.append('No valid data rows')
    coverage = 100*len(seen)/epochs if epochs else None
    if epoch_file:
        with epoch_file.open(encoding='utf-8-sig', newline='') as stream:
            expected = {(int(r['gps_week']),round(float(r['gps_tow_s']),1))
                        for r in csv.DictReader(stream)}
        unexpected = seen-expected
        if unexpected:
            errors.append(f'{len(unexpected)} epoch(s) are not present in the OBS epoch list')
        coverage = 100*len(seen & expected)/len(expected) if expected else 0
    return errors, count, coverage


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('files', nargs='+', type=Path)
    parser.add_argument('--week', type=int, default=0, help='Expected week; 0 disables this check')
    parser.add_argument('--epochs', type=int, default=0, help='Expected count; candidate coverage only')
    parser.add_argument('--epoch-file', type=Path, help='Actual OBS epoch list for exact matching')
    args = parser.parse_args()
    status = 0
    for path in args.files:
        try:
            errors, count, coverage = check(path,args.week,args.epochs,args.epoch_file)
        except (OSError, ValueError, csv.Error) as exc:
            errors, count, coverage = [str(exc)], 0, None
        print(path.name, ':', count, 'rows;', 'FORMAT OK' if not errors else 'FAILED')
        for e in errors[:20]: print(' ERROR:',e)
        if coverage is not None:
            print(f' OBS/count coverage: {coverage:.2f}% (ranking needs >=90%; no accuracy check)')
        status = max(status, bool(errors))
    return int(status)


if __name__ == '__main__':
    raise SystemExit(main())
