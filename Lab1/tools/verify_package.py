#!/usr/bin/env python3
"""Reproduce historical numbers, exercise the GPS solver and validate CSV failures."""
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import tempfile

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import WLSE_spp as spp
from tools.check_submission import check


def main():
    arrays, epochs = spp.read_inputs(ROOT/'example')
    rows, omitted = spp.run(arrays, epochs, 'historical')
    expected = np.genfromtxt(ROOT/'example/reference_results/estimated_positions.csv',
                             delimiter=',', names=True)
    reference = np.column_stack([expected[k] for k in ['ECEF_X_m','ECEF_Y_m','ECEF_Z_m']])
    actual = np.array([r['position'] for r in rows])
    max_diff = float(np.max(np.abs(actual-reference)))
    assert len(rows)==186 and not omitted and max_diff < 1e-5, max_diff
    metrics = spp.reference_metrics(rows,ROOT/'example/NAV-HPPOSECEF.csv',epochs)
    assert abs(metrics['rmse_2d_m']-37.9107584162)<1e-6
    assert abs(metrics['rmse_3d_m']-51.0428443069)<1e-6
    print('PASS historical: 186 epochs; max ECEF delta',max_diff)
    # Independent synthetic GPS geometry with a known clock and Sagnac term.
    truth=np.array([-2418200.,5385800.,2405700.]); bias=120000.
    sp=np.array([[15600000,7540000,20140000],[-18760000,2750000,18610000],
                 [17610000,14630000,13480000],[19170000,610000,-18390000],
                 [-17800000,-12300000,15400000],[1200000,24500000,8900000],
                 [-22000000,11000000,-9300000],[21000000,-12400000,-8800000]],dtype=float)
    ranges=np.linalg.norm(sp-truth,axis=1)
    sag=7.2921151467e-5/299792458.*(sp[:,0]*truth[1]-sp[:,1]*truth[0])
    z=np.zeros(len(sp))
    pos,clock=spp.gps_solution(sp,ranges+sag+bias,z,z,z,z,np.zeros(3))
    assert np.linalg.norm(pos-truth)<1e-3 and abs(clock-bias)<1e-3
    print('PASS synthetic GPS: position and clock recovered within 1 mm')
    # A missing epoch must be omitted, not filled with the previous solution.
    damaged={k:v.copy() for k,v in arrays.items()}
    damaged['pseudoranges_meas'][:,10]=np.nan
    reduced,omitted=spp.run(damaged,epochs,'historical')
    assert len(reduced)==185 and len(omitted)==1
    assert not any(r['gps_tow_s']==epochs[10][1] for r in reduced)
    print('PASS missing epoch: omitted, not duplicated or forward-filled')
    out=ROOT/'outputs/verification';out.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(dir=out) as d:
        d=Path(d)
        spp.write_results(d,'G01_owndata.csv',rows,[],epochs,'historical')
        errors,n,cov=check(d/'G01_owndata.csv',epoch_file=ROOT/'example/epochs.csv')
        assert not errors and n==186 and cov==100
        original=(d/'G01_owndata.csv').read_text()
        cases={'duplicate':original+original.splitlines()[1]+'\n',
               'nan':original.replace(original.splitlines()[1].split(',')[2],'nan',1),
               'extra_column':original.replace(spp.HEADER[-1],spp.HEADER[-1]+',extra',1),
               'wrong_unit':original.replace(original.splitlines()[1].split(',')[2],'1.000',1)}
        for label,contents in cases.items():
            target=d/'G02_benchmark.csv';target.write_text(contents,encoding='utf-8')
            assert check(target)[0],label
        empty=d/'G03_benchmark.csv'
        empty.write_text(','.join(spp.HEADER)+'\n',encoding='utf-8')
        assert check(empty)[0], 'Header-only file must be rejected'
        reversed_file=d/'G04_benchmark.csv'
        lines=original.splitlines()
        reversed_file.write_text(lines[0]+'\n'+'\n'.join(reversed(lines[1:]))+'\n',encoding='utf-8')
        assert check(reversed_file)[0], 'Unsorted times must be rejected'
        blank_file=d/'G05_owndata.csv'
        blank_file.write_text(original+'\n',encoding='utf-8')
        assert not check(blank_file)[0], 'Server-compatible blank trailing lines are accepted'
    print('PASS checker: valid file accepted; duplicates/NaN/header/unit/empty/order errors rejected')
    benchmark=ROOT/'benchmark/gps_csv'
    if benchmark.exists():
        provenance=json.loads((benchmark/'manifest.json').read_text(encoding='utf-8'))
        for name,digest in provenance['input_sha256'].items():
            assert hashlib.sha256((ROOT/'benchmark'/name).read_bytes()).hexdigest()==digest
        for item in provenance['files']:
            assert hashlib.sha256((benchmark/item['path']).read_bytes()).hexdigest()==item['sha256']
        assert hashlib.sha256((ROOT/'rinex2csv.m').read_bytes()).hexdigest()==provenance['converter_sha256']
        b,e=spp.read_inputs(benchmark)
        solved,missing=spp.run(b,e,'gps')
        assert len(e)==1534 and len(solved)>=.9*len(e)
        report=spp.write_results(out/'benchmark','G01_benchmark.csv',solved,missing,e,'gps')
        errors,n,cov=check(out/'benchmark/G01_benchmark.csv',2158,1534,benchmark/'epochs.csv')
        assert not errors and cov>=90
        print('PASS benchmark:',len(solved),'/',len(e),'actual OBS epochs, no reference trajectory used')
    else:
        raise AssertionError('Benchmark preprocessing has not yet been generated')
    manifest=ROOT/'docs/file_manifest.json'
    if manifest.exists():
        entries=json.loads(manifest.read_text(encoding='utf-8'))['files']
        for item in entries:
            p=ROOT/item['path']
            assert p.is_file(),str(p)
            assert hashlib.sha256(p.read_bytes()).hexdigest()==item['sha256'],str(p)
        print('PASS package hashes:',len(entries),'files')
    print('ALL CHECKS PASSED')


if __name__=='__main__':
    main()
