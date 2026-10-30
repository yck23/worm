#!/usr/bin/env python3
"""Continue incomplete disorder scans from preserved checkpoints."""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
import copy
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import h5py
import numpy as np
from run_disorder_ensemble import _write_json, _run_realization
from disorder_response_scan import point_rows, analyze

REPO=Path(__file__).resolve().parents[1]

def stamp():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def patch_checkpoint(source, destination, seconds):
    """Copy state/RNG/statistics; change only the time limit in both ALPS forms."""
    if destination.exists():
        raise ValueError(f'Refusing existing checkpoint: {destination}')
    shutil.copy2(source,destination)
    with h5py.File(destination,'r+') as f:
        parameters=f['parameters']
        assert parameters['checkpoint'][()].decode()=='state.clone.h5'
        assert parameters['outputfile'][()].decode()=='result.out.h5'
        parameters['runtimelimit'][()]=seconds
        keys=list(parameters.attrs['ini_keys'])
        values=np.array(parameters.attrs['ini_values'],copy=True)
        values[keys.index('runtimelimit')]=str(seconds)
        parameters.attrs.modify('ini_values',values)
        assert int(parameters['thermalization'][()])==50000


def missing_production(manifest_path, run):
    """Detect absent production; malformed existing output must be inspected."""
    output=Path(manifest_path).parent/run['output_file']
    if not output.exists():
        return True
    with h5py.File(output,'r') as f:
        names=['Equal_time_susceptibility','Winding_number_squared',
               'Winding_number_squared_total','Number_of_particles']
        return any(f'simulation/results/{name}/count' not in f or
                   int(f[f'simulation/results/{name}/count'][()]) <= 0 for name in names)


def prepare(scan_file,seconds,workers,only_missing_measurements=False):
    path=Path(scan_file).resolve(strict=True);root=path.parent
    state=json.loads(path.read_text())
    if state['status']=='running':raise ValueError('Refusing to modify a running scan')
    if seconds <= 0 or workers <= 0:raise ValueError('seconds and workers must be positive')
    backup=root/'scan.before_continuation.json'
    if backup.exists():
        raise ValueError('Continuation was already prepared; inspect existing state')
    assert sha(state['executable'])==state['executable_sha256']
    selected=[p for p in state['points'] if p['status']!='completed']
    if not selected:raise ValueError('No incomplete points to continue')
    continued_runs=0
    prepared=[]
    for point in selected:
        manifest_path=root/point['manifest']
        original=json.loads(manifest_path.read_text())
        manifest=copy.deepcopy(original)
        selected_ids=[]
        for run in manifest['runs']:
            if only_missing_measurements and not missing_production(manifest_path,run):
                if run['status']!='completed' or run.get('returncode')!=0:
                    raise ValueError(f'Inspect process failure for {point["id"]}/{run["id"]} before recovery')
                continue
            selected_ids.append(run['id'])
            continued_runs+=1
            source=manifest_path.parent/run['checkpoint_file']
            assert source.is_file()
            relative=Path(run['directory'])/'continuation01'
            case=manifest_path.parent/relative
            if case.exists():raise ValueError(f'Refusing existing continuation: {case}')
            case.mkdir();(case/'tmp').mkdir()
            with h5py.File(source) as f:
                saved_sweeps=int(f['simulation/realizations/0/clones/0/checkpoint/sweeps'][()])
                assert float(f['parameters/beta'][()])==float(state['config']['beta'])
                assert int(f['parameters/canonical'][()])==state['config']['particles'][str(point['L'])]
                assert int(f['parameters/Lx'][()])==point['L']
                potential=np.asarray(f['simulation/realizations/0/clones/0/checkpoint/hamiltonian/onsite_potential'][()],dtype='<f8')
                assert hashlib.sha256(potential.tobytes()).hexdigest()==run['potential_sha256']
            patch_checkpoint(source,case/'start.clone.h5',seconds)
            shutil.copy2(case/'start.clone.h5',case/'state.clone.h5')
            run['continuation_from']=dict(checkpoint=str(source.relative_to(root)),
                checkpoint_sha256=sha(source),sweeps=saved_sweeps,statistics='retained')
            run.update(directory=relative.as_posix(),input_file=(relative/'start.clone.h5').as_posix(),
                input_sha256=sha(case/'start.clone.h5'),output_file=(relative/'result.out.h5').as_posix(),
                checkpoint_file=(relative/'state.clone.h5').as_posix(),log_file=(relative/'run.log').as_posix(),
                status='prepared',returncode=None)
            for key in ['started_utc','finished_utc','error']:run.pop(key,None)
        if not selected_ids:raise ValueError(f'{point["id"]} needs inspection; no missing production found')
        manifest.pop('finished_utc',None)
        manifest.update(status='prepared',continuation=dict(runtime_seconds=seconds,prepared_utc=stamp(),
            runs=selected_ids,
            note='Same maps and RNG chains; prior state and measurements retained; original artifacts preserved'))
        prepared.append((point,manifest_path,original,manifest))
        print('Prepared checkpoint copies:',point['id'],flush=True)
    _write_json(backup,state)
    for point,manifest_path,original,manifest in prepared:
        _write_json(manifest_path.with_name('manifest.before_continuation.json'),original)
        _write_json(manifest_path,manifest)
        point.update(status='prepared')
        point.pop('error',None)
        point.pop('finished_utc',None)
    state.pop('finished_utc',None)
    state.pop('controller_pid',None)
    state.update(status='prepared_continuation',continuation=dict(runtime_seconds=seconds,workers=workers,
        points=[p['id'] for p in selected],prepared_utc=stamp()))
    _write_json(path,state)
    print(f'Prepared {continued_runs} warm continuations; original results/checkpoints retained.',flush=True)


def refresh(path,state):
    try:
        analyze(path)
        if state['config'].get('refresh_paired_notebook',False):return
        with (path.parent/'paired_notebook_refresh.log').open('a',encoding='utf-8') as log:
            result=subprocess.run([sys.executable,str(REPO/'scripts/analyze_disorder_response.py')],
                cwd=REPO,stdout=log,stderr=subprocess.STDOUT,timeout=600,
                creationflags=subprocess.CREATE_NO_WINDOW if os.name=='nt' else 0)
        if result.returncode:raise RuntimeError('Paired notebook refresh failed; see paired_notebook_refresh.log')
    except Exception as exc:
        state['analysis_error']=str(exc)
        _write_json(path,state)
        print('Analysis refresh issue:',exc,flush=True)


def run(scan_file):
    path=Path(scan_file).resolve(strict=True);root=path.parent
    state=json.loads(path.read_text())
    assert state['status']=='prepared_continuation'
    assert sha(state['executable'])==state['executable_sha256']
    jobs=[];manifests={}
    for point in state['points']:
        if point['id'] not in state['continuation']['points']:continue
        manifest_path=root/point['manifest'];manifest=json.loads(manifest_path.read_text())
        manifests[point['id']]=(point,manifest_path,manifest)
        for entry in manifest['runs']:
            if entry['status']=='completed':continue
            assert entry['status']=='prepared'
            assert not (manifest_path.parent/entry['output_file']).exists()
            assert sha(manifest_path.parent/entry['input_file'])==entry['input_sha256']
            jobs.append((point,manifest_path,manifest,entry))
    state.update(status='running',controller_pid=os.getpid(),continuation_started_utc=stamp())
    _write_json(path,state)
    workers=state['continuation']['workers'];next_job=0
    with ThreadPoolExecutor(max_workers=workers) as pool:
        active={}
        while active or next_job<len(jobs):
            while next_job<len(jobs) and len(active)<workers:
                point,mp,manifest,entry=jobs[next_job]
                entry['status']='running';manifest['status']='running';point['status']='running'
                _write_json(mp,manifest)
                active[pool.submit(_run_realization,mp.parent,dict(entry),Path(state['executable']))]=(point,mp,manifest,entry)
                next_job+=1
            _write_json(path,state)
            completed,_=wait(active,return_when=FIRST_COMPLETED)
            for future in completed:
                point,mp,manifest,entry=active.pop(future)
                try:entry.update(future.result())
                except Exception as exc:entry.update(status='failed',returncode=None,error=str(exc),finished_utc=stamp())
                _write_json(mp,manifest)
                print(point['id'],entry['id'],entry['status'],flush=True)
                if all(r['status'] in ['completed','failed'] for r in manifest['runs']):
                    manifest['status']='completed' if all(r['status']=='completed' for r in manifest['runs']) else 'failed'
                    manifest['finished_utc']=stamp();_write_json(mp,manifest)
                    try:
                        if manifest['status']!='completed':raise RuntimeError('A continuation process failed')
                        point_rows(root,point)
                        point.update(status='completed',finished_utc=stamp());point.pop('error',None)
                    except Exception as exc:point.update(status='failed',error=str(exc))
                    _write_json(path,state)
                    refresh(path,state)
                _write_json(path,state)
    state.update(status='completed' if all(p['status']=='completed' for p in state['points']) else 'failed',finished_utc=stamp())
    _write_json(path,state);refresh(path,state)
    return 0 if state['status']=='completed' else 1


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=['prepare','run'])
    parser.add_argument('--scan',type=Path,required=True)
    parser.add_argument('--seconds',type=int,default=7200)
    parser.add_argument('--workers',type=int,default=16)
    parser.add_argument('--only-missing-measurements',action='store_true',
        help='Keep valid production runs; continue only chains with absent production')
    args=parser.parse_args()
    if args.command=='prepare':prepare(args.scan,args.seconds,args.workers,args.only_missing_measurements);return 0
    return run(args.scan)


if __name__=='__main__':raise SystemExit(main())
