#!/usr/bin/env python3
"""Preserved-checkpoint drift audit; never promotes a point to equilibrated."""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
from datetime import datetime, timezone
import hashlib,json,os,shutil,subprocess,threading
from pathlib import Path
import h5py
import numpy as np
from run_disorder_ensemble import _write_json
REPO=Path(__file__).resolve().parents[1]
CP='simulation/realizations/0/clones/0/checkpoint'
OBS=['Equal_time_susceptibility','Winding_number_squared_total','Total_Energy','Number_of_particles']
def now():return datetime.now(timezone.utc).isoformat()
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def saved(p):
    with h5py.File(p) as f:return int(f[CP+'/sweeps'][()])
def patch(source,dest,thermalization,production,seconds):
    """Change only run controls; retain configuration, map, RNG and statistics."""
    if dest.exists():raise ValueError(f'Already exists: {dest}')
    shutil.copy2(source,dest)
    with h5py.File(dest,'r+') as f:
        p=f['parameters'];keys=list(p.attrs['ini_keys']);values=np.array(p.attrs['ini_values'],copy=True)
        for k,v in dict(thermalization=thermalization,sweeps=production,runtimelimit=seconds).items():
            p[k][()]=v;values[keys.index(k)]=str(v)
        p.attrs.modify('ini_values',values)
        assert p['checkpoint'][()].decode()=='state.clone.h5'
        assert p['outputfile'][()].decode()=='result.out.h5'
def controls(stage,threshold,window):return threshold,10 if stage=='discard' else window
def statistics(path):
    result={}
    with h5py.File(path) as f:
        for name in OBS:
            g=f['simulation/results/'+name]
            count=int(g['count'][()]);mean=np.asarray(g['mean/value'][()]);error=np.asarray(g['mean/error'][()])
            if count<=0 or not np.all(np.isfinite(mean)) or not np.all(np.isfinite(error)):
                raise ValueError('Missing/nonfinite '+name)
            result[name]=dict(count=count,mean=mean.tolist(),error=error.tolist())
    return result
def slots(main_plan):
    # Increase only after the entire original pipeline finished successfully.
    return 8 if json.loads(Path(main_plan).read_text())['status']=='completed' else 2

def prepare(manifest,root,main_plan,discard=200000,window=200000,seconds=3600):
    manifest=Path(manifest).resolve(strict=True);root=Path(root).resolve();m=json.loads(manifest.read_text())
    assert m['status']=='completed' and len(m['runs'])==8
    assert not root.exists(),'Refusing existing audit'
    assert discard>0 and window>=100 and seconds>0
    root.mkdir(parents=True)
    state=dict(schema_version=1,kind='same_map_equilibration_drift_audit',status='prepared',created_utc=now(),
        source_manifest=str(manifest),source_manifest_sha256=sha(manifest),main_plan=str(Path(main_plan).resolve(strict=True)),
        executable=m['executable'],executable_sha256=sha(m['executable']),discard=discard,window=window,chunk_seconds=seconds,
        max_chunks_per_stage=48,interpretation='PROVISIONAL: absence of drift is not proof of equilibration; adjacent windows are correlated and independent starts remain desirable.',runs=[])
    for run in m['runs']:
        src=manifest.parent/run['checkpoint_file'];case=root/run['id'];case.mkdir()
        snapshot=case/'source.clone.h5';shutil.copy2(src,snapshot)
        assert sha(src)==sha(snapshot)
        with h5py.File(snapshot) as f:
            assert float(f['parameters/beta'][()])==.6 and int(f['parameters/Lx'][()])==128
            assert int(f['parameters/canonical'][()])==30000
            pot=np.asarray(f[CP+'/hamiltonian/onsite_potential'][()],dtype='<f8')
            assert hashlib.sha256(pot.tobytes()).hexdigest()==run['potential_sha256']
        state['runs'].append(dict(id=run['id'],status='prepared',source=str(src),snapshot=str(snapshot.relative_to(root)),
            source_sha256=sha(snapshot),source_sweeps=saved(snapshot),potential_sha256=run['potential_sha256'],stages=[]))
    _write_json(root/'audit.json',state)
    return root/'audit.json'

def run_audit(path):
    path=Path(path).resolve(strict=True);root=path.parent;state=json.loads(path.read_text())
    assert state['status']=='prepared','Inspect running/failed audit before restarting; never reset a partial window automatically'
    assert sha(state['executable'])==state['executable_sha256']
    lock=threading.RLock()
    def update():
        with lock:_write_json(path,state)
    def worker(entry):
        source=root/entry['snapshot'];assert sha(source)==entry['source_sha256']
        threshold=entry['source_sweeps']+state['discard'];entry['status']='running';update()
        try:
            for stage in ['discard','window_A','window_B']:
                thermal,production=controls(stage,threshold,state['window']);target=thermal+production
                record=dict(name=stage,status='running',target_sweeps=target,chunks=[])
                with lock:entry['stages'].append(record);update()
                for index in range(state['max_chunks_per_stage']):
                    case=root/entry['id']/stage/f'chunk{index:04d}';case.mkdir(parents=True);(case/'tmp').mkdir()
                    start=case/'start.clone.h5';patch(source,start,thermal,production,state['chunk_seconds'])
                    shutil.copy2(start,case/'state.clone.h5')
                    reset=stage!='discard' and index==0
                    command=[state['executable'],str(start)]+(['--reset-statistics'] if reset else [])
                    chunk=dict(index=index,status='running',started_utc=now(),reset_statistics=reset,
                        checkpoint=str((case/'state.clone.h5').relative_to(root)),output=str((case/'result.out.h5').relative_to(root)),
                        input_sha256=sha(start),source_sweeps=saved(source))
                    env=dict(os.environ,OMP_NUM_THREADS='1',TMPDIR=str(case/'tmp'),TMP=str(case/'tmp'),TEMP=str(case/'tmp'))
                    env['PATH']=str(Path(state['executable']).parent)+os.pathsep+env.get('PATH','')
                    with (case/'run.log').open('w',encoding='utf-8') as log:
                        proc=subprocess.Popen(command,cwd=case,env=env,stdout=log,stderr=subprocess.STDOUT,
                            creationflags=subprocess.CREATE_NO_WINDOW if os.name=='nt' else 0)
                        with lock:chunk['pid']=proc.pid;record['chunks'].append(chunk);update()
                        rc=proc.wait()
                    source=case/'state.clone.h5';end=saved(source)
                    with lock:chunk.update(status='completed' if rc==0 else 'failed',returncode=rc,finished_utc=now(),saved_sweeps=end);update()
                    if rc:raise RuntimeError(f'{stage} chunk {index} exit {rc}')
                    # Reset rewinds the counter to thermalization once, so compare against that baseline.
                    baseline=thermal if reset else chunk['source_sweeps']
                    if end<=baseline:raise RuntimeError('No closed-excursion progress; inspect log')
                    if end>=target:
                        if stage!='discard':
                            stats=statistics(case/'result.out.h5')
                            assert stats['Total_Energy']['count']>=state['window']
                            assert stats['Equal_time_susceptibility']['count']>=state['window']//100
                            assert abs(stats['Number_of_particles']['mean']-30000)<1e-6
                            record['statistics']=stats
                        with lock:record.update(status='completed',finished_utc=now(),final_checkpoint=str(source.relative_to(root)));update()
                        break
                else:raise RuntimeError(f'{stage} reached chunk guard; preserve state and inspect before extending')
            entry.update(status='completed',finished_utc=now());update()
        except Exception as exc:
            with lock:entry.update(status='failed',error=str(exc),finished_utc=now());update()
            raise
    state.update(status='running',controller_pid=os.getpid(),started_utc=now());update()
    with ThreadPoolExecutor(max_workers=8) as pool:
        active={};pending=iter(state['runs']);exhausted=False
        while active or not exhausted:
            limit=slots(state['main_plan'])
            state['worker_limit']=limit
            while not exhausted and len(active)<limit:
                entry=next(pending,None)
                if entry is None:exhausted=True;break
                active[pool.submit(worker,entry)]=entry['id']
            update()
            if active:
                finished,_=wait(active,timeout=30,return_when=FIRST_COMPLETED)
                for f in finished:
                    ident=active.pop(f)
                    try:f.result();print(ident,'audit windows complete',flush=True)
                    except Exception as exc:print(ident,'FAILED',exc,flush=True)
    state.update(status='completed' if all(r['status']=='completed' for r in state['runs']) else 'failed',finished_utc=now());update()
    if state['status']=='completed':
        comparisons=[]
        for r in state['runs']:
            a,b=r['stages'][1]['statistics'],r['stages'][2]['statistics']
            comparisons.append(dict(id=r['id'],window_A=a,window_B=b,
                B_minus_A={k:np.asarray(b[k]['mean']).item()-np.asarray(a[k]['mean']).item() for k in OBS}))
        _write_json(root/'window_comparison.json',dict(interpretation=state['interpretation'],maps=comparisons))
    return 0 if state['status']=='completed' else 1

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('command',choices=['prepare','run','status'])
    p.add_argument('--audit',type=Path,required=True);p.add_argument('--manifest',type=Path);p.add_argument('--main-plan',type=Path)
    a=p.parse_args()
    if a.command=='prepare':print(prepare(a.manifest,a.audit.parent,a.main_plan));return 0
    if a.command=='run':return run_audit(a.audit)
    s=json.loads(a.audit.read_text());print(s['status'],'workers',s.get('worker_limit'))
    for r in s['runs']:
        print(r['id'],r['status'],[(t['name'],t['status'],t['chunks'][-1] if t['chunks'] else None) for t in r['stages']])
    return 0
if __name__=='__main__':raise SystemExit(main())
