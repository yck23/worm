#!/usr/bin/env python3
"""Keep one existing QMC process and run prepared sweeps on the spare cores.

A ready/commit handoff retains a Windows process handle before the former
controller is stopped. The QMC process is never stopped or restarted.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import threading
import time

from run_disorder_ensemble import _write_json
from run_disorder_queue import REPO, read, local_path, scan_reference, process_alive, execute_stage, status
from disorder_response_scan import point_rows, analyze


def now():
    return datetime.now(timezone.utc).isoformat()


class ProcessHandle:
    def __init__(self, pid, expected_creation=None):
        import ctypes
        from ctypes import wintypes
        self.ctypes=ctypes
        self.kernel=ctypes.WinDLL('kernel32',use_last_error=True)
        k=self.kernel
        k.OpenProcess.argtypes=[wintypes.DWORD,wintypes.BOOL,wintypes.DWORD]
        k.OpenProcess.restype=wintypes.HANDLE
        k.GetExitCodeProcess.argtypes=[wintypes.HANDLE,ctypes.POINTER(wintypes.DWORD)]
        k.GetExitCodeProcess.restype=wintypes.BOOL
        k.GetProcessTimes.argtypes=[wintypes.HANDLE,*([ctypes.POINTER(wintypes.FILETIME)]*4)]
        k.GetProcessTimes.restype=wintypes.BOOL
        k.WaitForSingleObject.argtypes=[wintypes.HANDLE,wintypes.DWORD]
        k.WaitForSingleObject.restype=wintypes.DWORD
        k.CloseHandle.argtypes=[wintypes.HANDLE]
        k.CloseHandle.restype=wintypes.BOOL
        self.handle=k.OpenProcess(0x100000|0x1000,False,int(pid))
        if not self.handle:raise ctypes.WinError(ctypes.get_last_error())
        try:
            times=[wintypes.FILETIME() for _ in range(4)]
            if not k.GetProcessTimes(self.handle,*[ctypes.byref(t) for t in times]):
                raise ctypes.WinError(ctypes.get_last_error())
            self.creation=(times[0].dwHighDateTime<<32)|times[0].dwLowDateTime
            if expected_creation is not None and self.creation!=expected_creation:
                raise ValueError('PID creation time changed; refusing process adoption')
        except BaseException:
            self.close();raise

    def poll(self):
        from ctypes import wintypes
        code=wintypes.DWORD()
        if not self.kernel.GetExitCodeProcess(self.handle,self.ctypes.byref(code)):
            raise self.ctypes.WinError(self.ctypes.get_last_error())
        return None if code.value==259 else code.value

    def wait_tick(self, milliseconds=1000):
        result=self.kernel.WaitForSingleObject(self.handle,milliseconds)
        if result not in (0,258):raise self.ctypes.WinError(self.ctypes.get_last_error())

    def close(self):
        if self.handle:
            self.kernel.CloseHandle(self.handle);self.handle=None


def prepare(plan_file, old_queue, pid, scans, ensemble_workers=7):
    plan_file=Path(plan_file).resolve()
    if plan_file.exists():raise ValueError('Plan already exists')
    old_queue=Path(old_queue).resolve(strict=True);old=read(old_queue)
    if old['status']!='recovering_dependency' or not process_alive(old['controller_pid']):
        raise ValueError('Expected the live recovery controller')
    dependency=local_path(old['dependency']);scan=read(dependency)
    candidates=[]
    for point in scan['points']:
        mp=dependency.parent/point['manifest'];m=read(mp)
        for run in m['runs']:
            if run['status']=='running':candidates.append((point,mp,m,run))
            elif run['status']!='completed':raise ValueError('Expected all other dependency runs completed')
    if len(candidates)!=1:raise ValueError('Exactly one running QMC chain is supported')
    point,mp,manifest,run=candidates[0]
    if not 1<=ensemble_workers<=7:raise ValueError('Reserve a core for the adopted process')
    handle=ProcessHandle(pid)
    try:
        if handle.poll() is not None:raise ValueError('QMC process already exited')
        creation=handle.creation
    finally:handle.close()
    changes=[];stages=[]
    for path in scans:
        path=Path(path).resolve(strict=True);state=read(path)
        if state['status']!='prepared':raise ValueError('Only never-started sweeps can be scheduled')
        if state['config']['workers']//state['config']['realizations']!=2:
            raise ValueError('Expected two concurrent ensembles')
        if hashlib.sha256(Path(state['executable']).read_bytes()).hexdigest()!=state['executable_sha256']:
            raise ValueError('Executable changed')
        for entry in state['points']:
            target=path.parent/entry['manifest'];m=read(target)
            if m['status']!='prepared' or any(r['status']!='prepared' for r in m['runs']):
                raise ValueError('A future ensemble has already started')
            backup=target.with_name('manifest.before_parallel_20261102.json')
            if backup.exists():raise ValueError('Worker handoff already prepared')
            changes.append((target,backup,m))
        stages.append(dict(beta=state['config']['beta'],scan=scan_reference(path),status='queued'))
    if [s['scan'] for s in stages]!=[s['scan'] for s in old['stages']]:
        raise ValueError('Prepared stages must match the old queue order')
    for target,backup,m in changes:
        _write_json(backup,m)
        m['workers']=ensemble_workers
        m['execution_override']=dict(reason='Reserve one active recovery process; two ensembles use at most 14 workers',
                                     previous_workers=read(backup)['workers'],changed_utc=now())
        _write_json(target,m)
    plan_file.parent.mkdir(parents=True,exist_ok=True)
    plan=dict(schema_version=1,kind='disorder_scan_queue',status='prepared',created_utc=now(),
        old_queue=scan_reference(old_queue),old_controller_pid=old['controller_pid'],
        dependency=scan_reference(dependency),stages=stages,
        ensemble_workers=ensemble_workers,max_simultaneous_qmc=2*ensemble_workers+1,
        adopted=dict(pid=pid,creation_filetime=creation,point=point['id'],run=run['id'],
                     manifest=scan_reference(mp),input_file=run['input_file'],status='running'))
    _write_json(plan_file,plan)
    shutil.copyfile(__file__,plan_file.parent/'run_parallel_disorder_tail.py')
    print('Prepared parallel handoff; QMC process remains untouched.',flush=True)


def finalize_adopted(plan, code):
    path=local_path(plan['dependency']);scan=read(path)
    point=next(p for p in scan['points'] if p['id']==plan['adopted']['point'])
    mp=path.parent/point['manifest'];manifest=read(mp)
    run=next(r for r in manifest['runs'] if r['id']==plan['adopted']['run'])
    if run['status']!='running':raise ValueError('Adopted run metadata was modified by another owner')
    success=code==0 and (mp.parent/run['output_file']).is_file()
    run.update(status='completed' if success else 'failed',returncode=code,finished_utc=now(),
               adopted_process_pid=plan['adopted']['pid'])
    if not success:run['error']='Adopted QMC process failed or its output is missing'
    manifest.update(status='completed' if all(r['status']=='completed' for r in manifest['runs']) else 'failed',finished_utc=now())
    _write_json(mp,manifest)
    try:
        if manifest['status']!='completed':raise RuntimeError(run.get('error','Incomplete ensemble'))
        point_rows(path.parent,point)
        point.update(status='completed',finished_utc=now());point.pop('error',None)
    except Exception as exc:point.update(status='failed',error=str(exc),finished_utc=now())
    scan.update(status='completed' if all(p['status']=='completed' for p in scan['points']) else 'failed',
                controller_pid=os.getpid(),finished_utc=now())
    _write_json(path,scan)
    try:analyze(path)
    except Exception as exc:
        scan['analysis_error']=str(exc);_write_json(path,scan)
        print('Recovery analysis refresh failed:',exc,flush=True)
    return 0 if scan['status']=='completed' else 1


def run_sweeps(plan, save):
    codes=[]
    for stage in plan['stages']:
        path=local_path(stage['scan'])
        if read(path)['status']!='prepared':raise ValueError('Refusing to restart a scan')
        stage.update(status='running',started_utc=now());save()
        print(f'{now()} Starting beta={stage["beta"]} alongside final recovery',flush=True)
        code=execute_stage(path);state=read(path)
        good=code==0 and state['status']=='completed'
        stage.update(status='completed' if good else 'failed',returncode=code,finished_utc=now());save()
        codes.append(0 if good else 1)
        if state['status'] not in ('completed','failed'):
            raise RuntimeError('Sweep controller interrupted; stop before further launches')
    return max(codes,default=0)


def run(plan_file):
    plan_file=Path(plan_file).resolve(strict=True);lock=plan_file.with_suffix('.lock')
    with lock.open('x') as f:f.write(str(os.getpid()))
    handle=None;plan=None
    try:
        plan=read(plan_file)
        if plan['status']!='prepared':raise ValueError('Only a new plan can be run')
        handle=ProcessHandle(plan['adopted']['pid'],plan['adopted']['creation_filetime'])
        mutex=threading.Lock()
        def save():
            with mutex:_write_json(plan_file,plan)
        plan.update(status='awaiting_handoff',controller_pid=os.getpid(),adoption_ready_utc=now());save()
        deadline=time.monotonic()+180
        while not plan_file.with_name('handoff.commit').exists():
            if time.monotonic()>deadline:raise TimeoutError('Handoff not committed; existing controller retains ownership')
            time.sleep(.25)
        old=read(local_path(plan['old_queue']))
        if process_alive(plan['old_controller_pid']) or old['status']!='superseded':
            raise ValueError('Former controller has not released ownership')
        dependency_path=local_path(plan['dependency'])
        dependency=read(dependency_path)
        if dependency.get('controller_pid')!=plan['old_controller_pid']:
            raise ValueError('Dependency ownership changed unexpectedly')
        backup=dependency_path.with_name('scan.before_parallel_handoff_20261102.json')
        if backup.exists():raise ValueError('Dependency handoff already recorded')
        _write_json(backup,dependency)
        dependency.update(controller_pid=os.getpid(),parallel_plan=str(plan_file),owner_handoff_utc=now())
        _write_json(dependency_path,dependency)
        plan.update(status='running',started_utc=now());save()
        def adopt():
            last=0
            while handle.poll() is None:
                handle.wait_tick()
                if time.monotonic()-last>=15:
                    plan['last_checked_utc']=now();save();last=time.monotonic()
            code=handle.poll()
            result=finalize_adopted(plan,code)
            plan['adopted'].update(status='completed' if result==0 else 'failed',returncode=code,finished_utc=now());save()
            return result
        with ThreadPoolExecutor(max_workers=2) as pool:
            recovery=pool.submit(adopt)
            sweeps=pool.submit(run_sweeps,plan,save)
            codes=[recovery.result(),sweeps.result()]
        plan.update(status='completed' if not any(codes) else 'completed_with_failures',finished_utc=now());save()
        return int(any(codes))
    except Exception as exc:
        if plan is not None:
            plan.update(status='failed',error=str(exc),finished_utc=now());_write_json(plan_file,plan)
        print('Parallel handoff failed:',exc,flush=True)
        return 1
    finally:
        if handle:handle.close()
        lock.unlink()


def main():
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='command',required=True)
    setup=sub.add_parser('prepare');setup.add_argument('--plan',type=Path,required=True)
    setup.add_argument('--old-queue',type=Path,required=True);setup.add_argument('--pid',type=int,required=True)
    setup.add_argument('--scan',type=Path,action='append',required=True)
    for name in ['run','status']:sub.add_parser(name).add_argument('--plan',type=Path,required=True)
    args=parser.parse_args()
    if args.command=='prepare':prepare(args.plan,args.old_queue,args.pid,args.scan)
    elif args.command=='run':return run(args.plan)
    else:status(args.plan)
    return 0


if __name__=='__main__':raise SystemExit(main())
