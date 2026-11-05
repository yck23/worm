#!/usr/bin/env python3
"""Run prepared disorder scans in order after an existing scan or queue completes.

The detached controller only reads the dependency while it is running. It never
restarts a scan or modifies that scan's inputs, checkpoints, or measurements.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

from run_disorder_ensemble import _write_json

REPO = Path(__file__).resolve().parents[1]


def now():
    return datetime.now(timezone.utc).isoformat()


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def local_path(value):
    path = Path(value)
    return path if path.is_absolute() else REPO / path


def scan_reference(path):
    return str(Path(path).resolve(strict=True).relative_to(REPO))


def process_alive(pid):
    """Check a controller without signaling it (os.kill is unsafe on Windows)."""
    if not pid:
        return False
    if os.name == 'nt':
        import ctypes
        from ctypes import wintypes
        kernel = ctypes.WinDLL('kernel32', use_last_error=True)
        kernel.OpenProcess.argtypes = [wintypes.DWORD, wintypes.BOOL, wintypes.DWORD]
        kernel.OpenProcess.restype = wintypes.HANDLE
        kernel.GetExitCodeProcess.argtypes = [wintypes.HANDLE, ctypes.POINTER(wintypes.DWORD)]
        kernel.GetExitCodeProcess.restype = wintypes.BOOL
        kernel.CloseHandle.argtypes = [wintypes.HANDLE]
        kernel.CloseHandle.restype = wintypes.BOOL
        handle = kernel.OpenProcess(0x1000, False, int(pid))
        if not handle:
            return ctypes.get_last_error() == 5  # Access denied: do not assume dead.
        try:
            code = wintypes.DWORD()
            if not kernel.GetExitCodeProcess(handle, ctypes.byref(code)):
                raise ctypes.WinError(ctypes.get_last_error())
            return code.value == 259  # STILL_ACTIVE
        finally:
            kernel.CloseHandle(handle)
    try:
        os.kill(int(pid), 0)
        return True
    except ProcessLookupError:
        return False
    except PermissionError:
        return True


def prepare(queue_path, dependency, scans, poll_seconds=60,
            recover_missing_dependency=False, recovery_seconds=43200, recovery_workers=2):
    queue_path = Path(queue_path).resolve()
    if queue_path.exists():
        raise ValueError('Queue already exists; use a new queue file for a new sequence')
    if not 1 <= poll_seconds <= 60:
        raise ValueError('poll_seconds must be between 1 and 60')
    dependency = scan_reference(dependency)
    dependency_state = read(local_path(dependency))
    allowed=('prepared','waiting_for_dependency','recovering_dependency','running','completed')
    if dependency_state['status'] not in allowed and not (recover_missing_dependency and dependency_state['status']=='failed'):
        raise ValueError('Dependency must be prepared, waiting, running, completed, or recoverable failed')
    if recover_missing_dependency:
        if dependency_state.get('kind')!='disorder_response_scan':
            raise ValueError('Checkpoint recovery requires a scan dependency')
        if recovery_seconds <= 0 or recovery_workers <= 0:
            raise ValueError('Recovery seconds and workers must be positive')
    stages = []
    for scan in scans:
        reference = scan_reference(scan)
        if reference == dependency or reference in [s['scan'] for s in stages]:
            raise ValueError('Duplicate scan in queue')
        state = read(local_path(reference))
        if state['status'] != 'prepared':
            raise ValueError(f'Stage is not newly prepared: {reference}')
        actual_hash = hashlib.sha256(Path(state['executable']).read_bytes()).hexdigest()
        if actual_hash != state['executable_sha256']:
            raise ValueError(f'Executable changed: {reference}')
        stages.append(dict(beta=state['config']['beta'], scan=reference, status='queued'))
    if not stages:
        raise ValueError('At least one follow-on scan is required')
    queue_path.parent.mkdir(parents=True, exist_ok=True)
    state = dict(schema_version=1, kind='disorder_scan_queue', status='prepared',
                 created_utc=now(), dependency=dependency, poll_seconds=poll_seconds,
                 stages=stages)
    if recover_missing_dependency:
        state['dependency_recovery']=dict(runtime_seconds=recovery_seconds,workers=recovery_workers,
            only_missing_measurements=True,max_attempts=1)
    _write_json(queue_path, state)
    shutil.copyfile(__file__, queue_path.parent / 'run_disorder_queue.py')
    print(f'Prepared {len(stages)} sequential scans: {queue_path}', flush=True)


def execute_stage(scan_path):
    root = Path(scan_path).parent
    command = [sys.executable, str(REPO/'scripts/disorder_response_scan.py'),
               'run', '--scan', str(scan_path)]
    with (root/'controller.stdout.log').open('a', encoding='utf-8') as stdout, \
         (root/'controller.stderr.log').open('a', encoding='utf-8') as stderr:
        return subprocess.run(command, cwd=REPO, stdout=stdout, stderr=stderr,
                              creationflags=getattr(subprocess, 'CREATE_NO_WINDOW', 0)).returncode


def recover_dependency(scan_path, options):
    from resume_disorder_response import prepare as prepare_continuation, run as run_continuation
    prepare_continuation(scan_path, options['runtime_seconds'], options['workers'],
                         only_missing_measurements=True)
    return run_continuation(scan_path)


def run(queue_path):
    queue_path = Path(queue_path).resolve(strict=True)
    lock_path = queue_path.with_suffix('.lock')
    # An exclusive lock plus the prepared-state check prevent duplicate launches.
    with lock_path.open('x', encoding='utf-8') as lock:
        lock.write(str(os.getpid()))
    state = None
    try:
        state = read(queue_path)
        if state['status'] != 'prepared':
            raise ValueError('Only a newly prepared queue can be launched')
        state.update(status='waiting_for_dependency', controller_pid=os.getpid(), started_utc=now())
        _write_json(queue_path, state)
        absent_checks = 0
        while True:
            dependency = read(local_path(state['dependency']))
            state.update(last_checked_utc=now(), dependency_status=dependency['status'])
            _write_json(queue_path, state)
            if dependency['status'] == 'completed':
                entries=dependency['stages'] if dependency.get('kind')=='disorder_scan_queue' else dependency['points']
                if not all(p['status'] == 'completed' for p in entries):
                    raise RuntimeError('Dependency claims completion but has incomplete points or stages')
                break
            if dependency['status']=='failed' and state.get('dependency_recovery') and not state.get('recovery_attempted'):
                # Final analysis may still be running after the failed scan is written.
                # Wait for its controller to exit before copying any checkpoint.
                if process_alive(dependency.get('controller_pid')):
                    time.sleep(state['poll_seconds'])
                    continue
                state.update(status='recovering_dependency',recovery_attempted=True,recovery_started_utc=now())
                _write_json(queue_path,state)
                print(f'{now()} Continuing missing-production dependency chains from checkpoints',flush=True)
                code=recover_dependency(local_path(state['dependency']),state['dependency_recovery'])
                state.update(recovery_returncode=code,recovery_finished_utc=now(),status='waiting_for_dependency')
                _write_json(queue_path,state)
                if code:
                    raise RuntimeError('Checkpoint recovery did not complete; inspect continuation logs')
                continue
            if dependency['status'] not in ('prepared', 'waiting_for_dependency', 'recovering_dependency', 'running'):
                raise RuntimeError(f'Dependency stopped with status {dependency["status"]}')
            if dependency['status'] in ('running', 'waiting_for_dependency', 'recovering_dependency'):
                absent_checks = 0 if process_alive(dependency.get('controller_pid')) else absent_checks + 1
                if absent_checks >= 3:
                    raise RuntimeError('Dependency controller absent for three consecutive checks; inspect its logs')
            time.sleep(state['poll_seconds'])
        state.update(status='running', dependency_finished_utc=dependency.get('finished_utc', now()))
        _write_json(queue_path, state)
        for stage in state['stages']:
            scan_path = local_path(stage['scan'])
            if read(scan_path)['status'] != 'prepared':
                raise RuntimeError(f'Refusing to restart existing scan: {scan_path}')
            stage.update(status='running', started_utc=now())
            state['current_beta'] = stage['beta']
            _write_json(queue_path, state)
            print(f'{now()} Starting beta={stage["beta"]:.2f}: {scan_path}', flush=True)
            code = execute_stage(scan_path)
            result = read(scan_path)
            successful = code == 0 and result['status'] == 'completed'
            stage.update(status='completed' if successful else 'failed', returncode=code, finished_utc=now())
            _write_json(queue_path, state)
            if result['status'] not in ('completed', 'failed'):
                raise RuntimeError(f'Scan controller interrupted ({result["status"]}); inspect before starting another scan')
            print(f'{now()} beta={stage["beta"]:.2f}: {stage["status"]}', flush=True)
        state.update(status='completed' if all(s['status']=='completed' for s in state['stages'])
                     else 'completed_with_failures', finished_utc=now())
        state.pop('current_beta', None)
        _write_json(queue_path, state)
        return 0 if state['status'] == 'completed' else 1
    except Exception as exc:
        if state is not None and state.get('controller_pid') == os.getpid():
            state.update(status='failed', finished_utc=now(), error=str(exc))
            _write_json(queue_path, state)
        print(f'Queue stopped: {exc}', file=sys.stderr, flush=True)
        return 1
    finally:
        lock_path.unlink()


def status(queue_path):
    state = read(queue_path)
    print(f'Queue: {state["status"]}; controller PID: {state.get("controller_pid", "not launched")}')
    print(f'Last dependency check: {state.get("last_checked_utc", "not checked")}')
    for label, path in [('dependency', state['dependency']), *[(f'beta={s["beta"]:.2f}', s['scan']) for s in state['stages']]]:
        scan = read(local_path(path))
        is_queue=scan.get('kind')=='disorder_scan_queue'
        entries=scan['stages'] if is_queue else scan['points']
        count = sum(p['status']=='completed' for p in entries)
        unit='stages' if is_queue else 'points'
        print(f'{label}: {scan["status"]}; {count}/{len(entries)} {unit}; {path}')
    if state.get('error'):
        print('Error: '+state['error'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    setup = commands.add_parser('prepare')
    setup.add_argument('--queue', type=Path, required=True)
    setup.add_argument('--after', type=Path, required=True)
    setup.add_argument('--scan', type=Path, action='append', required=True)
    setup.add_argument('--poll-seconds', type=int, default=60)
    setup.add_argument('--recover-missing-dependency',action='store_true')
    setup.add_argument('--recovery-seconds',type=int,default=43200)
    setup.add_argument('--recovery-workers',type=int,default=2)
    for name in ('run', 'status'):
        commands.add_parser(name).add_argument('--queue', type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'prepare':
        prepare(args.queue, args.after, args.scan, args.poll_seconds, args.recover_missing_dependency, args.recovery_seconds, args.recovery_workers)
    elif args.command == 'run':
        return run(args.queue)
    else:
        status(args.queue)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
