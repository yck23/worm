"""Collect isolated continuation histories for winding and equal-time susceptibility.

Examples (from the repository root):
  python scripts/collect_response_measurements.py --points l128t030 l128t045 --prepare
  python scripts/collect_response_measurements.py --points l128t030 l128t045 --run --seconds 600 --jobs 2
  python scripts/collect_response_measurements.py --points l128t030 l128t045 --verify

Original checkpoints are read only. Re-running --run continues the SAME rank
histories; it never counts those continuations as independent new replicas.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import time
import h5py
import numpy as np

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / 'results' / 'response_measurements'
BASE = 'simulation/realizations/0/clones/0'
EXE = REPO / 'build-repro' / 'worm-square-response' / 'qmc_worm_mpi.exe'
DLL = REPO / 'dist' / 'windows-square'
MPI = Path(r'C:\Program Files\Microsoft MPI\Bin\mpiexec.exe')


def timestamp():
    return datetime.now(timezone.utc).isoformat()


def write_json(path, value):
    temporary = path.with_suffix(path.suffix + '.new')
    temporary.write_text(json.dumps(value, indent=2) + '\n', encoding='utf-8')
    temporary.replace(path)


def set_parameter(handle, key, value):
    parameters = handle['parameters']
    if key in parameters:
        parameters[key][()] = value
    else:
        parameters.create_dataset(key, data=value, dtype=parameters['Nmeasure2'].dtype)
        parameters[key].attrs['description'] = 'Response collection parameter'
        parameters[key].attrs['defnumber'] = 1 + max(int(ds.attrs.get('defnumber', 0))
                                                  for ds in parameters.values())
    keys = [v.decode() if isinstance(v, bytes) else str(v) for v in parameters.attrs['ini_keys']]
    values = [v.decode() if isinstance(v, bytes) else str(v) for v in parameters.attrs['ini_values']]
    if key in keys:
        values[keys.index(key)] = str(value)
    else:
        keys.append(key)
        values.append(str(value))
    parameters.attrs['ini_keys'] = np.asarray(keys, dtype=h5py.string_dtype())
    parameters.attrs['ini_values'] = np.asarray(values, dtype=h5py.string_dtype())


def checkpoint_paths(directory, ranks):
    return [directory / ('response.clone.h5' + (f'.{rank}' if rank else ''))
            for rank in range(ranks)]


def prepare(name, ranks, interval):
    if not re.fullmatch(r'l(32|64|128)t0\d{2}', name):
        raise ValueError(f'Invalid scan point: {name}')
    directory = ROOT / 'sampling' / name
    if directory.exists():
        metadata = json.loads((directory / 'manifest.json').read_text())
        if (metadata['ranks'], metadata['interval']) != (ranks, interval):
            raise ValueError('Existing point has a different rank count or response interval')
        return directory
    source_dir = REPO / 'results' / 'cleanRbRDM' / name
    if not source_dir.is_dir():
        source_dir = REPO / 'results' / 'clean_Rb87' / name
    sources = list(source_dir.glob('*.clone.h5'))
    if len(sources) != 1:
        raise ValueError(f'Expected one source checkpoint family in {source_dir}')
    sources = [Path(str(sources[0]) + (f'.{rank}' if rank else '')) for rank in range(ranks)]
    provenance = []
    for rank, source in enumerate(sources):
        stat = source.stat()
        with h5py.File(source, 'r') as f:
            q = f['parameters']
            if int(q['canonical'][()]) < 0 or int(f[BASE + '/checkpoint/sweeps'][()]) <= int(q['thermalization'][()]):
                raise ValueError(f'Source must be an already thermalized canonical checkpoint: {source}')
            if name != f'l{int(q["Lx"][()])}t{round(100*float(q["beta"][()])):03d}':
                raise ValueError('Source physical parameters disagree with requested point')
        provenance.append(dict(rank=rank, source=str(source.relative_to(REPO)),
                               bytes=stat.st_size, mtime_ns=stat.st_mtime_ns))
    directory.mkdir(parents=True)
    for rank, (source, target) in enumerate(zip(sources, checkpoint_paths(directory, ranks))):
        with h5py.File(source, 'r') as src, h5py.File(target, 'x') as dst:
            for key, value in src.attrs.items():
                dst.attrs[key] = value
            src.copy('parameters', dst)
            clone = dst.require_group(BASE)
            src.copy(BASE + '/checkpoint', clone, name='checkpoint')
            clone.create_group('measurements')
            for key, value in dict(Nresponse_measure=interval, Nresponse_max_blocks=1000000,
                                   Nmeasure=1, Nmeasure2=100, sweeps=1000000000,
                                   runtimelimit=600, Tmin=1, Tmax=15,
                                   checkpoint='response.clone.h5', outputfile='run01.h5').items():
                set_parameter(dst, key, value)
    write_json(directory / 'manifest.json', dict(point=name, ranks=ranks, interval=interval,
        Nmeasure=1, Nmeasure2=100, created_utc=timestamp(), sources=provenance,
        independent_new_seeds=False, provenance='Warm continuation of copied rank states and RNG streams; source archives unchanged.',
        runs=[], executable=str(EXE.relative_to(REPO))))
    with h5py.File(checkpoint_paths(directory, ranks)[0], 'r') as f:
        lines = ['# Effective parameters; restart reads the HDF5 checkpoint.']
        for key, ds in f['parameters'].items():
            value = ds[()]
            value = json.dumps(value.decode()) if isinstance(value, bytes) else str(np.asarray(value).item())
            lines.append(f'{key} = {value}')
        (directory / 'parameters.ini').write_text('\n'.join(lines) + '\n', encoding='utf-8')
    print('Prepared', name, ranks, 'ranks', flush=True)
    return directory


def check_save_warnings(directory, number):
    """MPI may return success after failing to save a rank checkpoint."""
    for stream in ('stdout', 'stderr'):
        path = directory / f'run{number:02d}.{stream}.log'
        for line in path.read_text(encoding='utf-8', errors='replace').splitlines():
            if re.search(r'warning.*(?:saving checkpoint|failed to save checkpoint)', line, re.I):
                raise RuntimeError(f'Checkpoint-save warning in {path.name}: {line.strip()}')


def observable_count(group, name):
    # ALPS omits accumulator groups with zero samples when saving an archive.
    return int(group[name + '/count'][()]) if group is not None and name in group else 0


def export_history(directory, metadata, merged_output):
    """Validate rank histories against the merged output before publishing them."""
    payload = {}
    checks = []
    observables = [('Superfluid_susceptibility', 'chi_s'), ('Winding_number_squared_sum', 'w2')]
    totals = {obs: dict(count=0, total=0.) for obs, _ in observables}
    physical = None
    for rank, path in enumerate(checkpoint_paths(directory, metadata['ranks'])):
        with h5py.File(path, 'r') as f:
            q = f['parameters']
            group = f[BASE + '/checkpoint/response']
            values = {}
            for name in ['chi_s', 'w2', 'density', 'end_sweep', 'sweep_count']:
                ds = group[name]
                # A saved empty std::vector has an HDF5 null dataspace, not shape (0,).
                values[name] = np.empty(0, dtype=ds.dtype) if ds.shape is None else np.asarray(ds[()])
                if values[name].ndim != 1 or not np.isfinite(values[name]).all():
                    raise ValueError(f'Invalid rank {rank} response vector {name}')
                payload[f'rank{rank}_{name}'] = values[name]
            n = len(values['chi_s'])
            if any(len(value) != n for value in values.values()):
                raise ValueError(f'Inconsistent response vector lengths on rank {rank}')
            interval = int(group['interval'][()])
            partial = int(group['partial_count'][()])
            if interval != metadata['interval'] or not 0 <= partial < interval:
                raise ValueError(f'Invalid response interval or partial count on rank {rank}')
            if not np.all(values['sweep_count'] == interval):
                raise ValueError(f'Unexpected response block lengths on rank {rank}')
            measurements = f[BASE + '/measurements']
            for obs, field in observables:
                count = observable_count(measurements, obs)
                if count != n:
                    raise ValueError(f'{obs} count/history mismatch on rank {rank}: {count} != {n}')
                if n:
                    observed = float(measurements[obs + '/mean/value'][()])
                    if not np.isclose(observed, values[field].mean(), rtol=1e-10, atol=1e-12):
                        raise ValueError(f'{obs} mean/history mismatch on rank {rank}')
                totals[obs]['count'] += count
                totals[obs]['total'] += float(values[field].sum())
            current = dict(L=int(q['Lx'][()]), beta=float(q['beta'][()]),
                           N=int(q['canonical'][()]), nmax=int(q['nmax'][()]), interval=interval)
            if physical is not None and current != physical:
                raise ValueError(f'Physical parameters differ between rank checkpoints: rank {rank}')
            physical = current
            checks.append(dict(rank=rank, blocks=n, partial_count=partial,
                               sweeps=int(f[BASE + '/checkpoint/sweeps'][()])))
    merged_checks = []
    with h5py.File(merged_output, 'r') as f:
        q = f['parameters']
        current = dict(L=int(q['Lx'][()]), beta=float(q['beta'][()]),
                       N=int(q['canonical'][()]), nmax=int(q['nmax'][()]),
                       interval=int(q['Nresponse_measure'][()]))
        if current != physical:
            raise ValueError('Merged output physical parameters differ from rank checkpoints')
        results = f.get('simulation/results')
        for obs, _ in observables:
            expected = totals[obs]
            count = observable_count(results, obs)
            if count != expected['count']:
                raise ValueError(f'{obs} merged/checkpoint count mismatch: {count} != {expected["count"]}')
            observed = float(results[obs + '/mean/value'][()]) if count else None
            reference = expected['total'] / count if count else None
            if count and not np.isclose(observed, reference, rtol=1e-10, atol=1e-12):
                raise ValueError(f'{obs} merged/checkpoint mean mismatch: {observed} != {reference}')
            merged_checks.append(dict(observable=obs, count=count, merged_mean=observed,
                                      rank_weighted_mean=reference))
    payload.update({name: np.asarray(value) for name, value in physical.items()})
    temporary = directory / 'history.npz.new'
    with temporary.open('wb') as stream:
        np.savez_compressed(stream, **payload)
    temporary.replace(directory / 'history.npz')
    write_json(directory / 'history_checks.json', dict(time_utc=timestamp(), ranks=checks,
        merged_output=merged_output.name, merged_checks=merged_checks, status='verified',
        history_includes='All completed response blocks since initial reset; partial blocks only in checkpoints.'))
    return checks


def verify_point(name):
    """Re-export a completed point without starting or modifying simulations."""
    if not re.fullmatch(r'l(32|64|128)t0\d{2}', name):
        raise ValueError(f'Invalid scan point: {name}')
    directory = ROOT / 'sampling' / name
    if (directory / 'RUNNING.lock').exists():
        raise RuntimeError(f'Cannot verify a point with an active run lock: {name}')
    metadata = json.loads((directory / 'manifest.json').read_text())
    if not metadata['runs']:
        raise ValueError(f'No completed run is recorded for {name}')
    latest = metadata['runs'][-1]
    check_save_warnings(directory, latest['number'])
    checks = export_history(directory, metadata, directory / latest['output'])
    print('Verified', name, 'rank blocks', [c['blocks'] for c in checks], flush=True)
    return checks


def run_point(name, ranks, interval, seconds):
    directory = prepare(name, ranks, interval)
    lock = directory / 'RUNNING.lock'
    with lock.open('x') as f:
        f.write(str(os.getpid()))
    try:
        metadata = json.loads((directory / 'manifest.json').read_text())
        number = len(metadata['runs']) + 1
        output = f'run{number:02d}.h5'
        if (directory / output).exists():
            raise FileExistsError(directory / output)
        for checkpoint in checkpoint_paths(directory, ranks):
            with h5py.File(checkpoint, 'r+') as f:
                set_parameter(f, 'runtimelimit', seconds)
                set_parameter(f, 'outputfile', output)
        env = os.environ.copy()
        env['PATH'] = str(DLL) + os.pathsep + env.get('PATH', '')
        (directory / 'tmp').mkdir(exist_ok=True)
        env['TMPDIR'] = (directory / 'tmp').as_posix() + '/'
        command = [str(MPI), '-n', str(ranks), str(EXE), 'response.clone.h5']
        if number == 1:
            command.append('--reset-statistics')
        started = timestamp()
        with (directory / f'run{number:02d}.stdout.log').open('wb') as stdout, (directory / f'run{number:02d}.stderr.log').open('wb') as stderr:
            process = subprocess.Popen(command, cwd=directory, env=env, stdout=stdout, stderr=stderr,
                creationflags=getattr(subprocess, 'CREATE_NO_WINDOW', 0))
            write_json(directory / 'active.json', dict(pid=process.pid, run=number, started_utc=started, runtime_seconds=seconds))
            print('Started', name, 'run', number, 'seconds', seconds, 'pid', process.pid, flush=True)
            while True:
                try:
                    code = process.wait(timeout=30)
                    break
                except subprocess.TimeoutExpired:
                    print('Running', name, 'run', number, flush=True)
        if code:
            raise RuntimeError(f'{name} exited {code}; inspect run{number:02d}.stderr.log')
        check_save_warnings(directory, number)
        checks = export_history(directory, metadata, directory / output)
        metadata['runs'].append(dict(number=number, started_utc=started, completed_utc=timestamp(),
             runtime_seconds=seconds, output=output, executable_sha256=hashlib.sha256(EXE.read_bytes()).hexdigest(),
             rank_blocks=[c['blocks'] for c in checks]))
        write_json(directory / 'manifest.json', metadata)
        write_json(directory / 'active.json', dict(state='complete', completed_utc=timestamp(), run=number))
        print('Completed', name, 'rank blocks', [c['blocks'] for c in checks], flush=True)
    finally:
        lock.unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--points', nargs='+', required=True)
    parser.add_argument('--ranks', type=int, default=8)
    parser.add_argument('--interval', type=int, default=10)
    parser.add_argument('--seconds', type=int, default=600)
    parser.add_argument('--jobs', type=int, default=1)
    parser.add_argument('--prepare', action='store_true')
    parser.add_argument('--run', action='store_true')
    parser.add_argument('--verify', action='store_true', help='Validate/export completed histories without running simulations')
    args = parser.parse_args()
    if len(set(args.points)) != len(args.points) or not 1 <= args.ranks <= 24 or args.jobs * args.ranks > 24:
        parser.error('Points must be unique; use 1-24 ranks and at most 24 total simultaneous ranks.')
    if args.interval < 1 or args.seconds < 1 or args.jobs < 1:
        parser.error('Interval, runtime and jobs must be positive.')
    if args.verify and (args.prepare or args.run):
        parser.error('--verify cannot be combined with --prepare or --run.')
    if args.prepare:
        for point in args.points:
            prepare(point, args.ranks, args.interval)
    if args.run:
        if not EXE.is_file() or not MPI.is_file():
            raise FileNotFoundError('Build the response executable and install Microsoft MPI first.')
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            futures = [pool.submit(run_point, p, args.ranks, args.interval, args.seconds) for p in args.points]
            for future in futures:
                future.result()
    if args.verify:
        for point in args.points:
            verify_point(point)
    if not args.prepare and not args.run and not args.verify:
        parser.error('Choose --prepare and/or --run, or --verify.')


if __name__ == '__main__':
    main()
