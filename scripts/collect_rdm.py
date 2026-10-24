"""Collect additional equal-time density matrices under results/cleanRbRDM.

Initialize with --init, run the bounded queue with --run, inspect with --status.
Configurations/RNG states are copied from clean_Rb87; all measurements are new.
"""
from __future__ import annotations
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import h5py
import numpy as np

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / 'results' / 'cleanRbRDM'
BASE = 'simulation/realizations/0/clones/0'
POINTS = ('l128t035', 'l128t039', 'l128t045', 'l128t052',
          'l128t060', 'l128t070', 'l64t060', 'l64t070')
RANKS = 8
INTERVAL = 1000
BLOCK_SECONDS = 900
PILOT_SECONDS = 180
BLOCKS = 3
MPI = Path(r'C:\Program Files\Microsoft MPI\Bin\mpiexec.exe')
PACKAGE = REPO / 'dist' / 'windows-square'
EXE = PACKAGE / 'qmc_worm_canonical_mpi.exe'


def now():
    return datetime.now(timezone.utc).isoformat()


def owned(path):
    path = Path(path).resolve()
    if not path.is_relative_to(ROOT.resolve()):
        raise ValueError(f'Refusing to write outside collection folder: {path}')
    return path


def write_json(path, data):
    path = owned(path)
    temporary = path.with_suffix(path.suffix + '.new')
    temporary.write_text(json.dumps(data, indent=2) + '\n', encoding='utf-8')
    temporary.replace(path)


def event(kind, **details):
    record = dict(time_utc=now(), event=kind, **details)
    with (ROOT / 'status.jsonl').open('a', encoding='utf-8') as stream:
        stream.write(json.dumps(record) + '\n')
    print(json.dumps(record), flush=True)


def checkpoint_paths(directory):
    base = directory / 'rdm.clone.h5'
    return [Path(str(base) + ('' if rank == 0 else f'.{rank}')) for rank in range(RANKS)]


def initialize():
    if ROOT.exists():
        raise FileExistsError(f'Collection folder already exists: {ROOT}; inspect --status')
    if not EXE.is_file() or not MPI.is_file():
        raise FileNotFoundError('Required canonical MPI executable or Microsoft MPI launcher is absent')
    plan = []
    for name in POINTS:
        source_dir = REPO / 'results' / 'clean_Rb87' / name
        masters = list(source_dir.glob('*.clone.h5'))
        if len(masters) != 1:
            raise ValueError(f'Expected one checkpoint base in {source_dir}')
        master = masters[0]
        sources = [Path(str(master) + ('' if rank == 0 else f'.{rank}')) for rank in range(RANKS)]
        provenance = []
        for rank, path in enumerate(sources):
            with h5py.File(path, 'r') as handle:
                p = handle['parameters']
                if int(handle[BASE + '/checkpoint/sweeps'][()]) <= int(p['thermalization'][()]):
                    raise ValueError(f'Rank {rank} has not passed its thermalization threshold: {path}')
                L, beta, particles = int(p['Lx'][()]), float(p['beta'][()]), int(p['canonical'][()])
                if name != f'l{L}t{round(100 * beta):03d}' or int(p['Ly'][()]) != L:
                    raise ValueError(f'Checkpoint parameters disagree with folder: {path}')
                if particles <= 0:
                    raise ValueError(f'Expected exact-N canonical checkpoint: {path}')
            stat = path.stat()
            provenance.append(dict(rank=rank, path=str(path), size=stat.st_size,
                                   mtime_ns=stat.st_mtime_ns))
        plan.append(dict(folder=name, L=L, beta=beta, canonical=particles,
                         sources=provenance, prepared=False, completed=[]))
    ROOT.mkdir()
    manifest = dict(created_utc=now(), processes=RANKS, Nmeasure=1, Nmeasure2=INTERVAL,
                    block_seconds=BLOCK_SECONDS, production_blocks=BLOCKS,
                    initial_pilot_seconds=PILOT_SECONDS,
                    executable=str(EXE), executable_sha256=hashlib.sha256(EXE.read_bytes()).hexdigest(),
                    independent_new_seeds=False,
                    provenance='Continuation of copied rank 0..7 configurations and RNG states; new statistics.',
                    points=plan)
    write_json(ROOT / 'manifest.json', manifest)
    for point in plan:
        (ROOT / point['folder']).mkdir()
    event('initialized', points=list(POINTS), ranks=RANKS, Nmeasure2=INTERVAL)


def set_parameter(handle, name, value):
    parameters = handle['parameters']
    parameters[name][()] = value
    # ALPS reparses serialized raw INI values after loading the typed dictionary.
    # Keep both representations synchronized, or restoration reverts our edits.
    keys = [item.decode('utf-8') if isinstance(item, bytes) else str(item)
            for item in parameters.attrs['ini_keys']]
    if name in keys:
        values = [item.decode('utf-8') if isinstance(item, bytes) else str(item)
                  for item in parameters.attrs['ini_values']]
        values[keys.index(name)] = str(value)
        parameters.attrs.modify('ini_values', np.asarray(values, dtype=object))


def write_ini(checkpoint, target):
    lines = ['# Effective parameters of the copied checkpoint; restore uses HDF5, not this INI.',
             '# Statistics reset before each production; physical parameters and RNG state preserved.']
    with h5py.File(checkpoint, 'r') as handle:
        for name, ds in handle['parameters'].items():
            if not isinstance(ds, h5py.Dataset):
                continue
            value = ds[()]
            if isinstance(value, bytes):
                value = json.dumps(value.decode('utf-8'))
            elif isinstance(value, str):
                value = json.dumps(value)
            else:
                value = str(np.asarray(value).item())
            lines.append(f'{name} = {value}')
    owned(target).write_text('\n'.join(lines) + '\n', encoding='utf-8')


def prepare_point(point, manifest):
    directory = ROOT / point['folder']
    targets = checkpoint_paths(directory)
    if point['prepared']:
        if not all(path.is_file() for path in targets):
            raise RuntimeError(f'Prepared checkpoint missing for {point["folder"]}')
        return
    for source, target in zip(point['sources'], targets):
        if target.exists():
            raise FileExistsError(f'Unfinished preparation requires inspection: {target}')
        source_path = Path(source['path'])
        stat = source_path.stat()
        if stat.st_mtime_ns != source['mtime_ns'] or stat.st_size != source['size']:
            raise RuntimeError(f'Source checkpoint changed since initialization: {source_path}')
        # Repack only live checkpoint data; do not copy old measurement groups/free space.
        # The source archive is read-only. Group attributes and dataset types are retained.
        temporary = owned(target.with_name(target.name + '.preparing'))
        with h5py.File(source_path, 'r') as src, h5py.File(temporary, 'x') as dst:
            for key, value in src.attrs.items():
                dst.attrs[key] = value
            src.copy('parameters', dst)
            clone = dst.require_group(BASE)
            src.copy(BASE + '/checkpoint', clone, name='checkpoint')
            clone.create_group('measurements')
            set_parameter(dst, 'Nmeasure', 1)
            set_parameter(dst, 'Nmeasure2', INTERVAL)
            set_parameter(dst, 'runtimelimit', BLOCK_SECONDS)
            set_parameter(dst, 'Tmin', 1)
            set_parameter(dst, 'Tmax', 30)
            set_parameter(dst, 'checkpoint', 'rdm.clone.h5')
            set_parameter(dst, 'outputfile', 'rdmprod01.h5')
        temporary.replace(target)
        event('rank_prepared', folder=point['folder'], rank=source['rank'], bytes=target.stat().st_size)
    point['prepared'] = True
    write_json(ROOT / 'manifest.json', manifest)
    write_ini(targets[0], directory / 'parameters.ini')


def arm_block(directory, output_name, seconds):
    if (directory / output_name).exists():
        raise FileExistsError(f'Refusing to overwrite an existing result: {directory / output_name}')
    for path in checkpoint_paths(directory):
        owned(path)
        with h5py.File(path, 'r+') as handle:
            if int(handle[BASE + '/checkpoint/sweeps'][()]) <= int(handle['parameters/thermalization'][()]):
                raise RuntimeError(f'Cannot reset before thermalization: {path}')
            # Prevent ALPS append-mode saving from preserving a previous observable.
            del handle[BASE + '/measurements']
            handle.create_group(BASE + '/measurements')
            handle[BASE + '/checkpoint/counter'][...] = 0
            handle[BASE + '/checkpoint/configuration/hist_densmat'][...] = 0
            set_parameter(handle, 'checkpoint', 'rdm.clone.h5')
            set_parameter(handle, 'Nmeasure', 1)
            set_parameter(handle, 'Nmeasure2', INTERVAL)
            set_parameter(handle, 'Tmin', 1)
            set_parameter(handle, 'Tmax', 30)
            set_parameter(handle, 'outputfile', output_name)
            set_parameter(handle, 'runtimelimit', seconds)
    write_ini(checkpoint_paths(directory)[0], directory / (output_name + '.ini'))


def validate_result(directory, filename, label):
    path = directory / filename
    rank_rows, rank_payload = [], {}
    for rank, checkpoint in enumerate(checkpoint_paths(directory)):
        with h5py.File(checkpoint, 'r') as handle:
            measurements = handle[BASE + '/measurements']
            diagonal = int(measurements['Winding_number_squared/count'][()])
            obs = measurements.get('Density_Matrix')
            count = int(obs['count'][()]) if obs else 0
            if count != diagonal // INTERVAL:
                raise RuntimeError(f'Fresh histogram count mismatch for rank {rank}: {count} vs {diagonal}')
            rank_rows.append(dict(rank=rank, diagonal_count=diagonal, histogram_count=count))
            if obs:
                rank_payload[f'rank{rank}_mean'] = obs['mean/value'][()]
                rank_payload[f'rank{rank}_count'] = np.asarray(count)
    with h5py.File(path, 'r') as handle:
        p = handle['parameters']
        if int(p['Nmeasure2'][()]) != INTERVAL:
            raise RuntimeError('Requested measurement interval was not applied')
        results = handle['simulation/results']
        obs = results.get('Density_Matrix')
        count = int(obs['count'][()]) if obs else 0
        if count != sum(row['histogram_count'] for row in rank_rows):
            raise RuntimeError('Merged G count does not match rank counts')
        fraction = float(np.isfinite(obs['mean/error'][()]).mean()) if obs else 0.
        metrics = dict(file=str(path), diagonal_count=int(results['Winding_number_squared/count'][()]),
                       histogram_count=count, finite_error_fraction=fraction,
                       min_histograms_per_rank=min(row['histogram_count'] for row in rank_rows),
                       max_histograms_per_rank=max(row['histogram_count'] for row in rank_rows),
                       ranks=rank_rows)
    write_json(directory / f'{label}.metrics.json', metrics)
    np.savez_compressed(directory / f'{label}.rank_means.npz', **rank_payload)
    return metrics


class StopRequested(Exception):
    pass


def check_stop():
    if (ROOT / 'STOP').exists():
        raise StopRequested('Stop requested; current block has finished, no next block started')


def run_block(point, manifest, label, seconds, filename):
    check_stop()
    if label in point['completed']:
        return
    directory = ROOT / point['folder']
    arm_block(directory, filename, seconds)
    env = os.environ.copy()
    env['PATH'] = str(PACKAGE) + os.pathsep + env.get('PATH', '')
    temporary = directory / 'tmp'
    temporary.mkdir(exist_ok=True)
    env['TMPDIR'] = temporary.as_posix() + '/'
    args = [str(MPI), '-n', str(RANKS), str(EXE), 'rdm.clone.h5', '--reset-statistics']
    with (directory / f'{label}.stdout.log').open('wb') as stdout, (directory / f'{label}.stderr.log').open('wb') as stderr:
        process = subprocess.Popen(args, cwd=directory, env=env, stdout=stdout, stderr=stderr,
                                   creationflags=getattr(subprocess, 'CREATE_NO_WINDOW', 0))
        write_json(ROOT / 'active.json', dict(controller_pid=os.getpid(), mpi_pid=process.pid,
                   folder=point['folder'], block=label, started_utc=now(), runtime_seconds=seconds))
        event('block_started', folder=point['folder'], block=label, mpi_pid=process.pid, runtime_seconds=seconds)
        # Heartbeats let the user inspect progress while native logs are buffered.
        started = time.monotonic()
        while True:
            try:
                code = process.wait(timeout=30)
                break
            except subprocess.TimeoutExpired:
                event('heartbeat', folder=point['folder'], block=label, elapsed_seconds=round(time.monotonic()-started))
        if code:
            raise RuntimeError(f'{point["folder"]} {label}: MPI exited {code}; inspect stderr/stdout')
    metrics = validate_result(directory, filename, label)
    point['completed'].append(label)
    write_json(ROOT / 'manifest.json', manifest)
    event('block_complete', folder=point['folder'], block=label,
          histogram_count=metrics['histogram_count'], finite_error_fraction=metrics['finite_error_fraction'],
          min_histograms_per_rank=metrics['min_histograms_per_rank'])


def run_queue():
    if not (ROOT / 'manifest.json').is_file():
        raise RuntimeError('Initialize the collection with --init first')
    lock = owned(ROOT / 'queue.lock')
    with lock.open('x', encoding='utf-8') as stream:
        stream.write(str(os.getpid()))
    try:
        existing = subprocess.run(['powershell.exe', '-NoProfile', '-Command',
                    "@(Get-Process -Name 'qmc_worm*','mpiexec' -ErrorAction SilentlyContinue).Count"],
                    capture_output=True, text=True, check=True,
                    creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
        if existing.stdout.strip() != '0':
            raise RuntimeError('Another simulation is already running; leave it undisturbed and retry later')
        manifest = json.loads((ROOT / 'manifest.json').read_text(encoding='utf-8'))
        event('queue_started', controller_pid=os.getpid(), ranks=RANKS)
        for i, point in enumerate(manifest['points']):
            check_stop()
            prepare_point(point, manifest)
            if i == 0:
                run_block(point, manifest, 'pilot', PILOT_SECONDS, 'rdm_pilot.h5')
            for block in range(1, BLOCKS + 1):
                run_block(point, manifest, f'prod{block:02d}', BLOCK_SECONDS, f'rdmprod{block:02d}.h5')
            # The collation records queued/missing points too, so progress remains visible.
            os.environ['MPLBACKEND'] = 'Agg'
            from analyze_gr import collate_gr, export_gr
            export_gr(collate_gr(ROOT))
        write_json(ROOT / 'active.json', dict(state='complete', completed_utc=now()))
        event('queue_complete')
    except StopRequested as exc:
        write_json(ROOT / 'active.json', dict(state='stopped', stopped_utc=now(), reason=str(exc)))
        event('queue_stopped', reason=str(exc))
    except BaseException as exc:
        write_json(ROOT / 'active.json', dict(state='failed', failed_utc=now(), error=str(exc)))
        event('queue_failed', error=str(exc))
        raise
    finally:
        lock.unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument('--init', action='store_true')
    action.add_argument('--run', action='store_true')
    action.add_argument('--status', action='store_true')
    action.add_argument('--stop', action='store_true', help='Finish the current block, then stop the queue')
    args = parser.parse_args()
    if args.init:
        initialize()
    elif args.run:
        run_queue()
    elif args.stop:
        owned(ROOT / 'STOP').write_text('Stop requested at ' + now() + '\n', encoding='utf-8')
        print('Stop requested: the current block will finish, and no next block will start.')
    else:
        for name in ['active.json', 'manifest.json']:
            p = ROOT / name
            if p.exists():
                print(name + '\n' + p.read_text(encoding='utf-8'))
        log = ROOT / 'status.jsonl'
        if log.exists():
            print('\nRecent events:\n' + '\n'.join(log.read_text(encoding='utf-8').splitlines()[-8:]))


if __name__ == '__main__':
    main()
