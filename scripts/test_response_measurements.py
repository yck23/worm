"""Tiny canonical integration regressions for joint scalar response measurements.

Run with a Python containing NumPy/h5py; miniforge3/python.exe is available here.
This checks bookkeeping and reproducibility, not physical equilibration.
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
import time
import traceback
import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
BASE = 'simulation/realizations/0/clones/0'
RESP = BASE + '/checkpoint/response'
OBS = ('Superfluid_susceptibility', 'Winding_number_squared_sum')
DEFAULTS = dict(model='BoseHubbard', runtimelimit=60, sweeps=2000, thermalization=17,
    Lx=3, Ly=3, Lz=1, pbcx=1, pbcy=1, pbcz=0, t_hop=1., U_on=.7, V_nn=0.,
    mu=-2.4, beta=.3, nmax=5, canonical=13, canonical_window=.5,
    initial_occupancy=1, seed=337711, Nmeasure=1, Nmeasure2=10,
    Nresponse_measure=10, Nresponse_max_blocks=10000, Ntest=1000000,
    Nsave=1000000, C_worm=2., p_insertworm=1., p_moveworm=.4,
    p_insertkink=.3, p_deletekink=.2, p_glueworm=.1,
    checkpoint='state.h5', outputfile='result.h5')


def strings(items):
    return [x.decode() if isinstance(x, bytes) else str(x) for x in items]


def patch_parameters(path, changes=None, remove=()):
    """Synchronize ALPS's typed parameters and serialized raw INI values."""
    with h5py.File(path, 'r+') as f:
        group = f['parameters']
        raw = dict(zip(strings(group.attrs['ini_keys']), strings(group.attrs['ini_values'])))
        for name in remove:
            if name in group:
                del group[name]
            raw.pop(name, None)
        for name, value in (changes or {}).items():
            if name in group:
                group[name][()] = value
            else:
                highest = max(int(ds.attrs.get('defnumber', -1)) for ds in group.values())
                ds = group.create_dataset(name, data=np.uint32(value))
                ds.attrs['defnumber'] = np.int32(highest + 1)
                ds.attrs['description'] = 'Response regression migration'
            raw[name] = str(value)
        dtype = h5py.string_dtype('utf-8')
        for name, values in [('ini_keys', list(raw)), ('ini_values', list(raw.values()))]:
            del group.attrs[name]
            group.attrs.create(name, np.asarray(values, dtype=dtype), dtype=dtype)


def mean(f, name, checkpoint=False):
    prefix = BASE + '/measurements' if checkpoint else 'simulation/results'
    return np.asarray(f[f'{prefix}/{name}/mean/value'][()])


def snapshot(path):
    with h5py.File(path, 'r') as f:
        return {name: (np.empty(0, dtype=ds.dtype) if ds.shape is None else ds[()])
                for name, ds in f[RESP].items()}


def identical_response(left, right):
    a, b = snapshot(left), snapshot(right)
    assert a.keys() == b.keys()
    for name in a:
        np.testing.assert_array_equal(a[name], b[name], err_msg='response/' + name)


def identical_configuration(left, right):
    with h5py.File(left, 'r') as a, h5py.File(right, 'r') as b:
        for name in ('random', 'configuration/operator_string', 'configuration/mWinding',
                     'configuration/nrvertex', 'configuration/nr_of_particles',
                     'configuration/state', 'configuration/hist_densmat'):
            key = f'{BASE}/checkpoint/{name}'
            np.testing.assert_array_equal(a[key][()], b[key][()], err_msg=name)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', type=Path,
                        default=ROOT / 'build-repro/worm-square-response/qmc_worm.exe')
    parser.add_argument('--output', type=Path,
                        default=ROOT / 'results/response_measurements/validation')
    args = parser.parse_args()
    executable = args.executable.resolve()
    folder = args.output.resolve() / datetime.now(timezone.utc).strftime('run_%Y%m%dT%H%M%S_%fZ')
    folder.mkdir(parents=True, exist_ok=False)
    report = dict(schema_version=1, executable=str(executable),
        executable_sha256=hashlib.sha256(executable.read_bytes()).hexdigest(),
        started_utc=datetime.now(timezone.utc).isoformat(), checks=[], status='running')
    started = time.monotonic()

    def passed(name, **metrics):
        report['checks'].append(dict(name=name, passed=True, **metrics))
        print('PASS ' + name, flush=True)

    def fresh(name, **overrides):
        case = folder / name
        case.mkdir()
        params = DEFAULTS | overrides
        (case / 'input.ini').write_text('\n'.join(
            f'{key} = {json.dumps(value) if isinstance(value, str) else value}'
            for key, value in params.items()) + '\n', encoding='utf-8')
        return case

    def copy_state(name, source):
        case = folder / name
        case.mkdir()
        shutil.copy2(source / 'state.h5', case / 'state.h5')
        return case

    def run(case, restore=False, reset=False, expect=None):
        tmp = case / 'tmp'
        tmp.mkdir(exist_ok=True)
        env = os.environ.copy()
        env['PATH'] = str(ROOT / 'dist/windows-square') + os.pathsep + env.get('PATH', '')
        env['TMPDIR'] = tmp.as_posix() + '/'
        command = [str(executable), 'state.h5' if restore else 'input.ini']
        if reset:
            command.append('--reset-statistics')
        log = case / f'run_{len(list(case.glob("run_*.log"))) + 1:02d}.log'
        with log.open('wb') as stream:
            result = subprocess.run(command, cwd=case, env=env, stdout=stream,
                stderr=subprocess.STDOUT, timeout=90,
                creationflags=getattr(subprocess, 'CREATE_NO_WINDOW', 0))
        text = log.read_text(encoding='utf-8', errors='replace')
        if expect:
            assert result.returncode != 0 and expect in text, text[-3000:]
            return
        assert result.returncode == 0, text[-4000:]
        assert 'WARNING: Failed to save checkpoint' not in text, text[-4000:]
        assert (case / 'state.h5').is_file() and (case / 'result.h5').is_file()

    try:
        aligned = fresh('aligned')
        run(aligned)
        with h5py.File(aligned / 'result.h5', 'r') as f:
            g = mean(f, 'Density_Matrix')
            chi = float(mean(f, OBS[0]))
            w2 = float(mean(f, OBS[1]))
            expected_chi = float(g.sum() - g[0] + 13 / 9)
            expected_w2 = float(mean(f, 'Winding_number_squared').sum())
            np.testing.assert_allclose(chi, expected_chi, rtol=2e-13, atol=1e-13)
            np.testing.assert_allclose(w2, expected_w2, rtol=2e-13, atol=1e-13)
        h = snapshot(aligned / 'state.h5')
        assert len(h['chi_s']) == 200 and int(h['partial_count']) == 0
        np.testing.assert_allclose(h['chi_s'].mean(), chi, rtol=2e-13)
        np.testing.assert_allclose(h['w2'].mean(), w2, rtol=2e-13)
        np.testing.assert_allclose(h['density'], 13 / 9, rtol=2e-15)
        np.testing.assert_array_equal(h['sweep_count'], np.full(200, 10))
        np.testing.assert_array_equal(h['end_sweep'], 17 + np.arange(1, 201) * 10)
        passed('aligned_scalar_vector_identity', blocks=200, chi_s=chi, w2sum=w2,
               chi_difference=abs(chi-expected_chi), w2_difference=abs(w2-expected_w2))

        disabled = fresh('default_off', Nresponse_measure=0)
        ini = disabled / 'input.ini'
        ini.write_text('\n'.join(x for x in ini.read_text().splitlines()
                                 if not x.startswith('Nresponse_')) + '\n')
        run(disabled)
        with h5py.File(disabled / 'state.h5', 'r') as f:
            assert RESP not in f
            assert all(name not in f[BASE + '/measurements'] for name in OBS)
        identical_configuration(aligned / 'state.h5', disabled / 'state.h5')
        passed('default_off_preserves_trajectory')

        legacy = copy_state('legacy_off', disabled)
        patch_parameters(legacy / 'state.h5', {'sweeps': 2010},
                         remove=('Nresponse_measure', 'Nresponse_max_blocks'))
        run(legacy, restore=True)
        with h5py.File(legacy / 'state.h5', 'r') as f:
            assert RESP not in f
            assert int(f[BASE + '/checkpoint/sweeps'][()]) == 2027
        passed('legacy_checkpoint_without_new_parameters')

        migrated = copy_state('legacy_enabled', disabled)
        patch_parameters(migrated / 'state.h5', remove=('Nresponse_measure', 'Nresponse_max_blocks'))
        patch_parameters(migrated / 'state.h5',
                         {'sweeps': 2067, 'Nresponse_measure': 31, 'Nresponse_max_blocks': 10000})
        run(migrated, restore=True)
        h = snapshot(migrated / 'state.h5')
        assert len(h['chi_s']) == 2 and int(h['partial_count']) == 5
        np.testing.assert_array_equal(h['end_sweep'], [2048, 2079])
        passed('explicit_legacy_activation', blocks=2, partial_sweeps=5)

        full = fresh('uninterrupted', sweeps=1003, Nmeasure2=29, Nresponse_measure=31)
        run(full)
        split = fresh('restart', sweeps=397, Nmeasure2=29, Nresponse_measure=31)
        run(split)
        assert int(snapshot(split / 'state.h5')['partial_count']) == 25
        shutil.copy2(split / 'state.h5', split / 'segment1.h5')
        patch_parameters(split / 'state.h5', {'sweeps': 1003})
        run(split, restore=True)
        identical_response(full / 'state.h5', split / 'state.h5')
        identical_configuration(full / 'state.h5', split / 'state.h5')
        h = snapshot(split / 'state.h5')
        assert len(h['chi_s']) == 32 and int(h['partial_count']) == 11
        passed('deterministic_partial_block_restart', initial_partial_sweeps=25,
               final_blocks=32, final_partial_sweeps=11, bitwise_identical=True)

        with h5py.File(full / 'state.h5', 'r') as f:
            total = int(np.sum(h['sweep_count'])) + int(h['partial_count'])
            w2_all = (float(np.dot(h['w2'], h['sweep_count'])) + float(h['w2sum'])) / total
            density_all = (float(np.dot(h['density'], h['sweep_count'])) + float(h['density_sum'])) / total
            np.testing.assert_allclose(w2_all, mean(f, 'Winding_number_squared', True).sum(),
                                       rtol=2e-13, atol=1e-13)
            np.testing.assert_allclose(density_all, float(mean(f, 'Number_of_particles', True)) / 9,
                                       rtol=2e-13, atol=1e-13)
        passed('partial_winding_density_accounting', represented_sweeps=total)

        reset = copy_state('reset', full)
        patch_parameters(reset / 'state.h5', {'sweeps': 120})
        run(reset, restore=True, reset=True)
        reference = copy_state('reset_empty_reference', full)
        patch_parameters(reference / 'state.h5', {'sweeps': 120})
        with h5py.File(reference / 'state.h5', 'r+') as f:
            del f[RESP]
            for name in OBS:
                del f[f'{BASE}/measurements/{name}']
        run(reference, restore=True, reset=True)
        identical_response(reset / 'state.h5', reference / 'state.h5')
        identical_configuration(reset / 'state.h5', reference / 'state.h5')
        h = snapshot(reset / 'state.h5')
        assert len(h['chi_s']) == 3 and int(h['partial_count']) == 27
        np.testing.assert_array_equal(h['end_sweep'], [48, 79, 110])
        with h5py.File(reset / 'state.h5', 'r') as f:
            assert int(f[f'{BASE}/measurements/{OBS[0]}/count'][()]) == 3
            assert int(f[f'{BASE}/measurements/Number_of_particles/count'][()]) == 120
        passed('reset_discards_old_history_and_partial_sums', blocks=3, partial_sweeps=27)

        burnin = fresh('burnin_only', thermalization=100000000, sweeps=1000, runtimelimit=1)
        run(burnin)
        run(burnin, restore=True)
        h = snapshot(burnin / 'state.h5')
        assert len(h['chi_s']) == 0 and int(h['partial_count']) == 0
        assert all(float(h[key]) == 0 for key in ('offsite_sum', 'density_sum', 'w2sum'))
        with h5py.File(burnin / 'state.h5', 'r') as f:
            assert int(f[BASE + '/checkpoint/sweeps'][()]) < 100000000
        passed('burnin_empty_history_save_restore')

        bad_interval = copy_state('changed_interval_guard', full)
        patch_parameters(bad_interval / 'state.h5', {'Nresponse_measure': 32})
        run(bad_interval, restore=True, expect='Unsupported response schema or changed response interval')
        bad_disable = copy_state('disabled_history_guard', full)
        patch_parameters(bad_disable / 'state.h5', {'Nresponse_measure': 0})
        run(bad_disable, restore=True, expect='refusing to discard it silently')
        bad_limit = fresh('history_limit_guard', sweeps=30, Nresponse_max_blocks=2)
        run(bad_limit, expect='Response history limit reached')
        bad_boundary = fresh('boundary_guard', pbcx=0)
        run(bad_boundary, expect='requires BoseHubbard on a periodic, single-basis lattice')
        passed('invalid_response_configuration_guards', cases=4)
        report['status'] = 'passed'
    except Exception:
        report['status'] = 'failed'
        report['traceback'] = traceback.format_exc()
        raise
    finally:
        report['elapsed_seconds'] = round(time.monotonic() - started, 3)
        report['finished_utc'] = datetime.now(timezone.utc).isoformat()
        target = folder / 'report.json'
        target.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
        print('REPORT ' + str(target), flush=True)


if __name__ == '__main__':
    main()
