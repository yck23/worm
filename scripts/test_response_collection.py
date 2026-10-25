"""Synthetic HDF5 regressions for response-history export; launches no simulations."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import tempfile
import h5py
import numpy as np
import collect_response_measurements as collector


def fixture(directory, all_empty=False):
    metadata = dict(ranks=2, interval=10)
    for rank, target in enumerate(collector.checkpoint_paths(directory, 2)):
        n = 0 if rank or all_empty else 2
        with h5py.File(target, 'w') as f:
            q = f.create_group('parameters')
            for key, value in dict(Lx=3, beta=.3, canonical=13, nmax=5, Nresponse_measure=10).items():
                q[key] = value
            response = f.create_group(collector.BASE + '/checkpoint/response')
            response['interval'] = 10
            response['partial_count'] = 0
            f[collector.BASE + '/checkpoint/sweeps'] = 17 + n * 10
            values = dict(chi_s=[2., 4.], w2=[1., 3.], density=[13/9, 13/9],
                          end_sweep=[27, 37], sweep_count=[10, 10])
            for key, value in values.items():
                if n:
                    response[key] = value
                else:
                    response.create_dataset(key, dtype=np.asarray(value).dtype)
            measurements = f.create_group(collector.BASE + '/measurements')
            if n:
                for name, value in [('Superfluid_susceptibility', 3.), ('Winding_number_squared_sum', 2.)]:
                    measurements[name + '/count'] = n
                    measurements[name + '/mean/value'] = value
    merged = directory / 'run01.h5'
    with h5py.File(merged, 'w') as f:
        q = f.create_group('parameters')
        for key, value in dict(Lx=3, beta=.3, canonical=13, nmax=5, Nresponse_measure=10).items():
            q[key] = value
        if not all_empty:
            for name, value in [('Superfluid_susceptibility', 3.), ('Winding_number_squared_sum', 2.)]:
                f['simulation/results/' + name + '/count'] = 2
                f['simulation/results/' + name + '/mean/value'] = value
    for stream in ('stdout', 'stderr'):
        (directory / f'run01.{stream}.log').write_text('Completed normally.\n')
    return metadata, merged


def expect_failure(action, phrase):
    try:
        action()
    except (ValueError, RuntimeError) as exc:
        assert phrase in str(exc), str(exc)
    else:
        raise AssertionError('Expected rejection: ' + phrase)


def main():
    checks = []
    with tempfile.TemporaryDirectory(prefix='worm-response-export-') as temporary:
        directory = Path(temporary)
        metadata, merged = fixture(directory)
        collector.check_save_warnings(directory, 1)
        collector.export_history(directory, metadata, merged)
        with np.load(directory / 'history.npz', allow_pickle=False) as archive:
            assert archive['rank1_chi_s'].shape == (0,)
            np.testing.assert_array_equal(archive['rank0_chi_s'], [2., 4.])
        checks.append('mixed_populated_and_null_empty_rank_passes')
        good_export = hashlib.sha256((directory / 'history.npz').read_bytes()).hexdigest()

        with h5py.File(merged, 'r+') as f:
            f['simulation/results/Superfluid_susceptibility/count'][()] = 3
        expect_failure(lambda: collector.export_history(directory, metadata, merged), 'count mismatch')
        assert hashlib.sha256((directory / 'history.npz').read_bytes()).hexdigest() == good_export
        checks.append('stale_or_inconsistent_merged_count_rejected_without_overwriting_export')

        metadata, merged = fixture(directory)
        with h5py.File(merged, 'r+') as f:
            f['simulation/results/Winding_number_squared_sum/mean/value'][()] = 2.1
        expect_failure(lambda: collector.export_history(directory, metadata, merged), 'mean mismatch')
        checks.append('inconsistent_merged_mean_rejected')

        metadata, merged = fixture(directory)
        with h5py.File(collector.checkpoint_paths(directory, 2)[0], 'r+') as f:
            del f[collector.BASE + '/measurements/Superfluid_susceptibility']
        expect_failure(lambda: collector.export_history(directory, metadata, merged), 'count/history mismatch')
        checks.append('missing_nonempty_rank_accumulator_rejected')

        metadata, merged = fixture(directory, all_empty=True)
        collector.export_history(directory, metadata, merged)
        result = json.loads((directory / 'history_checks.json').read_text())
        assert all(row['count'] == 0 and row['merged_mean'] is None for row in result['merged_checks'])
        checks.append('all_empty_histories_and_omitted_result_groups_pass')

        (directory / 'run01.stderr.log').write_text('# WARNING on rank 3 saving checkpoint: disk full\n')
        expect_failure(lambda: collector.check_save_warnings(directory, 1), 'Checkpoint-save warning')
        (directory / 'run01.stderr.log').write_text('')
        (directory / 'run01.stdout.log').write_text('# WARNING: Failed to save checkpoint: disk full\n')
        expect_failure(lambda: collector.check_save_warnings(directory, 1), 'Checkpoint-save warning')
        checks.append('mpi_and_serial_checkpoint_save_warnings_rejected')
    target = collector.ROOT / 'validation/collector_checks.json'
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(dict(status='passed', simulations_launched=0,
        test_script='scripts/test_response_collection.py',
        collector_sha256=hashlib.sha256(Path(collector.__file__).read_bytes()).hexdigest(),
        checks=checks), indent=2) + '\n')
    print(f'PASS {len(checks)} collector regression groups; {target}')


if __name__ == '__main__':
    main()
