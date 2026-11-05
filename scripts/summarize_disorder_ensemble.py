#!/usr/bin/env python3
"""Equal-weight averages over independent, equilibrated disorder realizations.

Usage: python scripts/summarize_disorder_ensemble.py path/to/manifest.json
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import h5py
import numpy as np

REQUIRED = {'Total_Energy', 'Number_of_particles', 'Number_of_particles_squared', 'Density_Distribution'}


def _json_value(value):
    array = np.asarray(value)
    if array.ndim:
        return [_json_value(v) for v in array]
    scalar = array.item()
    if isinstance(scalar, bytes):
        return scalar.decode('utf-8')
    if isinstance(scalar, float) and not np.isfinite(scalar):
        return None
    return scalar


def _statistics(means, errors, counts):
    samples = np.stack(means)
    count = len(samples)
    mean = samples.mean(axis=0)
    sem = samples.std(axis=0, ddof=1) / np.sqrt(count) if count > 1 else np.full_like(mean, np.nan)
    errors = np.stack(errors)
    valid = np.all(np.isfinite(errors) & (errors >= 0), axis=0)
    thermal = np.where(valid, np.sqrt(np.sum(errors**2, axis=0)) / count, np.nan)
    return dict(mean=_json_value(mean), standard_error=_json_value(sem),
                thermal_standard_error=_json_value(thermal), counts=counts,
                realization_means=_json_value(samples))


def summarize_manifest(path: str | Path) -> dict:
    path = Path(path).resolve(strict=True)
    manifest = json.loads(path.read_text(encoding='utf-8'))
    if manifest.get('schema_version') != 1 or manifest.get('kind') != 'independent_disorder_ensemble':
        raise ValueError('Expected an independent_disorder_ensemble schema-1 manifest')
    runs = manifest.get('runs', [])
    if not runs or len(runs) != manifest.get('realization_count'):
        raise ValueError('Manifest must contain every requested disorder realization')
    if manifest.get('status') != 'completed':
        raise ValueError('Ensemble is not completed; refusing a partial disorder average')
    for key in ['id', 'disorder_seed', 'output_file']:
        if len({run[key] for run in runs}) != len(runs):
            raise ValueError(f'Duplicate {key}; independent realizations and separate outputs are required')
    rows, signatures, warnings = [], [], []
    for run in runs:
        if run.get('status') != 'completed' or run.get('returncode') != 0:
            raise ValueError(f"{run['id']} did not complete successfully; no maps are silently omitted")
        output = (path.parent / run['output_file']).resolve()
        if not output.is_relative_to(path.parent):
            raise ValueError('Manifest output path leaves ensemble directory')
        with h5py.File(output, 'r') as data:
            if 'simulation/results' not in data:
                raise ValueError(f"{run['id']} has no production measurements")
            if 'hamiltonian/onsite_potential' not in data:
                raise ValueError(f"{run['id']} has no resolved potential")
            potential = np.asarray(data['hamiltonian/onsite_potential'][()], dtype='<f8')
            digest = hashlib.sha256(potential.tobytes(order='C')).hexdigest()
            if not np.all(np.isfinite(potential)) or digest != run.get('potential_sha256'):
                raise ValueError(f"{run['id']} potential hash does not match the manifest")
            signature = {}
            for key in ['hamiltonian/model', 'hamiltonian/lattice', 'hamiltonian/shape',
                        'hamiltonian/periodic', 'hamiltonian/t_hop', 'hamiltonian/U_on',
                        'hamiltonian/mu', 'hamiltonian/V_nn', 'hamiltonian/nmax',
                        'parameters/beta', 'parameters/canonical']:
                if key in data:
                    signature[key] = _json_value(data[key][()])
            signatures.append(signature)
            groups = data['simulation/results']
            if REQUIRED - set(groups):
                raise ValueError(f"{run['id']} lacks required production results: {sorted(REQUIRED-set(groups))}")
            row = {}
            for name, group in groups.items():
                if 'mean/value' not in group or 'count' not in group or int(group['count'][()]) <= 0:
                    if name in REQUIRED:
                        raise ValueError(f"{run['id']}/{name} has zero or missing measurements")
                    continue
                mean = np.asarray(group['mean/value'][()], dtype=float)
                error = np.asarray(group['mean/error'][()], dtype=float) if 'mean/error' in group else np.full_like(mean, np.nan)
                if not np.all(np.isfinite(mean)) or mean.shape != error.shape:
                    raise ValueError(f"{run['id']}/{name} has invalid mean/error data")
                if not np.all(np.isfinite(error) & (error >= 0)):
                    warnings.append(f"{run['id']}/{name}: thermal error unavailable for some components")
                row[name] = (mean, error, int(group['count'][()]))
            rows.append(row)
    if any(signature != signatures[0] for signature in signatures[1:]):
        raise ValueError('Different physical parameters or lattice geometry across realizations')
    common = set.intersection(*(set(row) for row in rows))
    omitted = sorted(set.union(*(set(row) for row in rows)) - common)
    if omitted:
        warnings.append('Observables missing from some realizations: ' + ', '.join(omitted))
    if len(runs) == 1:
        warnings.append('One realization cannot estimate the uncertainty of a disorder average')
    observables = {}
    for name in sorted(common):
        if len({row[name][0].shape for row in rows}) != 1:
            raise ValueError(f'Different {name} shapes across realizations')
        observables[name] = _statistics([row[name][0] for row in rows],
                                        [row[name][1] for row in rows], [row[name][2] for row in rows])
    variances = []
    for row in rows:
        number, _, count = row['Number_of_particles']
        number2, _, count2 = row['Number_of_particles_squared']
        if number.ndim or number2.ndim or count != count2:
            raise ValueError('Number moments must be scalars measured on the same samples')
        variance = float(number2-number**2)
        tolerance = 64*np.finfo(float).eps*max(1., abs(float(number2)), float(number**2))
        if variance < -tolerance:
            raise ValueError('Negative per-realization number variance; check moment data')
        variances.append(np.asarray(max(variance, 0.)))
    observables['Number_variance'] = _statistics(variances, [np.asarray(np.nan)]*len(runs),
                                                [row['Number_of_particles'][2] for row in rows])
    observables['Number_variance']['definition'] = 'mean over maps of (<N^2>_map - <N>_map^2)'
    observables['Number_variance']['thermal_error_note'] = 'Not propagated: covariance of number moments is unavailable'
    return dict(schema_version=1, manifest=str(path), realizations=len(runs),
                run_ids=[run['id'] for run in runs], weighting='equal_per_realization',
                standard_error_definition='sample SD of realization means / sqrt(realizations); includes thermal noise',
                thermal_standard_error_definition='sqrt(sum(per-map MC errors squared))/realizations; diagnostic component, do not add to standard_error',
                convergence_note='Nonempty production results do not establish equilibration; check thermal and disorder-sample convergence separately',
                warnings=warnings, omitted_observables=omitted, observables=observables)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('manifest', type=Path)
    parser.add_argument('--output', type=Path, help='Default: summary.json next to manifest')
    args = parser.parse_args(argv)
    try:
        result = summarize_manifest(args.manifest)
        output = args.output or args.manifest.parent/'summary.json'
        if output.resolve() == args.manifest.resolve():
            raise ValueError('Summary output must not overwrite the manifest')
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    except (ValueError, OSError, KeyError) as exc:
        parser.error(str(exc))
    print(f"Equal-weight average of {result['realizations']} disorder realizations")
    for name in ['Total_Energy', 'Number_of_particles', 'Number_variance']:
        obs = result['observables'][name]
        print(f"{name}: {obs['mean']:.8g}; ensemble SE={obs['standard_error']}")
    print(f'Saved {output}')
    if result['warnings']:
        print(f"{len(result['warnings'])} diagnostic note(s) in summary.json")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
