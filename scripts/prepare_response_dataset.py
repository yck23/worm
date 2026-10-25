"""Rebuild the 27-point final-production baseline from compact checkpoint audits.

No simulation or checkpoint access is performed. Run from any working directory:
    python scripts/prepare_response_dataset.py
The rank jackknife is conditional on independent, equilibrated rank streams.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil

import numpy as np

REPO = Path(__file__).resolve().parents[1]
RELATIVE_TARGET = 0.05
W2_ABSOLUTE_TARGET = 0.02
COARSE_HISTOGRAM_THRESHOLD = 20


def read_csv(path):
    with path.open(newline='', encoding='utf-8-sig') as stream:
        return list(csv.DictReader(stream))


def write_csv(path, rows):
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def truth(value):
    return str(value).lower() == 'true'


def point_key(row):
    return int(row['L']), round(float(row['beta']), 8)


def rank_jackknife(values, counts):
    """Delete one whole rank and recompute its count-weighted ratio estimate."""
    values = np.asarray(values, dtype=float)
    counts = np.asarray(counts, dtype=float)
    if len(values) < 2 or not np.isfinite(values).all() or not (counts > 0).all():
        raise ValueError('Need at least two finite ranks with positive counts')
    weighted_sum = np.dot(values, counts)
    total = counts.sum()
    mean = weighted_sum / total
    leave_one_out = (weighted_sum - values * counts) / (total - counts)
    R = len(values)
    se = np.sqrt((R - 1) / R * np.sum((leave_one_out - leave_one_out.mean()) ** 2))
    return float(mean), float(se), float(np.std(values, ddof=1))


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build(audit_dir, output_root):
    rank_input = audit_dir / 'checkpoint_rank_observables.csv'
    match_input = audit_dir / 'checkpoint_production_matches.csv'
    drift_input = audit_dir / 'drift_audit.csv'
    source_ranks = read_csv(rank_input)
    matches = {point_key(row): row for row in read_csv(match_input)}
    grouped = {}
    for row in source_ranks:
        grouped.setdefault(point_key(row), []).append(row)
    if len(grouped) != 27 or len(source_ranks) != 520 or set(grouped) != set(matches):
        raise ValueError('Expected 27 matched points and 520 rank records')

    data_dir = output_root / 'data'
    provenance_dir = output_root / 'provenance'
    data_dir.mkdir(parents=True, exist_ok=True)
    provenance_dir.mkdir(parents=True, exist_ok=True)
    clean_ranks, summaries, arrays = [], [], {}
    match_checks = ('g_matches_prod03', 'w_matches_prod03', 'all_metadata_match',
                    'all_thermalized', 'all_g_count_fresh', 'all_g_finite', 'all_w_finite')
    rank_checks = ('thermalized', 'parameters_match', 'g_count_fresh', 'g_finite', 'w_finite')

    for (L, beta), records in sorted(grouped.items()):
        records = sorted(records, key=lambda row: int(row['rank']))
        match = matches[L, beta]
        if not all(truth(match[name]) for name in match_checks):
            raise ValueError(f'Unmatched or invalid production at {(L, beta)}')
        rank_ids = [int(row['rank']) for row in records]
        if rank_ids != list(range(len(records))) or len(records) != int(match['ranks']):
            raise ValueError(f'Missing or duplicated ranks at {(L, beta)}')
        N = int(match['N'])
        density = N / L**2
        scale = L**1.75
        chi, winding, g_counts, w_counts = [], [], [], []
        for original in records:
            if not all(truth(original[name]) for name in rank_checks):
                raise ValueError(f'Invalid rank record at {(L, beta, original["rank"])}')
            if int(original['N']) != N or not np.isclose(float(original['mean_N']), N):
                raise ValueError('Canonical particle number disagrees across ranks')
            if original['prod03_path'] != match['prod03_path']:
                raise ValueError('Rank and merged-production provenance disagree')
            # The audit extracted chi directly from each full checkpoint G vector.
            # Its off-site contribution is retained explicitly in this compact dataset.
            offsite = float(original['chi']) - density
            rank_chi = density + offsite
            rank_winding = float(original['Wx2']) + float(original['Wy2'])
            if not np.isclose(rank_chi / scale, float(original['chi_scaled'])):
                raise ValueError('Susceptibility scaling disagrees with source audit')
            if not np.isclose(rank_winding, float(original['W2'])):
                raise ValueError('Winding components disagree with source audit')
            row = dict(original)
            row.update(L=L, beta=beta, T=1 / beta, rank=int(original['rank']), N=N,
                       density=density, nmax=5, production=3,
                       g_count=int(original['g_count']), w_count=int(original['w_count']),
                       G_offsite_sum=offsite, chi=rank_chi, chi_scaled=rank_chi / scale,
                       W2=rank_winding, equilibrium_status='conditional_unverified')
            clean_ranks.append(row)
            chi.append(rank_chi)
            winding.append(rank_winding)
            g_counts.append(row['g_count'])
            w_counts.append(row['w_count'])
        chi_mean, chi_se, chi_sd = rank_jackknife(chi, g_counts)
        W2_mean, W2_se, W2_sd = rank_jackknife(winding, w_counts)
        if sum(g_counts) != int(match['g_count_prod03']) or sum(w_counts) != int(match['w_count_prod03']):
            raise ValueError(f'Counts do not reproduce prod03 at {(L, beta)}')
        if not np.isclose(chi_mean, float(match['chi_prod03']), rtol=1e-10, atol=1e-12):
            raise ValueError(f'Susceptibility does not reproduce prod03 at {(L, beta)}')
        if not np.isclose(W2_mean, float(match['W2_prod03']), rtol=1e-10, atol=1e-12):
            raise ValueError(f'Winding does not reproduce prod03 at {(L, beta)}')
        chi_relative = chi_se / abs(chi_mean) if chi_mean else float('nan')
        W2_relative = W2_se / abs(W2_mean) if W2_mean else float('nan')
        summary = dict(
            L=L, beta=beta, T=1 / beta, N=N, density=density, nmax=5, production=3,
            source_group=match['source_group'], prod03_path=match['prod03_path'],
            ranks=len(records), jackknife_degrees_of_freedom=len(records) - 1,
            g_count=sum(g_counts), g_count_min=min(g_counts), g_count_max=max(g_counts),
            w_count=sum(w_counts), w_count_min=min(w_counts), w_count_max=max(w_counts),
            chi=chi_mean, chi_se=chi_se, chi_rank_sd=chi_sd, chi_relative_se=chi_relative,
            chi_scaled=chi_mean / scale, chi_scaled_se=chi_se / scale,
            W2=W2_mean, W2_se=W2_se, W2_rank_sd=W2_sd, W2_relative_se=W2_relative,
            relative_precision_target=RELATIVE_TARGET,
            chi_precision_5pct=chi_relative <= RELATIVE_TARGET,
            W2_precision_5pct=W2_relative <= RELATIVE_TARGET,
            W2_absolute_se_target=W2_ABSOLUTE_TARGET,
            W2_absolute_precision_met=W2_se <= W2_ABSOLUTE_TARGET,
            W2_precision_relative_or_absolute=W2_se <= max(
                RELATIVE_TARGET * abs(W2_mean), W2_ABSOLUTE_TARGET),
            coarse_histories=min(g_counts) < COARSE_HISTOGRAM_THRESHOLD,
            coarse_histogram_threshold=COARSE_HISTOGRAM_THRESHOLD,
            source_match_verified=True, configured_thermalization_passed=True,
            uncertainty_method='whole_rank_delete_one_jackknife',
            uncertainty_unit='one_standard_error',
            equilibrium_status='conditional_unverified',
            scalar_time_series_available=False, blocking_plateau_established=False)
        summaries.append(summary)
        tag = f'l{L}t{round(100 * beta):03d}'
        arrays[tag + '_rank'] = np.asarray(rank_ids)
        arrays[tag + '_chi'] = np.asarray(chi)
        arrays[tag + '_g_count'] = np.asarray(g_counts)
        arrays[tag + '_W2'] = np.asarray(winding)
        arrays[tag + '_w_count'] = np.asarray(w_counts)

    write_csv(data_dir / 'baseline_rank_observables.csv', clean_ranks)
    write_csv(data_dir / 'baseline_summary.csv', summaries)
    np.savez_compressed(data_dir / 'baseline_scalar_ranks.npz', **arrays)
    shutil.copy2(match_input, provenance_dir / 'source_production_matches.csv')
    shutil.copy2(drift_input, provenance_dir / 'drift_audit.csv')
    manifest = dict(
        dataset='final_prod03_baseline', points=len(summaries), rank_records=len(clean_ranks),
        inputs={str(path.relative_to(REPO)).replace('\\', '/'): sha256(path)
                for path in (rank_input, match_input, drift_input)},
        generator='scripts/prepare_response_dataset.py', generator_sha256=sha256(Path(__file__)),
        production=3, original_checkpoint_access_required=False,
        susceptibility='chi = N/L^2 + sum(G[r] for r != 0)',
        scaled_susceptibility='chi_scaled = chi / L^(7/4)',
        winding='W2 = Wx2 + Wy2',
        weights={'chi': 'complete G histogram counts', 'W2': 'diagonal winding measurement counts'},
        standard_error='Delete one complete rank; recompute weighted mean; sqrt((R-1)/R * sum((leave_one_out-mean_leave_one_out)^2)).',
        nmax=5, nmax_provenance='Every rank matched the selected prod03 parameters including nmax=5 in the checkpoint audit.',
        equilibrium_status='conditional_unverified',
        precision_targets={'relative_SE': RELATIVE_TARGET, 'optional_W2_absolute_SE': W2_ABSOLUTE_TARGET},
        limitations=[
            'Rank streams must be independent and equilibrated; configured thermalization is not evidence of equilibrium.',
            'Only prod03 enters these estimates. Do not attach these SEs to a pooled three-production mean.',
            'Whole-rank scalar resampling includes spatial covariance and covariance of the two winding components.',
            'No scalar time series or demonstrated blocking plateau is available in this baseline.',
            'The fewer-than-20-histograms flag is descriptive, not a statistical convergence test.',
            'Precision flags assess sampling SE only, not bias or equilibrium.',
            'Drift screening is retained separately; failure to reject drift is not proof of stationarity.'])
    (provenance_dir / 'baseline_manifest.json').write_text(
        json.dumps(manifest, indent=2) + '\n', encoding='utf-8')
    (data_dir / 'README.md').write_text(
        '# Final-production baseline data\n\n'
        '`baseline_summary.csv` contains all 27 (L, beta) points. '
        '`baseline_rank_observables.csv` retains all 520 rank summaries and source paths; '
        '`baseline_scalar_ranks.npz` contains the same scalar measurements and their separate counts. '
        'No checkpoint configurations are needed to inspect or reanalyse these files.\n\n'
        'Each point uses only its matched final `prod03` production. `chi` replaces G(0) '
        'with exact canonical density N/L^2. `chi_scaled` divides chi by L^(7/4). '
        '`W2` is Wx^2+Wy^2; it is not divided by two.\n\n'
        'The `_se` columns are one-standard-error estimates from deleting entire ranks '
        'and recomputing count-weighted means. Chi uses `g_count`; W2 uses `w_count`. '
        'The rank SD columns describe rank-to-rank spread and are not the plotted SEs. '
        'The jackknife degrees of freedom are recorded for context, not used to claim calibrated confidence intervals.\n\n'
        'All uncertainties are conditional on independent, equilibrated rank streams. '
        'Neither a scalar blocking plateau nor equilibrium has been established. '
        '`coarse_histories` marks points with fewer than 20 G histograms in at least one rank; '
        'this is a descriptive warning threshold.\n\n'
        'Precision flags compare SE/abs(mean) with 5%. The optional winding criterion '
        'uses SE <= max(0.05*abs(W2), 0.02), so an absolute SE target is available when '
        'W2 is small. These flags do not establish convergence.\n\n'
        'Rebuild with `python scripts/prepare_response_dataset.py`. The compact input audits '
        'under `../audit/` are required; hashes and formulas are recorded in '
        '`../provenance/baseline_manifest.json`. Source checkpoint paths are historical provenance.\n',
        encoding='utf-8')
    print(f'Prepared {len(summaries)} points and {len(clean_ranks)} ranks in {data_dir}')
    print('Chi relative-SE target:', sum(row['chi_precision_5pct'] for row in summaries), '/ 27')
    print('W2 relative-SE target:', sum(row['W2_precision_5pct'] for row in summaries), '/ 27')
    print('W2 relative-or-absolute target:', sum(row['W2_precision_relative_or_absolute'] for row in summaries), '/ 27')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--audit-dir', type=Path, default=REPO / 'results/response_measurements/audit')
    parser.add_argument('--output-root', type=Path, default=REPO / 'results/response_measurements')
    args = parser.parse_args()
    build(args.audit_dir.resolve(), args.output_root.resolve())


if __name__ == '__main__':
    main()
