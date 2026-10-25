"""Rank-separated analysis of joint winding and equal-time susceptibility histories.

NPZ scalars: L, beta, N, nmax, interval. For every rankK: chi_s, w2 (or
w2sum), density, end_sweep, sweep_count. Prefix each vector with rankK_.
Rebase reset windows onto a monotonic cumulative sweep coordinate. Histories
sharing an RNG lineage must not be supplied as separate independent ranks.
Outputs use [W2sum, chi_s/L**1.75]. Numerical diagnostics do not establish
equilibration, independence, cutoff convergence, or a transition temperature.
"""
from __future__ import annotations
import argparse
import csv
from dataclasses import dataclass
import glob
from pathlib import Path
import re
import numpy as np


def _portable_path(path):
    path = Path(path).resolve()
    root = Path(__file__).resolve().parents[1]
    try:
        return path.relative_to(root).as_posix()
    except ValueError:
        return path.as_posix()


@dataclass
class RankHistory:
    rank: int
    values: np.ndarray  # W2sum, chi_s/L**1.75
    counts: np.ndarray
    end_sweep: np.ndarray
    density: np.ndarray


def load_history(path):
    """Return metadata and validated chronological rank histories."""
    with np.load(path, allow_pickle=False) as data:
        meta = {key: np.asarray(data[key]).item()
                for key in ('L', 'beta', 'N', 'nmax', 'interval')}
        for key in ('L', 'N', 'nmax', 'interval'):
            if not np.isfinite(meta[key]) or meta[key] != int(meta[key]):
                raise ValueError(f'{key} must be a finite integer')
            meta[key] = int(meta[key])
        if (meta['L'] <= 0 or meta['beta'] <= 0 or not np.isfinite(meta['beta'])
                or meta['interval'] <= 0 or meta['nmax'] < 0 or meta['N'] < -1
                or meta['N'] > meta['nmax']*meta['L']**2):
            raise ValueError('Invalid physical or interval metadata')
        ranks = sorted(int(m.group(1)) for key in data.files
                       if (m := re.fullmatch(r'rank(\d+)_chi_s', key)))
        if not ranks:
            raise ValueError('No rank histories found')
        histories = []
        for rank in ranks:
            prefix = f'rank{rank}_'
            wkey = prefix + ('w2' if prefix+'w2' in data else 'w2sum')
            arrays = [np.asarray(data[key], dtype=float) for key in
                      (wkey, prefix+'chi_s', prefix+'density',
                       prefix+'end_sweep', prefix+'sweep_count')]
            if any(x.ndim != 1 for x in arrays):
                raise ValueError(f'Rank {rank}: histories must be one-dimensional')
            if not arrays[0].size or any(x.size != arrays[0].size for x in arrays):
                raise ValueError(f'Rank {rank}: empty or unequal-length histories')
            if not all(np.isfinite(x).all() for x in arrays):
                raise ValueError(f'Rank {rank}: nonfinite history')
            w2, chi, density, end, counts = arrays
            if (np.any(counts <= 0) or np.any(counts != np.floor(counts))
                    or np.any(end != np.floor(end)) or np.any(end < counts)):
                raise ValueError(f'Rank {rank}: invalid counts/endpoints')
            if np.any(counts != meta['interval']):
                raise ValueError(f'Rank {rank}: sweep_count differs from fixed interval metadata')
            if np.any(np.diff(end) <= 0) or np.any(end[1:]-counts[1:] < end[:-1]):
                raise ValueError(f'Rank {rank}: overlapping/nonmonotonic intervals')
            if np.any(w2 < -1e-12) or np.any(chi < -1e-12) or np.any(density < 0):
                raise ValueError(f'Rank {rank}: negative response/density')
            if meta['N'] >= 0 and not np.allclose(
                    density, meta['N']/meta['L']**2, rtol=1e-10, atol=1e-12):
                raise ValueError(f'Rank {rank}: canonical density identity failed')
            histories.append(RankHistory(rank, np.column_stack(
                (w2, chi/meta['L']**1.75)), counts, end, density))
    meta['source'] = _portable_path(path)
    meta['fixed_interval_schema'] = True
    return meta, histories


def weighted_cluster_jackknife(numerators, weights):
    """Weighted mean and delete-one-cluster covariance; weights are counts."""
    numerators, weights = np.asarray(numerators, float), np.asarray(weights, float)
    if (numerators.ndim != 2 or weights.shape != (len(numerators),)
            or np.any(weights <= 0) or not np.isfinite(numerators).all()
            or not np.isfinite(weights).all()):
        raise ValueError('Invalid cluster sums/counts')
    mean = numerators.sum(axis=0)/weights.sum()
    if len(weights) < 2:
        return mean, np.full((numerators.shape[1],)*2, np.nan)
    deleted = (numerators.sum(axis=0)-numerators)/(weights.sum()-weights)[:, None]
    centered = deleted-deleted.mean(axis=0)
    return mean, (len(weights)-1)/len(weights)*centered.T @ centered


def _rank_sums(histories):
    return (np.array([(h.values*h.counts[:, None]).sum(axis=0) for h in histories]),
            np.array([h.counts.sum() for h in histories]))


def _se(covariance):
    return np.sqrt(np.maximum(np.diag(covariance), 0))


def blocking_curve(histories):
    """Block within each rank, then combine independent-rank mean covariances."""
    rows, size = [], 1
    while min(len(h.counts)//size for h in histories) >= 2:
        sums, totals, covariances, numbers = [], [], [], []
        for h in histories:
            blocks = len(h.counts)//size
            stop = blocks*size
            counts = h.counts[:stop].reshape(blocks, size).sum(axis=1)
            nums = (h.values[:stop]*h.counts[:stop, None]).reshape(blocks, size, 2).sum(axis=1)
            _, covariance = weighted_cluster_jackknife(nums, counts)
            sums.append(nums.sum(axis=0))
            totals.append(counts.sum())
            covariances.append(covariance)
            numbers.append(blocks)
        weights = np.asarray(totals)/np.sum(totals)
        covariance = np.sum(weights[:, None, None]**2*np.asarray(covariances), axis=0)
        mean, errors = np.sum(sums, axis=0)/np.sum(totals), _se(covariance)
        rows.append(dict(block_intervals=size, min_blocks_per_rank=min(numbers),
                         total_blocks=sum(numbers), retained_sweeps=sum(totals),
                         W2sum=mean[0], scaled_chi=mean[1], W2sum_se=errors[0],
                         scaled_chi_se=errors[1], covariance=covariance[0, 1]))
        size *= 2
    return rows


def autocorrelation_diagnostics(histories):
    """Initial-positive paired-lag estimate from separate within-rank FFTs.

    Short histories can miss slow modes. Irregular cadence makes this diagnostic
    unavailable; it does not invalidate the whole-rank jackknife mean/error.
    """
    if any(not np.all(h.counts == h.counts[0])
           or np.any(np.diff(h.end_sweep) != h.counts[1:]) for h in histories):
        return np.full(2, np.nan), False
    if len({h.counts[0] for h in histories}) != 1:
        return np.full(2, np.nan), False
    maxlag = min(len(h.counts) for h in histories)//4
    if maxlag < 4:
        return np.full(2, np.nan), False
    numerator, pairs = np.zeros((maxlag+1, 2)), np.zeros(maxlag+1)
    for h in histories:
        values = h.values-h.values.mean(axis=0)
        n = len(values)
        length = 1 << (2*n-1).bit_length()
        spectrum = np.fft.rfft(values, n=length, axis=0)
        numerator += np.fft.irfft(spectrum.conj()*spectrum, n=length, axis=0)[:maxlag+1]
        pairs += n-np.arange(maxlag+1)
    covariance, tau = numerator/pairs[:, None], np.full(2, np.nan)
    for col in range(2):
        if covariance[0, col] <= 0:
            continue
        rho, total = covariance[:, col]/covariance[0, col], 0.0
        for lag in range(1, maxlag, 2):
            pair = rho[lag]+rho[lag+1]
            if pair <= 0:
                break
            total += pair
        tau[col] = max(0.5, 0.5+total)
    return tau, bool(np.isfinite(tau).all() and maxlag >= 5*np.max(tau))


def paired_half_drift(histories):
    """Late-minus-early difference; preserve pairing in a whole-rank jackknife."""
    early, late, ew, lw = [], [], [], []
    for h in histories:
        if len(h.counts) < 2:
            return np.full(2, np.nan), np.full((2, 2), np.nan)
        split = min(int(np.searchsorted(np.cumsum(h.counts), h.counts.sum()/2)+1),
                    len(h.counts)-1)
        early.append((h.values[:split]*h.counts[:split, None]).sum(axis=0))
        late.append((h.values[split:]*h.counts[split:, None]).sum(axis=0))
        ew.append(h.counts[:split].sum())
        lw.append(h.counts[split:].sum())
    early, late, ew, lw = map(np.asarray, (early, late, ew, lw))
    difference = late.sum(axis=0)/lw.sum()-early.sum(axis=0)/ew.sum()
    if len(histories) < 2:
        return difference, np.full((2, 2), np.nan)
    deleted = ((late.sum(axis=0)-late)/(lw.sum()-lw)[:, None]
               -(early.sum(axis=0)-early)/(ew.sum()-ew)[:, None])
    centered = deleted-deleted.mean(axis=0)
    return difference, (len(histories)-1)/len(histories)*centered.T @ centered


def analyze_history(path, relative_target=0.05, absolute_w2_target=0.02,
                    min_ranks=8, min_intervals=64):
    """Return summary, blocking rows, rank rows. No flag certifies equilibration."""
    meta, histories = load_history(path)
    sums, counts = _rank_sums(histories)
    mean, covariance = weighted_cluster_jackknife(sums, counts)
    errors, weights = _se(covariance), counts/counts.sum()
    curve = blocking_curve(histories)
    tau, tau_ok = autocorrelation_diagnostics(histories)
    drift, drift_cov = paired_half_drift(histories)
    drift_se = _se(drift_cov)
    drift_z = np.divide(np.abs(drift), drift_se, out=np.full(2, np.inf), where=drift_se > 0)
    drift_z[(drift == 0) & (drift_se == 0)] = 0
    drift_z[~np.isfinite(drift_se)] = np.nan
    plateau, blocking_se = None, np.full(2, np.nan)
    # Check the coarsest adequately populated triple, not an early apparent
    # plateau contradicted by later, still well-populated blocking levels.
    windows = [curve[i-2:i+1] for i in range(2, len(curve))
               if curve[i]['min_blocks_per_rank'] >= 4 and curve[i]['total_blocks'] >= 50]
    if tau_ok and windows:
        window = windows[-1]
        se_values = np.array([[r['W2sum_se'], r['scaled_chi_se']] for r in window])
        adequate = window[0]['block_intervals'] >= 4*np.max(tau)
        stable = (np.all(se_values > 0)
                  and np.all(se_values.max(axis=0)/se_values.min(axis=0) <= 1.25))
        if adequate and stable:
            plateau = window[-1]['block_intervals']
            blocking_se = se_values.max(axis=0)
    enough = (len(histories) >= min_ranks
              and min(len(h.counts) for h in histories) >= min_intervals
              and 1/np.sum(weights**2) >= 0.75*min_ranks)
    drift_ok = bool(np.isfinite(drift_z).all() and np.all(drift_z <= 3))
    targets = np.array([max(relative_target*abs(mean[0]), absolute_w2_target),
                        relative_target*abs(mean[1])])
    rank_precision = bool(np.isfinite(errors).all() and mean[1] > 0
                          and np.all(errors <= targets))
    precision = bool(rank_precision and np.isfinite(blocking_se).all()
                     and np.all(blocking_se <= targets))
    diagnostics = bool(enough and plateau is not None and drift_ok)
    flags = []
    if not enough:
        flags.append('too_few_ranks_or_intervals')
    if not tau_ok:
        flags.append('autocorrelation_window_or_cadence_inadequate')
    if plateau is None:
        flags.append('blocking_plateau_not_demonstrated')
    if not drift_ok:
        flags.append('paired_half_drift_requires_review')
    if not precision:
        flags.append('precision_target_not_met')
    summary = dict(meta, ranks=len(histories), effective_rank_count=1/np.sum(weights**2),
                   intervals=sum(len(h.counts) for h in histories),
                   min_intervals_per_rank=min(len(h.counts) for h in histories),
                   sweeps=int(counts.sum()), W2sum=mean[0], W2sum_se=errors[0],
                   scaled_chi=mean[1], scaled_chi_se=errors[1], covariance=covariance[0, 1],
                   W2sum_relative_se=errors[0]/abs(mean[0]) if mean[0] else np.nan,
                   scaled_chi_relative_se=errors[1]/abs(mean[1]) if mean[1] else np.nan,
                   W2sum_tau_intervals=tau[0], scaled_chi_tau_intervals=tau[1],
                   W2sum_half_drift=drift[0], W2sum_half_drift_se=drift_se[0],
                   W2sum_half_drift_z=drift_z[0], scaled_chi_half_drift=drift[1],
                   scaled_chi_half_drift_se=drift_se[1], scaled_chi_half_drift_z=drift_z[1],
                   plateau_block_intervals=plateau or 0, W2sum_blocking_se=blocking_se[0],
                   scaled_chi_blocking_se=blocking_se[1], rank_precision_target_met=rank_precision,
                   precision_target_met=precision,
                   statistical_diagnostics_met=diagnostics, ready_for_review=precision and diagnostics,
                   convergence_certified=False, method='whole_rank_delete_one_jackknife',
                   status='diagnostics_met_review_equilibration' if precision and diagnostics
                          else 'more_sampling_or_checks_required', flags=';'.join(flags))
    ranks = [dict(source=meta['source'], L=meta['L'], beta=meta['beta'], rank=h.rank,
                  intervals=len(h.counts), sweeps=int(count), W2sum=num[0]/count,
                  scaled_chi=num[1]/count, density=np.average(h.density, weights=h.counts))
             for h, num, count in zip(histories, sums, counts)]
    return summary, [dict(source=meta['source'], L=meta['L'], beta=meta['beta'], **r)
                     for r in curve], ranks


def self_test():
    """Meaningful seeded checks of joint covariance and temporal blocking."""
    rng, histories = np.random.default_rng(91625), []
    for rank in range(16):
        innovations = rng.normal(size=8192)
        series = np.empty_like(innovations)
        series[0] = innovations[0]/np.sqrt(1-0.9**2)
        for i in range(1, len(series)):
            series[i] = 0.9*series[i-1]+innovations[i]
        histories.append(RankHistory(rank, np.column_stack((series+20, 3*series+60)),
                                     np.ones(len(series)), np.arange(1, len(series)+1),
                                     np.ones(len(series))))
    _, covariance = weighted_cluster_jackknife(*_rank_sums(histories))
    assert np.allclose(covariance[1, 1], 9*covariance[0, 0])
    assert np.allclose(covariance[0, 1], 3*covariance[0, 0])
    curve = blocking_curve(histories)
    block128 = next(r for r in curve if r['block_intervals'] == 128)
    assert block128['W2sum_se'] > 3*curve[0]['W2sum_se'], 'AR(1) correlation lost'
    expected = np.sqrt((1/(1-0.9**2))*((1+0.9)/(1-0.9))/(16*8192))
    assert 0.65 < block128['W2sum_se']/expected < 1.35
    shifted = [RankHistory(h.rank, h.values+20*h.rank, h.counts, h.end_sweep, h.density)
               for h in histories]
    shifted128 = next(r for r in blocking_curve(shifted) if r['block_intervals'] == 128)
    assert np.isclose(block128['W2sum_se'], shifted128['W2sum_se']), 'Blocks crossed ranks'
    tau, adequate = autocorrelation_diagnostics(histories)
    assert adequate and np.all((tau > 6) & (tau < 14))
    print('Self-test passed: covariance/scaling, AR(1) blocking, separate rank histories.')


def _write_csv(path, rows):
    if not rows:
        return
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with Path(path).open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('paths', nargs='*', help='NPZ files or quoted glob patterns')
    parser.add_argument('--output-dir', type=Path, default=Path('results/response_measurements/analysis'))
    parser.add_argument('--relative-target', type=float, default=0.05)
    parser.add_argument('--absolute-w2-target', type=float, default=0.02)
    parser.add_argument('--min-ranks', type=int, default=8)
    parser.add_argument('--min-intervals', type=int, default=64)
    parser.add_argument('--self-test', action='store_true')
    args = parser.parse_args()
    if args.self_test:
        self_test()
        if not args.paths:
            return
    if args.relative_target <= 0 or args.absolute_w2_target < 0 or args.min_ranks < 2 or args.min_intervals < 2:
        parser.error('Invalid precision/minimum-history settings')
    patterns = args.paths or ['results/response_measurements/sampling/*/history.npz']
    files = sorted({Path(name) for pattern in patterns for name in glob.glob(pattern)})
    if not files:
        parser.error('No history NPZ files matched')
    summaries, curves, ranks = [], [], []
    for path in files:
        try:
            summary, curve, rank_rows = analyze_history(path, args.relative_target,
                args.absolute_w2_target, args.min_ranks, args.min_intervals)
            summaries.append(summary)
            curves.extend(curve)
            ranks.extend(rank_rows)
        except (ValueError, KeyError, OSError) as exc:
            summaries.append(dict(source=_portable_path(path), status='invalid_history', flags=str(exc),
                                  precision_target_met=False, statistical_diagnostics_met=False,
                                  ready_for_review=False, convergence_certified=False))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    _write_csv(args.output_dir/'summary.csv', summaries)
    _write_csv(args.output_dir/'blocking.csv', curves)
    _write_csv(args.output_dir/'ranks.csv', ranks)
    ready = sum(bool(row.get('ready_for_review')) for row in summaries)
    print(f'Analyzed {len(files)} points; {ready} meet numerical diagnostics, '
          f'pending independent equilibration review. Outputs: {args.output_dir}')


if __name__ == '__main__':
    main()
