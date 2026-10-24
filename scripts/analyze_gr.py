"""Collate equal-time Density_Matrix G(r) and audit its measurement statistics.

Run: python scripts/analyze_gr.py [results/clean_Rb87]
These functions support archived analyses; the current notebooks are self-contained.
Raw HDF5 files are opened read-only.
"""
from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path
import re

import h5py
import numpy as np

FOLDER_PATTERN = re.compile(r"l(?P<L>\d+)t(?P<beta>0\d{2})", re.I)


def _scalar(parameters, name, default=np.nan):
    return np.asarray(parameters[name]).item() if name in parameters else default


def _axis_groups(L):
    # The source stores dx + L * dy. Deduplicate the two points at r = L/2.
    for r in range(1, L // 2 + 1):
        indices = sorted({r, (-r) % L, L * r, L * ((-r) % L)})
        yield float(r), np.asarray(indices)


def _shell_groups(L):
    # Exact squared-distance shells, with minimum-image distances and r <= L/2.
    displacement = np.minimum(np.arange(L), L - np.arange(L))
    r2 = (displacement[:, None] ** 2 + displacement[None, :] ** 2).ravel()
    valid = np.flatnonzero((r2 > 0) & (4 * r2 <= L * L))
    for squared in np.unique(r2[valid]):
        yield float(np.sqrt(squared)), valid[r2[valid] == squared]


def collate_gr(data_dir, productions=(1, 2, 3)):
    """Return audited files, per-production curves, and descriptive pooled curves.

    Pooled G is weighted by complete histogram counts, never by HDF5 errors.
    Production SD is descriptive and does not assume independent files.
    Spatial errors cannot be propagated exactly without spatial covariance.
    """
    data_dir = Path(data_dir).resolve()
    audit, arrays, curves = [], {}, []
    folders = sorted(p for p in data_dir.iterdir()
                     if p.is_dir() and FOLDER_PATTERN.fullmatch(p.name))
    if not folders:
        raise ValueError(f"No l<L>t0xx directories in {data_dir}")
    geometry = {}
    for folder in folders:
        match = FOLDER_PATTERN.fullmatch(folder.name)
        L, beta = int(match['L']), int(match['beta']) / 100
        if L not in geometry:
            geometry[L] = {'axis': list(_axis_groups(L)),
                           'shell': list(_shell_groups(L))}
        for production in productions:
            files = sorted(folder.glob(f'*prod{production:02d}.h5'))
            row = dict(folder=folder.name, L=L, beta=beta, temperature=1 / beta,
                       production=production, file=str(files[0]) if len(files) == 1 else '',
                       status='missing_file', detail='', diagonal_count=0,
                       histogram_count=0, Nmeasure=0, Nmeasure2=0,
                       finite_error_fraction=np.nan, error_bin_levels=0,
                       repeated_array=False, include=False)
            audit.append(row)
            if len(files) != 1:
                row['status'] = 'missing_file' if not files else 'ambiguous_file'
                continue
            try:
                with h5py.File(files[0], 'r') as handle:
                    params = {k: v[()] for k, v in handle['parameters'].items()
                              if isinstance(v, h5py.Dataset)}
                    if (int(params['Lx']), int(params['Ly'])) != (L, L):
                        raise ValueError('Folder and HDF5 lattice sizes disagree')
                    if not np.isclose(float(params['beta']), beta):
                        raise ValueError('Folder and HDF5 beta disagree')
                    if not bool(params['pbcx']) or not bool(params['pbcy']):
                        raise ValueError('Minimum-image analysis requires periodic boundaries')
                    results = handle['simulation/results']
                    row.update(diagonal_count=int(results['Winding_number_squared/count'][()]),
                               Nmeasure=int(params['Nmeasure']), Nmeasure2=int(params['Nmeasure2']))
                    if 'Density_Matrix' not in results:
                        row.update(status='missing_G', detail='No stored Density_Matrix')
                        continue
                    obs = results['Density_Matrix']
                    mean = np.asarray(obs['mean/value'][()], dtype=float).reshape(-1)
                    error = (np.asarray(obs['mean/error'][()], dtype=float).reshape(-1)
                             if 'mean/error' in obs else np.full(mean.shape, np.nan))
                    count = int(obs['count'][()])
                    if mean.size != L * L or error.size != mean.size or count <= 0:
                        raise ValueError('Invalid Density_Matrix shape or count')
                    if not np.all(np.isfinite(mean)):
                        raise ValueError('Nonfinite Density_Matrix means')
                    row.update(histogram_count=count,
                               finite_error_fraction=float(np.isfinite(error).mean()),
                               error_bin_levels=int(obs['error_bins'].shape[0])
                               if 'error_bins' in obs else 0)
                    arrays[(folder.name, production)] = dict(
                        mean=mean.copy(), error=error.copy(), count=count,
                        tau=np.asarray(obs['tau'][()]).copy() if 'tau' in obs else np.array([]))
                    # Production logs use --reset-statistics. Under that protocol,
                    # complete histogram blocks cannot cover more sweeps than recorded.
                    maximum_sweeps = (row['diagonal_count'] + 1) * row['Nmeasure'] - 1
                    if count * row['Nmeasure2'] > maximum_sweeps:
                        row.update(status='stale_count_inconsistent',
                                   detail='Histogram count exceeds fresh production sweep budget')
                    else:
                        row.update(status='finite_errors' if np.isfinite(error).all()
                                   else 'error_unavailable', include=True)
            except (OSError, KeyError, ValueError, TypeError) as exc:
                row.update(status='invalid', detail=str(exc), include=False)

    for folder in folders:
        group = [row for row in audit if row['folder'] == folder.name
                 and (folder.name, row['production']) in arrays]
        for i, left in enumerate(group):
            for right in group[i + 1:]:
                a = arrays[(folder.name, left['production'])]['mean']
                b = arrays[(folder.name, right['production'])]['mean']
                if np.array_equal(a, b):
                    left['repeated_array'] = right['repeated_array'] = True
                    left['include'] = right['include'] = False
                    for item in (left, right):
                        if item['status'] != 'stale_count_inconsistent':
                            item.update(status='duplicate_G', detail='Identical G in another production')

    for row in audit:
        if not row['include']:
            continue
        mean = arrays[(row['folder'], row['production'])]['mean']
        for kind, groups in geometry[row['L']].items():
            for r, indices in groups:
                curves.append(dict(folder=row['folder'], L=row['L'], beta=row['beta'],
                                   temperature=row['temperature'], production=row['production'],
                                   kind=kind, r=r, n_displacements=len(indices),
                                   G=float(mean[indices].mean()),
                                   histogram_count=row['histogram_count'],
                                   finite_site_errors=row['status'] == 'finite_errors'))

    groups = defaultdict(list)
    for row in curves:
        groups[(row['folder'], row['kind'], row['r'])].append(row)
    summary = []
    for (folder, kind, r), group in groups.items():
        values = np.asarray([row['G'] for row in group])
        counts = np.asarray([row['histogram_count'] for row in group])
        first = group[0]
        summary.append(dict(folder=folder, L=first['L'], beta=first['beta'],
                            temperature=first['temperature'], kind=kind, r=r,
                            n_displacements=first['n_displacements'], n_productions=len(group),
                            stored_histogram_count_sum=int(counts.sum()),
                            G_mean=float(np.average(values, weights=counts)),
                            G_production_sd=float(np.std(values, ddof=1)) if len(group) > 1 else np.nan,
                            all_site_errors_finite=all(row['finite_site_errors'] for row in group)))
    summary.sort(key=lambda row: (row['L'], row['beta'], row['kind'], row['r']))
    audit.sort(key=lambda row: (row['L'], row['beta'], row['production']))
    return dict(data_dir=data_dir, audit=audit, arrays=arrays, curves=curves, summary=summary)


def plot_gr_overview(result, kind='axis', normalized=False):
    """Compare every usable production and its count-weighted mean at each L."""
    import matplotlib.pyplot as plt
    sizes = sorted({row['L'] for row in result['audit']})
    betas = sorted({row['beta'] for row in result['audit']})
    colors = {b: plt.get_cmap('viridis')(i / max(1, len(betas) - 1)) for i, b in enumerate(betas)}
    figure, axes = plt.subplots(1, len(sizes), figsize=(5 * len(sizes), 4.5),
                                squeeze=False, constrained_layout=True)
    for ax, L in zip(axes.ravel(), sizes):
        for beta in betas:
            selected = [row for row in result['summary']
                        if row['L'] == L and row['beta'] == beta and row['kind'] == kind]
            if not selected:
                continue
            r = np.asarray([row['r'] for row in selected])
            y = np.asarray([row['G_mean'] for row in selected])
            for production in sorted({row['production'] for row in result['audit']}):
                block = [row for row in result['curves'] if row['L'] == L
                         and row['beta'] == beta and row['kind'] == kind
                         and row['production'] == production]
                if block:
                    block_y = np.asarray([row['G'] for row in block])
                    if normalized:
                        block_y = block_y / block_y[0]
                    ax.loglog([row['r'] for row in block], block_y,
                              color=colors[beta], alpha=.18, linewidth=.8)
            if normalized:
                y = y / y[0]
            finite = all(row['all_site_errors_finite'] for row in selected)
            ax.loglog(r, y, '-' if finite else '--', color=colors[beta], linewidth=1.6,
                      label=rf'$\beta={beta:.2f}$' + ('' if finite else ' *'))
        if normalized:
            r = np.arange(1, L // 2 + 1)
            ax.loglog(r, r ** (-.25), ':', color='0.3', label=r'$r^{-1/4}$ shape reference')
        ax.set(title=f'L = {L}', xlabel='Minimum-image distance r',
               ylabel='G(r) / G(1)' if normalized else 'G(r)')
        ax.grid(alpha=.2, which='both')
        ax.legend(fontsize=8)
    figure.suptitle(f'G(r), {kind} average: thin = production; thick = histogram-count mean\n'
                   '* / dashed = stored site errors unavailable; suspect stale arrays excluded')
    return figure


def plot_gr_point(result, folder, kind='axis'):
    """Inspect production spread and a descriptive local log slope for one point."""
    import matplotlib.pyplot as plt
    summary = [row for row in result['summary'] if row['folder'] == folder and row['kind'] == kind]
    if not summary:
        raise ValueError(f'No usable {kind} G(r) for {folder}; inspect gr_audit')
    r = np.asarray([row['r'] for row in summary])
    y = np.asarray([row['G_mean'] for row in summary])
    figure, axes = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
    for production in sorted({row['production'] for row in result['audit']}):
        block = [row for row in result['curves'] if row['folder'] == folder
                 and row['kind'] == kind and row['production'] == production]
        if block:
            axes[0].loglog([row['r'] for row in block], [row['G'] for row in block],
                           alpha=.6, label=f'prod{production:02d}')
    axes[0].loglog(r, y, 'k-', linewidth=2, label='Histogram-count mean')
    axes[0].loglog(r, y[0] * (r / r[0]) ** (-.25), 'k:', label=r'$r^{-1/4}$ reference')
    axes[0].set(xlabel='Distance r', ylabel='G(r)')
    axes[0].legend(fontsize=8)
    # A secant over about a factor of two in r is less noisy than adjacent sites.
    lower = np.arange(1, int(r.max() // 2) + 1)
    upper = 2 * lower
    lo, hi = np.interp(lower, r, y), np.interp(upper, r, y)
    valid = (lo > 0) & (hi > 0)
    eta = -np.log(hi[valid] / lo[valid]) / np.log(2)
    axes[1].semilogx(np.sqrt(lower[valid] * upper[valid]), eta, 'o-', markersize=3)
    axes[1].axhline(.25, color='0.3', linestyle=':', label='1/4 reference')
    axes[1].set(xlabel=r'Geometric midpoint of [r, 2r]',
                ylabel=r'$-\log[G(2r)/G(r)]/\log 2$ (descriptive)')
    axes[1].legend(fontsize=8)
    for ax in axes:
        ax.grid(alpha=.2, which='both')
    figure.suptitle(f'{folder}, {kind}: production curves and local slope\n'
                   'No fitted exponent or propagated spatial uncertainty')
    return figure


def _write_csv(path, rows):
    if not rows:
        return
    with Path(path).open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def export_gr(result, output_dir=None):
    output_dir = Path(output_dir) if output_dir else result['data_dir'] / 'analysis' / 'gr'
    output_dir.mkdir(parents=True, exist_ok=True)
    for name in ('audit', 'curves', 'summary'):
        _write_csv(output_dir / f'gr_{name}.csv', result[name])
    payload = {}
    for (folder, production), data in result['arrays'].items():
        for key in ('mean', 'error', 'tau'):
            payload[f'{folder}_prod{production:02d}_{key}'] = data[key]
    np.savez_compressed(output_dir / 'density_matrix_arrays.npz', **payload)
    figures = {'gr_axis_all': plot_gr_overview(result),
               'gr_axis_normalized': plot_gr_overview(result, normalized=True),
               'gr_shell_all': plot_gr_overview(result, kind='shell')}
    import matplotlib.pyplot as plt
    for name, figure in figures.items():
        figure.savefig(output_dir / f'{name}.png', dpi=160, bbox_inches='tight')
        plt.close(figure)
    return output_dir


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('data_dir', nargs='?', default='results/clean_Rb87')
    parser.add_argument('--output-dir')
    args = parser.parse_args()
    result = collate_gr(args.data_dir)
    output = export_gr(result, args.output_dir)
    from collections import Counter
    print(dict(Counter(row['status'] for row in result['audit'])))
    print(f'Exported G(r) analysis to {output}')


if __name__ == '__main__':
    main()
