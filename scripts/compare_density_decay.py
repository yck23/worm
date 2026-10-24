"""Two-model density-matrix fits with explicitly approximate diagonal-error chi-square.

The available archives do not contain the full covariance needed for a calibrated
correlated-data goodness-of-fit test. All nominal p-values are conditional on
independent spatial directions, productions, and distances.
"""
from __future__ import annotations
import csv
import json
from pathlib import Path
import warnings

import numpy as np
from scipy.optimize import curve_fit, OptimizeWarning
from scipy.stats import chi2
from analyze_gr import _axis_groups, _shell_groups

MODEL_NAMES = ('algebraic', 'exponential')
MODEL_COLORS = {'algebraic': 'tab:blue', 'exponential': 'tab:red'}
SIZE_COLORS = {32: 'tab:green', 64: 'tab:blue', 128: 'tab:orange'}
ASSUMPTION = 'HDF5 marginal errors; spatial, production, and distance covariances neglected'


def spatial_profile(mean, error, L, kind='axis'):
    """Spatial mean and diagonal-only propagation of its stored marginal SEMs."""
    mean, error = np.asarray(mean).ravel(), np.asarray(error).ravel()
    if mean.size != L * L or error.size != mean.size:
        raise ValueError('Density-matrix dimensions disagree with L')
    if (error < 0).any():
        raise ValueError('Stored marginal errors must be nonnegative')
    if kind == 'axis':
        groups = list(_axis_groups(L))
    elif kind == 'shell':
        groups = list(_shell_groups(L))
    elif kind == 'positive_x':
        groups = [(float(r), np.array([r])) for r in range(1, L // 2 + 1)]
    else:
        raise ValueError('kind must be axis, shell, or positive_x')
    r = np.array([distance for distance, _ in groups])
    values = np.array([mean[indices].mean() for _, indices in groups])
    sigma = np.array([np.sqrt(np.sum(error[indices] ** 2)) / len(indices)
                      for _, indices in groups])
    return r, values, sigma


def pooled_profile(L, beta, entry, kind='axis'):
    arrays = [entry['result']['arrays'][(entry['folder'], p)] for p in (1, 2, 3)]
    counts = np.array([a['count'] for a in arrays], dtype=float)
    weights = counts / counts.sum()
    individual = [spatial_profile(a['mean'], a['error'], L, kind) for a in arrays]
    r = individual[0][0]
    if not all(np.array_equal(item[0], r) for item in individual):
        raise ValueError('Production distance grids differ')
    y = np.average(np.stack([item[1] for item in individual]), axis=0, weights=weights)
    sigma = np.sqrt(np.sum(np.stack([item[2] for item in individual]) ** 2 * weights[:, None] ** 2, axis=0))
    return dict(L=L, beta=beta, temperature=1 / beta, source=entry['source'],
                folder=entry['folder'], kind=kind, r=r, y=y, sigma=sigma,
                histogram_count_sum=int(counts.sum()), productions=3,
                sigma_assumption=ASSUMPTION)


def evaluate(model, r, log_amplitude_ref, decay, r_ref):
    x = np.log(np.asarray(r) / r_ref) if model == 'algebraic' else np.asarray(r) - r_ref
    return np.exp(log_amplitude_ref - decay * x)


def fit_model(r, y, sigma, model):
    """Minimize sum((G-f)/sigma)^2 in raw G space with two parameters.

    Internal amplitude is at a reference separation to reduce parameter correlation.
    Stored sigma is absolute: it is never rescaled by the residual or reduced chi2.
    """
    if model not in MODEL_NAMES:
        raise ValueError(model)
    r, y, sigma = [np.asarray(item, dtype=float) for item in (r, y, sigma)]
    if r.size <= 2 or not np.isfinite(sigma).all() or not (sigma > 0).all():
        raise ValueError('Fit requires >2 points with finite positive sigma')
    r_ref = float(np.exp(np.mean(np.log(r))))
    x = np.log(r / r_ref) if model == 'algebraic' else r - r_ref
    positive = y > 0
    if positive.sum() < 2:
        raise ValueError('Need at least two positive means for initialization')
    initial_slope, initial_intercept = np.polyfit(x[positive], np.log(y[positive]), 1)
    initial_decay = max(float(-initial_slope), 1e-6)

    def function(distance, log_a, decay):
        return evaluate(model, distance, log_a, decay, r_ref)

    def jacobian(distance, log_a, decay):
        values = function(distance, log_a, decay)
        xx = np.log(np.asarray(distance) / r_ref) if model == 'algebraic' else np.asarray(distance) - r_ref
        return np.column_stack((values, -values * xx))

    candidates = []
    warning_messages = []
    for multiplier in (1., .5, 2.):
        decay0 = initial_decay * multiplier
        log_a0 = float(np.mean(np.log(y[positive]) + decay0 * x[positive]))
        try:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter('always', OptimizeWarning)
                params, covariance = curve_fit(function, r, y, sigma=sigma,
                    p0=(log_a0, decay0), absolute_sigma=True,
                    bounds=([-np.inf, 0.], [np.inf, np.inf]), jac=jacobian,
                    method='trf', max_nfev=10000, ftol=1e-11, xtol=1e-11, gtol=1e-11)
            prediction = function(r, *params)
            statistic = float(np.sum(((y - prediction) / sigma) ** 2))
            if np.isfinite(statistic) and np.isfinite(params).all():
                candidates.append((statistic, params, covariance, prediction))
                warning_messages.extend(str(item.message) for item in caught)
        except (ValueError, RuntimeError, FloatingPointError):
            continue
    if not candidates:
        raise RuntimeError(f'{model} fit did not converge')
    statistic, params, covariance, prediction = min(candidates, key=lambda item: item[0])
    log_a, decay = map(float, params)
    nu = int(r.size - 2)
    boundary = decay <= 1e-8
    nominal_se = float(np.sqrt(covariance[1, 1])) if np.isfinite(covariance[1, 1]) and covariance[1, 1] >= 0 else np.nan
    log_A = log_a + decay * (np.log(r_ref) if model == 'algebraic' else r_ref)
    result = dict(model=model, fit_status='ok', n_points=int(r.size), n_parameters=2, dof=nu,
        r_min=float(r.min()), r_max=float(r.max()), r_reference=r_ref,
        amplitude_at_reference=float(np.exp(log_a)), log_amplitude_at_reference=log_a,
        amplitude_A=float(np.exp(log_A)), log_amplitude_A=float(log_A),
        decay_coefficient=decay, decay_nominal_se=nominal_se,
        eta=decay if model == 'algebraic' else np.nan,
        xi=1 / decay if model == 'exponential' and decay > 0 else np.nan,
        xi_nominal_se=nominal_se / decay ** 2 if model == 'exponential' and decay > 0 else np.nan,
        chi_squared_diag=statistic, reduced_chi_squared_diag=statistic / nu,
        p_nominal_diag=float(chi2.sf(statistic, nu)) if not boundary else np.nan,
        boundary_solution=boundary, formal_p_value_valid=False,
        covariance_condition=float(np.linalg.cond(covariance)),
        optimizer_warnings='; '.join(sorted(set(warning_messages))))
    return result, prediction


def analyze_sources(sources, kind='axis', lower_fraction=1/16, upper_fraction=1/4,
                    absolute_r_min=3, min_points=5, upper_checks=(1/8, 1/4, 1/3)):
    profiles, fit_rows, comparisons, residuals, windows, fit_cache = {}, [], [], [], [], {}
    for (L, beta), entry in sorted(sources.items()):
        profile = pooled_profile(L, beta, entry, kind)
        profiles[(L, beta)] = profile
        r, y, sigma = profile['r'], profile['y'], profile['sigma']
        lower = max(absolute_r_min, L * lower_fraction)
        base = dict(L=L, beta=beta, temperature=1/beta, folder=entry['folder'],
                    source=entry['source'], kind=kind, sigma_assumption=ASSUMPTION,
                    histogram_count_sum=profile['histogram_count_sum'])
        primary = {}
        for upper in sorted(set((*upper_checks, upper_fraction))):
            in_window = (r >= lower) & (r <= L * upper)
            usable = in_window & np.isfinite(y) & np.isfinite(sigma) & (sigma > 0)
            for model in MODEL_NAMES:
                if usable.sum() < min_points:
                    windows.append(dict(**base, upper_fraction=upper, model=model,
                        fit_status='insufficient_points', n_points=int(usable.sum()),
                        chi_squared_diag=np.nan, reduced_chi_squared_diag=np.nan,
                        p_nominal_diag=np.nan, dof=max(0, int(usable.sum())-2), eta=np.nan, xi=np.nan))
                    if np.isclose(upper, upper_fraction):
                        raise ValueError(f'Primary window has too few points: {L}, {beta}')
                    continue
                fitted, prediction = fit_model(r[usable], y[usable], sigma[usable], model)
                windows.append(dict(**base, upper_fraction=upper, **fitted))
                if np.isclose(upper, upper_fraction):
                    row = dict(**base, **fitted, n_excluded_sigma=int(in_window.sum()-usable.sum()))
                    fit_rows.append(row)
                    primary[model] = row
                    fit_cache[(L, beta, model)] = row
                    for distance, observed, se, expected in zip(r[usable], y[usable], sigma[usable], prediction):
                        residuals.append(dict(**base, model=model, r=float(distance), observed=float(observed),
                            sigma_diag=float(se), prediction=float(expected),
                            standardized_residual=float((observed-expected)/se),
                            chi_squared_contribution=float(((observed-expected)/se)**2)))
        alg, exp = primary['algebraic'], primary['exponential']
        delta = exp['chi_squared_diag'] - alg['chi_squared_diag']
        comparisons.append(dict(**base, n_points=alg['n_points'], dof=alg['dof'],
            r_min=alg['r_min'], r_max=alg['r_max'], eta=alg['eta'], xi=exp['xi'],
            chi2_algebraic_diag=alg['chi_squared_diag'], chi2_exponential_diag=exp['chi_squared_diag'],
            reduced_chi2_algebraic_diag=alg['reduced_chi_squared_diag'],
            reduced_chi2_exponential_diag=exp['reduced_chi_squared_diag'],
            p_algebraic_nominal=alg['p_nominal_diag'], p_exponential_nominal=exp['p_nominal_diag'],
            delta_chi2_exp_minus_alg_diag=delta,
            lower_diagonal_misfit_model='algebraic' if delta > 0 else 'exponential' if delta < 0 else 'tie',
            formal_p_value_valid=False))
    return dict(profiles=profiles, fits=fit_rows, comparisons=comparisons,
                residuals=residuals, windows=windows, fit_cache=fit_cache,
                settings=dict(kind=kind, lower_fraction=lower_fraction, upper_fraction=upper_fraction,
                              absolute_r_min=absolute_r_min, min_points=min_points,
                              upper_checks=list(upper_checks), sigma_assumption=ASSUMPTION))


def export_analysis(analysis, directory):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    for key, filename in [('fits','model_fits.csv'), ('comparisons','model_comparison.csv'),
                          ('residuals','model_fit_residuals.csv'), ('windows','model_fit_windows.csv')]:
        rows = analysis[key]
        fields = list(dict.fromkeys(k for row in rows for k in row))
        with (directory / filename).open('w', newline='', encoding='utf-8') as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader();writer.writerows(rows)
    payload = {}
    for (L, beta), profile in analysis['profiles'].items():
        prefix=f'l{L}t{round(beta*100):03d}'
        for name in ('r','y','sigma'):
            payload[prefix+'_'+name] = profile[name]
    np.savez_compressed(directory/'model_fit_profiles.npz', **payload)
    settings = dict(analysis['settings'], formal_p_values_calibrated=False,
        models={'algebraic':'A*r**(-eta)', 'exponential':'A*exp(-r/xi)'},
        fit_space='raw G, fixed propagated diagonal HDF5 errors, absolute_sigma=True',
        references=['https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html',
                    'https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.chi2.html'])
    (directory/'model_fit_settings.json').write_text(json.dumps(settings,indent=2)+'\n',encoding='utf-8')
    return directory


def plot_parameters(analysis):
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    panels = [(axes[0,0],'algebraic','eta',r'Algebraic exponent $\eta$'),
              (axes[0,1],'exponential','xi',r'Exponential length $\xi$'),
              (axes[1,0],'algebraic','amplitude_A','Algebraic amplitude A'),
              (axes[1,1],'exponential','amplitude_A','Exponential amplitude A')]
    for ax, model, field, ylabel in panels:
        for L in sorted({row['L'] for row in analysis['fits']}):
            rows = sorted((row for row in analysis['fits'] if row['L']==L and row['model']==model),key=lambda row:row['temperature'])
            ax.plot([row['temperature'] for row in rows],[row[field] for row in rows], 'o-',
                    color=SIZE_COLORS[L],markersize=4,label=f'L = {L}')
        ax.set(xlabel=r'Temperature $T=1/\beta$',ylabel=ylabel)
        if field=='xi':ax.set_yscale('log')
        ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.suptitle('Weighted fits in raw G space; parameter estimates under the diagonal-error approximation')
    return fig


def plot_chisquared(analysis):
    import matplotlib.pyplot as plt
    sizes = sorted({row['L'] for row in analysis['fits']})
    fig, axes = plt.subplots(2,len(sizes),figsize=(5*len(sizes),7),constrained_layout=True,squeeze=False)
    for col,L in enumerate(sizes):
        for model in MODEL_NAMES:
            rows=sorted((row for row in analysis['fits'] if row['L']==L and row['model']==model),key=lambda row:row['temperature'])
            for i,field in enumerate(('chi_squared_diag','reduced_chi_squared_diag')):
                axes[i,col].plot([row['temperature'] for row in rows],[row[field] for row in rows],
                                 'o-',color=MODEL_COLORS[model],markersize=4,label=model.capitalize())
        axes[0,col].axhline(rows[0]['dof'],color='0.5',linestyle=':',label='nominal dof')
        axes[1,col].axhline(1,color='0.5',linestyle=':',label='1 reference')
        for i in (0,1):
            axes[i,col].set(title=f'L = {L}',xlabel=r'Temperature $T=1/\beta$',
                ylabel=r'$\chi^2_{diag}$' if i==0 else r'Reduced $\chi^2_{diag}=\chi^2_{diag}/(n-2)$',yscale='log')
            axes[i,col].grid(alpha=.2,which='both');axes[i,col].legend(fontsize=8)
    fig.suptitle('Chi-squared model comparison: identical data, errors, fit window and two parameters for both models')
    return fig


def plot_pvalues(analysis, floor=1e-14):
    import matplotlib.pyplot as plt
    sizes=sorted({row['L'] for row in analysis['fits']})
    fig,axes=plt.subplots(1,len(sizes),figsize=(5*len(sizes),4.2),constrained_layout=True,squeeze=False)
    for ax,L in zip(axes.flat,sizes):
        for model in MODEL_NAMES:
            rows=sorted((row for row in analysis['fits'] if row['L']==L and row['model']==model),key=lambda row:row['temperature'])
            ax.plot([row['temperature'] for row in rows],[max(floor,row['p_nominal_diag']) for row in rows],
                    'o-',color=MODEL_COLORS[model],markersize=4,label=model.capitalize())
        ax.axhline(.05,color='0.5',linestyle=':',label='0.05 reference')
        ax.set(title=f'L = {L}',xlabel=r'Temperature $T=1/\beta$',ylabel='Nominal p (diagonal approximation)',yscale='log',ylim=(floor/2,1.3))
        ax.grid(alpha=.2,which='both');ax.legend(fontsize=8)
    fig.suptitle(f'Nominal goodness-of-fit p-values; calibration unverified; values below {floor:g} clipped only for display')
    return fig


def plot_delta(analysis):
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(9,5),constrained_layout=True)
    for L in sorted({row['L'] for row in analysis['comparisons']}):
        rows=sorted((row for row in analysis['comparisons'] if row['L']==L),key=lambda row:row['temperature'])
        ax.plot([row['temperature'] for row in rows],[row['delta_chi2_exp_minus_alg_diag'] for row in rows],
                'o-',color=SIZE_COLORS[L],label=f'L = {L}')
    ax.axhline(0,color='0.5',linestyle=':')
    ax.set(xlabel=r'Temperature $T=1/\beta$',ylabel=r'$\chi^2_{exp,diag}-\chi^2_{alg,diag}$',yscale='symlog',
           title='Positive: algebraic has lower weighted misfit; negative: exponential has lower weighted misfit')
    ax.grid(alpha=.2,which='both');ax.legend()
    return fig


def plot_gallery(analysis,L):
    import matplotlib.pyplot as plt
    points=sorted((point for point in analysis['profiles'] if point[0]==L),key=lambda p:1/p[1])
    nrows=int(np.ceil(len(points)/3))
    fig,axes=plt.subplots(nrows,3,figsize=(12,3.2*nrows),constrained_layout=True,squeeze=False)
    for ax,key in zip(axes.flat,points):
        profile=analysis['profiles'][key];r,y,sigma=profile['r'],profile['y'],profile['sigma']
        fitted=analysis['fit_cache'][(L,key[1],'algebraic')]
        ax.plot(r,y,'.',color='0.7',markersize=3,label='outside window')
        mask=(r>=fitted['r_min'])&(r<=fitted['r_max'])
        ax.errorbar(r[mask],y[mask],yerr=sigma[mask],fmt='k.',markersize=4,capsize=2,label=r'G +/- $\sigma_{diag}$')
        dense=np.linspace(fitted['r_min'],fitted['r_max'],200)
        labels=[]
        for model in MODEL_NAMES:
            fit=analysis['fit_cache'][(L,key[1],model)]
            ax.plot(dense,evaluate(model,dense,fit['log_amplitude_at_reference'],fit['decay_coefficient'],fit['r_reference']),
                    color=MODEL_COLORS[model],label=model.capitalize())
            labels.append(f"{model[:3]}: {fit['reduced_chi_squared_diag']:.2g}")
        ax.set(title=f"T={profile['temperature']:.3f}; beta={key[1]:.2f}",xlabel='r',ylabel='G(r)',yscale='log')
        ax.text(.97,.97,'reduced chi2_diag\n'+'; '.join(labels),transform=ax.transAxes,ha='right',va='top',fontsize=7)
        ax.grid(alpha=.15)
    for ax in list(axes.flat)[len(points):]:ax.set_visible(False)
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='outside lower center', ncol=4, fontsize=8)
    fig.suptitle(f'L = {L}: algebraic and exponential fits at every temperature; fits drawn only in the common window')
    return fig


def plot_window_sensitivity(analysis):
    import matplotlib.pyplot as plt
    sizes=sorted({row['L'] for row in analysis['windows']})
    uppers=sorted({row['upper_fraction'] for row in analysis['windows']})
    fig,axes=plt.subplots(1,len(sizes),figsize=(5*len(sizes),4.3),constrained_layout=True,squeeze=False)
    for ax,L in zip(axes.flat,sizes):
        for index,upper in enumerate(uppers):
            valid=[row for row in analysis['windows'] if row['L']==L and row['upper_fraction']==upper and row['fit_status']=='ok']
            points=sorted({row['beta'] for row in valid},key=lambda beta:1/beta)
            if not points:continue
            deltas=[]
            for beta in points:
                a=next(row for row in valid if row['beta']==beta and row['model']=='algebraic')
                e=next(row for row in valid if row['beta']==beta and row['model']=='exponential')
                deltas.append(e['reduced_chi_squared_diag']-a['reduced_chi_squared_diag'])
            ax.plot([1/b for b in points],deltas,'o-',color=f'C{index}',markersize=4,label=f'r_max/L={upper:.3g}')
        ax.axhline(0,color='0.5',linestyle=':')
        ax.set(title=f'L = {L}',xlabel=r'Temperature $T=1/\beta$',ylabel='Reduced chi2_exp - reduced chi2_alg',yscale='symlog')
        ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.suptitle('Fit-window sensitivity of relative diagonal weighted misfit; insufficient windows omitted')
    return fig
