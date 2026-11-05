#!/usr/bin/env python3
"""Prepare, run unattended, inspect, and plot a fixed-beta disorder scan."""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import contextlib
import csv
from datetime import datetime, timezone
import hashlib
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import threading

import h5py
import numpy as np
from run_disorder_ensemble import prepare_ensemble, run_ensemble, _write_json
from summarize_disorder_ensemble import summarize_manifest

REPO = Path(__file__).resolve().parents[1]
CHI = 'Equal_time_susceptibility'
WTOTAL = 'Winding_number_squared_total'
LOCK = threading.Lock()


def now():
    return datetime.now(timezone.utc).isoformat()


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def prepare(output_dir, config_path, executable=None):
    config = json.loads(Path(config_path).read_text(encoding='utf-8-sig'))
    root = Path(output_dir).resolve()
    if root.exists() and any(root.iterdir()):
        raise ValueError(f'Refusing populated scan directory: {root}')
    executable = Path(executable or REPO/'dist/windows-square-programmable/qmc_worm_canonical.exe').resolve(strict=True)
    if '_mpi' in executable.stem.lower():
        raise ValueError('Use a serial canonical executable')
    realizations, workers = int(config['realizations']), int(config['workers'])
    if realizations < 2 or workers < realizations or workers % realizations:
        raise ValueError('workers must be a positive multiple of realizations, and realizations >= 2')
    data_subdirectory = Path(config.get('data_subdirectory', ''))
    if data_subdirectory.is_absolute() or '..' in data_subdirectory.parts:
        raise ValueError('data_subdirectory must stay inside the scan output directory')
    root.mkdir(parents=True, exist_ok=True)
    data_root = root / data_subdirectory
    data_root.mkdir(parents=True, exist_ok=True)
    (root/'inputs').mkdir()
    state = dict(schema_version=1, kind='disorder_response_scan', status='preparing',
                 created_utc=now(), config=config, executable=str(executable),
                 executable_sha256=digest(executable), points=[])
    state['branch'] = subprocess.check_output(['git','branch','--show-current'], cwd=REPO,text=True).strip()
    state['head'] = subprocess.check_output(['git','rev-parse','HEAD'],cwd=REPO,text=True).strip()
    if state['branch'] != 'assess/site-programmable-potentials':
        raise ValueError('Run this test on assess/site-programmable-potentials')
    provenance=root/'provenance'
    provenance.mkdir()
    for name in ['scripts/disorder_response_scan.py','scripts/run_disorder_ensemble.py',
                 'scripts/generate_speckle.py','scripts/summarize_disorder_ensemble.py',
                 'scripts/analyze_disorder_response.py',
                 'src/worm.cpp','src/worm.hpp','src/worm.update.cpp','src/worm.output.cpp']:
        target=provenance/name
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(REPO/name,target)
    (provenance/'working_changes.diff').write_bytes(subprocess.check_output(['git','diff'],cwd=REPO))
    (provenance/'reference_notebook.ipynb').write_bytes(subprocess.check_output(
        ['git','show',config['reference']['git_ref']+':'+config['reference']['notebook']],cwd=REPO))
    for name in ['parameters','particle_numbers']:
        reference=config['reference'][name]
        (provenance/Path(reference).name).write_bytes(subprocess.check_output(
            ['git','show',config['reference']['git_ref']+':'+reference],cwd=REPO))
    state['reference_commit']=subprocess.check_output(
        ['git','rev-parse',config['reference']['git_ref']],cwd=REPO,text=True).strip()
    _write_json(root/'config.json',config)
    _write_json(root/'scan.json',state)
    for strength_index,strength in enumerate(config['strengths']):
        for L in config['sizes']:
            label=f'L{L}_V{strength:g}'.replace('.','p')
            params=dict(config['parameters'], Lx=str(L), Ly=str(L),
                        canonical=str(config['particles'][str(L)]), beta=str(config['beta']),
                        thermalization=str(config['thermalization']), sweeps=str(config['sweeps']),
                        runtimelimit=str(config['runtime_seconds'][str(L)]))
            ini=root/'inputs'/f'{label}.ini'
            ini.write_text('\n'.join(f'{key} = {value}' for key,value in params.items())+'\n',encoding='utf-8')
            manifest=prepare_ensemble(ini,data_root/label,realizations=realizations,workers=realizations,
                correlation_length=config['correlation_length'],mean_strength=strength,
                base_disorder_seed=config['base_disorder_seed']+L*1000,
                base_mc_seed=config['base_mc_seed']+L*1000+strength_index*realizations,
                canonical=True, executable=executable,prepare_only=True)
            state['points'].append(dict(id=label,L=L,strength=strength,
                manifest=str(manifest.relative_to(root)), status='prepared'))
            _write_json(root/'scan.json',state)
    state['status']='prepared'
    _write_json(root/'scan.json',state)
    make_notebook(root, [], [])
    print(f'Prepared {len(state["points"])} points / {len(state["points"])*realizations} simulations: {root}')
    return root/'scan.json'


def measurement(group):
    if int(group['count'][()]) <= 0:
        raise ValueError('No production measurements')
    value=np.asarray(group['mean/value'][()],dtype=float)
    if not np.all(np.isfinite(value)):
        raise ValueError('Nonfinite production mean')
    error=np.asarray(group['mean/error'][()],dtype=float) if 'mean/error' in group else np.full_like(value,np.nan)
    return value,error,int(group['count'][()])


def point_rows(root,point):
    manifest_path=root/point['manifest']
    # Includes completeness, potential hash, common physical-parameter checks.
    summary=summarize_manifest(manifest_path)
    manifest=json.loads(manifest_path.read_text())
    rows=[]
    for run in manifest['runs']:
        file=manifest_path.parent/run['output_file']
        with h5py.File(file,'r') as data:
            groups=data['simulation/results']
            chi,chi_error,chi_count=measurement(groups[CHI])
            w,w_error,w_count=measurement(groups['Winding_number_squared'])
            total,total_error,total_count=measurement(groups[WTOTAL])
            N=int(data['parameters/canonical'][()])
            measured_N=float(groups['Number_of_particles/mean/value'][()])
            if not np.isclose(measured_N,N,rtol=0,atol=1e-9) or w.shape != (2,):
                raise ValueError('Unexpected canonical number or winding shape')
            if not np.isclose(w.sum(),total,rtol=1e-10,atol=1e-10) or w_count != total_count:
                raise ValueError('Total winding must match component sum and sample counts')
            rows.append(dict(L=point['L'],V=point['strength'],beta=float(data['parameters/beta'][()]),
                N=N,realization=run['id'],disorder_seed=run['disorder_seed'],mc_seed=run['mc_seed'],
                chi_s=float(chi),chi_s_mc_se=float(chi_error),chi_count=chi_count,
                wx2=float(w[0]),wy2=float(w[1]),wtotal=float(total),
                wx2_mc_se=float(w_error[0]),wy2_mc_se=float(w_error[1]),
                wtotal_mc_se=float(total_error),winding_count=w_count,
                source_file=str(file.relative_to(root))))
    _write_json(manifest_path.parent/'summary.json',summary)
    mean=dict(L=point['L'],V=point['strength'],beta=rows[0]['beta'],N=rows[0]['N'],realizations=len(rows))
    for key in ['chi_s','wx2','wy2','wtotal']:
        values=np.array([r[key] for r in rows])
        mean[key]=float(values.mean())
        mean[key+'_se']=float(values.std(ddof=1)/np.sqrt(len(values)))
    mean['chi_scaled']=mean['chi_s']/point['L']**1.75
    mean['chi_scaled_se']=mean['chi_s_se']/point['L']**1.75
    return mean,rows


def write_csv(path,rows,fields):
    temporary=path.with_suffix('.tmp')
    with temporary.open('w',newline='',encoding='utf-8') as stream:
        writer=csv.DictWriter(stream,fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def _analyze_local(scan):
    scan=Path(scan).resolve()
    root=scan.parent
    state=json.loads(scan.read_text())
    means,records=[],[]
    for point in state['points']:
        if point['status']=='completed':
            mean,rows=point_rows(root,point)
            means.append(mean)
            records.extend(rows)
    out=root/'analysis'
    out.mkdir(exist_ok=True)
    fields=['L','V','beta','N','realizations','chi_s','chi_s_se','wx2','wx2_se','wy2','wy2_se',
            'wtotal','wtotal_se','chi_scaled','chi_scaled_se']
    write_csv(out/'summary.csv',means,fields)
    write_csv(out/'realizations.csv',records,list(records[0]) if records else
              ['L','V','beta','N','realization','chi_s','wx2','wy2','wtotal'])
    make_notebook(root,means,records)
    _write_json(out/'analysis_status.json',dict(updated_utc=now(),completed_points=len(means),
                requested_points=len(state['points']),scan_status=state['status'],
                weighting='equal realization means; SE = sample SD / sqrt(R); MC noise already included',
                incomplete_points=[p['id'] for p in state['points'] if p['status']!='completed']))
    print(f'Analysis refreshed: {len(means)}/{len(state["points"])} points',flush=True)


def analyze(scan):
    # The standalone beta dataset may also be a root in the paired exporter.
    # Serialize its local table/notebook writes with other exporters, then
    # release before spawning the paired exporter (which takes the same lock).
    from disorder_analysis_lock import AnalysisLock
    scan=Path(scan).resolve()
    lock=AnalysisLock(REPO/'results/.disorder_analysis.lock').acquire()
    try:
        _analyze_local(scan)
    finally:
        lock.release()
    root=scan.parent
    out=root/'analysis'
    out.mkdir(exist_ok=True)
    state=json.loads(scan.read_text())
    if state['config'].get('refresh_paired_notebook', False):
        try:
            with (root/'paired_notebook_refresh.log').open('a',encoding='utf-8') as log:
                refreshed = subprocess.run([sys.executable,str(REPO/'scripts/analyze_disorder_response.py')],
                    cwd=REPO,stdout=log,stderr=subprocess.STDOUT,timeout=600,
                    creationflags=subprocess.CREATE_NO_WINDOW if os.name=='nt' else 0)
            if refreshed.returncode:
                raise RuntimeError('Paired notebook refresh failed; see paired_notebook_refresh.log')
        except Exception as exc:
            _write_json(out/'paired_refresh_error.json',dict(error=str(exc),updated_utc=now()))
            print('Paired notebook refresh issue:',exc,flush=True)


def make_notebook(root,means,records):
    beta_label = f"{float(json.loads((root/'config.json').read_text())['beta']):g}"
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import base64
    cells=[]
    def markdown(text):
        cells.append(dict(cell_type='markdown',metadata={},source=text.splitlines(keepends=True)))
    def code(source):
        cells.append(dict(cell_type='code',metadata={},source=source.splitlines(keepends=True),
                          execution_count=None,outputs=[]))
    markdown(f'# Susceptibility and winding versus disorder, beta = {beta_label}\n\n'
        'L = 32, 64, 128; canonical N = 1875, 7500, 30000. '
        'V is the mean of the nonnegative speckle potential in hopping units; correlation length = 2 sites.\n\n'
        'Matches origin/main `results/cleanRbRDM/density_matrix_comparison.ipynb`: '
        '`chi_s = G.sum() - G[0,0] + N/L**2`, with no beta factor or extra site-count division. '
        'The new scalar accumulator sums off-site events before binning and uses the exact density. '
        'Total winding is also accumulated before estimating its thermal error.\n\n'
        'Each faint dot is an independently equilibrated disorder realization. Lines weight maps equally; '
        'error bars are sample SD of the whole-realization means divided by sqrt(R), including MC noise. '
        'Do not add marginal MC errors to those bars. Different strengths reuse the same unit-strength '
        'maps within each size; errors across strengths are correlated. '
        'Pending or failed points are omitted and listed in `analysis/analysis_status.json`. '
        'These tests retain the reference nmax=5; convergence is not established by successful execution.')
    code("from pathlib import Path\nimport csv\nimport numpy as np\nimport matplotlib.pyplot as plt\n\n"
        "root = Path.cwd()\nfor folder in (root, *root.parents):\n"
        "    if (folder / 'scan.json').is_file():\n        root = folder\n        break\n"
        "# If launched elsewhere, set root to the scan directory here.\n"
        "with (root / 'analysis/summary.csv').open(newline='') as stream:\n"
        "    rows = [{k: float(v) for k, v in row.items()} for row in csv.DictReader(stream)]\n"
        "with (root / 'analysis/realizations.csv').open(newline='') as stream:\n"
        "    runs = list(csv.DictReader(stream))\n"
        "import json\nsizes = json.loads((root / 'config.json').read_text())['sizes']\ncolours = ['tab:green', 'tab:blue', 'tab:orange']\n"
        "print(f'{len(rows)} completed points; {len(runs)} realization records')\n")
    plots=[('Equal-time susceptibility','chi_s',r'$\chi_s$','susceptibility_vs_disorder'),
           ('Scaled susceptibility','chi_scaled',r'$\chi_s/L^{7/4}$','susceptibility_scaled_vs_disorder'),
           ('Squared x winding','wx2',r'$\langle W_x^2\rangle$','winding_x_vs_disorder'),
           ('Squared y winding','wy2',r'$\langle W_y^2\rangle$','winding_y_vs_disorder'),
           ('Total squared winding','wtotal',r'$\langle W_x^2+W_y^2\rangle$','winding_total_vs_disorder')]
    for title,key,ylabel,filename in plots:
        markdown('## '+title)
        code(f"plt.figure(figsize=(7,4.5))\nfor L, colour in zip(sizes, colours):\n"
            f"    selected = sorted((r for r in rows if r['L'] == L), key=lambda r:r['V'])\n"
            f"    individual = [r for r in runs if int(r['L']) == L]\n"
            f"    values = [float(r['chi_s'])/L**1.75 if '{key}' == 'chi_scaled' else float(r['{key if key!='chi_scaled' else 'chi_s'}']) for r in individual]\n"
            f"    plt.scatter([float(r['V']) for r in individual], values, color=colour, alpha=.2, s=12)\n"
            f"    plt.errorbar([r['V'] for r in selected], [r['{key}'] for r in selected],\n"
            f"                 yerr=[r['{key}_se'] for r in selected], fmt='o-', capsize=4, color=colour, label=f'L = {{L}}')\n"
            f"plt.xlabel('Mean disorder strength V / t')\nplt.ylabel({ylabel!r})\n"
            f"plt.title('{title}; beta = {beta_label}')\nplt.legend()\nplt.grid(alpha=.2)\nplt.tight_layout()\n"
            f"plt.savefig(root / 'analysis/{filename}.png', dpi=180)\nplt.show()\n")
    markdown('## Susceptibility versus size at each disorder strength')
    code("plt.figure(figsize=(7,4.5))\nfor V in sorted({r['V'] for r in rows}):\n"
         "    selected = sorted((r for r in rows if r['V'] == V), key=lambda r:r['L'])\n"
         "    plt.errorbar([r['L'] for r in selected], [r['chi_s'] for r in selected],\n"
         "                 yerr=[r['chi_s_se'] for r in selected], fmt='o-', capsize=3, label=f'V/t = {V:g}')\n"
         "plt.xscale('log')\nplt.yscale('log')\nplt.xticks(sizes, sizes)\nplt.minorticks_off()\n"
         f"plt.xlabel('Linear size L')\nplt.ylabel(r'$\\chi_s$')\nplt.title('Susceptibility versus size; beta = {beta_label}')\n"
         "plt.legend(ncol=2, fontsize=8)\nplt.grid(alpha=.2)\nplt.tight_layout()\n"
         "plt.savefig(root / 'analysis/susceptibility_vs_size.png', dpi=180)\nplt.show()\n")
    markdown('## Susceptibility uncertainties in percent')
    code("fig, axes = plt.subplots(1,3,figsize=(12,3.6),constrained_layout=True)\n"
         "for L, colour, ax in zip(sizes,colours,axes):\n"
         "    selected = sorted((r for r in rows if r['L']==L),key=lambda r:r['V'])\n"
         "    x = np.array([r['V'] for r in selected])\n"
         "    relative = np.array([100*r['chi_s_se']/abs(r['chi_s']) for r in selected])\n"
         "    ax.errorbar(x, np.zeros(len(x)), yerr=relative, fmt='o', capsize=4, color=colour)\n"
         "    ax.axhline(0,color='.5',lw=.7)\n    ax.set_title(f'L = {L}')\n"
         "    ax.set_xlabel('Mean disorder strength V / t')\n    ax.set_ylabel('Relative one-SE error (%)')\n"
         "fig.suptitle('Whole-realization susceptibility uncertainty')\n"
         "plt.savefig(root/'analysis/susceptibility_relative_errors.png',dpi=180)\nplt.show()\n")
    notebook=dict(cells=cells,metadata=dict(kernelspec=dict(display_name='Python 3',language='python',name='python3'),
                  language_info=dict(name='python',version=sys.version.split()[0])),nbformat=4,nbformat_minor=4)
    if means:
        namespace={'__name__':'__main__','root':root}
        for count,cell in enumerate((c for c in cells if c['cell_type']=='code'),1):
            source=''.join(cell['source'])
            # Preserve the ordinary plotting source; only root discovery is resolved for this execution.
            if count==1:
                source=source.replace('root = Path.cwd()',f'root = Path({str(root)!r})')
            outputs=[]
            def show():
                for number in plt.get_fignums():
                    figure=plt.figure(number)
                    stream=io.BytesIO()
                    figure.savefig(stream,format='png',dpi=110,bbox_inches='tight')
                    outputs.append(dict(output_type='display_data',metadata={},
                        data={'image/png':base64.b64encode(stream.getvalue()).decode()}))
                    plt.close(figure)
            original_show=plt.show
            plt.show=show
            stdout=io.StringIO()
            try:
                with contextlib.redirect_stdout(stdout):
                    exec(compile(source,'disorder_response_comparison.ipynb','exec'),namespace)
            finally:
                plt.show=original_show
                plt.close('all')
            if stdout.getvalue():
                outputs.insert(0,dict(output_type='stream',name='stdout',text=stdout.getvalue().splitlines(keepends=True)))
            cell.update(execution_count=count,outputs=outputs)
    temporary=root/'disorder_response_comparison.ipynb.tmp'
    temporary.write_text(json.dumps(notebook,indent=1)+'\n',encoding='utf-8')
    temporary.replace(root/'disorder_response_comparison.ipynb')


def run(scan):
    scan=Path(scan).resolve(strict=True)
    state=json.loads(scan.read_text())
    if state['status']!='prepared':
        raise ValueError('Only a newly prepared scan can be launched; existing runs are preserved')
    if digest(state['executable'])!=state['executable_sha256']:
        raise ValueError('Executable changed after preparation; prepare a fresh scan')
    root=scan.parent
    state.update(status='running',controller_pid=os.getpid(),started_utc=now())
    _write_json(scan,state)
    def worker(point):
        with LOCK:
            point.update(status='running',started_utc=now())
            _write_json(scan,state)
        code=run_ensemble(root/point['manifest'])
        if code:
            raise RuntimeError(f'Ensemble failed: {point["id"]}; inspect its manifest and run logs')
        point_rows(root,point)
        return point
    slots=state['config']['workers']//state['config']['realizations']
    with ThreadPoolExecutor(max_workers=slots) as pool:
        active={pool.submit(worker,point):point for point in state['points']}
        for future in as_completed(active):
            point=active[future]
            with LOCK:
                try:
                    future.result()
                    point.update(status='completed',finished_utc=now())
                except Exception as exc:
                    point.update(status='failed',finished_utc=now(),error=str(exc))
                _write_json(scan,state)
            analyze(scan)
    state.update(status='completed' if all(p['status']=='completed' for p in state['points']) else 'failed',finished_utc=now())
    _write_json(scan,state)
    analyze(scan)
    return 0 if state['status']=='completed' else 1


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    commands=parser.add_subparsers(dest='command',required=True)
    setup=commands.add_parser('prepare')
    setup.add_argument('--output-dir',type=Path,required=True)
    setup.add_argument('--config',type=Path,default=REPO/'parameter_files/disorder_response_beta035.json')
    setup.add_argument('--executable',type=Path)
    for name in ['run','status','analyze']:
        commands.add_parser(name).add_argument('--scan',type=Path,required=True)
    args=parser.parse_args()
    if args.command=='prepare':
        prepare(args.output_dir,args.config,args.executable)
    elif args.command=='run':
        return run(args.scan)
    elif args.command=='analyze':
        analyze(args.scan)
    else:
        state=json.loads(args.scan.read_text())
        print(f"Scan: {state['status']}; controller PID: {state.get('controller_pid','not launched')}")
        for point in state['points']:
            path=args.scan.parent/point['manifest']
            manifest=json.loads(path.read_text())
            counts={name:sum(r['status']==name for r in manifest['runs'])
                    for name in ['prepared','running','completed','failed','skipped']}
            print(f"{point['id']}: {point['status']}; {counts}")
            if 'error' in point:
                print('  '+point['error'])
    return 0


if __name__=='__main__':
    raise SystemExit(main())
