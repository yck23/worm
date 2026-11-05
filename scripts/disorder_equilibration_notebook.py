"""Notebook cells exposing partial same-map equilibration audit separately."""

from pathlib import Path

import base64,io,json

from contextlib import redirect_stdout

import matplotlib.pyplot as plt



CODE=r'''

# Read the live audit snapshot; rerun this cell to see newly saved windows.

import json

import h5py

import numpy as np

import matplotlib.pyplot as plt

from datetime import datetime, timezone

audit_root=repo/'results/disorder_beta060_20261027/equilibration/L128_V9_20261102'

audit=json.loads((audit_root/'audit.json').read_text())

manifest_file=Path(audit['source_manifest'])

original_manifest=manifest_file.with_name('manifest.before_extended_thermalisation_20261102.json')
manifest=json.loads((original_manifest if original_manifest.exists() else manifest_file).read_text())

observables=['Equal_time_susceptibility','Winding_number_squared_total','Total_Energy']

labels=['Equal-time susceptibility','Total squared winding','Total energy']

original={}

for entry in manifest['runs']:

    with h5py.File(manifest_file.parent/entry['output_file']) as f:

        original[entry['id']]={name:dict(mean=float(f['simulation/results/'+name+'/mean/value'][()]),

                                      error=float(f['simulation/results/'+name+'/mean/error'][()])) for name in observables}

print('Audit snapshot UTC:',datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M'))

print('Map    extra warmup    window A    window B')

for entry in audit['runs']:

    stages={s['name']:s for s in entry['stages']}

    states=[stages.get(name,{}).get('status','queued') for name in ['discard','window_A','window_B']]

    discarded=0

    if 'discard' in stages:

        completed=[c for c in stages['discard']['chunks'] if c['status']=='completed']

        if completed:discarded=min(audit['discard'],max(0,completed[-1]['saved_sweeps']-entry['source_sweeps']))

    print(entry['id'],f'{discarded}/{audit["discard"]}',*states)

fig,axes=plt.subplots(3,1,figsize=(10,10),constrained_layout=True)

for ax,name,label in zip(axes,observables,labels):

    for series,offset,colour in [('Original',-.22,'0.55'),('window_A',0,'tab:blue'),('window_B',.22,'tab:orange')]:

        x,y,err=[],[],[]

        for i,entry in enumerate(audit['runs']):

            if series=='Original':value=original[entry['id']][name]

            else:

                stage=next((s for s in entry['stages'] if s['name']==series and s['status']=='completed'),None)

                if stage is None:continue

                value=stage['statistics'][name]

            x.append(i+offset);y.append(value['mean']);err.append(value['error'])

        if x:ax.errorbar(x,y,yerr=err,fmt='o',capsize=3,color=colour,label=series.replace('_',' '))

    ax.set_xticks(range(8),[e['id'] for e in audit['runs']]);ax.set_ylabel(label);ax.grid(alpha=.2);ax.legend()

axes[-1].set_xlabel('Same disorder map; absent fresh windows are still queued or running')

fig.suptitle('PROVISIONAL: beta=0.60, L=128, disorder 9 — original versus extended-warmup windows')

plt.show()

print('Window B minus A (same map; adjacent windows are correlated):')

for entry in audit['runs']:

    stages={s['name']:s for s in entry['stages'] if s['status']=='completed'}

    if all(k in stages for k in ['window_A','window_B']):

        print(entry['id'],{name:round(stages['window_B']['statistics'][name]['mean']-stages['window_A']['statistics'][name]['mean'],6) for name in observables})

print('No partial-map ensemble average replaces the original eight-map point. More evolution may be needed; completion is not proof of equilibration.')

'''



def notebook_cells(repo):

    root=repo/'results/disorder_beta060_20261027/equilibration/L128_V9_20261102'

    if not (root/'audit.json').exists():return []

    a=json.loads((root/'audit.json').read_text())

    completed=sum(e['status']=='completed' for e in a['runs'])

    intro=f'''## Extended thermalisation: beta=0.60, L=128, disorder 9



**PROVISIONAL — {completed}/8 maps have completed both fresh windows in this saved snapshot.**

Each map receives 200,000 additional discarded excursions, then two separate

200,000-excursion measurement windows with statistics reset between them.

These are continuations of the same maps, not new disorder points or independent

realizations. The main sweep panels retain the original eight-map result until

the full audit is assessed. Missing window markers below mean unfinished work.



Error bars below are each chain's Monte Carlo estimates, unlike the disorder

standard errors in the main panels. Adjacent windows are correlated. Inspect

energy as well as susceptibility and winding; stable susceptibility alone does

not establish equilibration. Rerun the cell to read current audit progress.

'''

    outputs=[];capture=io.StringIO();previous=plt.show

    def show():

        fig=plt.gcf();buf=io.BytesIO();fig.savefig(buf,format='png',dpi=115,bbox_inches='tight')

        fig.savefig(root/'window_comparison.png',dpi=150,bbox_inches='tight')

        outputs.append(dict(output_type='display_data',metadata={},data={'image/png':base64.b64encode(buf.getvalue()).decode()}));plt.close(fig)

    try:

        plt.show=show

        with redirect_stdout(capture):exec(compile(CODE,'equilibration notebook cell','exec'),dict(repo=repo,Path=Path))

    finally:plt.show=previous

    outputs.insert(0,dict(output_type='stream',name='stdout',text=capture.getvalue().splitlines(keepends=True)))

    return [dict(cell_type='markdown',metadata={},source=intro.splitlines(keepends=True)),

            dict(cell_type='code',metadata={},source=CODE.splitlines(keepends=True),execution_count=12,outputs=outputs)]



def progress_markdown(repo):
    from datetime import datetime,timezone
    path=repo/'results/disorder_beta060_20261027/equilibration/L128_V9_20261102/audit.json'
    if not path.exists():return ''
    a=json.loads(path.read_text())
    done=sum(r['status']=='completed' for r in a['runs'])
    lines=[f"### Extended warmup and production: beta=0.60, L=128, disorder 9 ? {done}/8 maps finished",
        '',f"Saved counters as of {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}. Running chunks may have advanced further.",
        '', '| Map | Extra discarded / 200,000 | Fresh A / 200,000 | Fresh B / 200,000 |',
        '|---|---:|---:|---:|']
    for r in a['runs']:
        stages={s['name']:s for s in r['stages']};values=[]
        for name in ['discard','window_A','window_B']:
            goal=a['discard'] if name=='discard' else a['window'];s=stages.get(name)
            if s is None:values.append('queued');continue
            if s['status']=='completed':values.append(f'{goal:,}');continue
            chunk=s['chunks'][-1] if s['chunks'] else None
            if chunk is None:count=0
            elif name!='discard' and chunk['reset_statistics'] and chunk['status']=='running':count=0
            else:
                base=r['source_sweeps']+(0 if name=='discard' else a['discard'])
                count=min(goal,max(0,chunk.get('saved_sweeps',chunk['source_sweeps'])-base))
            values.append(f'{count:,} ({s["status"]})')
        lines.append('| '+r['id']+' | '+' | '.join(values)+' |')
    lines+=['','The main disorder-9 point is replaced only after all eight maps finish. Original files are preserved; equilibration remains provisional.','']
    return '\n'.join(lines)
