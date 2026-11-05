#!/usr/bin/env python3
"""Refresh existing comparison notebooks as equilibration stages finish."""
import json,os,subprocess,sys,time
from pathlib import Path
REPO=Path(__file__).resolve().parents[1]
ROOT=REPO/'results/disorder_beta060_20261027/equilibration/L128_V9_20261102'
def signature(a):
    return a['status'],tuple((r['id'],r['status'],tuple((s['name'],s['status'],len(s['chunks']),s['chunks'][-1].get('saved_sweeps',s['chunks'][-1]['source_sweeps']) if s['chunks'] else 0) for s in r['stages'])) for r in a['runs'])
def main():
    (ROOT/'notebook_watcher_pid.txt').write_text(str(os.getpid()))
    last=None
    while True:
        try:
            a=json.loads((ROOT/'audit.json').read_text());current=signature(a)
        except (OSError,json.JSONDecodeError) as exc:
            print('Transient audit read error; retrying:',exc,flush=True);time.sleep(30);continue
        if current!=last:
            print('Refreshing notebook:',current,flush=True)
            with (ROOT/'notebook_refresh.log').open('a',encoding='utf-8') as log:
                result=subprocess.run([sys.executable,str(REPO/'scripts/analyze_disorder_response.py')],cwd=REPO,stdout=log,stderr=subprocess.STDOUT,
                    timeout=900,creationflags=subprocess.CREATE_NO_WINDOW if os.name=='nt' else 0)
            if result.returncode==0:last=current
            else:print('Exporter returned',result.returncode,'will retry',flush=True)
        if a['status'] in ['completed','failed'] and last==current:return 0
        time.sleep(30)
if __name__=='__main__':raise SystemExit(main())
