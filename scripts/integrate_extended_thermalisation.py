#!/usr/bin/env python3
"""Integrate all eight extended-warmup maps into the existing dataset once ready."""
import copy,json,os,shutil,subprocess,sys,time
from pathlib import Path
from run_disorder_ensemble import _write_json
from summarize_disorder_ensemble import summarize_manifest
from check_disorder_equilibration import sha,now
REPO=Path(__file__).resolve().parents[1]
ROOT=REPO/'results/disorder_beta060_20261027/equilibration/L128_V9_20261102'
def integrate(root=ROOT):
    a=json.loads((root/'audit.json').read_text())
    if a['status']!='completed' or len(a['runs'])!=8 or any(r['status']!='completed' for r in a['runs']):
        raise ValueError('Wait for all eight maps; no partial replacement')
    mp=Path(a['source_manifest']);original=json.loads(mp.read_text())
    backup=mp.with_name('manifest.before_extended_thermalisation_20261102.json')
    if original.get('extended_thermalisation',{}).get('audit')==str(root/'audit.json'):return mp
    assert sha(mp)==a['source_manifest_sha256'],'Source manifest changed; inspect before integration'
    new=copy.deepcopy(original)
    for run in new['runs']:
        entry=next(e for e in a['runs'] if e['id']==run['id'])
        stage=next(s for s in entry['stages'] if s['name']=='window_B')
        assert stage['status']=='completed'
        chunk=stage['chunks'][-1];src=(root/chunk['checkpoint']).parent
        relative=Path(run['id'])/'extended_thermalisation_20261102';dest=mp.parent/relative
        dest.mkdir(exist_ok=True)
        for name in ['start.clone.h5','state.clone.h5','result.out.h5','run.log']:
            target=dest/name
            if target.exists():assert sha(target)==sha(src/name)
            else:shutil.copy2(src/name,target)
        run.update(directory=relative.as_posix(),input_file=(relative/'start.clone.h5').as_posix(),
            input_sha256=sha(dest/'start.clone.h5'),checkpoint_file=(relative/'state.clone.h5').as_posix(),
            output_file=(relative/'result.out.h5').as_posix(),log_file=(relative/'run.log').as_posix(),
            status='completed',returncode=0,finished_utc=chunk['finished_utc'],
            extended_thermalisation=dict(audit=str(root/'audit.json'),discarded_excursions=a['discard'],
                production_window='B',production_excursions=a['window'],original_output=run['output_file']))
    new['extended_thermalisation']=dict(audit=str(root/'audit.json'),integrated_utc=now(),
        statistics='Fresh window B for all eight maps; original and window A retained separately',
        equilibration='unverified; completing additional warmup does not establish equilibration')
    candidate=mp.with_name('manifest.extended_thermalisation_candidate.json');_write_json(candidate,new)
    summarize_manifest(candidate)  # Full physical/map/completeness validation before changing the live manifest.
    if not backup.exists():_write_json(backup,original)
    _write_json(mp,new)
    _write_json(root/'integration.json',dict(status='integrated',finished_utc=now(),manifest=str(mp),backup=str(backup),provisional=True))
    return mp

def main():
    _write_json(ROOT/'integration.json',dict(status='waiting_for_all_eight_maps',controller_pid=os.getpid(),started_utc=now()))
    while True:
        try:a=json.loads((ROOT/'audit.json').read_text())
        except (OSError,json.JSONDecodeError) as exc:
            print('Transient audit read error; retrying:',exc,flush=True);time.sleep(30);continue
        if a['status']=='failed':raise RuntimeError('Thermalisation audit failed; original notebook dataset retained')
        if a['status']=='completed':break
        time.sleep(30)
    mp=integrate()
    scan=mp.parents[2]/'scan.json'
    for command in [[sys.executable,str(REPO/'scripts/disorder_response_scan.py'),'analyze','--scan',str(scan)],
                    [sys.executable,str(REPO/'scripts/analyze_disorder_response.py')]]:
        subprocess.run(command,cwd=REPO,check=True,creationflags=subprocess.CREATE_NO_WINDOW if os.name=='nt' else 0)
    print('Integrated all eight fresh windows into existing dataset and refreshed notebook.',flush=True)
    return 0
if __name__=='__main__':raise SystemExit(main())
