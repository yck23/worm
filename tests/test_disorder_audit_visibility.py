import sys,tempfile,json,unittest
from pathlib import Path
sys.path.insert(0,str(Path.cwd()/'scripts'))
from disorder_equilibration_notebook import progress_markdown
from watch_disorder_equilibration_notebook import signature
class ProgressTests(unittest.TestCase):
 def test_reset_counter_not_mistaken_for_production(self):
  with tempfile.TemporaryDirectory() as d:
   repo=Path(d);folder=repo/'results/disorder_beta060_20261027/equilibration/L128_V9_20261102';folder.mkdir(parents=True)
   a=dict(status='running',discard=200000,window=200000,runs=[dict(id='r0000',status='running',source_sweeps=80000,stages=[
    dict(name='discard',status='completed',chunks=[]),dict(name='window_A',status='running',chunks=[dict(reset_statistics=True,status='running',source_sweeps=280010)])])])
   (folder/'audit.json').write_text(json.dumps(a));text=progress_markdown(repo)
   self.assertIn('| r0000 | 200,000 | 0 (running) | queued |',text)
   c=a['runs'][0]['stages'][1]['chunks'][0];c.update(status='completed',saved_sweeps=300000)
   (folder/'audit.json').write_text(json.dumps(a));self.assertIn('20,000 (running)',progress_markdown(repo))
 def test_signature_changes_on_saved_progress(self):
  a=dict(status='running',runs=[dict(id='r0000',status='running',stages=[dict(name='discard',status='running',chunks=[dict(source_sweeps=80000)])])])
  old=signature(a);a['runs'][0]['stages'][0]['chunks'][0]['saved_sweeps']=120000
  self.assertNotEqual(old,signature(a))
if __name__=='__main__':unittest.main()
