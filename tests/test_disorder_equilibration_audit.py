import sys,tempfile,json,unittest
from pathlib import Path
sys.path.insert(0,str(Path.cwd()/'scripts'))
import check_disorder_equilibration as audit
import h5py,numpy as np
class AuditTests(unittest.TestCase):
 def test_controls_and_state_preserved(self):
  with tempfile.TemporaryDirectory() as d:
   a,b=Path(d)/'a.h5',Path(d)/'b.h5'
   with h5py.File(a,'w') as f:
    p=f.create_group('parameters')
    for k,v in dict(thermalization=50000,sweeps=1000000000,runtimelimit=18000).items():p[k]=v
    p['checkpoint']='state.clone.h5';p['outputfile']='result.out.h5'
    p.attrs['ini_keys']=['thermalization','sweeps','runtimelimit'];p.attrs['ini_values']=['50000','1000000000','18000']
    f[audit.CP+'/random']='unchanged rng';f[audit.CP+'/configuration/data']=np.arange(15);f[audit.CP+'/sweeps']=80000
   digest=audit.sha(a);audit.patch(a,b,280000,10,3600)
   self.assertEqual(audit.sha(a),digest)
   with h5py.File(a) as src,h5py.File(b) as dst:
    for key in [audit.CP+'/random',audit.CP+'/configuration/data',audit.CP+'/sweeps']:
     np.testing.assert_equal(src[key][()],dst[key][()])
    self.assertEqual(int(dst['parameters/thermalization'][()]),280000)
    self.assertEqual(list(dst['parameters'].attrs['ini_values']),['280000','10','3600'])
   with self.assertRaises(ValueError):audit.patch(a,b,280000,10,3600)
 def test_slot_limit_waits_for_entire_main_plan(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/'plan.json'
   for status in ['running','completed_with_failures','failed','completed']:
    p.write_text(json.dumps(dict(status=status)))
    self.assertEqual(audit.slots(p),8 if status=='completed' else 2)
 def test_production_targets_are_separate(self):
  self.assertEqual(audit.controls('discard',280000,200000),(280000,10))
  self.assertEqual(audit.controls('window_A',280000,200000),(280000,200000))
  self.assertEqual(audit.controls('window_B',280000,200000),(280000,200000))
if __name__=='__main__':unittest.main()
