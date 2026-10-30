"""Checkpoint recovery must keep completed chains and their data untouched."""
from pathlib import Path
import hashlib
import json
import sys
import tempfile
import unittest
from unittest import mock
import h5py
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import resume_disorder_response as resume


class ContinuationTests(unittest.TestCase):
    def setUp(self):
        tmp=tempfile.TemporaryDirectory();self.addCleanup(tmp.cleanup)
        self.root=Path(tmp.name)
        exe=self.root/'solver.exe';exe.write_bytes(b'test executable')
        self.mp=self.root/'point/manifest.json';self.mp.parent.mkdir()
        runs=[]
        for index in range(2):
            case=self.mp.parent/f'r{index:04}';case.mkdir()
            output=case/'result.out.h5'
            with h5py.File(output,'w') as f:
                if index==0:
                    for name in ['Equal_time_susceptibility','Winding_number_squared',
                                 'Winding_number_squared_total','Number_of_particles']:
                        f[f'simulation/results/{name}/count']=100
            checkpoint=case/'state.clone.h5'
            with h5py.File(checkpoint,'w') as f:
                f['parameters/checkpoint']=np.bytes_('state.clone.h5')
                f['parameters/outputfile']=np.bytes_('result.out.h5')
                for name,value in [('beta',.45),('canonical',5),('Lx',4),
                                   ('runtimelimit',60),('thermalization',50000)]:
                    f['parameters/'+name]=value
                f['parameters'].attrs['ini_keys']=np.array(['runtimelimit'],dtype=h5py.string_dtype())
                f['parameters'].attrs['ini_values']=np.array(['60'],dtype=h5py.string_dtype())
                base='simulation/realizations/0/clones/0/checkpoint/'
                f[base+'sweeps']=21054
                f[base+'random']=np.bytes_('retained RNG state')
                f[base+'hamiltonian/onsite_potential']=np.zeros(16)
                f[base+'configuration/state']=np.arange(16)
            relative=f'r{index:04}'
            runs.append(dict(id=relative,directory=relative,status='completed',returncode=0,
                checkpoint_file=relative+'/state.clone.h5',output_file=relative+'/result.out.h5',
                input_file=relative+'/parameters.ini',log_file=relative+'/run.log',
                potential_sha256=hashlib.sha256(np.zeros(16,dtype='<f8').tobytes()).hexdigest()))
        self.save(self.mp,dict(status='completed',runs=runs))
        self.scan=self.root/'scan.json'
        self.save(self.scan,dict(status='failed',executable=str(exe),executable_sha256=resume.sha(exe),
            config=dict(beta=.45,particles={'4':5},refresh_paired_notebook=True),
            points=[dict(id='L4_V8',L=4,strength=8,status='failed',manifest='point/manifest.json')]))

    def save(self,path,state):path.write_text(json.dumps(state))
    def read(self,path):return json.loads(path.read_text())

    def test_only_missing_chain_continues_and_retains_checkpoint(self):
        before=self.read(self.mp)
        existing_output=self.mp.parent/before['runs'][0]['output_file']
        original_hash=resume.sha(existing_output)
        resume.prepare(self.scan,43200,2,only_missing_measurements=True)
        prepared=self.read(self.mp)
        self.assertEqual(prepared['runs'][0],before['runs'][0])
        self.assertEqual(prepared['continuation']['runs'],['r0001'])
        job=prepared['runs'][1]
        self.assertEqual(job['status'],'prepared')
        start=self.mp.parent/job['input_file']
        source=self.mp.parent/before['runs'][1]['checkpoint_file']
        with h5py.File(source) as a,h5py.File(start) as b:
            self.assertEqual(int(b['parameters/runtimelimit'][()]),43200)
            self.assertEqual(list(b['parameters'].attrs['ini_values']),['43200'])
            base='simulation/realizations/0/clones/0/checkpoint/'
            for key in ['sweeps','random','configuration/state','hamiltonian/onsite_potential']:
                np.testing.assert_array_equal(a[base+key][()],b[base+key][()])
        with mock.patch.object(resume,'_run_realization',return_value=dict(status='completed',returncode=0)) as execute, \
             mock.patch.object(resume,'point_rows'),mock.patch.object(resume,'refresh'):
            self.assertEqual(resume.run(self.scan),0)
        self.assertEqual(execute.call_count,1)
        self.assertEqual(self.read(self.mp)['runs'][0],before['runs'][0])
        self.assertEqual(resume.sha(existing_output),original_hash)
        self.assertTrue((self.root/'scan.before_continuation.json').exists())
        self.assertTrue(self.mp.with_name('manifest.before_continuation.json').exists())

    def test_refuses_running_scan(self):
        state=self.read(self.scan);state['status']='running';self.save(self.scan,state)
        with self.assertRaisesRegex(ValueError,'running scan'):
            resume.prepare(self.scan,43200,2,True)
        self.assertFalse((self.root/'scan.before_continuation.json').exists())

    def test_refuses_repeated_preparation(self):
        resume.prepare(self.scan,43200,2,True)
        with self.assertRaisesRegex(ValueError,'already prepared'):
            resume.prepare(self.scan,43200,2,True)


if __name__=='__main__':unittest.main()
