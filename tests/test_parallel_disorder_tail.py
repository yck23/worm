from pathlib import Path
import hashlib,json,os,sys,tempfile,threading,unittest
from unittest import mock
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import run_parallel_disorder_tail as tail
import run_disorder_queue as queue


class TailTests(unittest.TestCase):
    def setUp(self):
        tmp=tempfile.TemporaryDirectory();self.addCleanup(tmp.cleanup)
        self.root=Path(tmp.name)
        self.addCleanup(mock.patch.stopall)
        mock.patch.object(queue,'REPO',self.root).start()
        self.exe=self.root/'solver.exe';self.exe.write_bytes(b'fixture')
        self.dep=self.root/'dependency/scan.json';self.dep.parent.mkdir()
        self.mp=self.dep.parent/'manifest.json'
        self.save(self.mp,dict(status='running',runs=[
            dict(id='r0',status='completed',returncode=0,output_file='keep.h5',note='unchanged'),
            dict(id='r1',status='running',input_file='start.clone.h5',output_file='adopt.h5')]))
        self.save(self.dep,dict(status='running',controller_pid=123,points=[dict(id='L128_V9',status='running',manifest='manifest.json')]))
        self.scans=[]
        for index,beta in enumerate([.52,.39]):
            path=self.root/f'stage{index}/scan.json';path.parent.mkdir()
            self.save(path,dict(status='prepared',executable=str(self.exe),
                executable_sha256=hashlib.sha256(self.exe.read_bytes()).hexdigest(),
                config=dict(workers=16,realizations=8,beta=beta),points=[dict(manifest='manifest.json')]))
            self.save(path.parent/'manifest.json',dict(status='prepared',workers=8,
                runs=[dict(id=f'r{i}',status='prepared',input_sha256=f'hash{i}') for i in range(8)]))
            self.scans.append(path)
        self.old=self.root/'old/queue.json';self.old.parent.mkdir()
        self.save(self.old,dict(status='recovering_dependency',controller_pid=123,
            dependency=queue.scan_reference(self.dep),stages=[dict(scan=queue.scan_reference(p)) for p in self.scans]))
        self.plan=self.root/'new/plan.json'

    def save(self,path,state):path.write_text(json.dumps(state))

    def prepared(self):
        handle=mock.Mock(creation=12345);handle.poll.return_value=None
        with mock.patch.object(tail,'ProcessHandle',return_value=handle),mock.patch.object(tail,'process_alive',return_value=True):
            tail.prepare(self.plan,self.old,456,self.scans)
        return queue.read(self.plan)

    def test_worker_limit_preserves_inputs_and_all_realizations(self):
        state=self.prepared()
        self.assertEqual(state['max_simultaneous_qmc'],15)
        for path in self.scans:
            m=queue.read(path.parent/'manifest.json')
            self.assertEqual(m['workers'],7)
            self.assertEqual(len(m['runs']),8)
            self.assertEqual([r['input_sha256'] for r in m['runs']],[f'hash{i}' for i in range(8)])
            self.assertTrue((path.parent/'manifest.before_parallel_20261102.json').exists())

    def test_refuses_previously_started_work_before_mutating_any_manifest(self):
        m=queue.read(self.scans[1].parent/'manifest.json');m['status']='running'
        self.save(self.scans[1].parent/'manifest.json',m)
        with mock.patch.object(tail,'ProcessHandle') as ph,mock.patch.object(tail,'process_alive',return_value=True):
            ph.return_value.poll.return_value=None
            with self.assertRaisesRegex(ValueError,'already started'):
                tail.prepare(self.plan,self.old,456,self.scans)
        self.assertEqual(queue.read(self.scans[0].parent/'manifest.json')['workers'],8)
        self.assertFalse(self.plan.exists())

    def test_adopted_exit_preserves_completed_entry(self):
        state=self.prepared();before=queue.read(self.mp)['runs'][0]
        (self.mp.parent/'adopt.h5').write_bytes(b'completed output')
        with mock.patch.object(tail,'point_rows') as validate,mock.patch.object(tail,'analyze'):
            self.assertEqual(tail.finalize_adopted(state,0),0)
            validate.assert_called_once()
        self.assertEqual(queue.read(self.mp)['runs'][0],before)
        self.assertEqual(queue.read(self.mp)['runs'][1]['returncode'],0)
        self.assertEqual(queue.read(self.dep)['status'],'completed')

    def test_sweeps_start_before_adopted_process_exits(self):
        state=self.prepared()
        old=queue.read(self.old);old['status']='superseded';self.save(self.old,old)
        self.plan.with_name('handoff.commit').write_text('go')
        running=threading.Event()
        handle=mock.Mock()
        handle.poll.side_effect=lambda: 0 if running.is_set() else None
        handle.wait_tick.side_effect=lambda: running.wait(.05)
        def sweeps(plan,save):
            running.set();return 0
        with mock.patch.object(tail,'ProcessHandle',return_value=handle), \
             mock.patch.object(tail,'process_alive',return_value=False), \
             mock.patch.object(tail,'run_sweeps',side_effect=sweeps), \
             mock.patch.object(tail,'finalize_adopted',return_value=0):
            self.assertEqual(tail.run(self.plan),0)
        self.assertTrue(running.is_set())
        self.assertEqual(queue.read(self.plan)['status'],'completed')

    @unittest.skipUnless(os.name=='nt','Windows process handle')
    def test_real_handle_checks_creation_and_keeps_current_process_alive(self):
        handle=tail.ProcessHandle(os.getpid())
        try:
            self.assertIsNone(handle.poll())
            self.assertGreater(handle.creation,0)
            handle.wait_tick(1)
            with self.assertRaisesRegex(ValueError,'creation time'):
                tail.ProcessHandle(os.getpid(),handle.creation+1)
        finally:handle.close()


if __name__=='__main__':unittest.main()
