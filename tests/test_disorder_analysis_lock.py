from pathlib import Path
import subprocess,sys,tempfile,unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from disorder_analysis_lock import AnalysisLock


class LockTests(unittest.TestCase):
    def test_another_process_cannot_enter_until_release(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'analysis.lock'
            lock=AnalysisLock(path).acquire()
            code="from disorder_analysis_lock import AnalysisLock; import sys; lock=AnalysisLock(sys.argv[1],timeout=.2).acquire(); lock.release()"
            env=__import__('os').environ.copy();env['PYTHONPATH']=str(Path(__file__).resolve().parents[1]/'scripts')
            try:
                blocked=subprocess.run([sys.executable,'-c',code,str(path)],capture_output=True,text=True,env=env)
                self.assertNotEqual(blocked.returncode,0)
                self.assertIn('TimeoutError',blocked.stderr)
            finally:lock.release()
            allowed=subprocess.run([sys.executable,'-c',code,str(path)],capture_output=True,text=True,env=env)
            self.assertEqual(allowed.returncode,0,allowed.stderr)


if __name__=='__main__':unittest.main()
