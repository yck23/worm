from pathlib import Path
import json,sys,tempfile,unittest
from unittest import mock
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import disorder_response_scan as scan
from disorder_analysis_lock import AnalysisLock


class RefreshLockTests(unittest.TestCase):
    def test_local_lock_is_released_before_paired_export(self):
        with tempfile.TemporaryDirectory() as folder:
            repo=Path(folder);path=repo/'scan.json';path.write_text(json.dumps(dict(config=dict(refresh_paired_notebook=True))))
            def exporter(*args,**kwargs):
                lock=AnalysisLock(repo/'results/.disorder_analysis.lock',timeout=.1).acquire();lock.release()
                return mock.Mock(returncode=0)
            with mock.patch.object(scan,'REPO',repo),mock.patch.object(scan,'_analyze_local') as local, \
                 mock.patch.object(scan.subprocess,'run',side_effect=exporter) as paired:
                scan.analyze(path)
                local.assert_called_once();paired.assert_called_once()

    def test_paired_failure_is_recorded_after_lock_release(self):
        with tempfile.TemporaryDirectory() as folder:
            repo=Path(folder);path=repo/'scan.json';path.write_text(json.dumps(dict(config=dict(refresh_paired_notebook=True))))
            with mock.patch.object(scan,'REPO',repo),mock.patch.object(scan,'_analyze_local'), \
                 mock.patch.object(scan.subprocess,'run',return_value=mock.Mock(returncode=1)):
                scan.analyze(path)
            error=json.loads((repo/'analysis/paired_refresh_error.json').read_text())
            self.assertIn('Paired notebook refresh failed',error['error'])

    def test_local_error_releases_lock(self):
        with tempfile.TemporaryDirectory() as folder:
            repo=Path(folder);path=repo/'scan.json';path.write_text('{}')
            with mock.patch.object(scan,'REPO',repo), \
                 mock.patch.object(scan,'_analyze_local',side_effect=RuntimeError('fixture failure')):
                with self.assertRaisesRegex(RuntimeError,'fixture failure'):scan.analyze(path)
            lock=AnalysisLock(repo/'results/.disorder_analysis.lock',timeout=.1).acquire();lock.release()


if __name__=='__main__':unittest.main()
