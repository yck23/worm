"""Exercise dependency ordering and failure handling without running QMC."""
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import run_disorder_queue as queue


class QueueTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.addCleanup(mock.patch.stopall)
        mock.patch.object(queue, 'REPO', self.root).start()
        self.dependency = self.root/'dependency.json'
        self.save(self.dependency, dict(status='running', controller_pid=123, points=[dict(status='running')]))
        self.scans = [self.root/'beta045.json', self.root/'beta052.json']
        for scan in self.scans:
            self.save(scan, dict(status='prepared', points=[dict(status='prepared')]))
        self.path = self.root/'queue.json'
        self.save(self.path, dict(status='prepared', dependency='dependency.json', poll_seconds=60,
            stages=[dict(beta=beta, scan=scan.name, status='queued')
                    for beta,scan in zip([.45,.52],self.scans)]))

    def save(self, path, state):
        path.write_text(json.dumps(state))

    def complete_dependency(self):
        self.save(self.dependency, dict(status='completed', points=[dict(status='completed')]))

    def test_wait_then_run_in_order(self):
        started = []
        def sleep(seconds):
            self.assertEqual(seconds, 60)
            self.assertEqual(started, [])
            self.complete_dependency()
        def execute(scan):
            self.assertEqual(queue.read(self.dependency)['status'], 'completed')
            started.append(scan)
            self.save(scan, dict(status='completed', points=[dict(status='completed')]))
            return 0
        with mock.patch.object(queue, 'process_alive', return_value=True), \
             mock.patch.object(queue.time, 'sleep', side_effect=sleep), \
             mock.patch.object(queue, 'execute_stage', side_effect=execute):
            self.assertEqual(queue.run(self.path), 0)
        self.assertEqual(started, self.scans)
        self.assertEqual(queue.read(self.path)['status'], 'completed')
        self.assertFalse(self.path.with_suffix('.lock').exists())

    def test_failed_dependency_starts_nothing(self):
        self.save(self.dependency, dict(status='failed', points=[]))
        with mock.patch.object(queue, 'execute_stage') as execute:
            self.assertEqual(queue.run(self.path), 1)
            execute.assert_not_called()
        self.assertIn('Dependency stopped', queue.read(self.path)['error'])

    def test_dead_dependency_requires_three_checks(self):
        with mock.patch.object(queue, 'process_alive', return_value=False), \
             mock.patch.object(queue.time, 'sleep') as sleep, \
             mock.patch.object(queue, 'execute_stage') as execute:
            self.assertEqual(queue.run(self.path), 1)
            self.assertEqual(sleep.call_count, 2)
            execute.assert_not_called()
        self.assertIn('three consecutive', queue.read(self.path)['error'])

    def test_normal_failed_scan_allows_next_independent_beta(self):
        self.complete_dependency()
        def execute(scan):
            failed = scan == self.scans[0]
            self.save(scan, dict(status='failed' if failed else 'completed', points=[]))
            return int(failed)
        with mock.patch.object(queue, 'execute_stage', side_effect=execute) as execute_mock:
            self.assertEqual(queue.run(self.path), 1)
            self.assertEqual(execute_mock.call_count, 2)
        self.assertEqual(queue.read(self.path)['status'], 'completed_with_failures')

    def test_interrupted_scan_stops_before_next_beta(self):
        self.complete_dependency()
        def execute(scan):
            self.save(scan, dict(status='running', points=[]))
            return 1
        with mock.patch.object(queue, 'execute_stage', side_effect=execute) as execute_mock:
            self.assertEqual(queue.run(self.path), 1)
            self.assertEqual(execute_mock.call_count, 1)
        self.assertEqual(queue.read(self.scans[1])['status'], 'prepared')

    def test_wait_for_parent_queue_then_start(self):
        self.save(self.dependency, dict(kind='disorder_scan_queue',status='running',
            controller_pid=123,stages=[dict(status='running'),dict(status='queued')]))
        def sleep(seconds):
            self.save(self.dependency,dict(kind='disorder_scan_queue',status='completed',
                stages=[dict(status='completed'),dict(status='completed')]))
        def execute(scan):
            self.assertEqual(queue.read(self.dependency)['status'],'completed')
            self.save(scan,dict(status='completed',points=[]))
            return 0
        with mock.patch.object(queue,'process_alive',return_value=True), \
             mock.patch.object(queue.time,'sleep',side_effect=sleep) as wait, \
             mock.patch.object(queue,'execute_stage',side_effect=execute) as execute_mock:
            self.assertEqual(queue.run(self.path),0)
            self.assertEqual(wait.call_count,1)
            self.assertEqual(execute_mock.call_count,2)

    def test_parent_queue_must_complete_every_stage(self):
        self.save(self.dependency,dict(kind='disorder_scan_queue',status='completed',
            stages=[dict(status='completed'),dict(status='failed')]))
        with mock.patch.object(queue,'execute_stage') as execute:
            self.assertEqual(queue.run(self.path),1)
            execute.assert_not_called()
        self.assertIn('incomplete points or stages',queue.read(self.path)['error'])

    def test_failed_parent_queue_starts_nothing(self):
        self.save(self.dependency,dict(kind='disorder_scan_queue',status='completed_with_failures',
            stages=[dict(status='failed'),dict(status='completed')]))
        with mock.patch.object(queue,'execute_stage') as execute:
            self.assertEqual(queue.run(self.path),1)
            execute.assert_not_called()

    def test_exclusive_lock_preserves_active_queue(self):
        self.path.with_suffix('.lock').write_text('9999')
        with self.assertRaises(FileExistsError):
            queue.run(self.path)
        self.assertEqual(queue.read(self.path)['status'], 'prepared')
        self.assertEqual(self.path.with_suffix('.lock').read_text(), '9999')


if __name__ == '__main__':
    unittest.main()
