"""Synthetic regression tests; no solver is launched.

Run: python -m unittest discover -s scripts -p test_tune_mu.py -v
"""
from contextlib import redirect_stdout
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

import h5py

sys.path.insert(0, str(Path(__file__).resolve().parent))
import tune_mu as tuner


def write_result(path, *, mean=9.0, error=.02, variance=1.0, count=4096,
                 squared_count=None, canonical=-1, mu=0.0, beta=2.0):
    with h5py.File(path, "w") as data:
        for key, value in {"canonical": canonical, "mu": mu, "beta": beta}.items():
            data[f"parameters/{key}"] = value
        group = data.create_group("simulation/results")
        group["Number_of_particles/mean/value"] = mean
        group["Number_of_particles/mean/error"] = error
        group["Number_of_particles/count"] = count
        group["Number_of_particles_squared/mean/value"] = mean * mean + variance
        group["Number_of_particles_squared/count"] = count if squared_count is None else squared_count


def estimate(mean, *, error=.01):
    block = tuner.Block(mean, error, 1.0, 4096, 4096.0)
    return tuner.summarize([[block, block], [block, block]], 2, 9, .1)


class Fixture(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.ini = self.root / "input.ini"
        self.exe = self.root / "qmc_worm.exe"
        self.exe.write_bytes(b"synthetic executable; never launched")
        self.output = self.root / "output"
        self.original = (
            'model = "BoseHubbard"\nLx = 3\nLy = 3\nLz = 1\n'
            'pbcx = 1\npbcy = 1\npbcz = 0\nnmax = 4\n'
            't_hop = .42\nU_on = 3.7\nV_nn = .12\nbeta = 2\n'
            'canonical = 9\nmu = 0\nthermalization = 20\nseed = 7\n'
            'onsite_potential_enabled = 0\nonsite_potential_file = "missing.dat"\n'
            'outputfile = "old.out.h5"\ncheckpoint = "old.state.h5"\n'
        )
        self.ini.write_text(self.original, encoding="utf-8")

    def args(self, *extra):
        return tuner.parser().parse_args([
            str(self.ini), "--executable", str(self.exe), "--output-dir", str(self.output),
            "--seconds", "1", "--max-seconds", "4", "--max-iterations", "6",
            "--mu-step", ".5", *extra,
        ])

    def run_tuner(self, evaluator, *extra):
        with redirect_stdout(io.StringIO()):
            return tuner.tune(self.args(*extra), evaluator=evaluator)

    def state(self):
        return json.loads((self.output / "result.json").read_text(encoding="utf-8"))


class SettingsTests(Fixture):
    def test_target_inferred_from_canonical_and_overridable(self):
        _, config = tuner.settings(self.args())
        self.assertEqual(config["target"], 9)
        _, config = tuner.settings(self.args("--target", "8.5"))
        self.assertEqual(config["target"], 8.5)

    def test_grand_canonical_source_requires_explicit_target(self):
        self.ini.write_text(self.original.replace("canonical = 9", "canonical = -1"), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "target"):
            tuner.settings(self.args())
        _, config = tuner.settings(self.args("--target", "9"))
        self.assertEqual(config["target"], 9)

    def test_disabled_missing_map_is_ignored(self):
        _, config = tuner.settings(self.args())
        self.assertIsNone(config["potential"])

    def test_enabled_missing_map_rejected(self):
        self.ini.write_text(self.original.replace("onsite_potential_enabled = 0", "onsite_potential_enabled = 1"), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "map"):
            tuner.settings(self.args())

    def test_invalid_targets_and_canonical_executable_rejected(self):
        for target in ("0", "36", "nan", "inf"):
            with self.subTest(target=target), self.assertRaises(ValueError):
                tuner.settings(self.args("--target", target))
        canonical_exe = self.root / "qmc_worm_canonical.exe"
        canonical_exe.write_bytes(b"dummy")
        with self.assertRaisesRegex(ValueError, "grand-canonical"):
            tuner.settings(self.args("--executable", str(canonical_exe)))

    def test_dry_run_leaves_no_output(self):
        evaluator = mock.Mock(side_effect=AssertionError("must not evaluate"))
        self.assertEqual(self.run_tuner(evaluator, "--dry-run"), 0)
        evaluator.assert_not_called()
        self.assertFalse(self.output.exists())


class ReadingTests(Fixture):
    def test_reads_number_fluctuations_and_effective_count(self):
        path = self.root / "block.h5"
        write_result(path)
        block = tuner.read_block(path, 0, 2)
        self.assertEqual(block.mean, 9)
        self.assertAlmostEqual(block.variance, 1)
        self.assertAlmostEqual(block.effective_samples, 2500)

    def test_insufficient_or_noisy_samples_are_rejected(self):
        for changes in ({"error": 0}, {"variance": 0}, {"variance": -1},
                        {"count": 127}, {"error": .2}, {"mean": float("nan")}):
            with self.subTest(changes=changes):
                path = self.root / "block.h5"
                write_result(path, **changes)
                with self.assertRaises(tuner.InsufficientSampling):
                    tuner.read_block(path, 0, 2)

    def test_wrong_counts_ensemble_and_parameters_are_rejected(self):
        for changes in ({"squared_count": 4095}, {"canonical": 9}, {"mu": .1}, {"beta": 3}):
            with self.subTest(changes=changes):
                path = self.root / "block.h5"
                write_result(path, **changes)
                with self.assertRaises(ValueError):
                    tuner.read_block(path, 0, 2)

    def test_missing_measurements_requests_more_sampling(self):
        path = self.root / "block.h5"
        with h5py.File(path, "w") as data:
            for key, value in {"canonical": -1, "mu": 0, "beta": 2}.items():
                data[f"parameters/{key}"] = value
        with self.assertRaises(tuner.InsufficientSampling):
            tuner.read_block(path, 0, 2)

    def test_drift_and_disagreeing_starts_prevent_acceptance(self):
        def block(value):
            return tuner.Block(value, .01, 1, 4096, 4096)
        drifting = tuner.summarize([[block(8.8), block(9.2)]] * 2, 2, 9, 1)
        self.assertFalse(drifting["stationary"])
        self.assertFalse(drifting["passed"])
        disagreeing = tuner.summarize([[block(8.9)] * 2, [block(9.1)] * 2], 2, 9, 1)
        self.assertFalse(disagreeing["starts_agree"])
        self.assertFalse(disagreeing["passed"])
        self.assertFalse(estimate(9, error=.08)["passed"])


class TuningTests(Fixture):
    def test_monotonic_response_converges_and_requires_fresh_confirmation(self):
        calls = []

        def evaluator(root, params, config, iteration, mu, seconds):
            calls.append((iteration, mu))
            return estimate(9 + 2 * (mu - .65))

        self.assertEqual(self.run_tuner(evaluator), 0)
        state = self.state()
        self.assertEqual(state["status"], "converged")
        self.assertAlmostEqual(state["mu"], .65)
        self.assertGreaterEqual(len(calls), 3)
        self.assertEqual(calls[-1][1], calls[-2][1])
        self.assertNotEqual(calls[-1][0], calls[-2][0])
        self.assertEqual(state["history"][-1]["phase"], "confirmation")
        self.assertEqual(state["history"][-2]["phase"], "tuning")
        tuned = tuner.read_ini(self.output / "tuned_parameters.ini")
        original = tuner.read_ini(self.ini)
        for key in ("Lx", "Ly", "Lz", "pbcx", "pbcy", "nmax", "t_hop", "U_on", "V_nn", "beta"):
            self.assertEqual(tuned[key], original[key], key)
        self.assertEqual(tuned["canonical"], "-1")
        self.assertEqual(self.ini.read_text(encoding="utf-8"), self.original)
        self.assertEqual((self.output / "source.ini").read_bytes(), self.ini.read_bytes())

    def test_failed_confirmation_does_not_create_tuned_parameters(self):
        responses = iter((estimate(9), estimate(9.5)))
        self.assertEqual(self.run_tuner(lambda *args: next(responses), "--max-iterations", "2"), 2)
        state = self.state()
        self.assertEqual(state["status"], "not_converged")
        self.assertEqual(state["history"][-1]["phase"], "confirmation")
        self.assertFalse((self.output / "tuned_parameters.ini").exists())

    def test_insufficient_sampling_increases_duration_without_moving_mu(self):
        calls = []

        def evaluator(root, params, config, iteration, mu, seconds):
            calls.append((mu, seconds))
            raise tuner.InsufficientSampling("number fluctuations unresolved")

        self.assertEqual(self.run_tuner(evaluator), 2)
        self.assertEqual(calls, [(0, 1), (0, 2), (0, 4)])
        state = self.state()
        self.assertEqual(state["status"], "not_converged")
        self.assertTrue(all("unresolved" in row["issue"] for row in state["history"]))
        self.assertFalse((self.output / "tuned_parameters.ini").exists())

    def test_noisy_target_estimate_is_not_accepted(self):
        self.assertEqual(self.run_tuner(lambda *args: estimate(9, error=.08), "--max-seconds", "1"), 2)
        self.assertEqual(self.state()["status"], "not_converged")

    def test_solver_failure_is_recorded_and_propagated(self):
        def evaluator(*args):
            raise RuntimeError("Solver failed (7); see production.log")

        with self.assertRaisesRegex(RuntimeError, "Solver failed"):
            self.run_tuner(evaluator)
        state = self.state()
        self.assertEqual(state["status"], "failed")
        self.assertIn("production.log", state["reason"])
        self.assertTrue((self.output / "history.csv").is_file())
        self.assertFalse((self.output / "tuned_parameters.ini").exists())

    def test_existing_results_are_never_overwritten(self):
        self.output.mkdir()
        marker = self.output / "result.json"
        marker.write_text("existing valuable results", encoding="utf-8")
        evaluator = mock.Mock()
        with self.assertRaisesRegex(ValueError, "empty output directory"):
            self.run_tuner(evaluator)
        evaluator.assert_not_called()
        self.assertEqual(marker.read_text(encoding="utf-8"), "existing valuable results")

    def test_enabled_potential_is_snapshotted_without_changing_values(self):
        potential = self.root / "potential.dat"
        potential.write_text("\n".join(str(i / 10) for i in range(9)), encoding="utf-8")
        source = self.original.replace("onsite_potential_enabled = 0", "onsite_potential_enabled = 1")
        self.ini.write_text(source.replace('"missing.dat"', '"potential.dat"'), encoding="utf-8")
        observed_paths = []

        def evaluator(root, params, config, iteration, mu, seconds):
            observed_paths.append(Path(tuner.unquote(params["onsite_potential_file"])))
            return estimate(9)

        self.assertEqual(self.run_tuner(evaluator), 0)
        self.assertTrue(all(path == self.output / "potential.dat" for path in observed_paths))
        self.assertEqual((self.output / "potential.dat").read_bytes(), potential.read_bytes())
        tuned = tuner.read_ini(self.output / "tuned_parameters.ini")
        self.assertEqual(Path(tuner.unquote(tuned["onsite_potential_file"])), self.output / "potential.dat")

    def test_evaluation_discards_burn_in_resets_blocks_and_uses_fresh_seeds(self):
        params, config = tuner.settings(self.args("--blocks", "2"))
        self.output.mkdir()
        launched = []

        def fake_launch(exe, arguments, folder, seconds, label):
            values = tuner.read_ini(folder / "parameters.ini")
            launched.append((folder, list(arguments), label, values))
            (folder / "state.h5").write_bytes(b"synthetic checkpoint")
            write_result(folder / "result.h5", mean=100 if label.startswith("burn") else 9,
                         mu=float(values["mu"]), beta=float(values["beta"]))

        with mock.patch.object(tuner, "launch", side_effect=fake_launch):
            first = tuner.evaluate(self.output, params, config, 1, .65, 1)
            second = tuner.evaluate(self.output, params, config, 2, .65, 1)
        self.assertEqual(first["mean_number"], 9)
        self.assertEqual(second["mean_number"], 9)
        self.assertEqual(len(launched), 12)
        starts, seeds = set(), set()
        for folder, arguments, label, values in launched:
            self.assertEqual(values["canonical"], "-1")
            self.assertEqual(values["beta"], params["beta"])
            self.assertEqual(values["U_on"], params["U_on"])
            if label.startswith("production"):
                self.assertEqual(arguments, ["state.h5", "--reset-statistics"])
            else:
                self.assertEqual(arguments, ["parameters.ini"])
                seeds.add(values["seed"])
                starts.add(values["initial_particle_number"])
        self.assertEqual(len(seeds), 4)
        self.assertEqual(len(starts), 2)
        self.assertEqual(len(list(self.output.glob("trial_*/chain_*/production_*.h5"))), 8)


class LaunchTests(Fixture):
    def test_zero_exit_does_not_hide_stale_or_missing_checkpoint(self):
        for scenario in ("stale", "warning", "missing", "nonzero"):
            with self.subTest(scenario=scenario):
                folder = self.root / scenario
                folder.mkdir()
                if scenario == "stale":
                    (folder / "state.h5").write_bytes(b"old checkpoint")
                (folder / "result.h5").write_bytes(b"result")

                def fake_run(command, **kwargs):
                    if scenario == "warning":
                        (folder / "state.h5").write_bytes(b"saved but corrupt")
                        kwargs["stdout"].write("WARNING: Failed to save checkpoint\n")
                    return mock.Mock(returncode=7 if scenario == "nonzero" else 0)

                with mock.patch.object(tuner.subprocess, "run", side_effect=fake_run):
                    with redirect_stdout(io.StringIO()), self.assertRaises(RuntimeError):
                        tuner.launch(self.exe, ["state.h5"], folder, 1, "production")

    def test_updated_checkpoint_is_accepted(self):
        (self.root / "state.h5").write_bytes(b"old checkpoint")

        def fake_run(command, **kwargs):
            (self.root / "state.h5").write_bytes(b"fresh checkpoint")
            (self.root / "result.h5").write_bytes(b"fresh result")
            return mock.Mock(returncode=0)

        with mock.patch.object(tuner.subprocess, "run", side_effect=fake_run) as run:
            with redirect_stdout(io.StringIO()):
                tuner.launch(self.exe, ["state.h5", "--reset-statistics"], self.root, 1, "production")
        run.assert_called_once()

if __name__ == "__main__":
    unittest.main()