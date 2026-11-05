"""Regression checks for quenched-disorder jobs and statistical reductions."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock
import warnings
import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import run_disorder_ensemble as launcher
from summarize_disorder_ensemble import summarize_manifest


class Fixture(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.ini = self.root / "source.ini"
        self.ini.write_text(
            'model = "BoseHubbard"\nLx = 8\nLy = 8\nLz = 1\n'
            'pbcx = 1\npbcy = 1\npbcz = 0\nt_hop = 1\nU_on = 4\nV_nn = 0\n'
            'mu = 0.7\nbeta = 1\nnmax = 4\ncanonical = -1\ncanonical_window = 0.5\n'
            'seed = 123\nsweeps = 100\nthermalization = 10\n'
            'onsite_potential_file = "unused.dat"\nonsite_potential_enabled = 0\n'
            'outputfile = "source.out.h5"\ncheckpoint = "source.clone.h5"\n',
            encoding="utf-8")

    def prepare(self, name="ensemble", **changes):
        options = dict(realizations=3, workers=2, correlation_length=1,
                       mean_strength=0.7, base_disorder_seed=41,
                       base_mc_seed=5001, prepare_only=True)
        options.update(changes)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            return Path(launcher.prepare_ensemble(self.ini, self.root / name, **options))

    @staticmethod
    def read(path):
        return json.loads(Path(path).read_text(encoding="utf-8"))

    @staticmethod
    def write(path, manifest):
        Path(path).write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


class LauncherTests(Fixture):
    def test_maps_and_seeds_are_independent_of_workers(self):
        one = self.prepare("one", workers=1, realizations=4)
        four = self.prepare("four", workers=4, realizations=4)
        left, right = self.read(one), self.read(four)
        for a, b in zip(left["runs"], right["runs"]):
            for key in ("id", "index", "disorder_seed", "mc_seed", "potential_sha256"):
                self.assertEqual(a[key], b[key])
            np.testing.assert_array_equal(np.loadtxt(one.parent / a["potential_file"]),
                                          np.loadtxt(four.parent / b["potential_file"]))
        for key in ("disorder_seed", "mc_seed", "potential_sha256"):
            self.assertEqual(len({run[key] for run in left["runs"]}), 4)

    def test_realizations_have_separate_inputs_and_outputs(self):
        original = self.ini.read_bytes()
        path = self.prepare(realizations=4)
        manifest = self.read(path)
        for key in ("directory", "input_file", "potential_file", "output_file",
                    "checkpoint_file", "log_file"):
            self.assertEqual(len({run[key] for run in manifest["runs"]}), 4, key)
        for run in manifest["runs"]:
            directory = (path.parent / run["directory"]).resolve()
            params = launcher.read_ini(path.parent / run["input_file"])
            self.assertEqual(int(params["seed"]), run["mc_seed"])
            self.assertIn(params["onsite_potential_enabled"].lower(), ("1", "true"))
            for a, b in (("outputfile", "output_file"), ("checkpoint", "checkpoint_file")):
                self.assertEqual(Path(params[a].strip('"')).name, Path(run[b]).name)
            for key in ("input_file", "potential_file", "output_file", "checkpoint_file"):
                self.assertEqual((path.parent / run[key]).resolve().parent, directory)
            values = np.loadtxt(path.parent / run["potential_file"])
            self.assertAlmostEqual(float(values.mean()), 0.7, places=13)
            self.assertEqual(run["potential_sha256"],
                             hashlib.sha256(np.asarray(values, dtype="<f8").tobytes()).hexdigest())
        self.assertEqual(self.ini.read_bytes(), original)

    def test_existing_ensemble_is_not_overwritten(self):
        path = self.prepare()
        snapshot = lambda: {p.relative_to(path.parent): p.read_bytes()
                            for p in path.parent.rglob("*") if p.is_file()}
        before = snapshot()
        with self.assertRaises((ValueError, FileExistsError, RuntimeError)):
            self.prepare()
        self.assertEqual(before, snapshot())

    def test_existing_user_file_is_preserved(self):
        directory = self.root / "occupied"
        directory.mkdir()
        sentinel = directory / "keep.txt"
        sentinel.write_text("prior evidence")
        with self.assertRaises((ValueError, FileExistsError, RuntimeError)):
            self.prepare("occupied")
        self.assertEqual(sentinel.read_text(), "prior evidence")
        self.assertEqual(list(directory.iterdir()), [sentinel])

    def test_canonical_flag_must_match_integer_particle_number(self):
        original = self.ini.read_text()
        for index, (number, flag) in enumerate([("-1", True), ("64", False), ("1.5", True)]):
            with self.subTest(number=number, flag=flag):
                self.ini.write_text(original.replace("canonical = -1", f"canonical = {number}"))
                with self.assertRaises((ValueError, RuntimeError)):
                    self.prepare(f"invalid_{index}", canonical=flag)
        self.ini.write_text(original.replace("canonical = -1", "canonical = 64"))
        path = self.prepare("canonical", canonical=True)
        self.assertTrue(self.read(path)["canonical"])
        for run in self.read(path)["runs"]:
            params = launcher.read_ini(path.parent / run["input_file"])
            self.assertEqual(int(params["canonical"]), 64)
            self.assertEqual(float(params["canonical_window"]), 0.5)

    def test_nonzero_worker_exit_is_failure_and_preserves_evidence(self):
        path = self.prepare(workers=1, executable=Path(sys.executable))
        def fail(command, **kwargs):
            message = "deliberate worker failure\n"
            if hasattr(kwargs.get("stdout"), "write"):
                kwargs["stdout"].write(message)
            return subprocess.CompletedProcess(command, 7, stdout=message)
        with mock.patch.object(launcher.subprocess, "run", side_effect=fail):
            self.assertNotEqual(launcher.run_ensemble(path), 0)
        manifest = self.read(path)
        self.assertEqual(manifest["status"], "failed")
        failed = [run for run in manifest["runs"] if run["status"] == "failed"]
        self.assertTrue(failed)
        self.assertEqual(failed[0]["returncode"], 7)
        self.assertIn("deliberate worker failure", (path.parent / failed[0]["log_file"]).read_text())
        self.assertFalse(any(run["status"] == "completed" for run in manifest["runs"]))
        for run in manifest["runs"]:
            self.assertTrue((path.parent / run["input_file"]).is_file())
            self.assertTrue((path.parent / run["potential_file"]).is_file())

    def test_existing_result_stops_all_workers_before_overwrite(self):
        path = self.prepare(workers=2, executable=Path(sys.executable))
        run = self.read(path)["runs"][-1]
        output = path.parent / run["output_file"]
        output.write_bytes(b"previous result")
        before = path.read_bytes()
        with mock.patch.object(launcher.subprocess, "run") as process:
            with self.assertRaises((ValueError, RuntimeError)):
                launcher.run_ensemble(path)
        process.assert_not_called()
        self.assertEqual(output.read_bytes(), b"previous result")
        self.assertEqual(path.read_bytes(), before)

    def test_zero_exit_without_result_is_failure(self):
        path = self.prepare(realizations=1, workers=1, executable=Path(sys.executable))
        with mock.patch.object(launcher.subprocess, "run",
                               return_value=subprocess.CompletedProcess([], 0, stdout="")):
            self.assertNotEqual(launcher.run_ensemble(path), 0)
        manifest = self.read(path)
        self.assertEqual(manifest["status"], "failed")
        self.assertEqual(manifest["runs"][0]["status"], "failed")


class SummaryTests(Fixture):
    def completed(self, name="ensemble", *, realizations=2):
        path = self.prepare(name, realizations=realizations)
        manifest = self.read(path)
        for index, run in enumerate(manifest["runs"]):
            n, error, count = 1 + 2 * index, 0.2 * (index + 1), [100, 1000][index % 2]
            results = {"Total_Energy": (n, error),
                       "Number_of_particles": (n, error),
                       "Number_of_particles_squared": (n * n + 1, error),
                       "Density_Distribution": (np.full(64, n / 64), np.full(64, error / 64)),
                       "Winding_number_squared": ([n, n + 1], [error, 2 * error])}
            with h5py.File(path.parent / run["output_file"], "w") as data:
                data["hamiltonian/onsite_potential"] = np.loadtxt(path.parent / run["potential_file"])
                for observable, (mean, uncertainty) in results.items():
                    group = data.require_group(f"simulation/results/{observable}")
                    group["mean/value"] = np.asarray(mean, dtype=float)
                    group["mean/error"] = uncertainty
                    group["count"] = count
            run.update(status="completed", returncode=0)
        manifest["status"] = "completed"
        self.write(path, manifest)
        return path

    def test_equal_weighting_and_no_thermal_error_double_counting(self):
        result = summarize_manifest(self.completed())
        self.assertEqual(result["realizations"], 2)
        self.assertEqual(result["weighting"], "equal_per_realization")
        energy = result["observables"]["Total_Energy"]
        self.assertEqual(energy["mean"], 2)
        self.assertNotAlmostEqual(energy["mean"], (100 * 1 + 1000 * 3) / 1100)
        self.assertAlmostEqual(energy["standard_error"], 1)
        self.assertAlmostEqual(energy["thermal_standard_error"], np.hypot(0.2, 0.4) / 2)
        self.assertEqual(energy["counts"], [100, 1000])
        self.assertEqual(energy["realization_means"], [1, 3])
        self.assertNotAlmostEqual(energy["standard_error"], np.sqrt(1 + 0.05))
        winding = result["observables"]["Winding_number_squared"]
        np.testing.assert_allclose(winding["mean"], [2, 3])
        np.testing.assert_allclose(winding["standard_error"], [1, 1])
        np.testing.assert_allclose(winding["thermal_standard_error"], [np.sqrt(0.05), np.sqrt(0.2)])

    def test_nonlinear_variance_is_computed_before_disorder_average(self):
        variance = summarize_manifest(self.completed())["observables"]["Number_variance"]
        self.assertEqual(variance["mean"], 1)
        self.assertEqual(variance["realization_means"], [1, 1])
        self.assertEqual(variance["standard_error"], 0)
        self.assertIsNone(variance["thermal_standard_error"])
        # Subtracting squared disorder mean would instead give 6 - 2**2 = 2.
        self.assertNotEqual(variance["mean"], 2)

    def test_single_realization_has_no_disorder_sem(self):
        energy = summarize_manifest(self.completed(realizations=1))["observables"]["Total_Energy"]
        self.assertEqual(energy["mean"], 1)
        self.assertIsNone(energy["standard_error"])
        self.assertAlmostEqual(energy["thermal_standard_error"], 0.2)

    def test_nonfinite_thermal_error_components_are_null(self):
        path = self.completed()
        with h5py.File(path.parent / self.read(path)["runs"][0]["output_file"], "a") as data:
            data["simulation/results/Total_Energy/mean/error"][...] = np.nan
            data["simulation/results/Winding_number_squared/mean/error"][...] = [0.2, np.nan]
        result = summarize_manifest(path)["observables"]
        self.assertIsNone(result["Total_Energy"]["thermal_standard_error"])
        self.assertEqual(result["Total_Energy"]["mean"], 2)
        self.assertAlmostEqual(result["Winding_number_squared"]["thermal_standard_error"][0], np.sqrt(0.05))
        self.assertIsNone(result["Winding_number_squared"]["thermal_standard_error"][1])
        json.dumps(result, allow_nan=False)

    def test_failed_or_unfinished_realizations_are_not_omitted(self):
        path = self.completed()
        original = self.read(path)
        for status, code in [("failed", 9), ("running", None), ("prepared", None)]:
            with self.subTest(status=status):
                manifest = json.loads(json.dumps(original))
                manifest["runs"][1].update(status=status, returncode=code)
                self.write(path, manifest)
                with self.assertRaises((ValueError, RuntimeError)):
                    summarize_manifest(path)

    def test_missing_result_is_not_omitted(self):
        path = self.completed()
        (path.parent / self.read(path)["runs"][1]["output_file"]).unlink()
        with self.assertRaises((ValueError, RuntimeError, FileNotFoundError)):
            summarize_manifest(path)

    def test_nonfinite_means_and_zero_counts_are_refused(self):
        for index, (dataset, value) in enumerate([("mean/value", np.nan), ("count", 0)]):
            with self.subTest(dataset=dataset):
                path = self.completed(f"invalid_{index}")
                with h5py.File(path.parent / self.read(path)["runs"][1]["output_file"], "a") as data:
                    data[f"simulation/results/Total_Energy/{dataset}"][...] = value
                with self.assertRaises((ValueError, RuntimeError)):
                    summarize_manifest(path)

    def test_changed_hamiltonian_map_is_refused(self):
        path = self.completed()
        with h5py.File(path.parent / self.read(path)["runs"][1]["output_file"], "a") as data:
            data["hamiltonian/onsite_potential"][0] += 1
        with self.assertRaises((ValueError, RuntimeError)):
            summarize_manifest(path)


if __name__ == "__main__":
    unittest.main()
