"""Numerical and interchange checks; run: python -m unittest discover -s tests."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import warnings

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from generate_speckle import generate_speckle, write_potential


class SpeckleTests(unittest.TestCase):
    def test_rectangular_map_mean_strength_and_seed(self):
        options = dict(nx=48, ny=32, dx=0.5, dy=0.75, correlation_length=1.5, mean_strength=2.7)
        first = generate_speckle(seed=7, **options)
        second = generate_speckle(seed=7, **options)
        different = generate_speckle(seed=8, **options)
        self.assertEqual(first.shape, (32, 48))
        self.assertTrue(np.all(first >= 0))
        self.assertAlmostEqual(float(first.mean()), 2.7, places=13)
        np.testing.assert_array_equal(first, second)
        self.assertFalse(np.array_equal(first, different))
        rescaled = generate_speckle(seed=7, **dict(options, mean_strength=5.4))
        np.testing.assert_allclose(rescaled, 2 * first, rtol=2e-15)

    def test_global_random_state_unchanged(self):
        np.random.seed(912)
        before = np.random.get_state()
        generate_speckle(24, correlation_length=2, seed=42)
        after = np.random.get_state()
        self.assertEqual(before[0], after[0])
        np.testing.assert_array_equal(before[1], after[1])
        self.assertEqual(before[2:], after[2:])

    def test_subtract_mean_preserves_fluctuations(self):
        positive = generate_speckle(32, correlation_length=2, mean_strength=3.0, seed=1)
        centered = generate_speckle(32, correlation_length=2, mean_strength=3.0, seed=1, subtract_mean=True)
        np.testing.assert_allclose(centered, positive - positive.mean(), atol=1e-14)
        self.assertAlmostEqual(float(centered.mean()), 0.0, places=13)
        self.assertAlmostEqual(float(centered.std()), float(positive.std()), places=13)
        self.assertLess(float(centered.min()), 0.0)

    def test_zero_strength(self):
        for subtract in [False, True]:
            actual = generate_speckle(16, correlation_length=2, mean_strength=0, subtract_mean=subtract)
            np.testing.assert_array_equal(actual, np.zeros((16, 16)))

    def test_length_units_rescale_without_changing_map(self):
        first = generate_speckle(48, 40, dx=0.5, dy=0.75, correlation_length=1.5, seed=5)
        rescaled = generate_speckle(48, 40, dx=5, dy=7.5, correlation_length=15, seed=5)
        np.testing.assert_allclose(first, rescaled, atol=5e-15, rtol=5e-15)

    def test_intensity_correlation_matches_requested_length(self):
        # Ensemble averaging tests the physical length convention without
        # demanding a particular fitted length for each finite realization.
        for ell in [2.0, 4.0]:
            correlations = np.zeros(3)
            for seed in range(40):
                values = generate_speckle(96, correlation_length=ell, seed=seed)
                fluctuations = values - values.mean()
                for index, lag in enumerate([0, int(ell), int(2 * ell)]):
                    correlations[index] += np.mean(fluctuations * np.roll(fluctuations, lag, axis=1))
            correlations /= correlations[0]
            np.testing.assert_allclose(correlations, np.exp(-0.5 * np.arange(3)**2), atol=0.035, rtol=0)

    def test_finite_mode_correlation(self):
        # For a relatively small box, compare against the exact phase-ensemble
        # prediction, including the negative background from fixed amplitudes.
        size, ell = 24, 3.0
        frequency = np.fft.fftfreq(size)
        power = np.exp(-4 * np.pi**2 * ell**2 * (frequency[:, None]**2 + frequency[None, :]**2))
        field_corr = np.fft.ifft2(power).real
        field_corr /= field_corr[0, 0]
        q = np.sum(power**2) / np.sum(power)**2
        expected = (field_corr[0]**2 - q) / (1 - q)
        actual = np.zeros(size)
        for seed in range(250):
            values = generate_speckle(size, correlation_length=ell, seed=seed, subtract_mean=True)
            for lag in range(size):
                actual[lag] += np.mean(values * np.roll(values, lag, axis=1))
        actual /= actual[0]
        np.testing.assert_allclose(actual, expected, atol=0.045, rtol=0)
        self.assertLess(actual[size // 2], 0)

    def test_anisotropic_spacing_gives_isotropic_physical_length(self):
        correlation_x, correlation_y = [], []
        for seed in range(30):
            values = generate_speckle(96, 48, dx=0.5, dy=1, correlation_length=2, seed=seed)
            centered = values - values.mean()
            variance = np.mean(centered**2)
            correlation_x.append(np.mean(centered * np.roll(centered, 4, axis=1)) / variance)
            correlation_y.append(np.mean(centered * np.roll(centered, 2, axis=0)) / variance)
        self.assertAlmostEqual(float(np.mean(correlation_x)), np.exp(-0.5), delta=0.04)
        self.assertAlmostEqual(float(np.mean(correlation_y)), np.exp(-0.5), delta=0.04)

    def test_bad_parameters_fail(self):
        for overrides in [
            dict(nx=0), dict(nx=2.5), dict(nx=True), dict(ny=1), dict(dx=0), dict(dy=-1),
            dict(correlation_length=0), dict(correlation_length=float("nan")),
            dict(mean_strength=-1), dict(mean_strength=float("inf")), dict(seed=-1), dict(seed=1.2),
        ]:
            with self.subTest(overrides=overrides), self.assertRaises(ValueError):
                generate_speckle(**dict(dict(nx=24, correlation_length=2), **overrides))

    def test_sampling_warnings(self):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            generate_speckle(16, dx=1, correlation_length=0.2, seed=0)
            generate_speckle(16, dx=1, correlation_length=4, seed=0)
        self.assertTrue(any("shorter than a grid spacing" in str(item.message) for item in caught))
        self.assertTrue(any("six correlation lengths" in str(item.message) for item in caught))

    def test_export_has_x_fastest_and_roundtrips(self):
        values = np.array([[0.1, 1.0, 2.0], [10.0, 11.0, -12.25]])
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "new_directory" / "potential.dat"
            write_potential(path, values, metadata={"label": "arbitrary programmable map"})
            np.testing.assert_array_equal(np.loadtxt(path), [0.1, 1, 2, 10, 11, -12.25])
            self.assertTrue(path.read_bytes().startswith(b"#"))
            self.assertIn("site=x+nx*y", path.read_text())

    def test_export_rejects_invalid_maps(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "potential.dat"
            for values in [np.ones(4), np.empty((0, 2)), [[float("nan")]], [[1 + 1j]]]:
                with self.subTest(values=values), self.assertRaises(ValueError):
                    write_potential(path, values)
            self.assertFalse(path.exists())

    def test_cli_metadata_and_energy_values_agree(self):
        script = Path(__file__).resolve().parents[1] / "scripts" / "generate_speckle.py"
        with tempfile.TemporaryDirectory() as directory:
            potential_path = Path(directory) / "maps" / "potential.dat"
            metadata_path = Path(directory) / "metadata" / "potential.json"
            completed = subprocess.run(
                [sys.executable, str(script), str(potential_path), "--nx", "32", "--ny", "24",
                 "--correlation-length", "2", "--mean-strength", "1.7", "--seed", "31",
                 "--subtract-mean", "--metadata-json", str(metadata_path)],
                capture_output=True, text=True, check=True,
            )
            self.assertIn("32 x 24", completed.stdout)
            values = np.loadtxt(potential_path).reshape(24, 32)
            metadata = json.loads(metadata_path.read_text())
            self.assertEqual(metadata["seed"], 31)
            self.assertEqual(metadata["dx"], 1)
            self.assertEqual(metadata["dy"], 1)
            self.assertEqual(metadata["correlation_length"], 2)
            self.assertAlmostEqual(metadata["mean"], float(values.mean()), places=14)
            self.assertAlmostEqual(metadata["standard_deviation"], float(values.std()), places=14)
            np.testing.assert_array_equal(values, generate_speckle(32, 24, correlation_length=2, mean_strength=1.7, seed=31, subtract_mean=True))

    def test_cli_rejects_metadata_overwriting_potential(self):
        script = Path(__file__).resolve().parents[1] / "scripts" / "generate_speckle.py"
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "potential.dat"
            completed = subprocess.run(
                [sys.executable, str(script), str(path), "--nx", "24", "--correlation-length", "2", "--metadata-json", str(path)],
                capture_output=True, text=True,
            )
            self.assertNotEqual(completed.returncode, 0)
            self.assertIn("must be different files", completed.stderr)
            self.assertFalse(path.exists())


if __name__ == "__main__":
    unittest.main()
