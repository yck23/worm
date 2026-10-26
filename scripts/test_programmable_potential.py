#!/usr/bin/env python3
"""End-to-end programmable-potential checks against a small exact Hamiltonian.

Requires numpy, h5py and the built square-lattice executable. All generated
inputs, logs and output are confined to a new subdirectory of --work-dir.
"""
from __future__ import annotations

import argparse
import itertools
import json
import os
from pathlib import Path
import subprocess
import tempfile

import h5py
import numpy as np


def exact_observables(potential, *, hopping=0.35, mu=0.25, beta=1.3,
                      nmax=1, interaction=1.7, canonical=-1, v_nn=0.0):
    """Full dense diagonalization of a 2x2 open Bose-Hubbard lattice."""
    basis = np.array([n for n in itertools.product(range(nmax + 1), repeat=4)
                      if canonical < 0 or sum(n) == canonical])
    lookup = {tuple(n): i for i, n in enumerate(basis)}
    bonds = [(0, 1), (0, 2), (1, 3), (2, 3)]
    diagonal = ((potential - mu) * basis +
                interaction * basis * (basis - 1) / 2).sum(axis=1)
    diagonal += sum(v_nn * basis[:, i] * basis[:, j] for i, j in bonds)
    matrix = np.diag(diagonal)
    for col, n in enumerate(basis):
        for i, j in bonds:
            for dest, source in [(i, j), (j, i)]:
                if n[source] > 0 and n[dest] < nmax:
                    new = n.copy()
                    new[dest] += 1
                    new[source] -= 1
                    matrix[lookup[tuple(new)], col] -= hopping * np.sqrt(
                        (n[dest] + 1) * n[source])
    energies, eigenvectors = np.linalg.eigh(matrix)
    weights = np.exp(-beta * (energies - energies.min()))
    weights /= weights.sum()
    probabilities = (eigenvectors ** 2) @ weights
    return {"Density_Distribution": probabilities @ basis,
            "Total_Energy": weights @ energies,
            "Number_of_particles": probabilities @ basis.sum(axis=1)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--canonical-executable", type=Path)
    parser.add_argument("--uniform-executable", type=Path)
    parser.add_argument("--mpi-executable", type=Path)
    parser.add_argument("--canonical-mpi-executable", type=Path)
    parser.add_argument("--mpiexec", default="mpiexec")
    parser.add_argument("--work-dir", type=Path, default=Path("build/potential-validation"))
    parser.add_argument("--sweeps", type=int, default=250000)
    parser.add_argument("--only-disabled", action="store_true",
                        help="Run only potential-disable, restore and override checks")
    args = parser.parse_args()
    for name in ["executable", "canonical_executable", "uniform_executable",
                 "mpi_executable", "canonical_mpi_executable"]:
        value = getattr(args, name)
        if value:
            setattr(args, name, value.resolve(strict=True))
    args.work_dir.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix="check-", dir=args.work_dir.resolve()))
    env = os.environ.copy()
    env["TMPDIR"] = work.as_posix() + "/"
    env["PATH"] = str(args.executable.parent) + os.pathsep + env.get("PATH", "")
    report = []

    def run(executable, inputs, name, ranks=1, expected_success=True):
        command = [str(executable), *map(str, inputs)]
        if ranks > 1:
            command = [args.mpiexec, "-n", str(ranks), *command]
        result = subprocess.run(command, cwd=work, env=env, text=True,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                timeout=180)
        (work / (name + ".log")).write_text(result.stdout, encoding="utf-8")
        if (result.returncode == 0) != expected_success:
            raise AssertionError(f"{name}: unexpected exit {result.returncode}; see {work / (name + '.log')}")
        return result.stdout

    def create_case(name, potential, *, executable=None, ranks=1, **overrides):
        values = dict(model='"BoseHubbard"', runtimelimit=120, sweeps=args.sweeps,
                      thermalization=20000, Lx=2, Ly=2, Lz=1, pbcx=0, pbcy=0,
                      pbcz=0, t_hop=0.35, U_on=1.7, V_nn=0.0, mu=0.25, beta=1.3,
                      nmax=1, E_off=1.0, canonical=-1, initial_occupancy=0,
                      seed=83192, Ntest=10000, Nsave=100000000,
                      Nmeasure=1, Nmeasure2=1, C_worm=2.0,
                      p_insertworm=1.0, p_moveworm=0.4, p_insertkink=0.3,
                      p_deletekink=0.2, p_glueworm=0.1,
                      outputfile=f'"{name}.out.h5"', checkpoint=f'"{name}.clone.h5"')
        if potential is not None:
            np.savetxt(work / (name + ".dat"), potential, fmt="%.17g",
                       header="site = x + Lx*y; simulation energy units")
            values["onsite_potential_file"] = f'"{name}.dat"'
        values.update(overrides)
        ini = work / (name + ".ini")
        ini.write_text("\n".join(f"{k} = {v}" for k, v in values.items()) + "\n", encoding="utf-8")
        run(executable or args.executable, [ini], name, ranks=ranks)
        return values

    def compare_exact(name, potential, values):
        expected = exact_observables(np.asarray(potential), hopping=values["t_hop"],
                                    mu=values["mu"], beta=values["beta"],
                                    nmax=values["nmax"], interaction=values["U_on"],
                                    canonical=values["canonical"], v_nn=values["V_nn"])
        row = {"case": name, "observables": {}}
        with h5py.File(work / (name + ".out.h5")) as data:
            np.testing.assert_array_equal(data["hamiltonian/onsite_potential"][()], potential)
            results = data["simulation/results"]
            for observable, exact in expected.items():
                result = results[observable]
                mean = result["mean/value"][()]
                error = result["mean/error"][()]
                assert np.all(np.isfinite(mean)) and np.all(np.isfinite(error))
                # Conservative regression threshold, not a precision benchmark.
                tolerance = 7 * np.asarray(error) + (0.06 if observable == "Total_Energy" else 0.025)
                if np.any(np.abs(mean - exact) > tolerance):
                    raise AssertionError(f"{name}/{observable}: {mean} +/- {error}; exact {exact}")
                row["observables"][observable] = {"mean": np.asarray(mean).tolist(),
                    "error": np.asarray(error).tolist(), "exact": np.asarray(exact).tolist()}
            density_sum = np.sum(results["Density_Distribution/mean/value"][()])
            mean_n = results["Number_of_particles/mean/value"][()]
            np.testing.assert_allclose(density_sum, mean_n, atol=1e-10)
        report.append(row)
        print(f"PASS {name}: density, number and energy agree with exact diagonalization", flush=True)

    def assert_disabled(name, *, checkpoint=False, extension=""):
        filename = name + (".clone.h5" if checkpoint else ".out.h5") + extension
        prefix = ("simulation/realizations/0/clones/0/checkpoint/hamiltonian"
                  if checkpoint else "hamiltonian")
        with h5py.File(work / filename) as data:
            np.testing.assert_array_equal(data[f"{prefix}/onsite_potential"][()], np.zeros(4))
            assert data[f"{prefix}/onsite_potential_enabled"][()] == 0

    def finish():
        (work / "report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        print(f"All checks passed. Evidence: {work}")

    potential = np.array([-0.4, 0.2, 0.7, 1.1])
    # The off switch must work both with an existing nonzero map and when the
    # configured source cannot be opened. Both realize the zero-potential H.
    missing_source = '"does-not-exist-disabled.dat"'
    for name, profile, changes in [
        ("disabled_existing", potential, {}),
        ("disabled_missing", None, {"onsite_potential_file": missing_source}),
    ]:
        values = create_case(name, profile, onsite_potential_enabled=0, **changes)
        compare_exact(name, np.zeros(4), values)
        assert_disabled(name)
        assert_disabled(name, checkpoint=True)
    assert not (work / "does-not-exist-disabled.dat").exists()

    (work / "disabled_existing.dat").rename(work / "disabled_existing.removed.dat")
    run(args.executable, ["disabled_existing.clone.h5", "--reset-statistics"], "disabled_restore_missing")
    assert_disabled("disabled_existing")
    assert_disabled("disabled_existing", checkpoint=True)
    for suffix, override in [
        ("enable", "--onsite_potential_enabled=1"),
        ("disable_again", "--onsite_potential_enabled=0"),
        ("file", "--onsite_potential_file=disabled_missing.dat"),
    ]:
        log = run(args.executable, ["disabled_existing.clone.h5", override],
                  "disabled_reject_restore_" + suffix, expected_success=False)
        key = "onsite_potential_file" if suffix == "file" else "onsite_potential_enabled"
        assert f"Do not override {key} on restore" in log, log
    print("PASS disabled potential survives missing-source restart and rejects restore overrides", flush=True)

    if args.canonical_executable:
        values = create_case("canonical_disabled", potential, executable=args.canonical_executable,
                             canonical=2, canonical_window=0.5, onsite_potential_enabled=0)
        compare_exact("canonical_disabled", np.zeros(4), values)
        assert_disabled("canonical_disabled")
        print("PASS canonical disabled potential gives the zero-map fixed-N Hamiltonian", flush=True)

    if args.uniform_executable:
        values = create_case("uniform_disabled", None, executable=args.uniform_executable,
                             onsite_potential_enabled=0, onsite_potential_file=missing_source)
        compare_exact("uniform_disabled", np.zeros(4), values)
        assert_disabled("uniform_disabled")
        print("PASS uniform build ignores a configured missing potential when disabled", flush=True)

    for executable, name, canonical in [
        (args.mpi_executable, "mpi_disabled", -1),
        (args.canonical_mpi_executable, "canonical_mpi_disabled", 2),
    ]:
        if executable:
            create_case(name, None, executable=executable, ranks=2,
                        sweeps=10000, thermalization=1000, canonical=canonical,
                        canonical_window=0.5, Tmin=1, Tmax=1,
                        onsite_potential_enabled=0, onsite_potential_file=missing_source)
            assert_disabled(name)
            for extension in ["", ".1"]:
                assert_disabled(name, checkpoint=True, extension=extension)
            print(f"PASS {name}: two ranks skip the missing file and store a disabled zero map", flush=True)

    if args.only_disabled:
        finish()
        return

    for name, profile, changes in [
        ("nonuniform", potential, {}),
        ("atomic", potential, {"t_hop": 0.0, "nmax": 3}),
        ("zero", np.zeros(4), {}),
        ("constant", np.full(4, 0.5), {}),
        ("shifted_mu", np.zeros(4), {"mu": -0.25}),
        ("interacting", potential, {"nmax": 2, "V_nn": 0.2}),
    ]:
        values = create_case(name, profile, **changes)
        compare_exact(name, profile, values)

    # Both absent and subsequently changed source files must have no effect.
    original = work / "nonuniform.dat"
    original.rename(work / "nonuniform.removed.dat")
    for suffix in ["missing", "changed"]:
        if suffix == "changed":
            original.write_text("999\n" * 4, encoding="utf-8")
        run(args.executable, ["nonuniform.clone.h5", "--reset-statistics"], "restore_" + suffix)
        with h5py.File(work / "nonuniform.out.h5") as data:
            np.testing.assert_array_equal(data["hamiltonian/onsite_potential"][()], potential)
    print("PASS checkpoint restart with missing and changed source maps", flush=True)

    for suffix, override in [
        ("equals", ["--onsite_potential_file=zero.dat"]),
        ("separate", ["--onsite_potential_file", "zero.dat"]),
    ]:
        log = run(args.executable, ["nonuniform.clone.h5", *override],
                  "reject_restore_override_" + suffix, expected_success=False)
        if suffix == "equals":
            assert "Do not override onsite_potential_file on restore" in log
        else:
            # ALPS expects --key=value; a separated filename may be rejected
            # by its INI parser before the simulator sees the restore request.
            assert ("Do not override onsite_potential_file on restore" in log
                    or "Cannot read INI file" in log)
    print("PASS checkpoint map override rejection", flush=True)
    checkpoint = work / "nonuniform.clone.h5"
    pristine_checkpoint = checkpoint.read_bytes()
    stored_hamiltonian = "simulation/realizations/0/clones/0/checkpoint/hamiltonian"
    try:
        for key, corrupt_value in [("schema_version", 999), ("mu", 99.0)]:
            checkpoint.write_bytes(pristine_checkpoint)
            with h5py.File(checkpoint, "r+") as data:
                data[f"{stored_hamiltonian}/{key}"][...] = corrupt_value
            log = run(args.executable, [checkpoint], "reject_checkpoint_" + key,
                      expected_success=False)
            assert "mismatch" in log or "does not match" in log
    finally:
        checkpoint.write_bytes(pristine_checkpoint)
    print("PASS checkpoint metadata mismatch rejection", flush=True)

    for name, content in [("short", "0\n"), ("long", "0\n" * 5),
                          ("nonfinite", "0 0 nan 0"), ("token", "0 0 broken 0")]:
        bad = work / (name + ".dat")
        bad.write_text(content, encoding="utf-8")
        log = run(args.executable, ["zero.ini", f"--onsite_potential_file={bad.as_posix()}"],
                  "reject_" + name, expected_success=False)
        diagnostic = {"short": "values; expected", "long": "Too many values",
                      "nonfinite": "finite real numbers", "token": "Invalid potential value"}[name]
        assert diagnostic in log, log
    print("PASS malformed map rejection", flush=True)

    # Rectangular periodic geometry tests mapping and winding availability.
    periodic = np.linspace(0.0, 0.3, 12)
    create_case("periodic", periodic, Lx=4, Ly=3, pbcx=1, pbcy=1,
                sweeps=10000, thermalization=1000)
    with h5py.File(work / "periodic.out.h5") as data:
        np.testing.assert_array_equal(data["hamiltonian/onsite_potential"][()], periodic)
        assert data["simulation/results/Winding_number_squared/count"][()] > 0
        assert len(data["simulation/results/Density_Distribution/mean/value"]) == 12
    print("PASS rectangular periodic map and winding output", flush=True)

    if args.canonical_executable:
        for name, profile in [("canonical", potential), ("canonical_constant", np.full(4, 0.5))]:
            values = create_case(name, profile, executable=args.canonical_executable,
                                 canonical=2, canonical_window=0.5)
            compare_exact(name, profile, values)
        (work / "canonical.dat").rename(work / "canonical.removed.dat")
        run(args.canonical_executable, ["canonical.clone.h5", "--reset-statistics"], "canonical_restore")
        with h5py.File(work / "canonical.out.h5") as data:
            np.testing.assert_array_equal(data["hamiltonian/onsite_potential"][()], potential)
            assert data["simulation/results/Number_of_particles/mean/value"][()] == 2
        print("PASS canonical restart preserves map and exact N", flush=True)

    for executable, name, canonical in [(args.mpi_executable, "mpi", -1),
                                        (args.canonical_mpi_executable, "canonical_mpi", 2)]:
        if executable:
            create_case(name, potential, executable=executable, ranks=2,
                        sweeps=10000, thermalization=1000, canonical=canonical,
                        canonical_window=0.5, Tmin=1, Tmax=1)
            with h5py.File(work / (name + ".out.h5")) as data:
                np.testing.assert_array_equal(data["hamiltonian/onsite_potential"][()], potential)
                assert data["simulation/results/Density_Distribution/count"][()] > 0
            (work / (name + ".dat")).rename(work / (name + ".removed.dat"))
            run(executable, [name + ".clone.h5", "--reset-statistics"], name + "_restore", ranks=2)
            for extension in ["", ".1"]:
                with h5py.File(work / (name + ".clone.h5" + extension)) as data:
                    np.testing.assert_array_equal(data["simulation/realizations/0/clones/0/checkpoint/hamiltonian/onsite_potential"][()], potential)
            print(f"PASS {name} fresh and restart, two ranks, identical map", flush=True)
            rank_checkpoint = work / (name + ".clone.h5.1")
            with h5py.File(rank_checkpoint, "r+") as data:
                stored = data[f"{stored_hamiltonian}/onsite_potential"]
                stored[0] = float(stored[0]) + 0.1
            log = run(executable, [name + ".clone.h5"], name + "_reject_mixed_maps",
                      ranks=2, expected_success=False)
            assert "MPI ranks resolved different on-site potentials" in log
            print(f"PASS {name} rejects inconsistent rank maps", flush=True)

    if args.uniform_executable:
        create_case("uniform", None, executable=args.uniform_executable)
        run(args.uniform_executable, ["zero.ini"], "uniform_reject_map", expected_success=False)
        with h5py.File(work / "uniform.out.h5") as data:
            for observable, exact in exact_observables(np.zeros(4)).items():
                mean = data[f"simulation/results/{observable}/mean/value"][()]
                err = data[f"simulation/results/{observable}/mean/error"][()]
                assert np.all(np.abs(mean-exact) < 7*np.asarray(err) + 0.06)
        print("PASS uniform build baseline and rejection of potential files", flush=True)

    finish()


if __name__ == "__main__":
    main()
