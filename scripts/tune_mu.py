#!/usr/bin/env python3
"""Tune grand-canonical mu to a requested mean N using the packaged serial solver.

Uses discarded burn-in, separate production blocks, independent low/high starts,
and fresh confirmation. This does not measure finite-N F(N+1)-F(N).
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys

import h5py
import numpy as np

try:
    from .run_disorder_ensemble import read_ini
except ImportError:
    from run_disorder_ensemble import read_ini

REPO = Path(__file__).resolve().parents[1]


@dataclass
class Block:
    mean: float
    error: float
    variance: float
    count: int
    effective_samples: float


class InsufficientSampling(ValueError):
    """Longer fresh evaluation may fix this measurement failure."""


def scalar(handle, key):
    return float(np.asarray(handle[key][()]).reshape(-1)[0])


def read_block(path, mu, beta):
    with h5py.File(path, "r") as data:
        if scalar(data, "parameters/canonical") != -1:
            raise ValueError(f"Expected grand-canonical output: {path}")
        for key, expected in (("mu", mu), ("beta", beta)):
            if not math.isclose(scalar(data, "parameters/" + key), expected,
                                rel_tol=1e-12, abs_tol=1e-12):
                raise ValueError(f"Output {key} differs from requested value: {path}")
        if "simulation/results/Number_of_particles" not in data:
            raise InsufficientSampling("No measurements; increase run duration.")
        group = data["simulation/results"]
        mean = scalar(group, "Number_of_particles/mean/value")
        error = scalar(group, "Number_of_particles/mean/error")
        second = scalar(group, "Number_of_particles_squared/mean/value")
        count = int(scalar(group, "Number_of_particles/count"))
        if count != int(scalar(group, "Number_of_particles_squared/count")):
            raise ValueError("N and N^2 have different measurement counts.")
    variance = second - mean * mean
    if not all(math.isfinite(x) for x in (mean, error, variance)) or error <= 0:
        raise InsufficientSampling("Missing/zero number error; number mixing is unresolved.")
    if count < 128 or variance <= 1e-12:
        raise InsufficientSampling("Too few samples or no resolved number fluctuations.")
    effective = min(float(count), variance / error**2)
    if effective < 50:
        raise InsufficientSampling("Fewer than 50 estimated independent number samples per block.")
    return Block(mean, error, variance, count, effective)


def summarize(chains, beta, target, tolerance):
    """Conservative errors, explicit checks for drift and start disagreement."""
    means, errors, variances = [], [], []
    stationary = True
    for blocks in chains:
        values = np.array([b.mean for b in blocks])
        # Adjacent blocks can correlate: do not reduce their errors by sqrt(B).
        se = max(float(np.mean([b.error for b in blocks])),
                 float(np.std(values, ddof=1) / math.sqrt(len(blocks))))
        means.append(float(values.mean()))
        errors.append(se)
        # Exclude between-block/start drift from the fluctuation response.
        variances.append(float(np.mean([b.variance for b in blocks])))
        stationary &= abs(blocks[-1].mean - blocks[0].mean) <= 3 * math.hypot(
            blocks[-1].error, blocks[0].error)
    mean = float(np.mean(means))
    error = max(math.hypot(*errors) / 2, abs(means[1] - means[0]) / 2)
    starts_agree = abs(means[1] - means[0]) <= 3 * math.hypot(*errors)
    response = beta * float(np.mean(variances))
    if response <= 0 or not math.isfinite(response):
        raise InsufficientSampling("No resolved positive number response.")
    residual = mean - target
    passed = bool(stationary and starts_agree and abs(residual) + 2 * error <= tolerance)
    return dict(mean_number=mean, number_se=error, residual=residual,
                variance=float(np.mean(variances)), response=response,
                mu_se=error / response, number_error_bound=abs(residual) + 2 * error,
                stationary=bool(stationary), starts_agree=bool(starts_agree),
                passed=passed, chain_means=means, chain_errors=errors,
                samples=sum(b.count for chain in chains for b in chain),
                min_effective_samples=min(b.effective_samples for chain in chains for b in chain))


def next_mu(mu, estimate, low, high, step):
    """Safeguarded response step, geometric search until a bracket exists."""
    residual, response = estimate["residual"], estimate["response"]
    direction = -1 if residual > 0 else 1
    correction = -residual / response if response > 1e-12 else direction * step
    if low is not None and high is not None:
        width = high - low
        return min(high - .1 * width, max(low + .1 * width, mu + correction))
    return mu + direction * min(step, max(.1 * step, abs(correction)))


def unquote(value):
    return value.strip().strip('"').strip("'")


def write_ini(path, values):
    Path(path).write_text("\n".join(f"{k} = {v}" for k, v in values.items()) + "\n",
                          encoding="utf-8")


def write_json(path, value):
    path = Path(path)
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(path)


def file_hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def settings(args):
    source = args.parameter_file.resolve(strict=True)
    params = read_ini(source)
    if unquote(params.get("model", '"BoseHubbard"')) != "BoseHubbard":
        raise ValueError("This helper supports the square-lattice BoseHubbard model.")
    lx, ly = int(params["Lx"]), int(params["Ly"])
    if min(lx, ly) < 3 or int(params.get("Lz", "1")) != 1:
        raise ValueError("Use a 2D square lattice with Lx, Ly >= 3 and Lz = 1.")
    capacity = lx * ly * int(params.get("nmax", "10"))
    target = args.target if args.target is not None else int(params.get("canonical", "-1"))
    if not math.isfinite(target) or not 0 < target < capacity:
        raise ValueError("Supply --target N (or canonical=N in the INI), strictly between 0 and sites*nmax.")
    beta = float(params["beta"])
    mu = args.mu_start if args.mu_start is not None else float(params.get("mu", ".7"))
    tolerance = args.number_tolerance if args.number_tolerance is not None else max(.1, .005 * target)
    step = args.mu_step if args.mu_step is not None else (max(.25, 1 / beta) if beta > 0 else 1)
    for name, val in (("beta", beta), ("number tolerance", tolerance), ("mu step", step)):
        if not math.isfinite(val) or val <= 0:
            raise ValueError(f"{name} must be finite and positive.")
    if not math.isfinite(mu):
        raise ValueError("Starting mu must be finite.")
    if args.seconds < 1 or args.max_seconds < args.seconds or args.max_iterations < 2:
        raise ValueError("Require seconds >= 1, max-seconds >= seconds, max-iterations >= 2.")
    if args.blocks < 2 or args.burn_in_blocks < 1:
        raise ValueError("Require at least two production blocks and one discarded burn-in block.")
    if not 10 <= args.sweeps <= 2000000000 or not 0 <= args.seed < 2000000000:
        raise ValueError("Invalid sweep limit or seed.")
    if args.thermalization is not None:
        params["thermalization"] = str(args.thermalization)
    if int(params.get("thermalization", "10000")) < 0:
        raise ValueError("thermalization must be nonnegative.")
    enabled = unquote(params.get("onsite_potential_enabled", "1")).lower()
    if enabled not in ("0", "1", "false", "true"):
        raise ValueError("onsite_potential_enabled must be 0, 1, false, or true.")
    map_name = unquote(params.get("onsite_potential_file", '""'))
    potential = None
    if enabled in ("1", "true") and map_name:
        if args.potential_file:
            potential = args.potential_file.resolve(strict=True)
        else:
            candidate = Path(map_name)
            if candidate.is_absolute():
                potential = candidate.resolve(strict=True)
            else:
                candidates = {p.resolve() for p in (Path.cwd()/candidate, source.parent/candidate,
                                                   REPO/candidate) if p.is_file()}
                if len(candidates) != 1:
                    raise ValueError("Cannot uniquely locate enabled map; supply --potential-file PATH.")
                potential = candidates.pop()
    elif args.potential_file:
        raise ValueError("--potential-file resolves an enabled INI map; enable/configure it there first.")
    package = "windows-square-programmable" if potential else "windows-square"
    exe = (args.executable or REPO/"dist"/package/"qmc_worm.exe").resolve(strict=True)
    if "canonical" in exe.stem.lower() or "mpi" in exe.stem.lower():
        raise ValueError("Use a serial grand-canonical executable, qmc_worm.exe.")
    spread = max(1, int(math.ceil(.2 * target)))
    starts = [max(0, int(target) - spread), min(capacity, int(math.ceil(target)) + spread)]
    return params, dict(source=str(source), executable=str(exe), potential=str(potential) if potential else None,
                        target=target, beta=beta, mu_start=mu, tolerance=tolerance, mu_step=step,
                        starts=starts, seconds=args.seconds, max_seconds=args.max_seconds,
                        max_iterations=args.max_iterations, blocks=args.blocks, burn_in_blocks=args.burn_in_blocks,
                        sweeps=args.sweeps, seed=args.seed)


def launch(exe, arguments, folder, seconds, label):
    env = os.environ.copy()
    env["PATH"] = str(exe.parent) + os.pathsep + env.get("PATH", "")
    env["TMPDIR"] = folder.as_posix() + "/"
    print(f"  {folder.name}: {label}", flush=True)
    checkpoint = folder / "state.h5"
    previous_checkpoint = file_hash(checkpoint) if checkpoint.is_file() else None
    log_path = folder / f"{label}.log"
    with log_path.open("w", encoding="utf-8") as log:
        # The solver checks its time limit between complete worm updates.
        completed = subprocess.run([str(exe), *arguments], cwd=folder, env=env,
                                   stdout=log, stderr=subprocess.STDOUT,
                                   timeout=max(120, 4 * seconds),
                                   creationflags=subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0)
    if completed.returncode:
        raise RuntimeError(f"Solver failed ({completed.returncode}); see {folder / (label + '.log')}")
    if "WARNING: Failed to save checkpoint" in log_path.read_text(encoding="utf-8", errors="replace"):
        raise RuntimeError(f"Solver could not save a fresh checkpoint; see {log_path}")
    if not checkpoint.is_file() or not (folder/"result.h5").is_file():
        raise RuntimeError(f"Solver did not save its checkpoint/result in {folder}")
    if previous_checkpoint is not None and file_hash(checkpoint) == previous_checkpoint:
        raise RuntimeError(f"Solver left a stale checkpoint; see {log_path}")


def evaluate(root, params, config, iteration, mu, seconds):
    chains = []
    trial = root / f"trial_{iteration:02d}"
    trial.mkdir()
    for chain, start in enumerate(config["starts"]):
        folder = trial / f"chain_{chain + 1}"
        folder.mkdir()
        seed = (config["seed"] + 1009 * iteration + 97 * chain) % 2000000000
        values = dict(params, canonical="-1", initial_particle_number=str(start),
                      mu=format(mu, ".17g"), seed=str(seed), runtimelimit=str(seconds),
                      sweeps=str(config["sweeps"]), outputfile='"result.h5"', checkpoint='"state.h5"')
        values.pop("SEED", None)
        write_ini(folder/"parameters.ini", values)
        exe = Path(config["executable"])
        for burn in range(config["burn_in_blocks"]):
            arguments = ["parameters.ini"] if burn == 0 else ["state.h5"]
            launch(exe, arguments, folder, seconds, f"burn_in_{burn + 1:02d}")
            shutil.copyfile(folder/"result.h5", folder/f"burn_in_{burn + 1:02d}.h5")
        with h5py.File(folder/"result.h5", "r") as data:
            if "simulation/results/Number_of_particles" not in data:
                raise InsufficientSampling("Burn-in ended before thermalization completed; increase seconds.")
        blocks = []
        for block in range(config["blocks"]):
            label = f"production_{block + 1:02d}"
            launch(exe, ["state.h5", "--reset-statistics"], folder, seconds, label)
            saved = folder/f"{label}.h5"
            shutil.copyfile(folder/"result.h5", saved)
            blocks.append(read_block(saved, mu, config["beta"]))
        chains.append(blocks)
    write_json(trial/"blocks.json", [[asdict(b) for b in chain] for chain in chains])
    return summarize(chains, config["beta"], config["target"], config["tolerance"])


def save_history(root, state):
    write_json(root/"result.json", state)
    columns = ["iteration", "phase", "mu", "seconds", "mean_number", "number_se", "residual",
               "mu_se", "response", "stationary", "starts_agree", "passed", "issue"]
    with (root/"history.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(state["history"])


def tune(args, evaluator=evaluate):
    params, config = settings(args)
    if args.dry_run:
        print(json.dumps(config, indent=2))
        return 0
    name = "mu_" + datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S_%f")
    root = (args.output_dir or REPO/"results"/name).resolve()
    if root.exists() and any(root.iterdir()):
        raise ValueError(f"Choose an empty output directory: {root}")
    root.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(config["source"], root/"source.ini")
    config["source_sha256"] = file_hash(root/"source.ini")
    config["executable_sha256"] = file_hash(config["executable"])
    if config["potential"]:
        shutil.copyfile(config["potential"], root/"potential.dat")
        config["potential_sha256"] = file_hash(root/"potential.dat")
        params["onsite_potential_file"] = '"' + (root/"potential.dat").as_posix() + '"'
    state = dict(status="running", definition="grand-canonical mu at <N>=target; not F(N+1)-F(N)",
                 uncertainty="Approximate local statistical SE; includes start scatter, not systematic equilibration error.",
                 config=config, history=[])
    save_history(root, state)
    print(f"Target <N>={config['target']:g}; number tolerance={config['tolerance']:g}\nResults: {root}", flush=True)
    mu, seconds, step = config["mu_start"], config["seconds"], config["mu_step"]
    low = high = None
    confirming = False
    try:
        for iteration in range(1, config["max_iterations"] + 1):
            phase = "confirmation" if confirming else "tuning"
            print(f"Trial {iteration}: {phase}, mu={mu:.10g}, {seconds}s per block", flush=True)
            row = dict(iteration=iteration, phase=phase, mu=mu, seconds=seconds)
            try:
                estimate = evaluator(root, params, config, iteration, mu, seconds)
                row.update(estimate)
            except InsufficientSampling as exc:
                row["issue"] = str(exc)
                estimate = None
            state["history"].append(row)
            save_history(root, state)
            if estimate:
                print(f"  <N>={estimate['mean_number']:.8g} +/- {estimate['number_se']:.3g}; "
                      f"residual={estimate['residual']:+.4g}; mu SE~{estimate['mu_se']:.3g}", flush=True)
            else:
                print("  " + row["issue"], flush=True)
            if estimate and estimate["passed"]:
                if confirming:
                    state.update(status="converged", mu=mu, mu_se=estimate["mu_se"],
                                 mean_number=estimate["mean_number"], number_se=estimate["number_se"],
                                 residual=estimate["residual"],
                                 mu_correction_estimate=-estimate["residual"] / estimate["response"])
                    tuned = dict(params, mu=format(mu, ".17g"), canonical="-1",
                                 initial_particle_number=str(round(config["target"])),
                                 outputfile='"tuned.out.h5"', checkpoint='"tuned.clone.h5"')
                    write_ini(root/"tuned_parameters.ini", tuned)
                    break
                confirming = True
                continue
            confirming = False
            if (estimate is None or not estimate["stationary"] or not estimate["starts_agree"]
                    or abs(estimate["residual"]) <= 2 * estimate["number_se"]
                    or 2 * estimate["number_se"] >= config["tolerance"]):
                if seconds >= config["max_seconds"]:
                    state["reason"] = ("Sampling or precision remains insufficient at max-seconds. "
                                       "Increase --seconds/--max-seconds, --burn-in-blocks or --sweeps; inspect history.")
                    break
                seconds = min(config["max_seconds"], seconds * 2)
                continue
            if estimate["residual"] < 0:
                low = mu if low is None else max(low, mu)
            else:
                high = mu if high is None else min(high, mu)
            if low is not None and high is not None and low >= high:
                state["reason"] = "Conflicting resolved estimates; increase equilibration and sampling."
                break
            proposed = next_mu(mu, estimate, low, high, step)
            if not math.isfinite(proposed) or proposed == mu:
                state["reason"] = "No resolved mu step; inspect number plateau and sampling."
                break
            mu = proposed
            if low is None or high is None:
                step *= 2
        if state["status"] != "converged":
            state["status"] = "not_converged"
            state.setdefault("reason", "Trial limit reached before a fresh confirmation passed.")
        save_history(root, state)
    except (Exception, KeyboardInterrupt) as exc:
        state.update(status="interrupted" if isinstance(exc, KeyboardInterrupt) else "failed", reason=str(exc))
        save_history(root, state)
        raise
    if state["status"] == "converged":
        print(f"mu = {state['mu']:.10g} (approximate statistical SE {state['mu_se']:.3g})\n"
              f"Confirmed <N> = {state['mean_number']:.8g} +/- {state['number_se']:.3g}\n"
              f"Parameters: {root / 'tuned_parameters.ini'}", flush=True)
        return 0
    print(f"NOT CONVERGED: {state['reason']}\nHistory: {root / 'history.csv'}", flush=True)
    return 2


def parser():
    result = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    result.add_argument("parameter_file", type=Path, help="Existing canonical or grand-canonical INI (never modified)")
    result.add_argument("--target", type=float, help="Target mean total N; defaults to canonical=N in the INI")
    result.add_argument("--number-tolerance", type=float, help="Absolute N tolerance; default max(0.1, 0.005*target)")
    result.add_argument("--mu-start", type=float, help="Initial mu; defaults to INI mu")
    result.add_argument("--mu-step", type=float, help="Initial search step; default max(0.25, 1/beta)")
    result.add_argument("--seconds", type=int, default=10, help="Seconds per solver call (default 10)")
    result.add_argument("--max-seconds", type=int, default=120, help="Maximum adaptive seconds per call (default 120)")
    result.add_argument("--max-iterations", type=int, default=12, help="Trial budget including confirmation (default 12)")
    result.add_argument("--blocks", type=int, default=3, help="Retained production blocks per chain (default 3)")
    result.add_argument("--burn-in-blocks", type=int, default=1, help="Discarded solver calls per chain (default 1)")
    result.add_argument("--sweeps", type=int, default=1000000000, help="Sweep cap per call; wall time normally limits it")
    result.add_argument("--thermalization", type=int, help="Override INI technical thermalization threshold")
    result.add_argument("--seed", type=int, default=731031, help="Base seed; each trial/start receives a different seed")
    result.add_argument("--potential-file", type=Path, help="Resolve enabled INI map if its relative path is ambiguous")
    result.add_argument("--executable", type=Path, help="Custom serial grand-canonical solver; otherwise select packaged build")
    result.add_argument("--output-dir", type=Path, help="Empty output directory; default results/mu_TIMESTAMP")
    result.add_argument("--dry-run", action="store_true", help="Validate/print settings without running or writing")
    return result


def main():
    args = parser().parse_args()
    try:
        return tune(args)
    except KeyboardInterrupt:
        print("Interrupted; completed trials remain in the output directory.", file=sys.stderr)
        return 130
    except (ValueError, OSError, RuntimeError, KeyError, subprocess.TimeoutExpired) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
