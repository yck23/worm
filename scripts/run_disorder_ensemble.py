#!/usr/bin/env python3
"""Prepare and run independent quenched-disorder realizations.

Each worker launches one serial square-lattice executable with its own fresh
configuration, speckle map, Monte Carlo seed, output and checkpoint. Workers
limit simultaneous processes; this does not pin processes to particular cores.
Physical parameters and run lengths come from the supplied flat INI. This
launcher deliberately sets onsite_potential_enabled=1 and replaces the map,
seed, outputfile and checkpoint for every realization. It does not pool samples.
Use --prepare-only to generate reviewable inputs without starting simulations.
"""
from __future__ import annotations

import argparse
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from datetime import datetime, timezone
import hashlib
import json
import operator
import os
from pathlib import Path
import re
import subprocess
from typing import Any

import numpy as np

try:
    from .generate_speckle import CORRELATION_CONVENTION, generate_speckle, write_potential
except ImportError:
    from generate_speckle import CORRELATION_CONVENTION, generate_speckle, write_potential


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _without_comment(line: str) -> str:
    quote = None
    escaped = False
    for index, char in enumerate(line):
        if escaped:
            escaped = False
        elif char == "\\" and quote is not None:
            escaped = True
        elif quote is not None:
            if char == quote:
                quote = None
        elif char in ('"', "'"):
            quote = char
        elif char in ("#", ";"):
            return line[:index].strip()
    if quote is not None:
        raise ValueError("unterminated quoted INI value")
    return line.strip()


def read_ini(path: str | Path) -> dict[str, str]:
    """Read flat key=value assignments, preserving RHS syntax and physical values.

    Blank lines and #/; comments are accepted; duplicate keys and sectioned INI
    files are rejected to avoid silently changing ambiguous simulation inputs.
    """
    values: dict[str, str] = {}
    for number, raw in enumerate(Path(path).read_text(encoding="utf-8-sig").splitlines(), 1):
        line = _without_comment(raw)
        if not line:
            continue
        match = re.fullmatch(r"([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(.+)", line)
        if not match:
            raise ValueError(f"{path}:{number}: expected a flat key=value assignment")
        key, value = match.group(1), match.group(2).strip()
        if key in values:
            raise ValueError(f"{path}:{number}: duplicate parameter {key}")
        values[key] = value
    return values


def _integer(name: str, value: Any, minimum: int = 0) -> int:
    try:
        result = operator.index(value)
    except TypeError as exc:
        raise ValueError(f"{name} must be an integer >= {minimum}") from exc
    if isinstance(value, (bool, np.bool_)) or result < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return result


def _ini_integer(parameters: dict[str, str], key: str, default: int | None = None) -> int:
    if key not in parameters:
        if default is not None:
            return default
        raise ValueError(f"Parameter file must explicitly define {key}")
    value = parameters[key]
    if not re.fullmatch(r"[+-]?\d+", value):
        raise ValueError(f"INI parameter {key} must be an integer")
    return int(value)


def _positive(name: str, value: float, *, zero: bool = False) -> float:
    value = float(value)
    if not np.isfinite(value) or (value < 0 if zero else value <= 0):
        raise ValueError(f"{name} must be finite and {'nonnegative' if zero else 'positive'}")
    return value


def _write_json(path: Path, value: dict) -> None:
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(path)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare_ensemble(
    parameter_file: str | Path,
    output_dir: str | Path,
    *,
    realizations: int,
    workers: int,
    correlation_length: float,
    mean_strength: float = 1.0,
    dx: float = 1.0,
    dy: float | None = None,
    base_disorder_seed: int = 42,
    base_mc_seed: int = 1000042,
    canonical: bool = False,
    executable: str | Path | None = None,
    subtract_mean: bool = False,
    prepare_only: bool = False,
) -> Path:
    """Prepare all fresh inputs and return manifest.json; never starts a process.

    prepare_only permits preparing inputs before an executable has been built.
    Per-realization seeds are base+index, independent of worker scheduling.
    potential_sha256 hashes contiguous little-endian float64 potential values;
    potential_file_sha256 additionally protects the exact text file and metadata.
    An existing empty output directory is accepted; populated ones are refused.
    """
    source = Path(parameter_file).resolve(strict=True)
    parameters = read_ini(source)
    if parameters.get("model", "").strip('"\'') != "BoseHubbard":
        raise ValueError('Parameter file must explicitly set model = "BoseHubbard"')
    nx, ny = _ini_integer(parameters, "Lx"), _ini_integer(parameters, "Ly")
    if nx < 2 or ny < 2 or _ini_integer(parameters, "Lz", 1) != 1:
        raise ValueError("Speckle ensembles require a square lattice with Lx,Ly >= 2 and Lz=1")
    particle_number = _ini_integer(parameters, "canonical", -1)
    if particle_number < -1:
        raise ValueError("canonical must be -1 (grand canonical) or a nonnegative integer N")
    if bool(canonical) != (particle_number >= 0):
        raise ValueError("--canonical must match canonical=N>=0 in the INI; omit it for canonical=-1")
    realizations = _integer("realizations", realizations, 1)
    workers = _integer("workers", workers, 1)
    base_disorder_seed = _integer("base_disorder_seed", base_disorder_seed)
    base_mc_seed = _integer("base_mc_seed", base_mc_seed)
    if base_mc_seed + realizations - 1 > 2**32 - 1:
        raise ValueError("Monte Carlo seeds must fit unsigned 32-bit integers")
    ell = _positive("correlation_length", correlation_length)
    strength = _positive("mean_strength", mean_strength, zero=True)
    dx = _positive("dx", dx)
    dy = dx if dy is None else _positive("dy", dy)
    if executable is None:
        suffix = "_canonical" if canonical else ""
        extension = ".exe" if os.name == "nt" else ""
        executable = Path(__file__).resolve().parents[1] / "dist" / "windows-square-programmable" / f"qmc_worm{suffix}{extension}"
    executable = Path(executable).resolve()
    if "_mpi" in executable.stem.lower():
        raise ValueError("Use a serial executable for independent disorder realizations, not an MPI executable")
    if not prepare_only and not executable.is_file():
        raise ValueError(f"Serial executable not found: {executable}; build it or use --prepare-only")
    root = Path(output_dir).resolve()
    if root.exists() and (not root.is_dir() or any(root.iterdir())):
        raise ValueError(f"Refusing to overwrite populated ensemble output directory: {root}")
    root.mkdir(parents=True, exist_ok=True)
    (root / "source_parameters.ini").write_bytes(source.read_bytes())
    manifest: dict[str, Any] = {
        "schema_version": 1,
        "kind": "independent_disorder_ensemble",
        "created_utc": _utc_now(),
        "status": "preparing",
        "parameter_source": str(source),
        "parameter_source_sha256": _sha256(source),
        "executable": str(executable),
        "canonical": bool(canonical),
        "workers": workers,
        "realization_count": realizations,
        "shared_parameters": parameters,
        "overridden_parameters": ["seed", "onsite_potential_enabled", "onsite_potential_file", "outputfile", "checkpoint"],
        "disorder": {"nx": nx, "ny": ny, "dx": dx, "dy": dy, "correlation_length": ell,
                     "mean_strength": strength, "subtract_mean": bool(subtract_mean),
                     "base_disorder_seed": base_disorder_seed, "base_mc_seed": base_mc_seed,
                     "correlation_convention": CORRELATION_CONVENTION},
        "runs": [],
    }
    manifest_path = root / "manifest.json"
    _write_json(manifest_path, manifest)
    try:
        for index in range(realizations):
            name = f"r{index:04d}"
            case = root / name
            case.mkdir()
            (case / "tmp").mkdir()
            disorder_seed, mc_seed = base_disorder_seed + index, base_mc_seed + index
            potential = generate_speckle(nx, ny, dx=dx, dy=dy, correlation_length=ell,
                                         mean_strength=strength, seed=disorder_seed, subtract_mean=subtract_mean)
            metadata = dict(manifest["disorder"], seed=disorder_seed, mc_seed=mc_seed,
                            numpy_version=np.__version__, mean=float(potential.mean()),
                            standard_deviation=float(potential.std()), minimum=float(potential.min()),
                            maximum=float(potential.max()), ordering="site=x+nx*y; potential[y,x]")
            write_potential(case / "potential.dat", potential, metadata=metadata)
            _write_json(case / "potential.json", metadata)
            values = dict(parameters, seed=str(mc_seed), onsite_potential_enabled="1",
                          onsite_potential_file='"' + (case / "potential.dat").as_posix() + '"',
                          outputfile='"result.out.h5"', checkpoint='"state.clone.h5"')
            (case / "parameters.ini").write_text(
                "# Independent disorder realization; physical parameters preserved from source.\n"
                + "\n".join(f"{key} = {value}" for key, value in values.items()) + "\n", encoding="utf-8")
            (case / "run.log").touch()
            manifest["runs"].append({
                "id": name, "index": index, "disorder_seed": disorder_seed, "mc_seed": mc_seed,
                "directory": name, "input_file": f"{name}/parameters.ini",
                "potential_file": f"{name}/potential.dat", "potential_metadata_file": f"{name}/potential.json",
                "output_file": f"{name}/result.out.h5", "checkpoint_file": f"{name}/state.clone.h5",
                "log_file": f"{name}/run.log",
                "potential_sha256": hashlib.sha256(np.ascontiguousarray(potential, dtype="<f8").tobytes()).hexdigest(),
                "potential_file_sha256": _sha256(case / "potential.dat"),
                "input_sha256": _sha256(case / "parameters.ini"),
                "status": "prepared", "returncode": None,
            })
            _write_json(manifest_path, manifest)
        manifest["status"] = "prepared"
        _write_json(manifest_path, manifest)
    except Exception as exc:
        manifest["status"] = "failed"
        manifest["error"] = f"Preparation failed: {exc}"
        _write_json(manifest_path, manifest)
        raise
    return manifest_path


def _run_realization(root: Path, run: dict, executable: Path) -> dict:
    case = root / run["directory"]
    env = os.environ.copy()
    env["PATH"] = str(executable.parent) + os.pathsep + env.get("PATH", "")
    env["TMPDIR"] = (case / "tmp").as_posix() + "/"
    env["OMP_NUM_THREADS"] = "1"
    result: dict[str, Any] = {"started_utc": _utc_now(), "returncode": None}
    try:
        with (root / run["log_file"]).open("w", encoding="utf-8") as log:
            completed = subprocess.run(
                [str(executable), str(root / run["input_file"])], cwd=case, env=env,
                stdout=log, stderr=subprocess.STDOUT, check=False,
                creationflags=subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0,
            )
        result["returncode"] = completed.returncode
        if completed.returncode != 0:
            raise RuntimeError(f"Simulation exited with code {completed.returncode}; see {run['log_file']}")
        if not (root / run["output_file"]).is_file():
            raise RuntimeError("Simulation exited successfully but did not produce result.out.h5")
        result["status"] = "completed"
    except Exception as exc:
        result.update(status="failed", error=str(exc))
    result["finished_utc"] = _utc_now()
    return result


def run_ensemble(manifest_path: str | Path) -> int:
    """Run prepared cases with bounded concurrency; stop queuing after a failure.

    Completed means the serial process succeeded and produced its output file;
    the separate summarizer checks HDF5 measurement counts and finite estimates.
    All manifest writes occur on this parent thread. Already launched jobs finish
    if another fails; unlaunched jobs become skipped. Failed jobs are never pooled.
    """
    path = Path(manifest_path).resolve(strict=True)
    root = path.parent
    manifest = json.loads(path.read_text(encoding="utf-8"))
    if manifest.get("schema_version") != 1 or manifest.get("kind") != "independent_disorder_ensemble":
        raise ValueError("Unsupported ensemble manifest")
    runs = manifest["runs"]
    if manifest.get("status") != "prepared" or len(runs) != manifest["realization_count"] or any(r["status"] != "prepared" for r in runs):
        raise ValueError("Only a fully prepared, never-started ensemble can be run; use a fresh directory")
    executable = Path(manifest["executable"])
    if not executable.is_file():
        raise ValueError(f"Serial executable not found: {executable}")
    workers = _integer("workers", manifest["workers"], 1)
    # Refuse stale outputs or changed inputs before starting any process.
    for run in runs:
        for key in ["output_file", "checkpoint_file"]:
            if (root / run[key]).exists():
                raise ValueError(f"Refusing to overwrite existing run artifact: {root / run[key]}")
        if _sha256(root / run["potential_file"]) != run["potential_file_sha256"] or _sha256(root / run["input_file"]) != run["input_sha256"]:
            raise ValueError(f"Prepared inputs changed for {run['id']}; prepare a fresh ensemble")
    manifest.update(status="running", started_utc=_utc_now())
    _write_json(path, manifest)
    next_index, failed = 0, False
    with ThreadPoolExecutor(max_workers=min(workers, len(runs))) as pool:
        active = {}
        while active or (next_index < len(runs) and not failed):
            while not failed and next_index < len(runs) and len(active) < workers:
                run = runs[next_index]
                run["status"] = "running"
                active[pool.submit(_run_realization, root, dict(run), executable)] = run
                next_index += 1
            _write_json(path, manifest)
            completed, _ = wait(active, return_when=FIRST_COMPLETED)
            for future in completed:
                run = active.pop(future)
                try:
                    run.update(future.result())
                except Exception as exc:
                    run.update(status="failed", error=str(exc), finished_utc=_utc_now())
                failed = failed or run["status"] == "failed"
                print(f"{run['id']}: {run['status']}", flush=True)
            _write_json(path, manifest)
    for run in runs[next_index:]:
        run.update(status="skipped", error="An earlier realization failed; not launched")
    manifest.update(status="failed" if failed else "completed", finished_utc=_utc_now())
    _write_json(path, manifest)
    return 1 if failed else 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parameter-file", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True, help="New or empty ensemble directory")
    parser.add_argument("--realizations", type=int, required=True, help="Number of independent disorder realizations")
    parser.add_argument("--workers", type=int, default=1, help="Maximum simultaneous serial simulations; default: 1")
    parser.add_argument("--correlation-length", type=float, required=True)
    parser.add_argument("--mean-strength", type=float, default=1.0)
    parser.add_argument("--dx", type=float, default=1.0)
    parser.add_argument("--dy", type=float)
    parser.add_argument("--base-disorder-seed", type=int, default=42)
    parser.add_argument("--base-mc-seed", type=int, default=1000042)
    parser.add_argument("--subtract-mean", action="store_true")
    parser.add_argument("--canonical", action="store_true", help="Use the canonical executable; requires canonical=N>=0 in the INI")
    parser.add_argument("--executable", type=Path, help="Override packaged serial square programmable executable")
    parser.add_argument("--prepare-only", action="store_true", help="Write all inputs and manifest without launching simulations")
    args = parser.parse_args(argv)
    try:
        manifest = prepare_ensemble(**vars(args))
        print(f"Prepared ensemble: {manifest}", flush=True)
        if args.prepare_only:
            return 0
        status = run_ensemble(manifest)
        print(f"Ensemble {'completed' if status == 0 else 'failed'}; manifest: {manifest}")
        return status
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
