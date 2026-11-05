#!/usr/bin/env python3
"""Generate periodic 2D speckle in numerical length and energy units.

The requested correlation length ell is the Gaussian standard deviation of the
connected intensity correlation: C(r)/C(0) ~= exp[-r**2 / (2*ell**2)]. This holds
in the resolved, many-mode continuum limit; finite periodic grids give corrections.
No experimental waist, calibration factor, hbar, or Hz conversion is used.

Example (dx=dy=1 means ell is in lattice spacings)::

    python scripts/generate_speckle.py speckle.dat --nx 64 --ny 48 \
        --correlation-length 2.5 --mean-strength 0.8 --seed 42 \
        --metadata-json speckle.json

Output maps potential[y, x] to site = x + nx*y. Use onsite_potential_file in a
square-lattice PROGRAMMABLE_POTENTIAL build to load it.
"""

from __future__ import annotations

import argparse
import json
import operator
from pathlib import Path
from typing import Any
import warnings

import numpy as np


CORRELATION_CONVENTION = (
    "connected intensity C(r)/C(0) ~= exp(-r^2/(2*correlation_length^2)); "
    "well-resolved many-mode continuum limit"
)


def _positive_float(name: str, value: float) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a finite positive number") from exc
    if not np.isfinite(result) or result <= 0:
        raise ValueError(f"{name} must be a finite positive number")
    return result


def _grid_size(name: str, value: int) -> int:
    try:
        result = operator.index(value)
    except TypeError as exc:
        raise ValueError(f"{name} must be an integer >= 2") from exc
    if isinstance(value, (bool, np.bool_)) or result < 2:
        raise ValueError(f"{name} must be an integer >= 2")
    return result


def _parameters(nx, ny, dx, dy, correlation_length, mean_strength, seed):
    nx = _grid_size("nx", nx)
    ny = nx if ny is None else _grid_size("ny", ny)
    dx = _positive_float("dx", dx)
    dy = dx if dy is None else _positive_float("dy", dy)
    correlation_length = _positive_float("correlation_length", correlation_length)
    try:
        mean_strength = float(mean_strength)
    except (TypeError, ValueError) as exc:
        raise ValueError("mean_strength must be finite and nonnegative") from exc
    if not np.isfinite(mean_strength) or mean_strength < 0:
        raise ValueError("mean_strength must be finite and nonnegative")
    if seed is not None:
        try:
            checked_seed = operator.index(seed)
        except TypeError as exc:
            raise ValueError("seed must be a nonnegative integer or None") from exc
        if isinstance(seed, (bool, np.bool_)) or checked_seed < 0:
            raise ValueError("seed must be a nonnegative integer or None")
        seed = checked_seed
    return nx, ny, dx, dy, correlation_length, mean_strength, seed


def _sampling_messages(nx, ny, dx, dy, correlation_length):
    messages = []
    if correlation_length < max(dx, dy):
        messages.append(
            "Correlation length is shorter than a grid spacing: the requested "
            "continuum correlation is poorly resolved. Reduce dx/dy or increase "
            "the correlation length."
        )
    if min(nx * dx, ny * dy) < 6 * correlation_length:
        messages.append(
            "A box dimension is shorter than six correlation lengths: periodic "
            "images and too few Fourier modes can change the correlation and "
            "reduce the speckle contrast. Increase nx/ny for bulk speckle."
        )
    return messages


def generate_speckle(
    nx: int,
    ny: int | None = None,
    *,
    correlation_length: float,
    mean_strength: float = 1.0,
    dx: float = 1.0,
    dy: float | None = None,
    seed: int | None = None,
    subtract_mean: bool = False,
) -> np.ndarray:
    """Return a float64 potential of shape (ny, nx), with x varying fastest.

    nx/ny must each be >= 2; ny defaults to nx and dy defaults to dx. Coordinates
    are x = arange(nx)*dx and y = arange(ny)*dy, with FFT periods nx*dx, ny*dy.
    correlation_length, dx, and dy use the same arbitrary numerical length unit.
    mean_strength is the nonnegative spatial mean BEFORE optional subtraction,
    in the same energy units as the simulation's hopping/chemical potential.
    With subtraction, the returned mean is zero but fluctuations are unchanged.
    The standard deviation is not independently normalized; it approaches
    mean_strength in the many-mode limit.

    Fourier phases are independent uniform [0, 2*pi), as in the experimental
    construction. For amplitude filter A(f)=exp[-2*pi**2*ell**2*|f|**2], the
    power spectrum gives field correlation exp[-r**2/(4*ell**2)]. The connected
    intensity correlation approaches its square. More precisely, for discrete
    P=|A|**2, g(r)=sum(P*exp(2*pi*i*f.r))/sum(P), and q=sum(P**2)/sum(P)**2,
    the phase-ensemble result is C(r)/C(0)=(|g(r)|**2-q)/(1-q). Thus ell is a
    target spectral parameter, not a fitted length for every finite realization.

    A local NumPy Generator honors seed without changing NumPy's global RNG.
    Each realization is normalized to the requested mean to floating-point
    precision. No experimental calibration factor is needed.
    """
    nx, ny, dx, dy, ell, mean_strength, seed = _parameters(
        nx, ny, dx, dy, correlation_length, mean_strength, seed
    )
    for message in _sampling_messages(nx, ny, dx, dy, ell):
        warnings.warn(message, UserWarning, stacklevel=2)
    if mean_strength == 0:
        return np.zeros((ny, nx), dtype=np.float64)

    fx = np.fft.fftfreq(nx, d=dx)
    fy = np.fft.fftfreq(ny, d=dy)
    # ifft2 uses unshifted frequencies. Its normalization cancels below.
    amplitude = np.exp(-2 * np.pi**2 * ell**2 * (fy[:, None]**2 + fx[None, :]**2))
    rng = np.random.default_rng(seed)
    phase = rng.uniform(0.0, 2 * np.pi, size=(ny, nx))
    field = np.fft.ifft2(amplitude * np.exp(1j * phase))
    potential = np.abs(field)**2
    potential /= np.mean(potential)
    potential *= mean_strength
    if not np.all(np.isfinite(potential)):
        raise ValueError("mean_strength is too large for finite float64 intensities")
    if subtract_mean:
        potential -= np.mean(potential)
    return potential


def write_potential(
    path: str | Path,
    potential: np.ndarray,
    *,
    metadata: dict[str, Any] | None = None,
) -> None:
    """Export any finite 2D map as a dense onsite_potential_file.

    The map is indexed [y, x]; one scalar per line is emitted in site order
    x + nx*y. Comment lines carry shape/order and optional JSON metadata.
    Seventeen significant digits preserve float64 values on round-trip.
    Missing parent directories are created.
    """
    values = np.asarray(potential)
    if np.iscomplexobj(values):
        raise ValueError("potential must be real")
    values = np.asarray(values, dtype=np.float64)
    if values.ndim != 2 or not values.size or not np.all(np.isfinite(values)):
        raise ValueError("potential must be a nonempty finite 2D array")
    ny, nx = values.shape
    header = ["Dense onsite potential; one value per site", f"nx={nx} ny={ny}; site=x+nx*y; map[y,x]"]
    if metadata is not None:
        header.append(json.dumps(metadata, sort_keys=True, allow_nan=False))
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(path, values.ravel(order="C"), fmt="%.17g", header="\n".join(header), encoding="utf-8")


def _metadata(potential, *, nx, ny, dx, dy, correlation_length, mean_strength, seed, subtract_mean):
    return {
        "generator": "Gaussian-filtered independent random Fourier phases",
        "numpy_version": np.__version__,
        "nx": nx,
        "ny": ny,
        "dx": dx,
        "dy": dy,
        "period_x": nx * dx,
        "period_y": ny * dy,
        "correlation_length": correlation_length,
        "correlation_convention": CORRELATION_CONVENTION,
        "mean_strength_before_subtraction": mean_strength,
        "subtract_mean": subtract_mean,
        "seed": seed,
        "mean": float(np.mean(potential)),
        "standard_deviation": float(np.std(potential)),
        "minimum": float(np.min(potential)),
        "maximum": float(np.max(potential)),
        "ordering": "site=x+nx*y; potential[y,x]",
        "length_units": "same numerical units as dx and dy",
        "energy_units": "same as simulation t_hop and mu; no Hz conversion",
        "sampling_warnings": _sampling_messages(nx, ny, dx, dy, correlation_length),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("output", type=Path, help="Dense onsite potential file to write")
    parser.add_argument("--nx", type=int, required=True, help="Number of sites along x (simulation Lx)")
    parser.add_argument("--ny", type=int, help="Number of sites along y (simulation Ly); default: nx")
    parser.add_argument("--dx", type=float, default=1.0, help="Numerical grid spacing along x; default: 1")
    parser.add_argument("--dy", type=float, help="Numerical grid spacing along y; default: dx")
    parser.add_argument("--correlation-length", type=float, required=True, help="Gaussian connected-intensity correlation length in the same units as dx/dy")
    parser.add_argument("--mean-strength", type=float, default=1.0, help="Mean speckle energy before optional mean subtraction; default: 1")
    parser.add_argument("--seed", type=int, default=42, help="Local random generator seed; default: 42")
    parser.add_argument("--subtract-mean", action="store_true", help="Subtract the realized mean, retaining the original fluctuation strength")
    parser.add_argument("--metadata-json", type=Path, help="Also save generation parameters and realized statistics to this JSON file")
    args = parser.parse_args(argv)
    try:
        nx, ny, dx, dy, ell, strength, seed = _parameters(
            args.nx, args.ny, args.dx, args.dy, args.correlation_length, args.mean_strength, args.seed
        )
        if args.metadata_json is not None and args.output.resolve() == args.metadata_json.resolve():
            raise ValueError("output and metadata JSON must be different files")
        potential = generate_speckle(
            nx, ny, dx=dx, dy=dy, correlation_length=ell,
            mean_strength=strength, seed=seed, subtract_mean=args.subtract_mean,
        )
        metadata = _metadata(
            potential, nx=nx, ny=ny, dx=dx, dy=dy, correlation_length=ell,
            mean_strength=strength, seed=seed, subtract_mean=args.subtract_mean,
        )
        write_potential(args.output, potential, metadata=metadata)
        if args.metadata_json is not None:
            args.metadata_json.parent.mkdir(parents=True, exist_ok=True)
            args.metadata_json.write_text(json.dumps(metadata, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    except (ValueError, OSError) as exc:
        parser.error(str(exc))
    print(
        f"Wrote {nx} x {ny} sites to {args.output}; "
        f"mean={metadata['mean']:.8g}, std={metadata['standard_deviation']:.8g}, "
        f"correlation_length={ell:.8g}, seed={seed}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
