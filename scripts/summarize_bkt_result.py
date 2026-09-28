#!/usr/bin/env python3
"""Print the basic canonical or grand-canonical observables from one result."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import numpy as np


def scalar(dataset) -> float:
    return float(np.asarray(dataset[()]).reshape(-1)[0])


def parameter(handle: h5py.File, name: str):
    value = np.asarray(handle[f"/parameters/{name}"][()]).reshape(-1)[0]
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return value.item() if hasattr(value, "item") else value


def mean_and_error(results: h5py.Group, name: str):
    group = results[name]
    mean = np.asarray(group["mean/value"][()], dtype=float)
    error = np.asarray(group["mean/error"][()], dtype=float)
    count = int(np.asarray(group["count"][()]).reshape(-1)[0])
    tau = np.asarray(group["tau"][()], dtype=float)
    return mean, error, count, tau


def fmt(value) -> str:
    array = np.asarray(value).reshape(-1)
    if array.size == 1:
        return f"{array[0]:.12g}"
    return "[" + ", ".join(f"{entry:.12g}" for entry in array) + "]"


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Summarize particle number, energy, and winding from a WORM HDF5 result."
    )
    parser.add_argument("result", type=Path, help="Path to a .out.h5 result file")
    args = parser.parse_args()

    if not args.result.is_file():
        parser.error(f"result file not found: {args.result}")

    with h5py.File(args.result, "r") as handle:
        results = handle["/simulation/results"]
        lx = int(parameter(handle, "Lx"))
        ly = int(parameter(handle, "Ly"))
        sites = lx * ly
        beta = float(parameter(handle, "beta"))
        mu = float(parameter(handle, "mu"))
        canonical = int(parameter(handle, "canonical"))

        n_mean, n_error, count, n_tau = mean_and_error(results, "Number_of_particles")
        n2_mean, _, _, _ = mean_and_error(results, "Number_of_particles_squared")
        kinetic, kinetic_error, _, _ = mean_and_error(results, "Kinetic_Energy")
        potential, potential_error, _, _ = mean_and_error(results, "Potential_Energy")
        total, total_error, _, _ = mean_and_error(results, "Total_Energy")
        winding, winding_error, _, winding_tau = mean_and_error(
            results, "Winding_number_squared"
        )

        n_value = scalar(n_mean)
        n2_value = scalar(n2_mean)
        variance = max(0.0, n2_value - n_value * n_value)
        density = n_value / sites
        winding_sum = float(np.sum(winding))
        winding_sum_error = float(np.sqrt(np.sum(np.square(winding_error))))
        stiffness = winding_sum / (2.0 * beta)
        bkt_jump_ratio = np.pi * winding_sum / 4.0
        physical_potential = scalar(potential) + mu * n_value
        physical_total = scalar(total) + mu * n_value

        ensemble = "grand canonical (mu fixed, N measured)" if canonical < 0 else "canonical (N fixed)"
        print(f"file: {args.result.resolve()}")
        print(f"ensemble: {ensemble}")
        print(f"lattice: {lx} x {ly} = {sites} sites")
        print(f"beta: {beta:.12g}")
        print(f"mu: {mu:.12g}")
        print(f"samples: {count}")
        print(f"<N>: {n_value:.12g} +/- {scalar(n_error):.6g}")
        print(f"<N>/sites: {density:.12g}")
        print(f"<N^2>: {n2_value:.12g}")
        print(f"Var(N): {variance:.12g}")
        print(f"sqrt(Var(N)): {np.sqrt(variance):.12g}")
        print(f"tau_N: {fmt(n_tau)}")
        if canonical < 0:
            susceptibility = beta * variance / sites
            print(f"beta*Var(N)/sites: {susceptibility:.12g}  (number susceptibility per site)")
        print(f"Kinetic_Energy: {fmt(kinetic)} +/- {fmt(kinetic_error)}")
        print(f"Potential_Energy (interaction - mu*N): {fmt(potential)} +/- {fmt(potential_error)}")
        print(f"Total_Energy <H0-mu*N>: {fmt(total)} +/- {fmt(total_error)}")
        print(f"interaction energy without -mu*N: {physical_potential:.12g}")
        print(f"physical internal energy <H0>: {physical_total:.12g}")
        print("No error is printed for H0 because its covariance with N is not stored.")
        print(f"<Wx^2>, <Wy^2>: {fmt(winding)}")
        print(f"component errors: {fmt(winding_error)}")
        print(f"winding sum: {winding_sum:.12g} +/- {winding_sum_error:.6g} (quadrature only)")
        print(f"finite-size superfluid stiffness: {stiffness:.12g}")
        print(f"pi*(winding sum)/4: {bkt_jump_ratio:.12g}  (thermodynamic BKT jump reference = 1)")
        print(f"winding tau: {fmt(winding_tau)}")

        if canonical < 0 and variance < 0.25:
            print("WARNING: almost no particle-number variance; inspect number-sector mixing.")
        print("Convergence requires agreement across late blocks and independent starts; this summary alone is not a convergence test.")


if __name__ == "__main__":
    main()
