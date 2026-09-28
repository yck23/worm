#!/usr/bin/env python3
"""Check grand-canonical d<N>/dmu = beta Var(N) from two audit roots."""

from __future__ import annotations

import argparse
import math
import statistics
from pathlib import Path

import h5py
import numpy as np


def scalar(dataset):
    value = np.asarray(dataset[()]).reshape(-1)[0]
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return value.item() if hasattr(value, "item") else value


def parameter(handle: h5py.File, name: str):
    return scalar(handle[f"/parameters/{name}"])


def observable(handle: h5py.File, name: str):
    group = handle[f"/simulation/results/{name}"]
    mean = np.asarray(group["mean/value"][()], dtype=float).reshape(-1)
    count = int(scalar(group["count"]))
    return mean, count


def mean_sem(values):
    values = list(values)
    mean = statistics.fmean(values)
    if len(values) < 2:
        return mean, math.nan
    return mean, statistics.stdev(values) / math.sqrt(len(values))


def summarize(root: Path):
    family_summaries = []
    reference = None
    for family_dir in sorted(path for path in root.iterdir() if path.is_dir()):
        files = sorted(family_dir.glob("production_*.out.h5"))
        if not files:
            continue
        rows = []
        for path in files:
            with h5py.File(path, "r") as handle:
                point = {
                    "Lx": int(parameter(handle, "Lx")),
                    "Ly": int(parameter(handle, "Ly")),
                    "Lz": int(parameter(handle, "Lz")),
                    "pbcx": int(parameter(handle, "pbcx")),
                    "pbcy": int(parameter(handle, "pbcy")),
                    "pbcz": int(parameter(handle, "pbcz")),
                    "beta": float(parameter(handle, "beta")),
                    "mu": float(parameter(handle, "mu")),
                    "canonical": int(parameter(handle, "canonical")),
                    "model": str(parameter(handle, "model")),
                    "t_hop": float(parameter(handle, "t_hop")),
                    "U_on": float(parameter(handle, "U_on")),
                    "V_nn": float(parameter(handle, "V_nn")),
                    "nmax": int(parameter(handle, "nmax")),
                }
                if point["canonical"] >= 0:
                    raise ValueError(f"not grand canonical: {path}")
                if reference is None:
                    reference = point
                elif point != reference:
                    raise ValueError(f"parameter mismatch: {path}")
                number, count = observable(handle, "Number_of_particles")
                number2, count2 = observable(handle, "Number_of_particles_squared")
                winding, winding_count = observable(handle, "Winding_number_squared")
                if count <= 0 or count != count2 or count != winding_count:
                    raise ValueError(f"observable-count mismatch: {path}")
                rows.append(
                    {
                        "count": count,
                        "N": float(number[0]),
                        "N2": float(number2[0]),
                        "W": float(np.sum(winding)),
                    }
                )
        weights = np.asarray([row["count"] for row in rows], dtype=float)
        n_mean = float(np.average([row["N"] for row in rows], weights=weights))
        n2_mean = float(np.average([row["N2"] for row in rows], weights=weights))
        variance = n2_mean - n_mean**2
        tolerance = 1e-10 * max(1.0, abs(n2_mean), n_mean**2)
        if variance < -tolerance:
            raise ValueError(f"negative pooled variance in {family_dir}: {variance}")
        family_summaries.append(
            {
                "family": family_dir.name,
                "blocks": len(rows),
                "N": n_mean,
                "var_N": max(0.0, variance),
                "W": float(np.average([row["W"] for row in rows], weights=weights)),
            }
        )
    if reference is None or not family_summaries:
        raise RuntimeError(f"no family production blocks found below {root}")
    if len({entry["blocks"] for entry in family_summaries}) != 1:
        raise ValueError(f"families have unequal production-block counts below {root}")
    n_mean, n_sem = mean_sem(entry["N"] for entry in family_summaries)
    variance_mean, variance_sem = mean_sem(
        entry["var_N"] for entry in family_summaries
    )
    winding_mean, winding_sem = mean_sem(entry["W"] for entry in family_summaries)
    lx = reference["Lx"]
    ly = reference["Ly"]
    beta = reference["beta"]
    mu = reference["mu"]
    signature_keys = (
        "Lx", "Ly", "Lz", "pbcx", "pbcy", "pbcz", "beta", "canonical",
        "model", "t_hop", "U_on", "V_nn", "nmax",
    )
    return {
        "root": root,
        "L": lx,
        "sites": lx * ly,
        "beta": beta,
        "mu": mu,
        "signature": tuple(reference[name] for name in signature_keys),
        "families": len(family_summaries),
        "blocks_per_family": family_summaries[0]["blocks"],
        "N": n_mean,
        "N_sem": n_sem,
        "density": n_mean / (lx * ly),
        "var_N": variance_mean,
        "var_N_sem": variance_sem,
        "W": winding_mean,
        "W_sem": winding_sem,
        "family": family_summaries,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Compare two retained grand-canonical mu points with the fluctuation-response identity."
    )
    parser.add_argument("roots", nargs=2, type=Path, help="two roots containing family/production_*.out.h5")
    args = parser.parse_args()

    points = sorted((summarize(path.resolve()) for path in args.roots), key=lambda row: row["mu"])
    low, high = points
    if low["signature"] != high["signature"]:
        raise ValueError(
            "the two roots must have identical model, geometry, boundary, cutoff, "
            "and beta parameters; only mu may differ"
        )
    delta_mu = high["mu"] - low["mu"]
    if delta_mu == 0:
        raise ValueError("the two roots have the same mu")
    derivative = (high["N"] - low["N"]) / delta_mu
    derivative_sem = math.hypot(low["N_sem"], high["N_sem"]) / abs(delta_mu)
    fdt = low["beta"] * 0.5 * (low["var_N"] + high["var_N"])
    fdt_sem = low["beta"] * 0.5 * math.hypot(low["var_N_sem"], high["var_N_sem"])

    print("Grand-canonical fluctuation-response audit")
    print(f"L={low['L']}, beta={low['beta']:.12g}")
    for point in points:
        print(
            f"mu={point['mu']:.8g}: N={point['N']:.6f} +/- {point['N_sem']:.6f}, "
            f"density={point['density']:.8f}, Var(N)={point['var_N']:.6f} +/- "
            f"{point['var_N_sem']:.6f}, Wsum={point['W']:.6f} +/- {point['W_sem']:.6f}"
        )
        print(
            "  family N: "
            + ", ".join(f"{entry['family']}={entry['N']:.6f}" for entry in point["family"])
        )
    print(f"Delta<N>/Delta mu = {derivative:.6f} +/- {derivative_sem:.6f} (descriptive start SEM)")
    print(f"beta * endpoint-average Var(N) = {fdt:.6f} +/- {fdt_sem:.6f}")
    print(f"response/FDT ratio = {derivative / fdt:.6f}")
    print(f"relative difference = {abs(derivative - fdt) / abs(fdt):.4%}")
    print("The finite mu interval gives a secant response, so exact equality is not expected.")


if __name__ == "__main__":
    main()
