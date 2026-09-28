#!/usr/bin/env python3
"""Analyze independent-start grand-canonical BKT audit blocks."""

from __future__ import annotations

import argparse
import csv
import math
import re
from collections import defaultdict
from pathlib import Path

import h5py
import numpy as np


BLOCK_RE = re.compile(r"production_(\d+)\.out\.h5$")
W_UNIVERSAL_SUM = 4.0 / math.pi


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
    error = np.asarray(group["mean/error"][()], dtype=float).reshape(-1)
    tau = np.asarray(group["tau"][()], dtype=float).reshape(-1)
    count = int(scalar(group["count"]))
    return mean, error, tau, count


def mean_sem(values):
    array = np.asarray(values, dtype=float)
    if array.size == 0:
        return math.nan, math.nan
    if array.size == 1:
        return float(array[0]), math.nan
    return float(np.mean(array)), float(np.std(array, ddof=1) / math.sqrt(array.size))


def weighted_mean(rows, key):
    weights = np.asarray([row["count"] for row in rows], dtype=float)
    values = np.asarray([row[key] for row in rows], dtype=float)
    if np.any(weights <= 0):
        raise ValueError(f"non-positive measurement count while pooling {key}")
    return float(np.average(values, weights=weights))


def load_rows(roots: list[Path]):
    rows = []
    for root in roots:
        manifest_path = root / "audit_manifest.csv"
        if not manifest_path.is_file():
            raise FileNotFoundError(f"audit manifest not found: {manifest_path}")
        with manifest_path.open(newline="", encoding="utf-8-sig") as stream:
            jobs = list(csv.DictReader(stream))
        for job in jobs:
            job_dir = Path(job["JobDirectory"])
            files = sorted(
                job_dir.glob("production_*.out.h5"),
                key=lambda path: int(BLOCK_RE.search(path.name).group(1))
                if BLOCK_RE.search(path.name)
                else -1,
            )
            for path in files:
                match = BLOCK_RE.search(path.name)
                if not match:
                    continue
                with h5py.File(path, "r") as handle:
                    canonical = int(parameter(handle, "canonical"))
                    if canonical >= 0:
                        raise ValueError(f"not grand canonical: {path}")
                    lx = int(parameter(handle, "Lx"))
                    ly = int(parameter(handle, "Ly"))
                    beta = float(parameter(handle, "beta"))
                    mu = float(parameter(handle, "mu"))
                    expected_size = int(job["Size"])
                    expected_beta = float(job["Beta"])
                    expected_mu = float(job["Mu"])
                    expected_initial = int(job["InitialParticles"])
                    expected_seed = int(job["Seed"])
                    exact_parameters = {
                        "Lx": expected_size,
                        "Ly": expected_size,
                        "Lz": 1,
                        "pbcx": 1,
                        "pbcy": 1,
                        "pbcz": 0,
                        "nmax": 16,
                        "initial_particle_number": expected_initial,
                        "seed": expected_seed,
                    }
                    for name, expected in exact_parameters.items():
                        actual = int(parameter(handle, name))
                        if actual != expected:
                            raise ValueError(
                                f"parameter mismatch in {path}: {name}={actual}, "
                                f"manifest/runner expects {expected}"
                            )
                    float_parameters = {
                        "beta": expected_beta,
                        "mu": expected_mu,
                        "t_hop": 1.0,
                        "U_on": 0.152,
                        "V_nn": 0.0,
                    }
                    for name, expected in float_parameters.items():
                        actual = float(parameter(handle, name))
                        if not math.isclose(actual, expected, rel_tol=1e-12, abs_tol=1e-12):
                            raise ValueError(
                                f"parameter mismatch in {path}: {name}={actual}, "
                                f"manifest/runner expects {expected}"
                            )
                    if str(parameter(handle, "model")) != "BoseHubbard":
                        raise ValueError(f"not a BoseHubbard result: {path}")
                    n_mean, n_error, n_tau, count = observable(
                        handle, "Number_of_particles"
                    )
                    n2_mean, _, _, n2_count = observable(
                        handle, "Number_of_particles_squared"
                    )
                    winding, winding_error, winding_tau, winding_count = observable(
                        handle, "Winding_number_squared"
                    )
                    energy, energy_error, energy_tau, energy_count = observable(
                        handle, "Total_Energy"
                    )
                    if len(n_mean) != 1 or len(n2_mean) != 1:
                        raise ValueError(f"particle-number observable has wrong shape: {path}")
                    if len(winding) < 2:
                        raise ValueError(f"winding observable has fewer than two components: {path}")
                    if len({count, n2_count, winding_count, energy_count}) != 1:
                        raise ValueError(f"O(1) observable counts disagree: {path}")
                sites = lx * ly
                raw_variance = float(n2_mean[0] - n_mean[0] ** 2)
                variance_tolerance = 1e-10 * max(
                    1.0, abs(float(n2_mean[0])), float(n_mean[0] ** 2)
                )
                if raw_variance < -variance_tolerance:
                    raise ValueError(
                        f"unphysical negative particle-number variance {raw_variance} in {path}"
                    )
                variance = max(0.0, raw_variance)
                winding_sum = float(np.sum(winding))
                rows.append(
                    {
                        "root": str(root.resolve()),
                        "job": job["JobName"],
                        "start": job["Start"],
                        "initial_N": int(job["InitialParticles"]),
                        "block": int(match.group(1)),
                        "file": str(path.resolve()),
                        "L": lx,
                        "sites": sites,
                        "beta": beta,
                        "temperature_over_t": 1.0 / beta,
                        "mu": mu,
                        "count": count,
                        "N": float(n_mean[0]),
                        "N2": float(n2_mean[0]),
                        "N_internal_error": float(n_error[0]),
                        "N_tau": float(n_tau[0]),
                        "density": float(n_mean[0] / sites),
                        "var_N": variance,
                        "number_susceptibility_per_site": beta * variance / sites,
                        "Wx2": float(winding[0]),
                        "Wy2": float(winding[1]),
                        "W_sum": winding_sum,
                        "W_internal_error": float(
                            np.sqrt(np.sum(np.square(winding_error)))
                        ),
                        "W_tau_max": float(np.max(winding_tau)),
                        "bkt_jump_ratio": math.pi * winding_sum / 4.0,
                        "grand_energy": float(energy[0]),
                        "physical_energy": float(energy[0] + mu * n_mean[0]),
                        "energy_internal_error": float(energy_error[0]),
                        "energy_tau": float(energy_tau[0]),
                    }
                )
    if not rows:
        raise RuntimeError("no production_*.out.h5 files were found")
    return rows


def aggregate_rows(rows):
    groups = defaultdict(list)
    for row in rows:
        groups[(row["L"], row["beta"], row["mu"])].append(row)

    summaries = []
    for (size, beta, mu), point_rows in sorted(groups.items()):
        starts = defaultdict(list)
        for row in point_rows:
            # Keep repeated audit roots independent even when both use labels
            # such as low/mid/high.
            starts[(row["root"], row["job"])].append(row)

        family = {}
        for start, start_rows in starts.items():
            family_n = weighted_mean(start_rows, "N")
            family_n2 = weighted_mean(start_rows, "N2")
            family_variance = family_n2 - family_n**2
            variance_tolerance = 1e-10 * max(1.0, abs(family_n2), family_n**2)
            if family_variance < -variance_tolerance:
                raise ValueError(f"unphysical pooled N variance for family {start}")
            family[start] = {
                key: weighted_mean(start_rows, key)
                for key in (
                    "density",
                    "Wx2",
                    "Wy2",
                    "W_sum",
                    "grand_energy",
                    "physical_energy",
                )
            }
            family[start]["N"] = family_n
            family[start]["N2"] = family_n2
            family[start]["var_N"] = max(0.0, family_variance)

        start_names = sorted(family)
        n_mean, n_sem = mean_sem([family[name]["N"] for name in start_names])
        density_mean, density_sem = mean_sem(
            [family[name]["density"] for name in start_names]
        )
        variance_mean, variance_sem = mean_sem(
            [family[name]["var_N"] for name in start_names]
        )
        wx_mean, wx_sem = mean_sem([family[name]["Wx2"] for name in start_names])
        wy_mean, wy_sem = mean_sem([family[name]["Wy2"] for name in start_names])
        winding_mean, winding_sem = mean_sem(
            [family[name]["W_sum"] for name in start_names]
        )
        grand_energy_mean, grand_energy_sem = mean_sem(
            [family[name]["grand_energy"] for name in start_names]
        )
        physical_energy_mean, physical_energy_sem = mean_sem(
            [family[name]["physical_energy"] for name in start_names]
        )
        family_n = [family[name]["N"] for name in start_names]
        family_w = [family[name]["W_sum"] for name in start_names]
        block_numbers = sorted({row["block"] for row in point_rows})
        block_n = {
            block: float(np.mean([row["N"] for row in point_rows if row["block"] == block]))
            for block in block_numbers
        }
        block_w = {
            block: float(
                np.mean([row["W_sum"] for row in point_rows if row["block"] == block])
            )
            for block in block_numbers
        }
        if len(block_numbers) > 1:
            n_first_last_shift = block_n[block_numbers[-1]] - block_n[block_numbers[0]]
            w_first_last_shift = block_w[block_numbers[-1]] - block_w[block_numbers[0]]
        else:
            n_first_last_shift = math.nan
            w_first_last_shift = math.nan
        initial_values = [
            row["initial_N"]
            for row in point_rows
            if row["block"] == block_numbers[0]
        ]
        summaries.append(
            {
                "L": size,
                "beta": beta,
                "temperature_over_t": 1.0 / beta,
                "mu": mu,
                "families": len(start_names),
                "blocks": len(point_rows),
                "blocks_per_family_min": min(len(starts[name]) for name in starts),
                "samples_total": sum(row["count"] for row in point_rows),
                "N": n_mean,
                "N_start_sem": n_sem,
                "N_start_range": max(family_n) - min(family_n),
                "initial_N_range": max(initial_values) - min(initial_values),
                "N_first_last_shift": n_first_last_shift,
                "density": density_mean,
                "density_start_sem": density_sem,
                "var_N": variance_mean,
                "var_N_start_sem": variance_sem,
                "number_susceptibility_per_site": beta
                * variance_mean
                / (size * size),
                "Wx2": wx_mean,
                "Wx2_start_sem": wx_sem,
                "Wy2": wy_mean,
                "Wy2_start_sem": wy_sem,
                "xy_fractional_difference": abs(wx_mean - wy_mean)
                / max(0.5 * (wx_mean + wy_mean), np.finfo(float).tiny),
                "W_sum": winding_mean,
                "W_start_sem": winding_sem,
                "W_start_range": max(family_w) - min(family_w),
                "W_first_last_shift": w_first_last_shift,
                "bkt_jump_ratio": math.pi * winding_mean / 4.0,
                "grand_energy": grand_energy_mean,
                "grand_energy_start_sem": grand_energy_sem,
                "physical_energy": physical_energy_mean,
                "physical_energy_start_sem": physical_energy_sem,
                "N_tau_max": max(row["N_tau"] for row in point_rows),
                "W_tau_max": max(row["W_tau_max"] for row in point_rows),
                "N_effective_min": min(
                    row["count"] / (2.0 * max(row["N_tau"], 0.5))
                    for row in point_rows
                ),
                "W_effective_min": min(
                    row["count"] / (2.0 * max(row["W_tau_max"], 0.5))
                    for row in point_rows
                ),
            }
        )
    return summaries


def crossing_estimates(summaries):
    by_size = defaultdict(list)
    for row in summaries:
        by_size[row["L"]].append(row)
    crossings = {}
    for size, points in by_size.items():
        points.sort(key=lambda row: row["beta"])
        crossing = None
        for left, right in zip(points, points[1:]):
            y0 = left["W_sum"] - W_UNIVERSAL_SUM
            y1 = right["W_sum"] - W_UNIVERSAL_SUM
            if y0 == 0:
                crossing = left["beta"]
                break
            if y0 * y1 <= 0 and right["W_sum"] != left["W_sum"]:
                fraction = (W_UNIVERSAL_SUM - left["W_sum"]) / (
                    right["W_sum"] - left["W_sum"]
                )
                crossing = left["beta"] + fraction * (
                    right["beta"] - left["beta"]
                )
                break
        crossings[size] = crossing
    return crossings


def write_csv(path: Path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_markdown(path: Path, summaries, crossings):
    lines = [
        "# Grand-canonical BKT audit summary",
        "",
        "Uncertainties labelled start SEM are calculated from the low, central, and high",
        "independent-start family means. They are more conservative than the pooled",
        "ALPS accumulator error when particle-number autocorrelation is long, but three",
        "families give only a descriptive uncertainty rather than a precision confidence interval.",
        "",
        "| L | beta | T/t | density | chi N/site | Wx2 | Wy2 | W sum | start SEM | pi W/4 |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in summaries:
        lines.append(
            "| {L} | {beta:.6f} | {temperature_over_t:.4f} | "
            "{density:.5f} | {number_susceptibility_per_site:.4f} | "
            "{Wx2:.5f} | {Wy2:.5f} | {W_sum:.5f} | {W_start_sem:.5f} | "
            "{bkt_jump_ratio:.4f} |".format(**row)
        )
    lines.extend(
        [
            "",
            "## Convergence diagnostics",
            "",
            "The initial N span is the deliberately imposed low-to-high starting",
            "difference. The retained N span compares retained independent-family means.",
            "First-to-last shifts compare numbered block averages and are blank when",
            "only one production block exists.",
            "",
            "| L | beta | families | blocks/family | initial N span | retained N span | first-last N | first-last W | x-y fraction | min effective N |",
            "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in summaries:
        n_shift = (
            "--"
            if math.isnan(row["N_first_last_shift"])
            else f"{row['N_first_last_shift']:.2f}"
        )
        w_shift = (
            "--"
            if math.isnan(row["W_first_last_shift"])
            else f"{row['W_first_last_shift']:.4f}"
        )
        lines.append(
            f"| {row['L']} | {row['beta']:.6f} | {row['families']} | "
            f"{row['blocks_per_family_min']} | {row['initial_N_range']:.0f} | "
            f"{row['N_start_range']:.2f} | {n_shift} | {w_shift} | "
            f"{row['xy_fractional_difference']:.4f} | "
            f"{row['N_effective_min']:.0f} |"
        )
    lines.extend(
        [
            "",
            "The thermodynamic universal-jump reference is W sum = 4/pi = "
            f"{W_UNIVERSAL_SUM:.6f}. A raw finite-L crossing is only a diagnostic;",
            "a final transition temperature requires logarithmic BKT finite-size scaling.",
            "",
            "## Raw 4/pi crossing estimates",
            "",
        ]
    )
    for size, crossing in sorted(crossings.items()):
        if crossing is None:
            lines.append(f"- L={size}: scan does not bracket 4/pi.")
        else:
            lines.append(
                f"- L={size}: beta approximately {crossing:.6f}, "
                f"T/t approximately {1.0 / crossing:.4f}."
            )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_plot(path: Path, summaries):
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return False

    fig, axes = plt.subplots(2, 1, figsize=(7.2, 7.2), sharex=True)
    sizes = sorted({row["L"] for row in summaries})
    for size in sizes:
        points = sorted(
            (row for row in summaries if row["L"] == size),
            key=lambda row: row["beta"],
        )
        beta = [row["beta"] for row in points]
        winding = [row["W_sum"] for row in points]
        winding_error = [
            0.0 if math.isnan(row["W_start_sem"]) else row["W_start_sem"]
            for row in points
        ]
        density = [row["density"] for row in points]
        density_error = [
            0.0 if math.isnan(row["density_start_sem"]) else row["density_start_sem"]
            for row in points
        ]
        axes[0].errorbar(beta, winding, yerr=winding_error, marker="o", label=f"L={size}")
        axes[1].errorbar(beta, density, yerr=density_error, marker="o", label=f"L={size}")
    axes[0].axhline(W_UNIVERSAL_SUM, color="black", linestyle="--", linewidth=1, label="4/pi")
    axes[0].set_ylabel("Wx^2 + Wy^2")
    axes[0].legend()
    axes[0].grid(alpha=0.25)
    axes[1].set_xlabel("beta t")
    axes[1].set_ylabel("<N>/L^2")
    axes[1].grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Analyze production blocks made by run_gc_bkt_audit.ps1."
    )
    parser.add_argument("roots", nargs="+", type=Path, help="Audit result directories")
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()

    roots = [root.resolve() for root in args.roots]
    output_dir = (
        args.output_dir.resolve()
        if args.output_dir
        else roots[0] / "analysis"
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    rows = load_rows(roots)
    summaries = aggregate_rows(rows)
    crossings = crossing_estimates(summaries)
    write_csv(output_dir / "audit_blocks.csv", rows)
    write_csv(output_dir / "audit_points.csv", summaries)
    write_markdown(output_dir / "AUDIT_SUMMARY.md", summaries, crossings)
    plotted = write_plot(output_dir / "audit_summary.png", summaries)

    print(f"read {len(rows)} production blocks from {len(roots)} audit root(s)")
    print(f"aggregated {len(summaries)} (L,beta,mu) points")
    print(f"wrote {output_dir / 'audit_blocks.csv'}")
    print(f"wrote {output_dir / 'audit_points.csv'}")
    print(f"wrote {output_dir / 'AUDIT_SUMMARY.md'}")
    if plotted:
        print(f"wrote {output_dir / 'audit_summary.png'}")
    else:
        print("matplotlib is unavailable; skipped the PNG plot")


if __name__ == "__main__":
    main()
