#!/usr/bin/env python3
"""Report current grand-canonical audit means and checkpoint endpoint sectors."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import h5py
import numpy as np


def scalar(dataset):
    value = np.asarray(dataset[()]).reshape(-1)[0]
    return value.item() if hasattr(value, "item") else value


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Inspect non-retained Fresh/Resume state from run_gc_bkt_audit.ps1. "
            "Use analyze_gc_bkt_audit.py for retained production blocks."
        )
    )
    parser.add_argument("root", type=Path, help="Audit result directory")
    args = parser.parse_args()

    root = args.root.resolve()
    manifest = root / "audit_manifest.csv"
    if not manifest.is_file():
        parser.error(f"manifest not found: {manifest}")

    with manifest.open(newline="", encoding="utf-8-sig") as stream:
        jobs = list(csv.DictReader(stream))

    print(
        "beta  start  initial  checkpoint endpoints  cumulative <N>  "
        "density   Wsum    samples   tau_N"
    )
    for job in jobs:
        directory = Path(job["JobDirectory"])
        result_path = directory / "result.out.h5"
        checkpoint_base = directory / "checkpoint.clone.h5"
        processes = int(job["Processes"])
        endpoints = []
        for rank in range(processes):
            checkpoint = (
                checkpoint_base
                if rank == 0
                else Path(str(checkpoint_base) + f".{rank}")
            )
            with h5py.File(checkpoint, "r") as handle:
                endpoints.append(
                    int(
                        scalar(
                            handle[
                                "/simulation/realizations/0/clones/0/"
                                "checkpoint/configuration/nr_of_particles"
                            ]
                        )
                    )
                )

        with h5py.File(result_path, "r") as handle:
            results = handle["/simulation/results"]
            number = float(scalar(results["Number_of_particles/mean/value"]))
            winding = float(
                np.asarray(
                    results["Winding_number_squared/mean/value"][()], dtype=float
                ).sum()
            )
            count = int(scalar(results["Number_of_particles/count"]))
            tau = float(scalar(results["Number_of_particles/tau"]))
            lx = int(scalar(handle["/parameters/Lx"]))
            ly = int(scalar(handle["/parameters/Ly"]))

        beta = float(job["Beta"])
        start = job["Start"]
        initial = int(job["InitialParticles"])
        print(
            f"{beta:.2f}  {start:>5}  "
            f"{initial:7d}  {str(endpoints):>20}  "
            f"{number:14.2f}  {number/(lx*ly):7.4f}  {winding:6.4f}  "
            f"{count:8d}  {tau:7.1f}"
        )


if __name__ == "__main__":
    main()
