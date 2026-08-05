#!/usr/bin/env python3
"""Run every full quantitative validation documented by LBFAST."""

import argparse
import csv
import math
import shutil
import subprocess
import sys
import time
from pathlib import Path


def read_one_row(path):
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != 1:
        raise RuntimeError("expected one data row in {}".format(path))
    result = {}
    for key, value in rows[0].items():
        number = float(value)
        if not math.isfinite(number):
            raise RuntimeError("non-finite {} in {}".format(key, path))
        result[key] = number
    return result


def read_many_rows(path):
    with path.open(newline="") as stream:
        source = list(csv.DictReader(stream))
    if not source:
        raise RuntimeError("no data rows in {}".format(path))
    result = []
    for source_row in source:
        row = {key: float(value) for key, value in source_row.items()}
        if not all(math.isfinite(value) for value in row.values()):
            raise RuntimeError("non-finite value in {}".format(path))
        result.append(row)
    return result


def read_key_values(path):
    values = {}
    for line in path.read_text().splitlines():
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        values[key.strip()] = float(value)
    return values


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def check_laplace(directory):
    rows = read_many_rows(directory / "laplace_csf_sweep_summary.csv")
    fit = read_key_values(directory / "laplace_csf_fit.txt")
    max_error = max(abs(row["relative_error"]) for row in rows)
    slope_error = abs(fit["slope_error"])
    require(len(rows) == 3, "Laplace sweep did not contain three radii")
    require(max_error <= 0.03, "maximum pointwise Laplace error exceeds 3%")
    require(slope_error <= 0.03, "Laplace slope error exceeds 3%")
    return slope_error, 0.03, "slope error={:.3%}; maximum point error={:.3%}".format(
        slope_error, max_error)


def check_capillary(directory):
    row = read_one_row(directory / "capillary_wave_summary.csv")
    frequency_error = abs(row["frequency_error_percent"]) / 100.0
    require(frequency_error <= 0.12, "capillary-wave frequency error exceeds 12%")
    require(row["max_mass_error"] <= 1.0e-6, "phase-mass error exceeds 1e-6")
    require(abs(row["force_ratio"] - 1.0) <= 0.02,
            "CSF projection differs from unity by more than 2%")
    require(row["kinematic_relative_rms_error"] <= 0.02,
            "phase kinematic RMS error exceeds 2%")
    return frequency_error, 0.12, "frequency error={:.3%}; force ratio={:.6f}".format(
        frequency_error, row["force_ratio"])


def check_lamb(directory):
    row = read_one_row(directory / "lamb_validation_summary.csv")
    period_error = abs(row["period_error_percent"]) / 100.0
    require(row["peaks"] >= 2.0, "fewer than two Lamb peaks were detected")
    require(period_error <= 0.10, "Lamb period error exceeds 10%")
    return period_error, 0.10, "period error={:.3%}; peaks={:.0f}".format(
        period_error, row["peaks"])


def check_taylor_green(directory):
    row = read_one_row(directory / "taylorgreen_validation_summary.csv")
    viscosity_error = abs(row["viscosity_error_percent"]) / 100.0
    require(viscosity_error <= 0.03, "Taylor-Green viscosity error exceeds 3%")
    require(row["fit_log_energy_rms"] <= 2.0e-3,
            "Taylor-Green log-energy fit RMS exceeds 2e-3")
    return viscosity_error, 0.03, "viscosity error={:.3%}; fit RMS={:.3e}".format(
        viscosity_error, row["fit_log_energy_rms"])


def check_poiseuille(directory):
    row = read_one_row(directory / "poiseuille_force_validation_summary.csv")
    relative_l2 = row["relative_l2_error"]
    center_error = abs(row["center_velocity_error_percent"]) / 100.0
    require(relative_l2 <= 0.01, "Poiseuille relative L2 error exceeds 1%")
    require(center_error <= 0.01, "Poiseuille centerline error exceeds 1%")
    require(row["statistics_analytical_max_difference"] <= 1.0e-12,
            "Python and statistics.f90 analytical profiles disagree")
    return relative_l2, 0.01, "relative L2={:.3%}; center error={:.3%}".format(
        relative_l2, center_error)


CASES = (
    ("laplace_csf", "run_laplace_csf_sweep.py", "laplace_csf_sweep", check_laplace),
    ("capillary_wave", "analyse_capillary_wave.py", "capillary_wave_run", check_capillary),
    ("lamb_oscillation", "run_lamb_validation.py", "lamb_validation_run", check_lamb),
    ("taylor_green", "run_taylorgreen_validation.py", "taylorgreen_validation_run",
     check_taylor_green),
    ("poiseuille_force", "run_poiseuille_force_validation.py",
     "poiseuille_force_validation_run", check_poiseuille),
)

MPI_DECOMPOSITIONS = {
    2: {
        "laplace_csf": (1, 1, 2),
        "capillary_wave": (1, 2, 1),
        "lamb_oscillation": (1, 1, 2),
        "taylor_green": (1, 1, 2),
        "poiseuille_force": (1, 1, 2),
    },
    4: {
        "laplace_csf": (1, 1, 4),
        "capillary_wave": (2, 2, 1),
        "lamb_oscillation": (1, 2, 2),
        "taylor_green": (1, 1, 4),
        "poiseuille_force": (1, 1, 4),
    },
}


def run_driver(command, root, log_path):
    print("+ " + " ".join(command), flush=True)
    with log_path.open("w") as stream:
        return subprocess.run(command, cwd=str(root), stdout=stream,
                              stderr=subprocess.STDOUT, check=False).returncode


def write_results(path, results):
    fields = ("case", "status", "mpi_procs", "decomposition", "seconds",
              "metric", "tolerance", "detail")
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(results)


def main():
    tests = Path(__file__).resolve().parent
    root = tests.parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gpu-cc", default="80")
    parser.add_argument("--make-target", default=None,
                        help="defaults to nvfortran or nvfortran-mpi according to --mpi-procs")
    parser.add_argument("--mpi-procs", type=int, choices=(1, 2, 4), default=1,
                        help="run sequentially or on 2/4 MPI ranks with compatible per-case grids")
    parser.add_argument("--launcher", choices=("mpirun", "srun"), default="mpirun")
    parser.add_argument("--case", action="append", choices=[case[0] for case in CASES])
    parser.add_argument("--keep-results", action="store_true")
    parser.add_argument("--fail-fast", action="store_true")
    args = parser.parse_args()
    make_target = args.make_target or ("nvfortran-mpi" if args.mpi_procs > 1 else "nvfortran")
    if args.mpi_procs > 1 and "mpi" not in make_target.lower():
        raise SystemExit("--mpi-procs greater than one requires an MPI Make target")
    if args.mpi_procs == 1 and "mpi" in make_target.lower():
        raise SystemExit("an MPI Make target requires --mpi-procs 2 or 4")

    selected = [case for case in CASES if not args.case or case[0] in args.case]
    output = tests / "full_results"
    if output.exists() and not args.keep_results:
        shutil.rmtree(str(output))
    output.mkdir(parents=True, exist_ok=True)
    defines = root / "defines.h"
    binary = root / "main.x"
    original_defines = defines.read_text()
    results = []

    try:
        for name, driver, directory_name, checker in selected:
            started = time.monotonic()
            status, metric, tolerance, detail = "FAIL", float("nan"), float("nan"), ""
            case_output = output / directory_name
            case_output.mkdir(parents=True, exist_ok=True)
            driver_log = case_output / "driver.log"
            command = [sys.executable, str(root / driver),
                       "--output", str(case_output),
                       "--gpu-cc", args.gpu_cc,
                       "--make-target", make_target,
                       "--force"]
            decomposition = None
            if args.mpi_procs > 1:
                decomposition = MPI_DECOMPOSITIONS[args.mpi_procs][name]
                command.extend(["--mpi-procs", str(args.mpi_procs),
                                "--decomposition"] +
                               [str(value) for value in decomposition] +
                               ["--launcher", args.launcher])
            try:
                code = run_driver(command, root, driver_log)
                if code:
                    raise RuntimeError("driver exited with status {}; see {}".format(
                        code, driver_log))
                metric, tolerance, detail = checker(case_output)
                status = "PASS"
            except Exception as error:
                detail = str(error)
            results.append({
                "case": name,
                "status": status,
                "mpi_procs": str(args.mpi_procs),
                "decomposition": "x".join(str(value) for value in decomposition)
                                 if decomposition else "1x1x1",
                "seconds": "{:.3f}".format(time.monotonic() - started),
                "metric": "{:.12e}".format(metric) if math.isfinite(metric) else "",
                "tolerance": "{:.12e}".format(tolerance) if math.isfinite(tolerance) else "",
                "detail": detail,
            })
            write_results(output / "full_validation_results.csv", results)
            print("{:<5} {:<20} {}".format(status, name, detail), flush=True)
            if status == "FAIL" and args.fail_fast:
                break
    finally:
        defines.write_text(original_defines)
        subprocess.run(["make", "clean"], cwd=str(root), check=False,
                       stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
        if binary.exists():
            binary.unlink()
        write_results(output / "full_validation_results.csv", results)

    passed = sum(result["status"] == "PASS" for result in results)
    failed = len(results) - passed
    print("\n{} passed, {} failed".format(passed, failed))
    print("results={}".format(output / "full_validation_results.csv"))
    raise SystemExit(0 if failed == 0 and len(results) == len(selected) else 1)


if __name__ == "__main__":
    main()
