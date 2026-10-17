#!/usr/bin/env python3
"""Build, run, and analyse the LBFAST planar capillary-wave benchmark."""

# Benchmark: planar capillary wave on a periodic liquid slab in the two-component solver with
# CSF (docs/CAPILLARY_WAVE_VALIDATION.md).
# Build: unless --skip-build is given, defines.h is edited in place to enable LATTICE 27,
# HIGHORDER, TWOCOMPONENT, DENSRATIO, CSF, DOBENCHMARK, CAPILLARYWAVE and PRINTPHI and to
# disable the LAPLACE, LAMBTEST, POISEUILLE, TWOPOISEUILLE and TAYLORGREEN selectors; then
# "make clean" and "make <target> GPUCC=<cc>" are run. Other macros (PRC, STRPRC, MIXEDPRC, ...)
# are kept as found.
# Run: capillary_wave.inp (96 x 288 x 8, initial amplitude uwall = 0.96, slab half-thickness
# radius = 72, 80000 steps, diagnostics every 20 steps) is copied to and run in
# capillary_wave_run/.
# Metrics, written to capillary_wave_summary.csv:
#  - angular frequency from the mean spacing of the zero crossings of the cosine amplitude
#    eta_cos, against the finite-depth inviscid value written by the solver;
#  - maximum relative error of the total phase mass;
#  - kinematic check: least-squares gain and relative RMS difference between d eta/dt and the
#    interface-normal velocity mode;
#  - ratio of the numerical to the theoretical Fourier coefficient of the initial CSF force.
# This script applies no tolerance; tests/run_full_validation.py requires a frequency error
# <= 12%, a mass error <= 1e-6, a force ratio within 2% of one and a kinematic RMS error <= 2%.

import argparse
import csv
import math
import re
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Sequence, Tuple
from validation_mpi import add_mpi_arguments, solver_command, validate_local_tiles, validate_mpi_arguments


REQUIRED_DEFINES = {
    "LATTICE": "27",
    "HIGHORDER": None,
    "TWOCOMPONENT": None,
    "DENSRATIO": None,
    "CSF": None,
    "DOBENCHMARK": None,
    "CAPILLARYWAVE": None,
    "PRINTPHI": None,
}
DISABLED_CASES = ("LAPLACE", "LAMBTEST", "POISEUILLE", "TWOPOISEUILLE", "TAYLORGREEN")


def set_define(text: str, name: str, value: str = None) -> str:
    replacement = f"#define {name}" + (f" {value}" if value is not None else "")
    pattern = re.compile(
        rf"(?m)^[ \t]*#define[ \t]+(?:no)?{re.escape(name)}"
        rf"(?:[ \t]+[^\n]*)?[ \t]*$"
    )
    updated, count = pattern.subn(replacement, text, count=1)
    return updated if count else text.rstrip() + f"\n{replacement}\n"


def disable_define(text: str, name: str) -> str:
    """Rename an enabled '#define NAME' line without value to '#define noNAME'."""
    pattern = re.compile(rf"(?m)^[ \t]*#define[ \t]+{re.escape(name)}[ \t]*$")
    updated, _ = pattern.subn(f"#define no{name}", text, count=1)
    return updated


def configure_defines(path: Path) -> None:
    original = path.read_text()
    configured = original
    for name, value in REQUIRED_DEFINES.items():
        configured = set_define(configured, name, value)
    for name in DISABLED_CASES:
        configured = disable_define(configured, name)
    if configured != original:
        path.write_text(configured)
        print(f"Updated planar capillary-wave macros in {path}", flush=True)
    else:
        print(f"Planar capillary-wave macros already configured in {path}", flush=True)


def run_command(command: Sequence[str], cwd: Path) -> None:
    print("+ " + " ".join(command), flush=True)
    completed = subprocess.run(command, cwd=cwd, check=False)
    if completed.returncode != 0:
        raise SystemExit(f"Command failed: {' '.join(command)}")


def build_solver(root: Path, defines: Path, gpu_cc: str, target: str) -> None:
    configure_defines(defines)
    run_command(["make", "clean"], root)
    run_command(["make", target, f"GPUCC={gpu_cc}"], root)


def zero_crossings(times: Sequence[float], signal: Sequence[float]) -> List[float]:
    """Times of the sign changes of signal, by linear interpolation between samples."""
    roots: List[float] = []
    for index in range(1, len(signal)):
        y0, y1 = signal[index - 1], signal[index]
        if y0 == 0.0:
            roots.append(times[index - 1])
        elif y0 * y1 < 0.0:
            fraction = abs(y0) / (abs(y0) + abs(y1))
            roots.append(times[index - 1] + fraction * (times[index] - times[index - 1]))
    return roots


def read_rows(path: Path) -> List[Tuple[float, ...]]:
    """Read capillary_wave.dat (statistics.f90).

    Columns: step, eta_cos, eta_sin, phase mass, relative mass error, k, omega_theory and, when
    present, d eta_cos/dt and the interface-normal velocity mode.
    """
    rows: List[Tuple[float, ...]] = []
    for line in path.read_text().splitlines():
        fields = line.split()
        if not fields or fields[0].startswith("#") or len(fields) < 7:
            continue
        try:
            rows.append(tuple(float(value) for value in fields))
        except ValueError:
            continue
    if len(rows) < 10:
        raise RuntimeError(f"Not enough capillary-wave samples in {path}")
    return rows


def read_force_projection(path: Path) -> Dict[str, float]:
    """Read step, numerical and theoretical Fourier coefficient of the CSF force, and their ratio."""
    for line in path.read_text().splitlines():
        fields = line.split()
        if not fields or fields[0].startswith("#") or len(fields) < 4:
            continue
        return {
            "force_mode": float(fields[1]),
            "force_theory": float(fields[2]),
            "force_ratio": float(fields[3]),
        }
    raise RuntimeError(f"No force-projection row found in {path}")


def analyse(data_path: Path, force_path: Path) -> Dict[str, float]:
    rows = read_rows(data_path)
    times = [row[0] for row in rows]
    signal = [row[1] for row in rows]
    roots = zero_crossings(times, signal)
    if len(roots) < 3:
        raise RuntimeError("At least three zero crossings are required")
    # Consecutive zero crossings are half a period apart: omega = pi / mean half period.
    half_periods = [roots[i] - roots[i - 1] for i in range(1, len(roots))]
    mean_half_period = sum(half_periods) / len(half_periods)
    omega = math.pi / mean_half_period
    # omega_theory = sqrt(sigma k^3 tanh(k h) / (rho_1 + rho_2)), computed by the solver.
    omega_theory = rows[-1][6]
    result = {
        "samples": float(len(rows)),
        "zero_crossings": float(len(roots)),
        "omega_theory": omega_theory,
        "omega_numerical": omega,
        "frequency_error_percent": 100.0 * (omega / omega_theory - 1.0),
        "period_numerical": 2.0 * mean_half_period,
        "max_mass_error": max(abs(row[4]) for row in rows),
    }
    # Kinematic check without the first sample, whose rate is zero: gain of the least-squares
    # fit d eta/dt = gain * velocity mode, and ||rate - velocity|| / ||velocity||.
    if len(rows[0]) >= 9:
        eta_rates = [row[7] for row in rows[1:]]
        velocity_modes = [row[8] for row in rows[1:]]
        velocity_norm = sum(value * value for value in velocity_modes)
        if velocity_norm > 0.0:
            result["kinematic_gain_deta_dt_over_vn"] = sum(
                rate * velocity for rate, velocity in zip(eta_rates, velocity_modes)
            ) / velocity_norm
            result["kinematic_relative_rms_error"] = math.sqrt(
                sum((rate - velocity) ** 2 for rate, velocity in zip(eta_rates, velocity_modes))
                / velocity_norm
            )
    result.update(read_force_projection(force_path))
    return result


def write_summary(path: Path, result: Dict[str, float]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(result))
        writer.writeheader()
        writer.writerow(result)


def main() -> None:
    root_default = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=root_default)
    parser.add_argument("--binary", type=Path, default=Path("main.x"))
    parser.add_argument("--defines", type=Path, default=Path("defines.h"))
    parser.add_argument("--input", type=Path, default=Path("capillary_wave.inp"))
    parser.add_argument("--output", type=Path, default=Path("capillary_wave_run"))
    parser.add_argument("--gpu-cc", default="80", help="NVIDIA compute capability")
    parser.add_argument("--make-target", default="nvfortran")
    add_mpi_arguments(parser)
    parser.add_argument("--skip-build", action="store_true", help="reuse an existing executable")
    parser.add_argument("--skip-run", action="store_true", help="analyse files already in --output")
    parser.add_argument("--force", action="store_true", help="overwrite an existing completed run")
    args = parser.parse_args()
    validate_mpi_arguments(args)

    root = args.root.resolve()
    binary = (root / args.binary).resolve()
    defines = (root / args.defines).resolve()
    input_path = (root / args.input).resolve()
    output = (root / args.output).resolve()
    data_path = output / "capillary_wave.dat"
    force_path = output / "capillary_force_projection.dat"

    if not args.skip_build:
        if not defines.is_file():
            raise SystemExit(f"Definitions file not found: {defines}")
        build_solver(root, defines, args.gpu_cc, args.make_target)
    if not binary.is_file() and not args.skip_run:
        raise SystemExit(f"Binary not found: {binary}")
    if not input_path.is_file() and not args.skip_run:
        raise SystemExit(f"Input file not found: {input_path}")
    if not args.skip_run:
        validate_local_tiles(input_path, args.decomposition)

    output.mkdir(parents=True, exist_ok=True)
    if not args.skip_run:
        if data_path.exists() and not args.force:
            raise SystemExit(f"Completed output already exists: {data_path}; use --force or --skip-run")
        local_input = output / input_path.name
        shutil.copy2(input_path, local_input)
        log_path = output / "run.log"
        print(f"Running planar capillary wave in {output}", flush=True)
        with log_path.open("w") as log:
            completed = subprocess.run(
                solver_command(args, binary, local_input.name),
                cwd=output,
                stdout=log,
                stderr=subprocess.STDOUT,
                check=False,
            )
        if completed.returncode != 0:
            raise SystemExit(f"Simulation failed; see {log_path}")

    result = analyse(data_path, force_path)
    summary_path = output / "capillary_wave_summary.csv"
    write_summary(summary_path, result)
    for key, value in result.items():
        print(f"{key}={value:.12e}")
    print(f"summary={summary_path}")


if __name__ == "__main__":
    main()
