#!/usr/bin/env python3
"""Build, run, and analyse the single-component body-force validation."""

import argparse
import csv
import math
import re
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple
from validation_mpi import add_mpi_arguments, solver_command, validate_local_tiles, validate_mpi_arguments


REQUIRED_DEFINES = {
    "LATTICE": "27",
    "HIGHORDER": None,
    "PRC": "8",
    "STRPRC": "8",
    "DOBENCHMARK": None,
    "POISEUILLE": None,
    "POISEUILLESTARTREST": None,
}
DISABLED_DEFINES = (
    "MIXEDPRC",
    "TWOCOMPONENT",
    "DENSRATIO",
    "CSF",
    "TWOPOISEUILLE",
    "TAYLORGREEN",
    "LAMBTEST",
    "LAPLACE",
    "CAPILLARYWAVE",
)


def set_macro(text: str, name: str, value: Optional[str], enabled: bool) -> str:
    prefix = "" if enabled else "no"
    replacement = f"#define {prefix}{name}"
    if enabled and value is not None:
        replacement += f" {value}"
    pattern = re.compile(
        rf"(?m)^[ \t]*#define[ \t]+(?:no)?{re.escape(name)}"
        rf"(?:[ \t]+[^\n]*)?[ \t]*$"
    )
    updated, count = pattern.subn(replacement, text, count=1)
    if count:
        return updated
    return text.rstrip() + "\n" + replacement + "\n"


def configure_defines(path: Path) -> None:
    original = path.read_text()
    configured = original
    for name, value in REQUIRED_DEFINES.items():
        configured = set_macro(configured, name, value, True)
    for name in DISABLED_DEFINES:
        configured = set_macro(configured, name, None, False)
    if configured != original:
        path.write_text(configured)
        print(f"Updated single-component force-validation macros in {path}", flush=True)
    else:
        print(f"Force-validation macros already configured in {path}", flush=True)


def run_command(command: Sequence[str], cwd: Path) -> None:
    print("+ " + " ".join(command), flush=True)
    completed = subprocess.run(command, cwd=cwd, check=False)
    if completed.returncode != 0:
        raise SystemExit(f"Command failed: {' '.join(command)}")


def build_solver(root: Path, defines: Path, gpu_cc: str, target: str) -> None:
    configure_defines(defines)
    run_command(["make", "clean"], root)
    run_command(["make", target, f"GPUCC={gpu_cc}"], root)


def read_profile(path: Path) -> List[Tuple[float, float, float]]:
    rows: List[Tuple[float, float, float]] = []
    for line in path.read_text().splitlines():
        fields = line.split()
        if len(fields) < 3:
            continue
        try:
            rows.append((float(fields[0]), float(fields[1]), float(fields[2])))
        except ValueError:
            continue
    if len(rows) < 3:
        raise RuntimeError(f"Too few profile points in {path}")
    return rows


def read_input_value(text: str, name: str) -> float:
    match = re.search(rf"(?m)^\s*{re.escape(name)}\s*=\s*([^!,\s]+)", text)
    if not match:
        raise RuntimeError(f"Could not read {name!r} from the input")
    return float(match.group(1).replace("D", "E").replace("d", "e"))


def analyse(path: Path, input_path: Path) -> Dict[str, float]:
    rows = read_profile(path)
    input_text = input_path.read_text()
    lx = read_input_value(input_text, "lx")
    force_z = read_input_value(input_text, "fz")
    viscosity = read_input_value(input_text, "visc1")

    # The regularized boundary reconstruction used by this benchmark has its
    # hydrodynamic no-slip planes at the two solid-node coordinates.  This is
    # the non-BOUNCE_BACK convention also used by statistics.f90.
    wall_lower = 1.0
    wall_upper = lx
    center = 0.5 * (wall_lower + wall_upper)
    half_width = 0.5 * (wall_upper - wall_lower)

    analytical_rows: List[Tuple[float, float, float]] = []
    statistics_differences: List[float] = []
    for position, numerical, statistics_exact in rows:
        distance = position - center
        if abs(distance) <= half_width:
            exact = force_z * (half_width**2 - distance**2) / (2.0 * viscosity)
        else:
            exact = 0.0
        analytical_rows.append((position, numerical, exact))
        statistics_differences.append(statistics_exact - exact)

    residuals = [numerical - exact for _, numerical, exact in analytical_rows]
    exact_norm = math.sqrt(sum(exact * exact for _, _, exact in analytical_rows))
    if exact_norm == 0.0:
        raise RuntimeError("Analytical profile has zero norm")
    relative_l2 = math.sqrt(sum(value * value for value in residuals)) / exact_norm
    scale = max(abs(exact) for _, _, exact in analytical_rows)
    max_error = max(abs(value) for value in residuals)
    center_row = max(analytical_rows, key=lambda row: row[2])
    return {
        "domain_x": lx,
        "domain_y": read_input_value(input_text, "ly"),
        "domain_z": read_input_value(input_text, "lz"),
        "nsteps": read_input_value(input_text, "nsteps"),
        "force_z": force_z,
        "kinematic_viscosity": viscosity,
        "wall_lower": wall_lower,
        "wall_upper": wall_upper,
        "channel_center": center,
        "channel_half_width": half_width,
        "profile_points": float(len(rows)),
        "center_position": center_row[0],
        "center_velocity_theory": center_row[2],
        "center_velocity_numerical": center_row[1],
        "center_velocity_error_percent": 100.0 * (center_row[1] / center_row[2] - 1.0),
        "relative_l2_error": relative_l2,
        "relative_max_error": max_error / scale,
        "absolute_max_error": max_error,
        "statistics_analytical_max_difference": max(
            abs(value) for value in statistics_differences
        ),
    }


def write_csv(path: Path, result: Dict[str, float]) -> None:
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
    parser.add_argument("--input", type=Path, default=Path("poiseuille_force.inp"))
    parser.add_argument("--output", type=Path, default=Path("poiseuille_force_validation_run"))
    parser.add_argument("--gpu-cc", default="80", help="NVIDIA compute capability")
    parser.add_argument("--make-target", default="nvfortran")
    add_mpi_arguments(parser)
    parser.add_argument("--skip-build", action="store_true")
    parser.add_argument("--skip-run", action="store_true")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    validate_mpi_arguments(args)

    root = args.root.resolve()
    binary = (root / args.binary).resolve()
    defines = (root / args.defines).resolve()
    source_input = (root / args.input).resolve()
    output = (root / args.output).resolve()
    local_input = output / source_input.name
    log_path = output / "run.log"
    profile_path = output / "plot_poiseuille.dat"

    if not args.skip_build:
        if not defines.is_file():
            raise SystemExit(f"Definitions file not found: {defines}")
        build_solver(root, defines, args.gpu_cc, args.make_target)
    if not args.skip_run and not binary.is_file():
        raise SystemExit(f"Binary not found: {binary}")
    if not args.skip_run and not source_input.is_file():
        raise SystemExit(f"Input file not found: {source_input}")

    output.mkdir(parents=True, exist_ok=True)
    if not args.skip_run:
        if profile_path.exists() and not args.force:
            raise SystemExit(f"Output already exists: {profile_path}; use --force or --skip-run")
        shutil.copy2(source_input, local_input)
        validate_local_tiles(local_input, args.decomposition)
        print(f"Running body-force validation in {output}", flush=True)
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

    if not local_input.is_file() or not log_path.is_file() or not profile_path.is_file():
        raise SystemExit(f"Incomplete result set in {output}")
    result = analyse(profile_path, local_input)
    summary = output / "poiseuille_force_validation_summary.csv"
    write_csv(summary, result)
    for key, value in result.items():
        print(f"{key}={value:.12e}")
    print(f"summary={summary}")


if __name__ == "__main__":
    main()
