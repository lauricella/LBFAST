#!/usr/bin/env python3
"""Configure, build, run, and analyse the LBFAST Taylor-Green test."""

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
    "TAYLORGREEN": None,
}
DISABLED_DEFINES = (
    "MIXEDPRC",
    "TWOCOMPONENT",
    "DENSRATIO",
    "CSF",
    "LAMBTEST",
    "LAPLACE",
    "CAPILLARYWAVE",
    "POISEUILLE",
    "TWOPOISEUILLE",
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
        print(f"Updated single-component Taylor-Green macros in {path}", flush=True)
    else:
        print(f"Taylor-Green macros already configured in {path}", flush=True)


def run_command(command: Sequence[str], cwd: Path) -> None:
    print("+ " + " ".join(command), flush=True)
    completed = subprocess.run(command, cwd=cwd, check=False)
    if completed.returncode != 0:
        raise SystemExit(f"Command failed: {' '.join(command)}")


def build_solver(root: Path, defines: Path, gpu_cc: str, target: str) -> None:
    configure_defines(defines)
    run_command(["make", "clean"], root)
    run_command(["make", target, f"GPUCC={gpu_cc}"], root)


def read_input_value(text: str, name: str) -> float:
    match = re.search(rf"(?m)^\s*{re.escape(name)}\s*=\s*([^!,\s]+)", text)
    if not match:
        raise RuntimeError(f"Could not read {name!r} from the input")
    return float(match.group(1).replace("D", "E").replace("d", "e"))


def read_energy(path: Path) -> List[Tuple[float, float]]:
    rows: List[Tuple[float, float]] = []
    for line in path.read_text().splitlines():
        fields = line.split()
        if len(fields) < 3:
            continue
        try:
            rows.append((float(fields[0]), float(fields[2])))
        except ValueError:
            continue
    if len(rows) < 3:
        raise RuntimeError(f"Too few energy samples in {path}")
    return rows


def linear_fit(rows: Sequence[Tuple[float, float]]) -> Tuple[float, float]:
    n = float(len(rows))
    sx = sum(row[0] for row in rows)
    sy = sum(row[1] for row in rows)
    sxx = sum(row[0] * row[0] for row in rows)
    sxy = sum(row[0] * row[1] for row in rows)
    denominator = n * sxx - sx * sx
    if denominator == 0.0:
        raise RuntimeError("Degenerate Taylor-Green fit window")
    slope = (n * sxy - sx * sy) / denominator
    intercept = (sy - slope * sx) / n
    return slope, intercept


def analyse(data_path: Path, input_path: Path, max_log_decay: float) -> Dict[str, float]:
    input_text = input_path.read_text()
    size_x = read_input_value(input_text, "lx")
    viscosity = read_input_value(input_text, "visc1")
    rows = read_energy(data_path)
    fit_rows = [row for row in rows if row[1] >= max_log_decay]
    if len(fit_rows) < 3:
        raise RuntimeError("Too few samples in the requested fit window")
    slope, intercept = linear_fit(fit_rows)
    wave_number = 2.0 * math.pi / size_x
    viscosity_numerical = -slope / (6.0 * wave_number * wave_number)
    fitted = [slope * time + intercept for time, _ in fit_rows]
    rms = math.sqrt(
        sum((row[1] - prediction) ** 2 for row, prediction in zip(fit_rows, fitted))
        / len(fit_rows)
    )
    return {
        "domain_x": size_x,
        "domain_y": read_input_value(input_text, "ly"),
        "domain_z": read_input_value(input_text, "lz"),
        "nsteps": read_input_value(input_text, "nsteps"),
        "initial_velocity": read_input_value(input_text, "uwall"),
        "viscosity_theory": viscosity,
        "viscosity_numerical": viscosity_numerical,
        "viscosity_error_percent": 100.0 * (viscosity_numerical / viscosity - 1.0),
        "fit_samples": float(len(fit_rows)),
        "fit_last_step": fit_rows[-1][0],
        "fit_slope": slope,
        "fit_intercept": intercept,
        "fit_log_energy_rms": rms,
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
    parser.add_argument("--input", type=Path, default=Path("taylorgreen128.inp"))
    parser.add_argument("--output", type=Path, default=Path("taylorgreen_validation_run"))
    parser.add_argument("--gpu-cc", default="80", help="NVIDIA compute capability")
    parser.add_argument("--make-target", default="nvfortran")
    add_mpi_arguments(parser)
    parser.add_argument(
        "--max-log-decay",
        type=float,
        default=-0.25,
        help="fit samples satisfying log(E/E0) >= this value",
    )
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
    data_path = output / "taylorgreen.dat"

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
        if data_path.exists() and not args.force:
            raise SystemExit(f"Output already exists: {data_path}; use --force or --skip-run")
        shutil.copy2(source_input, local_input)
        validate_local_tiles(local_input, args.decomposition)
        print(f"Running Taylor-Green validation in {output}", flush=True)
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

    if not local_input.is_file() or not log_path.is_file() or not data_path.is_file():
        raise SystemExit(f"Incomplete result set in {output}")
    result = analyse(data_path, local_input, args.max_log_decay)
    summary = output / "taylorgreen_validation_summary.csv"
    write_csv(summary, result)
    for key, value in result.items():
        print(f"{key}={value:.12e}")
    print(f"summary={summary}")


if __name__ == "__main__":
    main()
