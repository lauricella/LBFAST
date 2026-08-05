#!/usr/bin/env python3
"""Configure, build, run, and summarize the LBFAST Lamb benchmark."""

import argparse
import csv
import re
import shutil
import subprocess
from pathlib import Path
from typing import Dict, Optional, Sequence


REQUIRED_DEFINES = {
    "LATTICE": "27",
    "HIGHORDER": None,
    "TWOCOMPONENT": None,
    "DENSRATIO": None,
    "CSF": None,
    "DOBENCHMARK": None,
    "PRC": "8",
    "STRPRC": "8",
    "LAMBTEST": None,
}
DISABLED_DEFINES = (
    "MIXEDPRC",
    "LAPLACE",
    "CAPILLARYWAVE",
    "POISEUILLE",
    "TWOPOISEUILLE",
    "TAYLORGREEN",
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
        print(f"Updated Lamb benchmark macros in {path}", flush=True)
    else:
        print(f"Lamb benchmark macros already configured in {path}", flush=True)


def run_command(command: Sequence[str], cwd: Path) -> None:
    print("+ " + " ".join(command), flush=True)
    completed = subprocess.run(command, cwd=cwd, check=False)
    if completed.returncode != 0:
        raise SystemExit(f"Command failed: {' '.join(command)}")


def build_solver(root: Path, defines: Path, gpu_cc: str, target: str) -> None:
    configure_defines(defines)
    run_command(["make", "clean"], root)
    run_command(["make", target, f"GPUCC={gpu_cc}"], root)


def replace_parameter(text: str, name: str, value: str) -> str:
    pattern = re.compile(rf"(?m)^(\s*{re.escape(name)}\s*=)[^!\n]*")
    updated, count = pattern.subn(rf"\g<1>{value}", text, count=1)
    if count != 1:
        raise RuntimeError(f"Could not replace {name!r} in the input")
    return updated


def parse_scalar(log: str, label: str) -> float:
    match = re.search(rf"(?m)^{re.escape(label)}\s*=*\s*([-+0-9.EeDd]+)", log)
    if not match:
        raise RuntimeError(f"Could not find {label!r} in run.log")
    return float(match.group(1).replace("D", "E").replace("d", "e"))


def parse_result(log_path: Path, input_path: Path) -> Dict[str, float]:
    log = log_path.read_text()
    input_text = input_path.read_text()

    def input_value(name: str) -> float:
        match = re.search(rf"(?m)^\s*{re.escape(name)}\s*=\s*([^!,\s]+)", input_text)
        if not match:
            raise RuntimeError(f"Could not read {name!r} from {input_path}")
        return float(match.group(1).replace("D", "E").replace("d", "e"))

    size = input_value("lx")
    radius = parse_scalar(log, "Equivalent radius")
    period_theory = parse_scalar(log, "Theoretical period")
    omega_theory = parse_scalar(log, "Theoretical omega")
    period_numerical = parse_scalar(log, "Mean numerical period T_num")
    omega_numerical = parse_scalar(log, "Angular frequency omega_num")
    peaks = parse_scalar(log, "Number of peaks found")
    width = input_value("width")
    return {
        "domain_x": size,
        "domain_y": input_value("ly"),
        "domain_z": input_value("lz"),
        "radius_equivalent": radius,
        "interface_width": width,
        "cahn_number": width / radius,
        "nsteps": input_value("nsteps"),
        "peaks": peaks,
        "period_theory": period_theory,
        "period_numerical": period_numerical,
        "period_error_percent": 100.0 * (period_numerical / period_theory - 1.0),
        "omega_theory": omega_theory,
        "omega_numerical": omega_numerical,
        "omega_error_percent": 100.0 * (omega_numerical / omega_theory - 1.0),
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
    parser.add_argument("--input", type=Path, default=Path("lamb.inp"))
    parser.add_argument("--output", type=Path, default=Path("lamb_validation_run"))
    parser.add_argument("--nsteps", type=int, default=None)
    parser.add_argument("--gpu-cc", default="80", help="NVIDIA compute capability")
    parser.add_argument("--make-target", default="nvfortran")
    parser.add_argument("--skip-build", action="store_true")
    parser.add_argument("--skip-run", action="store_true")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    root = args.root.resolve()
    binary = (root / args.binary).resolve()
    defines = (root / args.defines).resolve()
    source_input = (root / args.input).resolve()
    output = (root / args.output).resolve()
    local_input = output / source_input.name
    log_path = output / "run.log"
    data_path = output / "lamb.dat"

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
        if args.nsteps is not None:
            local_input.write_text(
                replace_parameter(local_input.read_text(), "nsteps", str(args.nsteps))
            )
        print(f"Running Lamb validation in {output}", flush=True)
        with log_path.open("w") as log:
            completed = subprocess.run(
                [str(binary), local_input.name],
                cwd=output,
                stdout=log,
                stderr=subprocess.STDOUT,
                check=False,
            )
        if completed.returncode != 0:
            raise SystemExit(f"Simulation failed; see {log_path}")

    if not local_input.is_file() or not log_path.is_file() or not data_path.is_file():
        raise SystemExit(f"Incomplete result set in {output}")
    result = parse_result(log_path, local_input)
    summary = output / "lamb_validation_summary.csv"
    write_csv(summary, result)
    for key, value in result.items():
        print(f"{key}={value:.12e}")
    print(f"summary={summary}")


if __name__ == "__main__":
    main()
