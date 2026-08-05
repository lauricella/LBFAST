#!/usr/bin/env python3
"""Run a CSF Laplace sweep and fit Delta p as a function of 1/R."""

import argparse
import csv
import re
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple


RADII = (16, 24, 32)
REQUIRED_DEFINES = {
    "LATTICE": "27",
    "HIGHORDER": None,
    "TWOCOMPONENT": None,
    "DENSRATIO": None,
    "CSF": None,
    "DOBENCHMARK": None,
    "LAPLACE": None,
}
DISABLED_CASES = ("LAMBTEST", "POISEUILLE", "TWOPOISEUILLE", "TAYLORGREEN")


def set_define(text: str, name: str, value: str = None) -> str:
    """Enable one top-level option while preserving the defines.h layout."""
    replacement = f"#define {name}" + (f" {value}" if value is not None else "")
    pattern = re.compile(
        rf"(?m)^[ \t]*#define[ \t]+(?:no)?{re.escape(name)}"
        rf"(?:[ \t]+[^\n]*)?[ \t]*$"
    )
    updated, count = pattern.subn(replacement, text, count=1)
    if count:
        return updated
    return text.rstrip() + f"\n{replacement}\n"


def disable_define(text: str, name: str) -> str:
    """Disable an active top-level benchmark selector if it is present."""
    pattern = re.compile(
        rf"(?m)^[ \t]*#define[ \t]+{re.escape(name)}[ \t]*$"
    )
    return pattern.sub(f"#define no{name}", text, count=1)


def configure_defines(path: Path) -> None:
    text = path.read_text()
    configured = text
    for name, value in REQUIRED_DEFINES.items():
        configured = set_define(configured, name, value)
    for name in DISABLED_CASES:
        configured = disable_define(configured, name)
    if configured != text:
        path.write_text(configured)
        print(f"Updated CSF Laplace macros in {path}", flush=True)
    else:
        print(f"CSF Laplace macros already configured in {path}", flush=True)


def build_solver(root: Path, defines: Path, gpu_cc: str, target: str) -> None:
    configure_defines(defines)
    commands = (["make", "clean"], ["make", target, f"GPUCC={gpu_cc}"])
    for command in commands:
        print("+ " + " ".join(command), flush=True)
        completed = subprocess.run(command, cwd=root, check=False)
        if completed.returncode != 0:
            raise SystemExit(f"Build command failed: {' '.join(command)}")


def replace_parameter(text: str, name: str, value: str) -> str:
    pattern = re.compile(rf"(?m)^(\s*{re.escape(name)}\s*=)[^!\n]*")
    updated, count = pattern.subn(rf"\g<1>{value}", text, count=1)
    if count != 1:
        raise RuntimeError(f"Could not replace {name!r} in the input template")
    return updated


def build_input(template: str, radius: int, nsteps: int) -> str:
    size = 4 * radius
    center = size / 2
    values = {
        "nsteps": str(nsteps),
        "lx": str(size),
        "ly": str(size),
        "lz": str(size),
        "iprobe": "1",
        "jprobe": str(size // 2),
        "kprobe": str(size // 2),
        "radius": str(radius),
        "width": "4.0e0",
        "center": f"{center:.1f}, {center:.1f}, {center:.1f}",
    }
    result = template
    for name, value in values.items():
        result = replace_parameter(result, name, value)
    return result


def read_plateau(path: Path, start_fraction: float) -> Dict[str, float]:
    rows: List[Tuple[int, float, float, float]] = []
    for line in path.read_text().splitlines():
        fields = line.split()
        if not fields or fields[0].startswith("#") or len(fields) < 4:
            continue
        try:
            rows.append((int(fields[0]), float(fields[1]), float(fields[2]), float(fields[3])))
        except ValueError:
            continue
    if not rows:
        raise RuntimeError(f"No numerical rows found in {path}")
    final_step = max(row[0] for row in rows)
    plateau = [row for row in rows if row[0] >= start_fraction * final_step]
    if len(plateau) < 3:
        raise RuntimeError(f"Too few plateau samples in {path}")

    def mean(column: int) -> float:
        return sum(row[column] for row in plateau) / len(plateau)

    return {
        "samples": float(len(plateau)),
        "sigma_eff": mean(1),
        "delta_p": mean(2),
        "radius_measured": mean(3),
    }


def fit_results(rows: List[Dict[str, float]], sigma: float) -> Dict[str, float]:
    x = [1.0 / row["radius"] for row in rows]
    y = [row["delta_p"] for row in rows]
    slope_zero = sum(a * b for a, b in zip(x, y)) / sum(a * a for a in x)
    xmean = sum(x) / len(x)
    ymean = sum(y) / len(y)
    slope = sum((a - xmean) * (b - ymean) for a, b in zip(x, y)) / sum(
        (a - xmean) ** 2 for a in x
    )
    intercept = ymean - slope * xmean
    return {
        "slope_through_origin": slope_zero,
        "slope_error": slope_zero / (2.0 * sigma) - 1.0,
        "sigma_from_slope": 0.5 * slope_zero,
        "slope_with_intercept": slope,
        "intercept": intercept,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, default=Path("main.x"))
    parser.add_argument("--defines", type=Path, default=Path("defines.h"))
    parser.add_argument("--template", type=Path, default=Path("laplace.inp"))
    parser.add_argument("--output", type=Path, default=Path("laplace_csf_sweep"))
    parser.add_argument("--nsteps", type=int, default=10_000)
    parser.add_argument("--plateau-start", type=float, default=0.5)
    parser.add_argument("--gpu-cc", default="80", help="NVIDIA compute capability")
    parser.add_argument("--make-target", default="nvfortran")
    parser.add_argument(
        "--skip-build", action="store_true", help="reuse an existing executable"
    )
    parser.add_argument("--force", action="store_true", help="rerun completed cases")
    args = parser.parse_args()

    root = Path.cwd()
    binary = (root / args.binary).resolve()
    defines = (root / args.defines).resolve()
    template_path = (root / args.template).resolve()
    output = (root / args.output).resolve()
    if not args.skip_build:
        if not defines.is_file():
            raise SystemExit(f"Definitions file not found: {defines}")
        build_solver(root, defines, args.gpu_cc, args.make_target)
    if not binary.is_file():
        raise SystemExit(f"Binary not found: {binary}")
    if not 0.0 <= args.plateau_start < 1.0:
        raise SystemExit("--plateau-start must be in [0,1)")

    template = template_path.read_text()
    sigma_match = re.search(r"(?m)^\s*sigma\s*=\s*([^!\s]+)", template)
    if not sigma_match:
        raise SystemExit("Could not read sigma from the input template")
    sigma = float(sigma_match.group(1).replace("d", "e").replace("D", "E"))
    output.mkdir(parents=True, exist_ok=True)
    results: List[Dict[str, float]] = []

    for radius in RADII:
        size = 4 * radius
        case_dir = output / f"R{radius}_N{size}_W4"
        case_dir.mkdir(parents=True, exist_ok=True)
        input_path = case_dir / "laplace.inp"
        data_path = case_dir / "laplace.dat"
        log_path = case_dir / "run.log"
        input_path.write_text(build_input(template, radius, args.nsteps))

        if args.force or not data_path.is_file():
            print(f"Running R={radius}, N={size}, W=4 ...", flush=True)
            with log_path.open("w") as log:
                completed = subprocess.run(
                    [str(binary), input_path.name],
                    cwd=case_dir,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    check=False,
                )
            if completed.returncode != 0:
                raise SystemExit(f"Case R={radius} failed; see {log_path}")
        else:
            print(f"Reusing completed R={radius}, N={size}, W=4", flush=True)

        row = read_plateau(data_path, args.plateau_start)
        row.update(
            {
                "radius": float(radius),
                "size": float(size),
                "delta_p_theory": 2.0 * sigma / radius,
            }
        )
        row["relative_error"] = row["delta_p"] / row["delta_p_theory"] - 1.0
        results.append(row)

    fit = fit_results(results, sigma)
    summary_path = output / "laplace_csf_sweep_summary.csv"
    with summary_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(results[0]))
        writer.writeheader()
        writer.writerows(results)

    fit_path = output / "laplace_csf_fit.txt"
    fit_path.write_text("".join(f"{key}={value:.12e}\n" for key, value in fit.items()))

    print("\nR      N       delta_p_num    delta_p_exact  error[%]   sigma_eff")
    for row in results:
        print(
            f"{row['radius']:4.0f} {row['size']:6.0f}  {row['delta_p']:.9e}  "
            f"{row['delta_p_theory']:.9e}  {100*row['relative_error']:8.4f}  "
            f"{row['sigma_eff']:.9e}"
        )
    print(f"\nThrough-origin slope : {fit['slope_through_origin']:.12e}")
    print(f"Theoretical slope    : {2*sigma:.12e}")
    print(f"Slope error          : {100*fit['slope_error']:.6f}%")
    print(f"Sigma from slope     : {fit['sigma_from_slope']:.12e}")
    print(f"Free-fit intercept   : {fit['intercept']:.12e}")
    print(f"Summary              : {summary_path}")


if __name__ == "__main__":
    main()
