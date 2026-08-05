#!/usr/bin/env python3
"""Compile and run short regression checks for documented LBFAST cases."""

import argparse
import csv
import math
import re
import shutil
import subprocess
import time
from pathlib import Path


CASE_MACROS = ("LAPLACE", "CAPILLARYWAVE", "LAMBTEST", "TAYLORGREEN",
               "POISEUILLE", "POISEUILLESTARTREST", "TWOPOISEUILLE")


def set_macro(text, name, value=None, enabled=True):
    replacement = "#define {}{}".format("" if enabled else "no", name)
    if enabled and value is not None:
        replacement += " " + value
    pattern = re.compile(
        r"(?m)^[ \t]*#define[ \t]+(?:no)?{}(?:[ \t]+[^\n]*)?[ \t]*$".format(
            re.escape(name)))
    updated, count = pattern.subn(replacement, text, count=1)
    return updated if count else text.rstrip() + "\n" + replacement + "\n"


def configure(original, selected):
    text = original
    required = {"LATTICE": "27", "HIGHORDER": None, "PRC": "8",
                "STRPRC": "8", "DOBENCHMARK": None}
    for name, value in required.items():
        text = set_macro(text, name, value)
    text = set_macro(text, "MIXEDPRC", enabled=False)
    text = set_macro(text, "USEGNUPLOT", enabled=False)
    for name in CASE_MACROS:
        text = set_macro(text, name, enabled=(name == selected))
    one_component = selected in ("TAYLORGREEN", "POISEUILLE")
    for name in ("TWOCOMPONENT", "DENSRATIO", "CSF"):
        text = set_macro(text, name, enabled=not one_component)
    text = set_macro(text, "PRINTPHI",
                     enabled=selected in ("LAPLACE", "CAPILLARYWAVE", "LAMBTEST"))
    if selected == "POISEUILLE":
        text = set_macro(text, "POISEUILLESTARTREST")
    return text


def rows(path, minimum_columns):
    result = []
    for line in path.read_text().splitlines():
        fields = line.split()
        if not fields or fields[0].startswith("#") or len(fields) < minimum_columns:
            continue
        try:
            values = [float(value.replace("D", "E").replace("d", "e"))
                      for value in fields]
        except ValueError:
            continue
        if not all(math.isfinite(value) for value in values):
            raise RuntimeError("non-finite value in " + path.name)
        result.append(values)
    if not result:
        raise RuntimeError("no numerical rows in " + path.name)
    return result


def check_laplace(directory):
    last = rows(directory / "laplace.dat", 4)[-1]
    if last[1] <= 0.0 or last[2] <= 0.0 or not 4.0 < last[3] < 12.0:
        raise RuntimeError("invalid pressure jump, surface tension, or radius")
    return "positive Laplace jump; sigma_eff={:.4g}".format(last[1]), last[1]


def check_capillary(directory):
    data = rows(directory / "capillary_wave.dat", 7)
    ratio = rows(directory / "capillary_force_projection.dat", 4)[0][3]
    mass_error = max(abs(row[4]) for row in data)
    if mass_error > 1.0e-4:
        raise RuntimeError("relative phase-mass error exceeds 1e-4")
    if not 0.5 < ratio < 1.5:
        raise RuntimeError("CSF projection ratio is outside [0.5, 1.5]")
    return "mass conserved; CSF projection={:.4g}".format(ratio), ratio


def check_lamb(directory):
    last = rows(directory / "lamb.dat", 7)[-1]
    if last[2] <= 0.0 or last[3] <= 0.0 or abs(last[1]) >= 1.0:
        raise RuntimeError("invalid droplet dimensions or deformation")
    return "finite interface; deformation={:.4g}".format(last[1]), last[1]


def check_taylor_green(directory):
    data = rows(directory / "taylorgreen.dat", 3)
    if len(data) < 3:
        raise RuntimeError("fewer than three Taylor-Green samples")
    if data[-1][2] >= data[0][2]:
        raise RuntimeError("kinetic energy did not decay")
    return "energy decays; log(E/E0)={:.4g}".format(data[-1][2]), data[-1][2]


def check_poiseuille(directory):
    data = rows(directory / "plot_poiseuille.dat", 3)
    velocity = [row[1] for row in data]
    maximum = max(velocity)
    if maximum <= 0.0:
        raise RuntimeError("body force did not generate positive flow")
    if abs(velocity[0]) > 1.0e-12 or abs(velocity[-1]) > 1.0e-12:
        raise RuntimeError("solid-wall velocities are not zero")
    symmetry = max(abs(velocity[i] - velocity[-1-i]) for i in range(len(velocity)))
    if symmetry > max(1.0e-10, 1.0e-5 * maximum):
        raise RuntimeError("velocity profile is not symmetric")
    return "force generates symmetric flow; umax={:.4g}".format(maximum), maximum


CASES = (
    ("laplace_csf", "LAPLACE", "laplace_smoke.inp", check_laplace),
    ("capillary_wave", "CAPILLARYWAVE", "capillary_wave_smoke.inp", check_capillary),
    ("lamb_oscillation", "LAMBTEST", "lamb_smoke.inp", check_lamb),
    ("taylor_green", "TAYLORGREEN", "taylorgreen_smoke.inp", check_taylor_green),
    ("poiseuille_force", "POISEUILLE", "poiseuille_force_smoke.inp", check_poiseuille),
)


def command(arguments, directory, log=None):
    print("+ " + " ".join(arguments), flush=True)
    if log is None:
        return subprocess.run(arguments, cwd=str(directory), check=False).returncode
    with log.open("w") as stream:
        return subprocess.run(arguments, cwd=str(directory), stdout=stream,
                              stderr=subprocess.STDOUT, check=False).returncode


def write_summary(path, results):
    fields = ("case", "status", "seconds", "metric", "detail")
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(results)


def main():
    tests = Path(__file__).resolve().parent
    root = tests.parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gpu-cc", default="80")
    parser.add_argument("--make-target", default="nvfortran")
    parser.add_argument("--case", action="append", choices=[case[0] for case in CASES])
    parser.add_argument("--keep-results", action="store_true")
    args = parser.parse_args()
    selected = [case for case in CASES if not args.case or case[0] in args.case]

    defines = root / "defines.h"
    binary = root / "main.x"
    output = tests / "results"
    if output.exists() and not args.keep_results:
        shutil.rmtree(str(output))
    output.mkdir(parents=True, exist_ok=True)
    original = defines.read_text()
    results = []
    try:
        for name, macro, input_name, checker in selected:
            start = time.monotonic()
            status, metric, detail = "FAIL", float("nan"), ""
            case_dir = output / name
            case_dir.mkdir(parents=True, exist_ok=True)
            local_input = case_dir / input_name
            shutil.copy2(str(tests / "inputs" / input_name), str(local_input))
            try:
                defines.write_text(configure(original, macro))
                if command(["make", "clean"], root):
                    raise RuntimeError("make clean failed")
                build_log = case_dir / "build.log"
                if command(["make", args.make_target, "GPUCC=" + args.gpu_cc],
                           root, build_log):
                    raise RuntimeError("build failed; see {}".format(build_log))
                log = case_dir / "run.log"
                code = command([str(binary), local_input.name], case_dir, log)
                if code:
                    raise RuntimeError("solver exited with status {}; see {}".format(code, log))
                if re.search(r"CUDA error|segmentation fault|floating invalid",
                             log.read_text(errors="replace"), re.I):
                    raise RuntimeError("runtime error marker found in run.log")
                detail, metric = checker(case_dir)
                status = "PASS"
            except Exception as error:
                detail = str(error)
            results.append({"case": name, "status": status,
                            "seconds": "{:.3f}".format(time.monotonic() - start),
                            "metric": "{:.12e}".format(metric) if math.isfinite(metric) else "",
                            "detail": detail})
            print("{:<5} {:<20} {}".format(status, name, detail), flush=True)
    finally:
        defines.write_text(original)
        command(["make", "clean"], root)
        if binary.exists():
            binary.unlink()
        write_summary(output / "smoke_test_results.csv", results)

    passed = sum(result["status"] == "PASS" for result in results)
    failed = len(results) - passed
    print("\n{} passed, {} failed".format(passed, failed))
    print("results={}".format(output / "smoke_test_results.csv"))
    raise SystemExit(0 if failed == 0 else 1)


if __name__ == "__main__":
    main()
