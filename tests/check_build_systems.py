#!/usr/bin/env python3
"""Check that the Makefile and CMakeLists.txt build the same LBFAST.

The current defines.h configuration is built with make, in the repository
root as the smoke suite does, and with CMake (macros from defines.h) in a
separate build directory. The check passes when

  1. both compile the same source files;
  2. every source and the link step get the same compiler and the same
     flags, apart from include paths and module, output and library options,
     and NVML builds use the same nvml.h and libnvidia-ml;
  3. with --input, both executables write byte-identical output files,
     including the final RAW fields and the double-precision restart files
     that the check enables in its copy of the input;
  4. a CMake configuration with LBFAST_USE_DEFINES_H=OFF and no macro given
     writes a lbfast_config.h identical to the user block of defines.h;
  5. with --input, the CMake executable has no stub library directory (CUDA
     lib/stubs) in its run-time path.

defines.h is not modified. The script runs make clean before and after its
build and removes main.x, so no executable inconsistent with defines.h is
left in the root.
"""

import argparse
import filecmp
import json
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path


# make target -> (LBFAST_MPI, LBFAST_NVML, CMAKE_BUILD_TYPE)
TARGETS = {
    "nvfortran": ("OFF", "OFF", "Release"),
    "nvfortran-mpi": ("ON", "OFF", "Release"),
    "nvfortran-nvml": ("OFF", "ON", "Release"),
    "nvfortran-nvml-mpi": ("ON", "ON", "Release"),
    "nvfortran-debug": ("OFF", "OFF", "Debug"),
    "nvfortran-mpi-debug": ("ON", "OFF", "Debug"),
}

# Options that only locate files: they differ between the two systems by
# construction and do not change the generated code.
DROP_WITH_ARGUMENT = {"-o", "-module", "-isystem", "-rpath", "-I", "-L"}
DROP_PREFIXES = ("-I", "-L", "-l", "-Wl,")
DROP_EXACT = {"-c", "-pthread"}
DROP_SUFFIXES = (".o", ".f90", ".c", ".so", ".a", ".x")


def run(arguments, cwd, log=None):
    print("+ " + " ".join(str(a) for a in arguments), flush=True)
    if log is None:
        return subprocess.run([str(a) for a in arguments], cwd=str(cwd)).returncode
    with Path(log).open("w") as stream:
        return subprocess.run([str(a) for a in arguments], cwd=str(cwd),
                              stdout=stream, stderr=subprocess.STDOUT).returncode


def real_compiler(name, cache={}):
    """Resolve an MPI wrapper (mpif90, mpicc) to the compiler it calls."""
    base = Path(name).name
    if base not in cache:
        resolved = base
        if base.startswith("mpi"):
            try:
                out = subprocess.run([base, "-show"], stdout=subprocess.PIPE,
                                     stderr=subprocess.DEVNULL,
                                     universal_newlines=True).stdout.split()
                if out:
                    resolved = Path(out[0]).name
            except OSError:
                pass
        cache[base] = resolved
    return cache[base]


def normalise(command, link=False):
    """Return (compiler, sorted significant flags) of a compiler command."""
    tokens = shlex.split(command)
    compiler = real_compiler(tokens[0])
    flags = []
    skip = False
    for token in tokens[1:]:
        if skip:
            skip = False
            continue
        if token in DROP_WITH_ARGUMENT:
            skip = True
            continue
        if token in DROP_EXACT or token.startswith(DROP_PREFIXES) \
                or token.endswith(DROP_SUFFIXES):
            continue
        if link and token.startswith("-D"):
            continue  # macros on the link line of the Makefile have no effect
        flags.append(token)
    return compiler, sorted(flags)


def source_of(command):
    for token in shlex.split(command):
        if token.endswith((".f90", ".c")):
            return Path(token).name
    return None


def cmake_cache(build, name):
    for line in (build / "CMakeCache.txt").read_text().splitlines():
        if line.startswith(name + ":"):
            return line.split("=", 1)[1]
    return ""


def make_variables(build, gpu_cc):
    """make variables matching a CMake build. NVML is left to the Makefile,
    which finds it in the HPC SDK as CMake does; the check verifies that both
    find the same files."""
    return ["GPUCC=" + gpu_cc]


def same_file(a, b):
    try:
        return Path(a).resolve() == Path(b).resolve()
    except OSError:
        return False


def make_commands(root, target, variables):
    """Compile commands per source and the link command of a make build."""
    text = subprocess.run(["make", "-n", target] + variables,
                          cwd=str(root), stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT, universal_newlines=True).stdout
    text = text.replace("\\\n", " ")
    compile_commands, link = {}, None
    for line in text.splitlines():
        words = line.split()
        if not words or real_compiler(words[0]) not in ("nvfortran", "nvc", "nvcc"):
            continue
        if " -c " in line:
            compile_commands[source_of(line)] = line
        elif " -o " in line:
            link = line
    return compile_commands, link


def cmake_commands(build):
    entries = json.loads((build / "compile_commands.json").read_text())
    compile_commands = {Path(e["file"]).name: e["command"] for e in entries}
    link = (build / "CMakeFiles" / "lbfast.dir" / "link.txt").read_text().strip()
    return compile_commands, link


def user_block(defines):
    lines, inside = [], False
    for line in defines.read_text().splitlines():
        if line.startswith("#ifndef LBFAST_CMAKE_CONFIG"):
            inside = True
        elif inside and line.startswith("#else"):
            break
        elif inside:
            lines.append(line)
    return lines


def generated_block(header):
    text = header.read_text()
    if text.startswith("/*"):
        text = text[text.index("*/") + 2:].lstrip("\n")
    return text.splitlines()


# Output files that hold timings, not results.
TIMING_FILES = {"time.dat"}


def with_full_output(text):
    """Input variant that also writes the final fields (RAW) and the
    double-precision restart files, so that the comparison is bitwise on
    the whole state and not only on the printed diagnostics."""
    steps = re.search(r"^\s*nsteps\s*=\s*(\d+)", text, re.I | re.M)
    replacements = [(r"lwriterestart\s*=\s*\.false\.", "lwriterestart=.true."),
                    (r"lraw\s*=\s*\.false\.", "lraw=.true.")]
    if steps:
        replacements.append((r"^(\s*)stamp\s*=\s*\d+", r"\g<1>stamp=" + steps.group(1)))
    for pattern, replacement in replacements:
        text = re.sub(pattern, replacement, text, flags=re.I | re.M)
    return text


def runpath(executable):
    """RPATH/RUNPATH entries of an executable (empty without readelf)."""
    try:
        text = subprocess.run(["readelf", "-d", str(executable)], stdout=subprocess.PIPE,
                              stderr=subprocess.DEVNULL, universal_newlines=True).stdout
    except OSError:
        return []
    paths = []
    for match in re.finditer(r"\((?:RPATH|RUNPATH)\).*\[(.*)\]", text):
        paths += [path for path in match.group(1).split(":") if path]
    return paths


def output_files(directory, skip):
    return sorted(p.relative_to(directory) for p in directory.rglob("*")
                  if p.is_file() and p.name not in skip)


def check(root, work, target, gpu_cc, input_file=None, jobs=8):
    """Run all checks; return a list of failure messages (empty: pass)."""
    if target not in TARGETS:
        return ["unsupported make target " + target]
    if shutil.which("cmake") is None:
        return ["cmake not found"]
    mpi, nvml, build_type = TARGETS[target]
    failures = []
    if work.exists():
        shutil.rmtree(str(work))
    work.mkdir(parents=True)
    cmake_dir = work / "cmake"
    options = ["-DGPUCC=" + gpu_cc, "-DLBFAST_MPI=" + mpi, "-DLBFAST_NVML=" + nvml,
               "-DCMAKE_BUILD_TYPE=" + build_type]

    # 1-2: sources and flags
    if run(["cmake", "-S", root, "-B", cmake_dir] + options, root, work / "cmake_configure.log"):
        return ["cmake configuration failed; see {}".format(work / "cmake_configure.log")]
    variables = make_variables(cmake_dir, gpu_cc)
    if run(["make", "clean"], root, work / "make_clean.log"):
        return ["make clean failed"]
    make_compile, make_link = make_commands(root, target, variables)
    if nvml == "ON":
        make_nvml = subprocess.run(
            ["make", "-s", "-f", "-", "print"] + variables, cwd=str(root),
            input="include Makefile\nprint:\n\t@echo $(NVML_INC) $(NVML_LIB)\n",
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
            universal_newlines=True).stdout.split()
        cmake_nvml = [cmake_cache(cmake_dir, "NVML_INCLUDE_DIR"),
                      str(Path(cmake_cache(cmake_dir, "NVML_LIBRARY")).parent)]
        if len(make_nvml) != 2 or not all(map(same_file, make_nvml, cmake_nvml)):
            failures.append("NVML: make uses {}, cmake {}".format(make_nvml, cmake_nvml))
    cmake_compile, cmake_link = cmake_commands(cmake_dir)
    if set(make_compile) != set(cmake_compile):
        failures.append("different sources: make only {}, cmake only {}".format(
            sorted(set(make_compile) - set(cmake_compile)),
            sorted(set(cmake_compile) - set(make_compile))))
    for source in sorted(set(make_compile) & set(cmake_compile)):
        mc, mf = normalise(make_compile[source])
        cc, cf = normalise(cmake_compile[source])
        if source.endswith(".c") and {mc, cc} == {"nvcc", "nvc"}:
            mc = cc  # serial make builds compile the C helpers with nvcc
        if (mc, mf) != (cc, cf):
            failures.append("{}: make {} {} | cmake {} {}".format(
                source, mc, " ".join(mf), cc, " ".join(cf)))
    if make_link is None:
        failures.append("no link command in make -n output")
    elif normalise(make_link, link=True) != normalise(cmake_link, link=True):
        failures.append("link: make {} | cmake {}".format(
            normalise(make_link, link=True), normalise(cmake_link, link=True)))

    # 4: lbfast_config.h written in cache mode equals the defines.h block
    cache_dir = work / "cmake_cache_mode"
    if run(["cmake", "-S", root, "-B", cache_dir, "-DLBFAST_USE_DEFINES_H=OFF"] + options,
           root, work / "cmake_cache_mode.log"):
        failures.append("cmake configuration with LBFAST_USE_DEFINES_H=OFF failed")
    elif generated_block(cache_dir / "lbfast_config.h") != user_block(root / "defines.h"):
        failures.append("lbfast_config.h differs from the user block of defines.h")

    # 3: identical output of the two executables
    if input_file is not None:
        binary = root / "main.x"
        try:
            if run(["make", "clean"], root, work / "make_clean.log"):
                raise RuntimeError("make clean failed")
            if run(["make", target] + variables, root, work / "make_build.log"):
                raise RuntimeError("make build failed; see {}".format(work / "make_build.log"))
            shutil.copy2(str(binary), str(work / "main_make.x"))
        except RuntimeError as error:
            return failures + [str(error)]
        finally:
            run(["make", "clean"], root, work / "make_clean.log")
            if binary.exists():
                binary.unlink()
        if run(["cmake", "--build", cmake_dir, "-j", str(jobs)], root, work / "cmake_build.log"):
            return failures + ["cmake build failed; see {}".format(work / "cmake_build.log")]
        runs = {"make": work / "main_make.x", "cmake": cmake_dir / "main.x"}
        stubs = [path for path in runpath(runs["cmake"]) if path.rstrip("/").endswith("stubs")]
        if stubs:
            failures.append("the CMake executable has a stub library directory in its "
                            "RPATH: {}".format(stubs))
        for name, executable in runs.items():
            directory = work / ("run_" + name)
            directory.mkdir()
            (directory / input_file.name).write_text(with_full_output(input_file.read_text()))
            prefix = (["mpirun", "-np", "1", "--bind-to", "none"] if mpi == "ON" else [])
            arguments = [executable] + (["1", "1", "1"] if mpi == "ON" else []) + [input_file.name]
            if run(prefix + arguments, directory, directory / "run.log"):
                return failures + ["{} executable failed; see {}".format(
                    name, directory / "run.log")]
        skip = {"run.log", input_file.name} | TIMING_FILES
        make_files = output_files(work / "run_make", skip)
        cmake_files = output_files(work / "run_cmake", skip)
        if not make_files:
            failures.append("the run wrote no output file to compare")
        if make_files != cmake_files:
            failures.append("different output files: {} vs {}".format(
                [str(p) for p in make_files], [str(p) for p in cmake_files]))
        for relative in sorted(set(make_files) & set(cmake_files)):
            if not filecmp.cmp(str(work / "run_make" / relative),
                               str(work / "run_cmake" / relative), shallow=False):
                failures.append("output {} differs".format(relative))
    return failures


def main():
    root = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--gpu-cc", default="80")
    parser.add_argument("--make-target", default="nvfortran", choices=sorted(TARGETS))
    parser.add_argument("--input", type=Path,
                        help="input file run with both executables (compiled for defines.h)")
    parser.add_argument("--work", type=Path, default=root / "tests" / "results" / "build_systems")
    parser.add_argument("--jobs", type=int, default=8)
    args = parser.parse_args()
    input_file = args.input.resolve() if args.input else None
    failures = check(root, args.work.resolve(), args.make_target, args.gpu_cc, input_file, args.jobs)
    for failure in failures:
        print("FAIL  " + failure)
    print("PASS  make and cmake agree" if not failures else
          "{} difference(s) between make and cmake".format(len(failures)))
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
