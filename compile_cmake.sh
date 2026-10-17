#!/usr/bin/env bash
# Build the executable matrix of compile.sh with CMake, without editing
# defines.h (docs/CMAKE.md).
#
# Each variant is a preset of CMakePresets.json, configured in
# build/<preset>; its executable is written to build/bin/main_<preset>.x.
# Macros that the presets do not set (DOBENCHMARK, test-case selectors, ...)
# take their value from defines.h when a build directory is first configured.
#
# usage: ./compile_cmake.sh [-p PARALLEL] [-j JOBS] [preset ...] [-- cmake options]
#   -p PARALLEL  variants built at the same time (default 4)
#   -j JOBS      compile jobs per variant (default 4)
#   preset ...   presets to build (default: all 24 matrix presets)
#   -- options   passed to every configure step, e.g. -- -DGPUCC=90
set -euo pipefail
cd "$(dirname "$0")"

parallel=4
jobs=4
presets=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    -p) parallel="$2"; shift 2 ;;
    -j) jobs="$2"; shift 2 ;;
    -h|--help) sed -n '2,15p' "$0"; exit 0 ;;
    --) shift; break ;;
    *) presets+=("$1"); shift ;;
  esac
done
extra=("$@")

if [[ ${#presets[@]} -eq 0 ]]; then
  # every visible configure preset except "default"
  mapfile -t presets < <(cmake --list-presets=configure | sed -n 's/^ *"\([^"]*\)".*/\1/p' | grep -v '^default$')
fi

mkdir -p build/logs build/bin
build_one() {
  local preset="$1"; shift
  local log="build/logs/${preset}.log"
  if cmake --preset "$preset" "$@" > "$log" 2>&1 &&
     cmake --build --preset "$preset" -j "$jobs" >> "$log" 2>&1; then
    echo "OK    $preset"
  else
    echo "FAIL  $preset (see $log)"
    return 1
  fi
}
export -f build_one
export jobs

start=$SECONDS
status=0
printf '%s\n' "${presets[@]}" |
  xargs -P "$parallel" -I{} bash -c 'build_one "$@"' _ {} "${extra[@]}" || status=$?
echo "${#presets[@]} variants in $((SECONDS - start)) s; executables in build/bin"
exit $status
