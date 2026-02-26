#!/bin/bash
#set compute capability
export MYCC=80

#set where is nvml.h
export NVML_INC=$NVHPC_HOME/Linux_x86_64/24.3/cuda/12.3/targets/x86_64-linux/include 

#set where is libnvidia-ml
export NVML_LIB=$NVHPC_HOME/Linux_x86_64/24.3/cuda/12.3/targets/x86_64-linux/lib/stubs

#!/usr/bin/env bash
set -euo pipefail

FILE="defines.h"

# Expected defines (regex, but printed as plain text in the warning)
expected_lines=(
"#define HIGHORDER"
"#define LATTICE 27"
"#define TWOCOMPONENT"
"#define DENSRATIO"
"#define INTERFACE_INCOMP"
"#define PRC 8"
"#define noMIXEDPRC"
"#define STRPRC 8"
)

# Corresponding regex patterns (allow flexible whitespace)
patterns=(
'^#define[[:space:]]+HIGHORDER[[:space:]]*$'
'^#define[[:space:]]+LATTICE[[:space:]]+27[[:space:]]*$'
'^#define[[:space:]]+TWOCOMPONENT[[:space:]]*$'
'^#define[[:space:]]+DENSRATIO[[:space:]]*$'
'^#define[[:space:]]+INTERFACE_INCOMP[[:space:]]*$'
'^#define[[:space:]]+PRC[[:space:]]+8[[:space:]]*$'
'^#define[[:space:]]+noMIXEDPRC[[:space:]]*$'
'^#define[[:space:]]+STRPRC[[:space:]]+8[[:space:]]*$'
)

if [[ ! -f "$FILE" ]]; then
  echo "ERROR: '$FILE' not found." >&2
  exit 2
fi

missing=()

for i in "${!patterns[@]}"; do
  if ! grep -Eq "${patterns[$i]}" "$FILE"; then
    missing+=("${expected_lines[$i]}")
  fi
done

if ((${#missing[@]} > 0)); then
  cat >&2 <<'EOF'
================================================================================
WARNING: defines.h does not match the expected build configuration.
The file must contain ALL of the following lines exactly (whitespace may vary):
EOF
  for line in "${expected_lines[@]}"; do
    echo "  $line" >&2
  done

  echo >&2
  echo "Missing or mismatched lines detected:" >&2
  for line in "${missing[@]}"; do
    echo "  - $line" >&2
  done

  cat >&2 <<'EOF'
Fix defines.h and re-run.
================================================================================
EOF
  exit 1
fi

echo "OK: defines.h matches the expected configuration."

echo "compiling with compute capability: $MYCC"

make clean
make nvfortran-nvml-mpi GPUCC=$MYCC NVML_INC=$NVML_INC NNVML_LIB=$NVML_LIB
mv main.x main_2c_27high_d.x
make clean
sed -i 's/^#define[[:space:]]\+HIGHORDER/#define noHIGHORDER/' defines.h
make nvfortran-nvml-mpi GPUCC=$MYCC NVML_INC=$NVML_INC NNVML_LIB=$NVML_LIB
mv main.x main_2c_27_d.x
make clean
sed -i 's/^#define[[:space:]]\+LATTICE[[:space:]]\+27/#define LATTICE 19/' defines.h
make nvfortran-nvml-mpi GPUCC=$MYCC NVML_INC=$NVML_INC NNVML_LIB=$NVML_LIB
mv main.x main_2c_19_d.x
make clean
sed -i 's/^#define[[:space:]]\+LATTICE[[:space:]]\+19/#define LATTICE 15/' defines.h
make nvfortran-nvml-mpi GPUCC=$MYCC NVML_INC=$NVML_INC NNVML_LIB=$NVML_LIB
mv main.x main_2c_15_d.x
make clean
sed -i \
-e 's/^#define[[:space:]]\+PRC[[:space:]]\+8/#define PRC 4/' \
-e 's/^#define[[:space:]]\+STRPRC[[:space:]]\+8/#define STRPRC 4/' \
defines.h
make nvfortran-nvml-mpi GPUCC=$MYCC NVML_INC=$NVML_INC NNVML_LIB=$NVML_LIB
mv main.x main_2c_15.x
make clean
sed -i 's/^#define[[:space:]]\+LATTICE[[:space:]]\+15/#define LATTICE 19/' defines.h
make nvfortran-nvml-mpi GPUCC=$MYCC NVML_INC=$NVML_INC NNVML_LIB=$NVML_LIB
mv main.x main_2c_19.x
make clean
sed -i 's/^#define[[:space:]]\+LATTICE[[:space:]]\+19/#define LATTICE 27/' defines.h
make nvfortran-nvml-mpi GPUCC=$MYCC NVML_INC=$NVML_INC NNVML_LIB=$NVML_LIB
mv main.x main_2c_27.x
make clean
sed -i 's/^#define[[:space:]]\+noHIGHORDER/#define HIGHORDER/' defines.h
make nvfortran-nvml-mpi GPUCC=$MYCC NVML_INC=$NVML_INC NNVML_LIB=$NVML_LIB
mv main.x main_2c_27high.x
make clean
sed -i \
-e 's/^#define[[:space:]]\+noMIXEDPRC/#define MIXEDPRC/' \
-e 's/^#define[[:space:]]\+PRC[[:space:]]\+4/#define PRC 8/' \
defines.h
make nvfortran-nvml-mpi GPUCC=$MYCC NVML_INC=$NVML_INC NNVML_LIB=$NVML_LIB
mv main.x main_2c_27high_sd.x
make clean
sed -i 's/^#define[[:space:]]\+HIGHORDER/#define noHIGHORDER/' defines.h
make nvfortran-nvml-mpi GPUCC=$MYCC NVML_INC=$NVML_INC NNVML_LIB=$NVML_LIB
mv main.x main_2c_27_sd.x
make clean
sed -i 's/^#define[[:space:]]\+LATTICE[[:space:]]\+27/#define LATTICE 19/' defines.h
make nvfortran-nvml-mpi GPUCC=$MYCC NVML_INC=$NVML_INC NNVML_LIB=$NVML_LIB
mv main.x main_2c_19_sd.x
make clean
sed -i 's/^#define[[:space:]]\+LATTICE[[:space:]]\+19/#define LATTICE 15/' defines.h
make nvfortran-nvml-mpi GPUCC=$MYCC NVML_INC=$NVML_INC NNVML_LIB=$NVML_LIB
mv main.x main_2c_15_sd.x
make clean
sed -i \
-e 's/^#define[[:space:]]\+TWOCOMPONENT/#define noTWOCOMPONENT/' \
-e 's/^#define[[:space:]]\+DENSRATIO/#define noDENSRATIO/' \
-e 's/^#define[[:space:]]\+INTERFACE_INCOMP/#define noINTERFACE_INCOMP/' \
defines.h
make nvfortran-nvml-mpi GPUCC=$MYCC NVML_INC=$NVML_INC NNVML_LIB=$NVML_LIB
mv main.x main_1c_15_sd.x
make clean
sed -i 's/^#define[[:space:]]\+LATTICE[[:space:]]\+15/#define LATTICE 19/' defines.h
make nvfortran-nvml-mpi GPUCC=$MYCC NVML_INC=$NVML_INC NNVML_LIB=$NVML_LIB
mv main.x main_1c_19_sd.x
make clean
sed -i 's/^#define[[:space:]]\+LATTICE[[:space:]]\+19/#define LATTICE 27/' defines.h
make nvfortran-nvml-mpi GPUCC=$MYCC NVML_INC=$NVML_INC NNVML_LIB=$NVML_LIB
mv main.x main_1c_27_sd.x
make clean
sed -i 's/^#define[[:space:]]\+noHIGHORDER/#define HIGHORDER/' defines.h
make nvfortran-nvml-mpi GPUCC=$MYCC NVML_INC=$NVML_INC NNVML_LIB=$NVML_LIB
mv main.x main_1c_27high_sd.x
make clean
sed -i \
-e 's/^#define[[:space:]]\+MIXEDPRC/#define noMIXEDPRC/' \
-e 's/^#define[[:space:]]\+PRC[[:space:]]\+8/#define PRC 4/' \
defines.h
make nvfortran-nvml-mpi GPUCC=$MYCC NVML_INC=$NVML_INC NNVML_LIB=$NVML_LIB
mv main.x main_1c_27high.x
make clean
sed -i 's/^#define[[:space:]]\+HIGHORDER/#define noHIGHORDER/' defines.h
make nvfortran-nvml-mpi GPUCC=$MYCC NVML_INC=$NVML_INC NNVML_LIB=$NVML_LIB
mv main.x main_1c_27.x
make clean
sed -i 's/^#define[[:space:]]\+LATTICE[[:space:]]\+27/#define LATTICE 19/' defines.h
make nvfortran-nvml-mpi GPUCC=$MYCC NVML_INC=$NVML_INC NNVML_LIB=$NVML_LIB
mv main.x main_1c_19.x
make clean
sed -i 's/^#define[[:space:]]\+LATTICE[[:space:]]\+19/#define LATTICE 15/' defines.h
make nvfortran-nvml-mpi GPUCC=$MYCC NVML_INC=$NVML_INC NNVML_LIB=$NVML_LIB
mv main.x main_1c_15.x
make clean
sed -i \
-e 's/^#define[[:space:]]\+STRPRC[[:space:]]\+4/#define STRPRC 8/' \
-e 's/^#define[[:space:]]\+PRC[[:space:]]\+4/#define PRC 8/' \
defines.h
make nvfortran-nvml-mpi GPUCC=$MYCC NVML_INC=$NVML_INC NNVML_LIB=$NVML_LIB
mv main.x main_1c_15_d.x
make clean
sed -i 's/^#define[[:space:]]\+LATTICE[[:space:]]\+15/#define LATTICE 19/' defines.h
make nvfortran-nvml-mpi GPUCC=$MYCC NVML_INC=$NVML_INC NNVML_LIB=$NVML_LIB
mv main.x main_1c_19_d.x
make clean
sed -i 's/^#define[[:space:]]\+LATTICE[[:space:]]\+19/#define LATTICE 27/' defines.h
make nvfortran-nvml-mpi GPUCC=$MYCC NVML_INC=$NVML_INC NNVML_LIB=$NVML_LIB
mv main.x main_1c_27_d.x
make clean
sed -i 's/^#define[[:space:]]\+noHIGHORDER/#define HIGHORDER/' defines.h
make nvfortran-nvml-mpi GPUCC=$MYCC NVML_INC=$NVML_INC NNVML_LIB=$NVML_LIB
mv main.x main_1c_27high_d.x
make clean
sed -i \
-e 's/^#define[[:space:]]\+noTWOCOMPONENT/#define TWOCOMPONENT/' \
-e 's/^#define[[:space:]]\+noDENSRATIO/#define DENSRATIO/' \
-e 's/^#define[[:space:]]\+noINTERFACE_INCOMP/#define INTERFACE_INCOMP/' \
defines.h
