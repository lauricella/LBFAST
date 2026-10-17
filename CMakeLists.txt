# CMake build for LBFAST, an alternative to the Makefile (docs/CMAKE.md).
#
# Both build systems read the source lists and compiler flags from
# lbfast_build.mk. The compile-time macros come from defines.h, as with make,
# unless LBFAST_USE_DEFINES_H=OFF. Then every macro of the user block of
# defines.h becomes a cache variable whose default is its value there, and
# CMake writes the chosen values to lbfast_config.h in the build directory,
# leaving defines.h untouched. This allows several configurations to be
# built side by side, each in its own build directory.
#
#   cmake -S . -B build -DGPUCC=80          # macros from defines.h
#   cmake --build build -j 8
#
#   cmake --preset 2c_27high_d              # one of the compile.sh variants
#   cmake --build --preset 2c_27high_d -j 8

cmake_minimum_required(VERSION 3.20)

# An in-source build would replace the LBFAST Makefile with a generated one.
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
  message(FATAL_ERROR
    "In-source builds would overwrite the LBFAST Makefile. Remove the "
    "CMakeCache.txt and CMakeFiles/ just created here and use a separate "
    "build directory, e.g.  cmake -S . -B build")
endif()

# Default to the NVIDIA HPC SDK compilers when neither the command line nor
# the environment (FC, CC) chooses others.
if(NOT DEFINED CMAKE_Fortran_COMPILER AND NOT DEFINED ENV{FC})
  set(CMAKE_Fortran_COMPILER nvfortran)
endif()
if(NOT DEFINED CMAKE_C_COMPILER AND NOT DEFINED ENV{CC})
  set(CMAKE_C_COMPILER nvc)
endif()

project(LBFAST VERSION 1.1.0 LANGUAGES C Fortran)

if(NOT CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
  message(FATAL_ERROR "LBFAST needs nvfortran from the NVIDIA HPC SDK; "
    "CMake found ${CMAKE_Fortran_COMPILER_ID} (${CMAKE_Fortran_COMPILER}). "
    "Select it with -DCMAKE_Fortran_COMPILER=nvfortran in a new build directory.")
endif()

# ---------------------------------------------------------------------------
# Options
# ---------------------------------------------------------------------------
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING
    "Release (flags of make nvfortran*) or Debug (make nvfortran*-debug)" FORCE)
endif()
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Release Debug)
if(NOT CMAKE_BUILD_TYPE MATCHES "^(Release|Debug)$")
  message(FATAL_ERROR "CMAKE_BUILD_TYPE must be Release or Debug, not '${CMAKE_BUILD_TYPE}'.")
endif()

set(GPUCC "all" CACHE STRING
  "GPU compute capability, passed as -gpu=cc<GPUCC> (e.g. 80, 90, all), as make GPUCC=...")
option(LBFAST_MPI "Build the MPI version (-DMPI), as make nvfortran-mpi" ON)
option(LBFAST_NVML "Enable NVML power monitoring (-D_NVML), as make nvfortran-nvml" OFF)
option(LBFAST_USE_DEFINES_H
  "Take the compile-time macros from defines.h, as make does; OFF takes them from cache variables" ON)
set(LBFAST_EXTRA_DEFINES "" CACHE STRING
  "Additional macros for every source, e.g. BOUNCE_BACK;GETPOWER or NAME=value")
set(LBFAST_EXE_NAME "main.x" CACHE STRING "Name of the executable")
set(LBFAST_OUTPUT_DIR "" CACHE PATH "Directory of the executable (default: the build directory)")
# compile_commands.json, read by tests/check_build_systems.py and by editors.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Only the LBFAST flags below are used: drop the CMake defaults for each
# build type (e.g. -fast for Release), which the Makefile does not use.
foreach(lang IN ITEMS C Fortran)
  foreach(cfg IN ITEMS RELEASE DEBUG RELWITHDEBINFO MINSIZEREL)
    set(CMAKE_${lang}_FLAGS_${cfg} "")
  endforeach()
endforeach()

# ---------------------------------------------------------------------------
# Source lists and flags shared with the Makefile
# ---------------------------------------------------------------------------
set(LBFAST_MK "${CMAKE_CURRENT_SOURCE_DIR}/lbfast_build.mk")
set(LBFAST_DEFINES_H "${CMAKE_CURRENT_SOURCE_DIR}/defines.h")
# Re-run CMake when either file changes (new sources, new flags, new macros).
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
  "${LBFAST_MK}" "${LBFAST_DEFINES_H}")

# Read the "NAME = value" assignments of lbfast_build.mk into CMake lists.
function(lbfast_read_mk file)
  file(READ "${file}" text)
  string(REGEX REPLACE "\\\\\n" " " text "${text}")
  string(REPLACE "\n" ";" lines "${text}")
  foreach(line IN LISTS lines)
    if(line MATCHES "^([A-Za-z0-9_]+)[ \t]*=(.*)$")
      set(name "${CMAKE_MATCH_1}")
      separate_arguments(value UNIX_COMMAND "${CMAKE_MATCH_2}")
      set(${name} "${value}" PARENT_SCOPE)
    endif()
  endforeach()
endfunction()

lbfast_read_mk("${LBFAST_MK}")
foreach(var IN ITEMS LBFAST_C_SOURCES LBFAST_NVML_C_SOURCES LBFAST_NVML_F_SOURCES
    LBFAST_F_SOURCES LBFAST_HEADERS LBFAST_FFLAGS LBFAST_CUDAFLAGS LBFAST_INFO
    LBFAST_LDFLAGS LBFAST_LINK_GPU_OPTS LBFAST_RELEASE_OPT LBFAST_RELEASE_LINK_OPT
    LBFAST_RELEASE_GPU_OPTS LBFAST_DEBUG_OPT LBFAST_DEBUG_LINK_OPT
    LBFAST_DEBUG_GPU_OPTS LBFAST_DEBUG_CHECKS)
  if(NOT DEFINED ${var})
    message(FATAL_ERROR "${var} is missing from lbfast_build.mk")
  endif()
endforeach()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(opt ${LBFAST_DEBUG_OPT})
  set(link_opt ${LBFAST_DEBUG_LINK_OPT})
  set(gpu_opts ${LBFAST_DEBUG_GPU_OPTS})
  set(checks ${LBFAST_DEBUG_CHECKS})
else()
  set(opt ${LBFAST_RELEASE_OPT})
  set(link_opt ${LBFAST_RELEASE_LINK_OPT})
  set(gpu_opts ${LBFAST_RELEASE_GPU_OPTS})
  set(checks "")
endif()

# ---------------------------------------------------------------------------
# Compile-time macros
# ---------------------------------------------------------------------------
# The user block of defines.h is the part between
# "#ifndef LBFAST_CMAKE_CONFIG" and "#else". Each "#define NAME" (ON),
# "#define noNAME" (OFF) or "#define NAME value" line defines one macro.
file(READ "${LBFAST_DEFINES_H}" defines_text)
string(REPLACE "\n" ";" defines_lines "${defines_text}")
set(block_lines "")
set(macro_names "")
set(in_block FALSE)
foreach(line IN LISTS defines_lines)
  if(line MATCHES "^#ifndef[ \t]+LBFAST_CMAKE_CONFIG")
    set(in_block TRUE)
    continue()
  endif()
  if(NOT in_block)
    continue()
  endif()
  if(line MATCHES "^#else")
    break()
  endif()
  list(APPEND block_lines "${line}")
  if(line MATCHES "^#define[ \t]+([A-Za-z_][A-Za-z0-9_]*)[ \t]*(.*)$")
    set(name "${CMAKE_MATCH_1}")
    string(STRIP "${CMAKE_MATCH_2}" value)
    if(value STREQUAL "")
      if(name MATCHES "^no([A-Z_].*)$")
        set(name "${CMAKE_MATCH_1}")
        set(default_${name} OFF)
      else()
        set(default_${name} ON)
      endif()
      set(kind_${name} switch)
    else()
      set(default_${name} "${value}")
      set(kind_${name} value)
    endif()
    list(APPEND macro_names "${name}")
  endif()
endforeach()
if(NOT macro_names)
  message(FATAL_ERROR "No user block (#ifndef LBFAST_CMAKE_CONFIG ... #else) found in defines.h")
endif()

set(LBFAST_CONFIG_HEADER "${CMAKE_CURRENT_BINARY_DIR}/lbfast_config.h")
set(config_definitions "")
if(LBFAST_USE_DEFINES_H)
  # A macro given on the command line would be silently ignored: refuse it.
  set(ignored "")
  foreach(name IN LISTS macro_names)
    if(DEFINED CACHE{${name}})
      list(APPEND ignored "${name}")
    endif()
  endforeach()
  if(ignored)
    list(JOIN ignored ", " ignored)
    message(FATAL_ERROR "Macros set as cache variables (${ignored}) are ignored "
      "while LBFAST_USE_DEFINES_H=ON, which takes them from defines.h. Add "
      "-DLBFAST_USE_DEFINES_H=OFF, or use a new build directory to go back to defines.h.")
  endif()
else()
  # One cache variable per macro, defaulting to defines.h; the header keeps
  # the layout of the user block with the chosen values.
  set(config "/* Generated by CMake from the user block of defines.h and the cache\n"
             "   variables of this build directory: change those, not this file. */\n")
  foreach(line IN LISTS block_lines)
    if(line MATCHES "^#define[ \t]+(no)?([A-Za-z_][A-Za-z0-9_]*)")
      set(name "${CMAKE_MATCH_2}")
      if(NOT DEFINED kind_${name})
        set(name "${CMAKE_MATCH_1}${CMAKE_MATCH_2}")
      endif()
      if(kind_${name} STREQUAL "switch")
        option(${name} "defines.h macro ${name} (OFF writes #define no${name})" ${default_${name}})
        if(${name})
          set(line "#define ${name}")
        else()
          set(line "#define no${name}")
        endif()
      else()
        set(${name} "${default_${name}}" CACHE STRING "defines.h macro ${name}")
        if("${${name}}" STREQUAL "")
          message(FATAL_ERROR "Macro ${name} needs a value.")
        endif()
        set(line "#define ${name} ${${name}}")
      endif()
    endif()
    string(APPEND config "${line}\n")
  endforeach()
  file(CONFIGURE OUTPUT "${LBFAST_CONFIG_HEADER}" CONTENT "${config}" @ONLY)
  list(APPEND config_definitions LBFAST_CMAKE_CONFIG)
endif()

# ---------------------------------------------------------------------------
# Executable
# ---------------------------------------------------------------------------
set(c_sources ${LBFAST_C_SOURCES})
set(f_sources ${LBFAST_F_SOURCES})
if(LBFAST_NVML)
  list(APPEND c_sources ${LBFAST_NVML_C_SOURCES})
  list(PREPEND f_sources ${LBFAST_NVML_F_SOURCES})
endif()

# Flags, as composed in the Makefile: every Fortran source gets the common
# flags, lb_cuda_*.f90 add LBFAST_CUDAFLAGS, the others add -DMPI.
set(fortran_flags ${opt} ${LBFAST_FFLAGS} "-gpu=cc${GPUCC},${gpu_opts}" ${LBFAST_INFO} ${checks})
set(header_deps "")
foreach(header IN LISTS LBFAST_HEADERS)
  list(APPEND header_deps "${CMAKE_CURRENT_SOURCE_DIR}/${header}")
endforeach()
if(NOT LBFAST_USE_DEFINES_H)
  list(APPEND header_deps "${LBFAST_CONFIG_HEADER}")
endif()
if(LBFAST_MPI)
  find_package(MPI REQUIRED COMPONENTS Fortran)
endif()

# Module dependencies. The CMake dependency scanner of the Makefile
# generators misses "use" statements whose continuation lines contain
# preprocessor directives, which LBFAST uses throughout, and ignores
# per-source macros such as -DMPI. Each Fortran source is therefore an
# object library, ordered by the module/use statements found below with a
# plain text search. The search ignores #ifdef, so it may add an unneeded
# ordering, never miss one. Sources that do not depend on each other compile
# in parallel, and an object is rebuilt when any source it depends on,
# directly or not, changes.
set(fortran_targets "")
foreach(src IN LISTS f_sources)
  get_filename_component(stem "${src}" NAME_WE)
  set(tgt "lbfast_${stem}")
  list(APPEND fortran_targets ${tgt})
  set(src_${tgt} "${src}")
  file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/${src}" lines
    REGEX "^[ \t]*([Mm][Oo][Dd][Uu][Ll][Ee]|[Uu][Ss][Ee])[ \t,:]")
  set(uses_${tgt} "")
  foreach(line IN LISTS lines)
    string(TOLOWER "${line}" line)
    if(line MATCHES "^[ \t]*module[ \t]+([a-z0-9_]+)[ \t]*(!.*)?$"
        AND NOT CMAKE_MATCH_1 STREQUAL "procedure")
      set(provider_${CMAKE_MATCH_1} "${tgt}")
    elseif(line MATCHES "^[ \t]*use[ \t]*(,[ \t]*non_intrinsic[ \t]*)?(::)?[ \t]*([a-z0-9_]+)")
      list(APPEND uses_${tgt} "${CMAKE_MATCH_3}")
    endif()
  endforeach()
endforeach()

set(index 0)
foreach(tgt IN LISTS fortran_targets)
  set(src "${src_${tgt}}")
  add_library(${tgt} OBJECT "${src}")
  set(deps "")
  foreach(module IN LISTS uses_${tgt})
    set(provider "${provider_${module}}")
    if(provider STREQUAL "" OR provider STREQUAL tgt)
      continue()   # intrinsic or external (mpi, openacc, cudafor, ...), or own module
    endif()
    list(FIND fortran_targets ${provider} provider_index)
    if(provider_index GREATER index)
      message(FATAL_ERROR "${src} uses module ${module}, defined in ${src_${provider}}, "
        "which comes after it in LBFAST_F_SOURCES (lbfast_build.mk). Fix the order: "
        "the Makefile compiles the sources in that order.")
    endif()
    add_dependencies(${tgt} ${provider})
    list(APPEND deps "${CMAKE_CURRENT_SOURCE_DIR}/${src_${provider}}" ${deps_${provider}})
  endforeach()
  list(REMOVE_DUPLICATES deps)
  set(deps_${tgt} ${deps})
  math(EXPR index "${index} + 1")

  set_target_properties(${tgt} PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules")
  target_include_directories(${tgt} PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_BINARY_DIR}")
  target_compile_options(${tgt} PRIVATE ${fortran_flags})
  target_compile_definitions(${tgt} PRIVATE ${config_definitions} ${LBFAST_EXTRA_DEFINES})
  if(src MATCHES "^lb_cuda_")
    target_compile_options(${tgt} PRIVATE ${LBFAST_CUDAFLAGS})
  elseif(LBFAST_MPI)
    target_compile_definitions(${tgt} PRIVATE MPI)
  endif()
  if(LBFAST_NVML)
    target_compile_definitions(${tgt} PRIVATE _NVML)
  endif()
  if(LBFAST_MPI)
    target_link_libraries(${tgt} PRIVATE MPI::MPI_Fortran)
  endif()
  set_property(SOURCE "${src}" APPEND PROPERTY OBJECT_DEPENDS ${header_deps} ${deps})
endforeach()

add_executable(lbfast ${c_sources})
target_link_libraries(lbfast PRIVATE ${fortran_targets})
set_target_properties(lbfast PROPERTIES
  OUTPUT_NAME "${LBFAST_EXE_NAME}"
  SUFFIX ""
  LINKER_LANGUAGE Fortran)
if(LBFAST_OUTPUT_DIR)
  set_target_properties(lbfast PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${LBFAST_OUTPUT_DIR}")
endif()
target_include_directories(lbfast PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_BINARY_DIR}")
target_compile_options(lbfast PRIVATE ${opt})
target_compile_definitions(lbfast PRIVATE ${config_definitions} ${LBFAST_EXTRA_DEFINES})
if(LBFAST_NVML)
  set_property(SOURCE ${LBFAST_NVML_C_SOURCES} APPEND PROPERTY COMPILE_DEFINITIONS _NVML)
  set_property(SOURCE ${LBFAST_NVML_C_SOURCES} APPEND PROPERTY OBJECT_DEPENDS ${header_deps})
endif()
target_link_options(lbfast PRIVATE
  ${link_opt} ${LBFAST_LDFLAGS} "-gpu=cc${GPUCC},${LBFAST_LINK_GPU_OPTS}" ${LBFAST_INFO} ${checks})
if(LBFAST_MPI)
  target_link_libraries(lbfast PRIVATE MPI::MPI_Fortran)
endif()

if(LBFAST_NVML)
  # Same inputs as make NVML_INC=... NVML_LIB=..., then the CUDA toolkits
  # shipped with the HPC SDK of the selected compiler.
  get_filename_component(compiler_dir "${CMAKE_Fortran_COMPILER}" REALPATH)
  get_filename_component(compiler_dir "${compiler_dir}" DIRECTORY)
  set(nvml_hints "")
  foreach(root IN ITEMS "$ENV{NVHPC_ROOT}" "${compiler_dir}/../..")
    if(root)
      file(GLOB cuda_dirs LIST_DIRECTORIES true "${root}/cuda/[0-9]*")
      list(SORT cuda_dirs COMPARE NATURAL ORDER DESCENDING)
      foreach(dir IN LISTS cuda_dirs)
        list(APPEND nvml_hints "${dir}/targets/x86_64-linux" "${dir}")
      endforeach()
    endif()
  endforeach()
  find_path(NVML_INCLUDE_DIR nvml.h
    HINTS "$ENV{NVML_INC}" ENV CUDA_HOME ${nvml_hints} PATH_SUFFIXES include)
  find_library(NVML_LIBRARY nvidia-ml
    HINTS "$ENV{NVML_LIB}" ENV CUDA_HOME ${nvml_hints}
    PATH_SUFFIXES lib/stubs lib64/stubs lib lib64)
  if(NOT NVML_INCLUDE_DIR OR NOT NVML_LIBRARY)
    message(FATAL_ERROR "NVML not found: set NVML_INCLUDE_DIR (directory of nvml.h) "
      "and NVML_LIBRARY (libnvidia-ml), or configure with -DLBFAST_NVML=OFF.")
  endif()
  set_property(SOURCE ${LBFAST_NVML_C_SOURCES} APPEND PROPERTY INCLUDE_DIRECTORIES "${NVML_INCLUDE_DIR}")
  # Link as the Makefile does, -L<dir> -lnvidia-ml. With the full path CMake
  # would add the directory, usually CUDA's lib/stubs, to the RPATH, and the
  # executable would then load the stub libcuda.so of that directory instead
  # of the driver (cuInit error 34).
  get_filename_component(NVML_LIBRARY_DIR "${NVML_LIBRARY}" DIRECTORY)
  target_link_options(lbfast PRIVATE "-L${NVML_LIBRARY_DIR}")
  target_link_libraries(lbfast PRIVATE nvidia-ml)
endif()

install(TARGETS lbfast RUNTIME DESTINATION bin)

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
if(LBFAST_USE_DEFINES_H)
  set(macro_source "defines.h")
else()
  set(macro_source "${LBFAST_CONFIG_HEADER}")
endif()
message(STATUS "LBFAST ${PROJECT_VERSION}: ${CMAKE_BUILD_TYPE}, GPUCC=${GPUCC}, "
  "MPI=${LBFAST_MPI}, NVML=${LBFAST_NVML}, executable ${LBFAST_EXE_NAME}")
message(STATUS "LBFAST macros from ${macro_source}")
if(LBFAST_EXTRA_DEFINES)
  message(STATUS "LBFAST extra macros: ${LBFAST_EXTRA_DEFINES}")
endif()
