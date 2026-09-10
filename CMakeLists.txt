cmake_minimum_required(VERSION 3.25)

if(POLICY CMP0206)
  cmake_policy(SET CMP0206 NEW)
endif()

if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_CURRENT_BINARY_DIR)
  message(
    FATAL_ERROR
    "In-source builds are unsupported; configure with a preset or pass a "
    "separate build directory with -B"
  )
endif()

project(
  ParallelAssemblyCpp
  VERSION 5.0.0
  DESCRIPTION "Compute molecular assembly indices and pathways"
  HOMEPAGE_URL "https://github.com/ELIFE-ASU/parallelassemblycpp"
  LANGUAGES CXX
)

include(GNUInstallDirs)

set(
  PARALLELASSEMBLYCPP_EXECUTABLE_PERMISSIONS
  OWNER_READ OWNER_WRITE OWNER_EXECUTE
  GROUP_READ GROUP_EXECUTE
  WORLD_READ WORLD_EXECUTE
)
set(
  PARALLELASSEMBLYCPP_FILE_PERMISSIONS
  OWNER_READ OWNER_WRITE
  GROUP_READ
  WORLD_READ
)
set(
  CMAKE_INSTALL_DEFAULT_DIRECTORY_PERMISSIONS
  ${PARALLELASSEMBLYCPP_EXECUTABLE_PERMISSIONS}
)

option(
  PARALLELASSEMBLYCPP_BUILD_TELEMETRY
  "Build the instrumented ParallelAssemblyCppTelemetry executable"
  OFF
)
option(
  PARALLELASSEMBLYCPP_X86_64_V3
  "Target x86-64-v3 with POPCNT"
  OFF
)
option(
  PARALLELASSEMBLYCPP_ENABLE_LTO
  "Enable link-time optimization"
  OFF
)
option(
  PARALLELASSEMBLYCPP_STRICT_WARNINGS
  "Treat high-signal compiler warnings as errors"
  ON
)
set(
  PARALLELASSEMBLYCPP_PGO_PHASE
  "OFF"
  CACHE STRING
  "PGO phase: OFF, GENERATE, or USE"
)
set_property(
  CACHE PARALLELASSEMBLYCPP_PGO_PHASE
  PROPERTY STRINGS OFF GENERATE USE
)
set(
  PARALLELASSEMBLYCPP_PGO_PROFILE_DIR
  "${CMAKE_BINARY_DIR}/pgo-profiles"
  CACHE PATH
  "Shared PGO profile directory"
)
option(
  PARALLELASSEMBLYCPP_FULL_REGRESSION_TESTS
  "Run the full regression manifest"
  OFF
)
option(
  PARALLELASSEMBLYCPP_BUILD_OPENMP
  "Build the experimental OpenMP executable"
  OFF
)
option(
  PARALLELASSEMBLYCPP_BUILD_MPI
  "Build the experimental MPI executable"
  OFF
)
option(
  PARALLELASSEMBLYCPP_BUILD_HYBRID
  "Build the experimental MPI/OpenMP executable"
  OFF
)

set(_PARALLELASSEMBLYCPP_PGO_PHASES OFF GENERATE USE)
if(NOT PARALLELASSEMBLYCPP_PGO_PHASE IN_LIST _PARALLELASSEMBLYCPP_PGO_PHASES)
  message(
    FATAL_ERROR
    "PARALLELASSEMBLYCPP_PGO_PHASE must be OFF, GENERATE, or USE; got "
    "${PARALLELASSEMBLYCPP_PGO_PHASE}"
  )
endif()

if(PARALLELASSEMBLYCPP_PGO_PHASE STREQUAL "OFF")
  set(PARALLELASSEMBLYCPP_PGO_PROFILE_DIR_ABSOLUTE "")
else()
  if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    message(FATAL_ERROR "PGO requires GCC")
  endif()
  if(PARALLELASSEMBLYCPP_PGO_PROFILE_DIR STREQUAL "")
    message(FATAL_ERROR "PARALLELASSEMBLYCPP_PGO_PROFILE_DIR must not be empty")
  endif()
  cmake_path(
    ABSOLUTE_PATH PARALLELASSEMBLYCPP_PGO_PROFILE_DIR
    BASE_DIRECTORY "${CMAKE_BINARY_DIR}"
    NORMALIZE
    OUTPUT_VARIABLE PARALLELASSEMBLYCPP_PGO_PROFILE_DIR_ABSOLUTE
  )
  set(
    PARALLELASSEMBLYCPP_PGO_DIRECTORY_MARKER
    "${PARALLELASSEMBLYCPP_PGO_PROFILE_DIR_ABSOLUTE}/.parallelassemblycpp-pgo-profile-directory"
  )

  if(PARALLELASSEMBLYCPP_PGO_PHASE STREQUAL "GENERATE")
    file(MAKE_DIRECTORY "${PARALLELASSEMBLYCPP_PGO_PROFILE_DIR_ABSOLUTE}")
    if(NOT EXISTS "${PARALLELASSEMBLYCPP_PGO_DIRECTORY_MARKER}")
      file(
        GLOB PARALLELASSEMBLYCPP_PGO_PROFILE_DIRECTORY_ENTRIES
        LIST_DIRECTORIES true
        "${PARALLELASSEMBLYCPP_PGO_PROFILE_DIR_ABSOLUTE}/*"
      )
      if(PARALLELASSEMBLYCPP_PGO_PROFILE_DIRECTORY_ENTRIES)
        message(
          FATAL_ERROR
          "PGO profile directory is non-empty but unmarked: "
          "${PARALLELASSEMBLYCPP_PGO_PROFILE_DIR_ABSOLUTE}"
        )
      endif()
      file(
        WRITE "${PARALLELASSEMBLYCPP_PGO_DIRECTORY_MARKER}"
        "Dedicated ParallelAssemblyCpp PGO profile directory.\n"
      )
    elseif(
      IS_DIRECTORY "${PARALLELASSEMBLYCPP_PGO_DIRECTORY_MARKER}"
      OR IS_SYMLINK "${PARALLELASSEMBLYCPP_PGO_DIRECTORY_MARKER}"
    )
      message(FATAL_ERROR "PGO profile marker must be a regular file")
    endif()
  else()
    if(
      NOT EXISTS "${PARALLELASSEMBLYCPP_PGO_DIRECTORY_MARKER}"
      OR IS_DIRECTORY "${PARALLELASSEMBLYCPP_PGO_DIRECTORY_MARKER}"
      OR IS_SYMLINK "${PARALLELASSEMBLYCPP_PGO_DIRECTORY_MARKER}"
    )
      message(
        FATAL_ERROR
        "PGO profile marker is missing or invalid; configure preset "
        "pgo-generate first"
      )
    endif()
    set(
      PARALLELASSEMBLYCPP_PGO_COMPLETION_FILE
      "${PARALLELASSEMBLYCPP_PGO_PROFILE_DIR_ABSOLUTE}/parallelassemblycpp-pgo-complete.json"
    )
    if(NOT EXISTS "${PARALLELASSEMBLYCPP_PGO_COMPLETION_FILE}")
      message(
        FATAL_ERROR
        "PGO training is incomplete; build preset pgo-train first"
      )
    endif()
    file(READ "${PARALLELASSEMBLYCPP_PGO_COMPLETION_FILE}" PARALLELASSEMBLYCPP_PGO_COMPLETION)
    string(
      JSON PARALLELASSEMBLYCPP_PGO_COMPLETION_SCHEMA
      ERROR_VARIABLE PARALLELASSEMBLYCPP_PGO_COMPLETION_ERROR
      GET "${PARALLELASSEMBLYCPP_PGO_COMPLETION}" schema_version
    )
    if(
      PARALLELASSEMBLYCPP_PGO_COMPLETION_ERROR
      OR NOT PARALLELASSEMBLYCPP_PGO_COMPLETION_SCHEMA EQUAL 2
    )
      message(
        FATAL_ERROR
        "Invalid PGO training completion record: "
        "${PARALLELASSEMBLYCPP_PGO_COMPLETION_FILE}"
      )
    endif()
    string(
      JSON PARALLELASSEMBLYCPP_PGO_RECORDED_WEIGHTS_SHA256
      ERROR_VARIABLE PARALLELASSEMBLYCPP_PGO_WEIGHTS_ERROR
      GET "${PARALLELASSEMBLYCPP_PGO_COMPLETION}" weights sha256
    )
    file(
      SHA256 "${PROJECT_SOURCE_DIR}/benchmarks/pgo-training.tsv"
      PARALLELASSEMBLYCPP_PGO_CURRENT_WEIGHTS_SHA256
    )
    if(
      PARALLELASSEMBLYCPP_PGO_WEIGHTS_ERROR
      OR NOT PARALLELASSEMBLYCPP_PGO_RECORDED_WEIGHTS_SHA256
        STREQUAL PARALLELASSEMBLYCPP_PGO_CURRENT_WEIGHTS_SHA256
    )
      message(
        FATAL_ERROR
        "PGO training weights changed; rebuild preset pgo-train"
      )
    endif()
    set(
      PARALLELASSEMBLYCPP_PGO_MANIFEST
      "${PROJECT_SOURCE_DIR}/benchmarks/cases.tsv"
    )
    string(
      JSON PARALLELASSEMBLYCPP_PGO_RECORDED_MANIFEST_SHA256
      ERROR_VARIABLE PARALLELASSEMBLYCPP_PGO_MANIFEST_ERROR
      GET "${PARALLELASSEMBLYCPP_PGO_COMPLETION}" corpus manifest sha256
    )
    file(
      SHA256 "${PARALLELASSEMBLYCPP_PGO_MANIFEST}"
      PARALLELASSEMBLYCPP_PGO_CURRENT_MANIFEST_SHA256
    )
    if(
      PARALLELASSEMBLYCPP_PGO_MANIFEST_ERROR
      OR NOT PARALLELASSEMBLYCPP_PGO_RECORDED_MANIFEST_SHA256
        STREQUAL PARALLELASSEMBLYCPP_PGO_CURRENT_MANIFEST_SHA256
    )
      message(
        FATAL_ERROR
        "Benchmark manifest changed; rebuild preset pgo-train"
      )
    endif()
    string(
      JSON PARALLELASSEMBLYCPP_PGO_INPUT_COUNT
      ERROR_VARIABLE PARALLELASSEMBLYCPP_PGO_INPUT_COUNT_ERROR
      LENGTH "${PARALLELASSEMBLYCPP_PGO_COMPLETION}" corpus inputs
    )
    if(
      PARALLELASSEMBLYCPP_PGO_INPUT_COUNT_ERROR
      OR PARALLELASSEMBLYCPP_PGO_INPUT_COUNT LESS 1
    )
      message(FATAL_ERROR "Invalid PGO corpus input fingerprints")
    endif()
    math(EXPR PARALLELASSEMBLYCPP_PGO_LAST_INPUT "${PARALLELASSEMBLYCPP_PGO_INPUT_COUNT} - 1")
    set(PARALLELASSEMBLYCPP_PGO_CORPUS_FILES "${PARALLELASSEMBLYCPP_PGO_MANIFEST}")
    foreach(PARALLELASSEMBLYCPP_PGO_INPUT_INDEX RANGE ${PARALLELASSEMBLYCPP_PGO_LAST_INPUT})
      string(
        JSON PARALLELASSEMBLYCPP_PGO_RECORDED_INPUT_PATH
        ERROR_VARIABLE PARALLELASSEMBLYCPP_PGO_INPUT_PATH_ERROR
        GET "${PARALLELASSEMBLYCPP_PGO_COMPLETION}"
        corpus inputs ${PARALLELASSEMBLYCPP_PGO_INPUT_INDEX} path
      )
      string(
        JSON PARALLELASSEMBLYCPP_PGO_RECORDED_INPUT_SHA256
        ERROR_VARIABLE PARALLELASSEMBLYCPP_PGO_INPUT_SHA256_ERROR
        GET "${PARALLELASSEMBLYCPP_PGO_COMPLETION}"
        corpus inputs ${PARALLELASSEMBLYCPP_PGO_INPUT_INDEX} sha256
      )
      if(
        PARALLELASSEMBLYCPP_PGO_INPUT_PATH_ERROR
        OR PARALLELASSEMBLYCPP_PGO_INPUT_SHA256_ERROR
      )
        message(FATAL_ERROR "Invalid PGO corpus input fingerprint")
      endif()
      cmake_path(
        ABSOLUTE_PATH PARALLELASSEMBLYCPP_PGO_RECORDED_INPUT_PATH
        BASE_DIRECTORY "${PROJECT_SOURCE_DIR}"
        NORMALIZE
        OUTPUT_VARIABLE PARALLELASSEMBLYCPP_PGO_CURRENT_INPUT
      )
      if(
        NOT EXISTS "${PARALLELASSEMBLYCPP_PGO_CURRENT_INPUT}"
        OR IS_DIRECTORY "${PARALLELASSEMBLYCPP_PGO_CURRENT_INPUT}"
      )
        message(
          FATAL_ERROR
          "PGO input is missing: ${PARALLELASSEMBLYCPP_PGO_CURRENT_INPUT}; rebuild "
          "preset pgo-train"
        )
      endif()
      file(
        SHA256 "${PARALLELASSEMBLYCPP_PGO_CURRENT_INPUT}"
        PARALLELASSEMBLYCPP_PGO_CURRENT_INPUT_SHA256
      )
      if(
        NOT PARALLELASSEMBLYCPP_PGO_RECORDED_INPUT_SHA256
          STREQUAL PARALLELASSEMBLYCPP_PGO_CURRENT_INPUT_SHA256
      )
        message(
          FATAL_ERROR
          "PGO input changed: ${PARALLELASSEMBLYCPP_PGO_CURRENT_INPUT}; rebuild "
          "preset pgo-train"
        )
      endif()
      list(APPEND PARALLELASSEMBLYCPP_PGO_CORPUS_FILES "${PARALLELASSEMBLYCPP_PGO_CURRENT_INPUT}")
    endforeach()
    file(
      GLOB_RECURSE PARALLELASSEMBLYCPP_PGO_PROFILE_DATA
      LIST_DIRECTORIES false
      "${PARALLELASSEMBLYCPP_PGO_PROFILE_DIR_ABSOLUTE}/*.gcda"
    )
    if(NOT PARALLELASSEMBLYCPP_PGO_PROFILE_DATA)
      message(
        FATAL_ERROR
        "PGO data is missing from ${PARALLELASSEMBLYCPP_PGO_PROFILE_DIR_ABSOLUTE}; "
        "build preset pgo-train first"
      )
    endif()
    # Reconfigure when training inputs change to reject partial or stale data.
    set_property(
      DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
      "${PARALLELASSEMBLYCPP_PGO_COMPLETION_FILE}"
      "${PROJECT_SOURCE_DIR}/benchmarks/pgo-training.tsv"
      ${PARALLELASSEMBLYCPP_PGO_CORPUS_FILES}
    )
  endif()
endif()

if(PARALLELASSEMBLYCPP_ENABLE_LTO)
  include(CheckIPOSupported)
  check_ipo_supported(
    RESULT PARALLELASSEMBLYCPP_LTO_SUPPORTED
    OUTPUT PARALLELASSEMBLYCPP_LTO_ERROR
    LANGUAGES CXX
  )
  if(NOT PARALLELASSEMBLYCPP_LTO_SUPPORTED)
    message(FATAL_ERROR "LTO is unavailable: ${PARALLELASSEMBLYCPP_LTO_ERROR}")
  endif()
endif()

function(parallelassemblycpp_configure_target target)
  target_compile_features(${target} PRIVATE cxx_std_20)
  set_target_properties(
    ${target}
    PROPERTIES
      CXX_EXTENSIONS OFF
      CXX_STANDARD_REQUIRED ON
      INSTALL_REMOVE_ENVIRONMENT_RPATH ON
  )

  if(MSVC)
    # These are language-conformance requirements, independent of diagnostics.
    target_compile_options(
      ${target}
      PRIVATE
        /permissive-
        /Zc:__cplusplus
        /Zc:preprocessor
    )
  endif()

  if(PARALLELASSEMBLYCPP_STRICT_WARNINGS)
    set_property(TARGET ${target} PROPERTY COMPILE_WARNING_AS_ERROR ON)

    if(MSVC)
      target_compile_options(
        ${target}
        PRIVATE
          /W4
          # Match the GNU/Clang policy: conversion warnings are handled in
          # focused audits, cache-line alignment intentionally adds padding,
          # and forced-inline failures are optimizer decisions rather than bugs.
          /wd4244
          /wd4267
          /wd4324
          /wd4714
      )
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
      target_compile_options(
        ${target}
        PRIVATE
          -Wall
          -Wextra
          -Wpedantic
          -Wcast-align
          -Wcast-qual
          -Wdouble-promotion
          -Wextra-semi
          -Wformat=2
          -Wimplicit-fallthrough
          -Wnon-virtual-dtor
          -Wnull-dereference
          -Wold-style-cast
          -Woverloaded-virtual
          -Wshadow
          -Wswitch-enum
          -Wundef
          -Wzero-as-null-pointer-constant
      )

      if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        target_compile_options(
          ${target}
          PRIVATE
            -Walloc-zero
            -Warray-bounds=2
            -Wcast-function-type
            -Wdate-time
            -Wduplicated-branches
            -Wduplicated-cond
            -Wformat-overflow=2
            -Wformat-signedness
            -Wformat-truncation=2
            -Wlogical-op
            -Wredundant-decls
            -Wstringop-overflow=4
            -Wstringop-truncation
            -Wtrampolines
            -Wuseless-cast
            -Wvla
        )
      endif()
    else()
      message(
        FATAL_ERROR
        "PARALLELASSEMBLYCPP_STRICT_WARNINGS supports GNU, Clang-compatible, and "
        "MSVC toolchains; disable it only for unsupported compiler diagnosis"
      )
    endif()
  endif()

  if(PARALLELASSEMBLYCPP_ENABLE_LTO)
    set_property(
      TARGET ${target}
      PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE
    )

    if(
      CMAKE_CXX_COMPILER_ID STREQUAL "GNU"
      AND PARALLELASSEMBLYCPP_PGO_PHASE STREQUAL "OFF"
    )
      # One-partition LTO outperformed GCC's default without PGO.
      set(
        _PARALLELASSEMBLYCPP_GNU_LTO_OPTIONS
        -flto-partition=one
        -fipa-pta
        -fno-semantic-interposition
        -fdevirtualize-at-ltrans
        -fno-plt
      )
      target_compile_options(
        ${target}
        PRIVATE ${_PARALLELASSEMBLYCPP_GNU_LTO_OPTIONS}
      )
      target_link_options(
        ${target}
        PRIVATE ${_PARALLELASSEMBLYCPP_GNU_LTO_OPTIONS}
      )
    endif()
  endif()

  if(WIN32)
    target_compile_definitions(${target} PRIVATE NOMINMAX WIN32_LEAN_AND_MEAN)
  endif()

  if(PARALLELASSEMBLYCPP_X86_64_V3)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
      include(CheckCXXCompilerFlag)
      check_cxx_compiler_flag("-march=x86-64-v3" PARALLELASSEMBLYCPP_HAS_X86_64_V3)
      check_cxx_compiler_flag("-mpopcnt" PARALLELASSEMBLYCPP_HAS_MPOPCNT)
      if(NOT PARALLELASSEMBLYCPP_HAS_X86_64_V3 OR NOT PARALLELASSEMBLYCPP_HAS_MPOPCNT)
        message(
          FATAL_ERROR
          "PARALLELASSEMBLYCPP_X86_64_V3 requires compiler support for "
          "-march=x86-64-v3 and -mpopcnt"
        )
      endif()
      target_compile_options(${target} PRIVATE -march=x86-64-v3 -mpopcnt)
      if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        # Allow more inlining for template-heavy search loops.
        target_compile_options(
          ${target}
          PRIVATE --param=inline-unit-growth=104
        )
        if(
          NOT PARALLELASSEMBLYCPP_ENABLE_LTO
          AND PARALLELASSEMBLYCPP_PGO_PHASE STREQUAL "OFF"
        )
          # GCC's versioned AVX2 loops bloat short-circuit mask operations.
          # Keep straight-line SLP vectorization while disabling loop
          # vectorization.
          target_compile_options(
            ${target}
            PRIVATE -fno-tree-loop-vectorize
          )
        endif()
      endif()
    else()
      message(
        FATAL_ERROR
        "PARALLELASSEMBLYCPP_X86_64_V3 is supported only by GCC-compatible and "
        "Clang-compatible compilers"
      )
    endif()
  endif()
endfunction()

function(parallelassemblycpp_configure_test_target target)
  parallelassemblycpp_configure_target(${target})

  # Release-style configurations define NDEBUG globally. Keep standard assert
  # checks active in test executables until they use an always-active helper.
  if(MSVC)
    target_compile_options(${target} PRIVATE /UNDEBUG)
  else()
    target_compile_options(${target} PRIVATE -UNDEBUG)
  endif()
endfunction()

function(parallelassemblycpp_configure_pgo_target target)
  if(PARALLELASSEMBLYCPP_PGO_PHASE STREQUAL "GENERATE")
    set(
      _PARALLELASSEMBLYCPP_PGO_OPTIONS
      "-fprofile-generate=${PARALLELASSEMBLYCPP_PGO_PROFILE_DIR_ABSOLUTE}"
      "-fprofile-prefix-path=${CMAKE_BINARY_DIR}"
    )
  elseif(PARALLELASSEMBLYCPP_PGO_PHASE STREQUAL "USE")
    set(
      _PARALLELASSEMBLYCPP_PGO_OPTIONS
      "-fprofile-use=${PARALLELASSEMBLYCPP_PGO_PROFILE_DIR_ABSOLUTE}"
      "-fprofile-prefix-path=${CMAKE_BINARY_DIR}"
      -fno-tracer
      -fno-reorder-blocks-and-partition
      -Werror=missing-profile
      -Werror=coverage-mismatch
    )
  else()
    return()
  endif()

  # Only the CLI is trained; instrumenting other targets creates profile gaps.
  target_compile_options(${target} PRIVATE ${_PARALLELASSEMBLYCPP_PGO_OPTIONS})
  target_link_options(${target} PRIVATE ${_PARALLELASSEMBLYCPP_PGO_OPTIONS})
endfunction()

function(parallelassemblycpp_add_solver target)
  cmake_parse_arguments(
    PARSE_ARGV 1
    PARALLELASSEMBLYCPP_SOLVER
    "MPI;OPENMP;TELEMETRY"
    ""
    ""
  )
  if(PARALLELASSEMBLYCPP_SOLVER_UNPARSED_ARGUMENTS)
    message(
      FATAL_ERROR
      "Unknown options for ${target}: ${PARALLELASSEMBLYCPP_SOLVER_UNPARSED_ARGUMENTS}"
    )
  endif()

  add_executable(${target} v5/main.cpp)
  set(PARALLELASSEMBLYCPP_SOLVER_DEFINITIONS)
  set(PARALLELASSEMBLYCPP_SOLVER_LIBRARIES)
  if(PARALLELASSEMBLYCPP_SOLVER_MPI)
    list(APPEND PARALLELASSEMBLYCPP_SOLVER_DEFINITIONS PARALLELASSEMBLYCPP_USE_MPI)
    list(APPEND PARALLELASSEMBLYCPP_SOLVER_LIBRARIES MPI::MPI_CXX)
  endif()
  if(PARALLELASSEMBLYCPP_SOLVER_OPENMP)
    list(APPEND PARALLELASSEMBLYCPP_SOLVER_DEFINITIONS PARALLELASSEMBLYCPP_USE_OPENMP)
    list(APPEND PARALLELASSEMBLYCPP_SOLVER_LIBRARIES OpenMP::OpenMP_CXX)
  endif()
  if(PARALLELASSEMBLYCPP_SOLVER_TELEMETRY)
    list(APPEND PARALLELASSEMBLYCPP_SOLVER_DEFINITIONS ASSEMBLY_ENABLE_TELEMETRY)
  endif()
  if(PARALLELASSEMBLYCPP_SOLVER_DEFINITIONS)
    target_compile_definitions(
      ${target}
      PRIVATE ${PARALLELASSEMBLYCPP_SOLVER_DEFINITIONS}
    )
  endif()
  if(PARALLELASSEMBLYCPP_SOLVER_LIBRARIES)
    target_link_libraries(${target} PRIVATE ${PARALLELASSEMBLYCPP_SOLVER_LIBRARIES})
  endif()
  parallelassemblycpp_configure_target(${target})
endfunction()

add_library(ParallelAssemblyCppLibrary STATIC v5/main.cpp)
add_library(ParallelAssemblyCpp::Library ALIAS ParallelAssemblyCppLibrary)
target_compile_definitions(
  ParallelAssemblyCppLibrary
  PRIVATE PARALLELASSEMBLYCPP_LIBRARY_BUILD PARALLELASSEMBLYCPP_NO_MAIN
)
target_include_directories(
  ParallelAssemblyCppLibrary
  PUBLIC
    "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/v5>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/parallelassemblycpp>"
)
target_compile_features(ParallelAssemblyCppLibrary PUBLIC cxx_std_20)
set_target_properties(
  ParallelAssemblyCppLibrary
  PROPERTIES
    EXPORT_NAME Library
    OUTPUT_NAME parallelassemblycpp
    POSITION_INDEPENDENT_CODE ON
)
parallelassemblycpp_configure_target(ParallelAssemblyCppLibrary)

parallelassemblycpp_add_solver(ParallelAssemblyCpp)
parallelassemblycpp_configure_pgo_target(ParallelAssemblyCpp)

if(PARALLELASSEMBLYCPP_BUILD_OPENMP OR PARALLELASSEMBLYCPP_BUILD_HYBRID)
  find_package(OpenMP REQUIRED COMPONENTS CXX)
endif()
if(PARALLELASSEMBLYCPP_BUILD_MPI OR PARALLELASSEMBLYCPP_BUILD_HYBRID)
  find_package(MPI REQUIRED COMPONENTS CXX)
endif()

if(PARALLELASSEMBLYCPP_BUILD_OPENMP)
  parallelassemblycpp_add_solver(ParallelAssemblyCppOMP OPENMP)
endif()

if(PARALLELASSEMBLYCPP_BUILD_MPI)
  parallelassemblycpp_add_solver(ParallelAssemblyCppMPI MPI)
endif()

if(PARALLELASSEMBLYCPP_BUILD_HYBRID)
  parallelassemblycpp_add_solver(ParallelAssemblyCppHybrid MPI OPENMP)
endif()

if(PARALLELASSEMBLYCPP_BUILD_TELEMETRY)
  parallelassemblycpp_add_solver(ParallelAssemblyCppTelemetry TELEMETRY)

  if(PARALLELASSEMBLYCPP_BUILD_OPENMP)
    parallelassemblycpp_add_solver(ParallelAssemblyCppOMPTelemetry OPENMP TELEMETRY)
  endif()

  if(PARALLELASSEMBLYCPP_BUILD_MPI)
    parallelassemblycpp_add_solver(ParallelAssemblyCppMPITelemetry MPI TELEMETRY)
  endif()

  if(PARALLELASSEMBLYCPP_BUILD_HYBRID)
    parallelassemblycpp_add_solver(ParallelAssemblyCppHybridTelemetry MPI OPENMP TELEMETRY)
  endif()
endif()

if(PARALLELASSEMBLYCPP_X86_64_V3 AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  # Internalize single-file targets while preserving public library symbols.
  foreach(
    PARALLELASSEMBLYCPP_WHOLE_PROGRAM_TARGET
    ParallelAssemblyCpp
    ParallelAssemblyCppLibrary
    ParallelAssemblyCppOMP
    ParallelAssemblyCppMPI
    ParallelAssemblyCppHybrid
    ParallelAssemblyCppTelemetry
    ParallelAssemblyCppOMPTelemetry
    ParallelAssemblyCppMPITelemetry
    ParallelAssemblyCppHybridTelemetry
  )
    if(TARGET ${PARALLELASSEMBLYCPP_WHOLE_PROGRAM_TARGET})
      target_compile_options(
        ${PARALLELASSEMBLYCPP_WHOLE_PROGRAM_TARGET}
        PRIVATE -fwhole-program
      )
    endif()
  endforeach()
endif()

install(
  TARGETS ParallelAssemblyCpp ParallelAssemblyCppLibrary
  EXPORT ParallelAssemblyCppTargets
  RUNTIME
    DESTINATION "${CMAKE_INSTALL_BINDIR}"
    PERMISSIONS ${PARALLELASSEMBLYCPP_EXECUTABLE_PERMISSIONS}
  ARCHIVE
    DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    PERMISSIONS ${PARALLELASSEMBLYCPP_FILE_PERMISSIONS}
  INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/parallelassemblycpp"
)
install(
  FILES v5/parallelassemblycpp.h
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/parallelassemblycpp"
  PERMISSIONS ${PARALLELASSEMBLYCPP_FILE_PERMISSIONS}
)
install(
  FILES README.md License.md
  DESTINATION "${CMAKE_INSTALL_DATAROOTDIR}/doc/parallelassemblycpp"
  PERMISSIONS ${PARALLELASSEMBLYCPP_FILE_PERMISSIONS}
)
install(
  EXPORT ParallelAssemblyCppTargets
  FILE ParallelAssemblyCppTargets.cmake
  NAMESPACE ParallelAssemblyCpp::
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/ParallelAssemblyCpp"
)
include(CMakePackageConfigHelpers)
configure_package_config_file(
  cmake/ParallelAssemblyCppConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/ParallelAssemblyCppConfig.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/ParallelAssemblyCpp"
)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/ParallelAssemblyCppConfigVersion.cmake"
  VERSION "${PROJECT_VERSION}"
  COMPATIBILITY SameMajorVersion
)
install(
  FILES
    "${CMAKE_CURRENT_BINARY_DIR}/ParallelAssemblyCppConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/ParallelAssemblyCppConfigVersion.cmake"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/ParallelAssemblyCpp"
)
if(UNIX AND NOT APPLE)
  # Remove absolute RPATHs injected by Conda compiler wrappers.
  if(IS_ABSOLUTE "${CMAKE_INSTALL_BINDIR}")
    set(
      PARALLELASSEMBLYCPP_INSTALLED_EXECUTABLE
      "\$ENV{DESTDIR}${CMAKE_INSTALL_BINDIR}/ParallelAssemblyCpp"
    )
  else()
    set(
      PARALLELASSEMBLYCPP_INSTALLED_EXECUTABLE
      "\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}/ParallelAssemblyCpp"
    )
  endif()
  install(
    CODE
      "file(RPATH_REMOVE FILE \"${PARALLELASSEMBLYCPP_INSTALLED_EXECUTABLE}\")"
  )
endif()

include(CTest)
if(PARALLELASSEMBLYCPP_FULL_REGRESSION_TESTS AND NOT BUILD_TESTING)
  message(
    FATAL_ERROR
    "PARALLELASSEMBLYCPP_FULL_REGRESSION_TESTS requires BUILD_TESTING=ON"
  )
endif()
if(BUILD_TESTING OR PARALLELASSEMBLYCPP_PGO_PHASE STREQUAL "GENERATE")
  find_package(Python3 3.10 REQUIRED COMPONENTS Interpreter)
endif()

if(PARALLELASSEMBLYCPP_PGO_PHASE STREQUAL "GENERATE")
  add_custom_target(
    pgo-train
    COMMAND
      "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/benchmarks/pgo_training.py"
      --executable "$<TARGET_FILE:ParallelAssemblyCpp>"
      --profile-dir "${PARALLELASSEMBLYCPP_PGO_PROFILE_DIR_ABSOLUTE}"
      --weights "${PROJECT_SOURCE_DIR}/benchmarks/pgo-training.tsv"
    DEPENDS ParallelAssemblyCpp
    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
    USES_TERMINAL
    VERBATIM
  )
endif()

function(parallelassemblycpp_add_unit_test name target source)
  add_executable(${target} ${source})
  parallelassemblycpp_configure_test_target(${target})
  add_test(NAME unit.${name} COMMAND ${target} ${ARGN})
  set_tests_properties(
    unit.${name}
    PROPERTIES
      LABELS "unit"
      TIMEOUT 300
  )
endfunction()

function(parallelassemblycpp_add_python_test name timeout)
  add_test(NAME ${name} COMMAND "${Python3_EXECUTABLE}" ${ARGN})
  set_tests_properties(
    ${name}
    PROPERTIES
      LABELS "tooling"
      TIMEOUT ${timeout}
      WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
  )
endfunction()

if(BUILD_TESTING)
  parallelassemblycpp_add_unit_test(
    active-word-mask
    parallelassemblycpp_active_word_mask_tests
    unitTests/activeWordMaskTester.cpp
  )
  parallelassemblycpp_add_unit_test(
    assembly-transposition-table
    parallelassemblycpp_transposition_table_tests
    unitTests/assemblyTranspositionTableTester.cpp
  )

  add_executable(
    parallelassemblycpp_tree_canon_tests
    unitTests/treeCanonTester.cpp
  )
  parallelassemblycpp_configure_test_target(parallelassemblycpp_tree_canon_tests)
  if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    add_test(
      NAME unit.tree-canonicalisation
      COMMAND
        sh -c "ulimit -s 64 && exec \"$1\"" _
        "$<TARGET_FILE:parallelassemblycpp_tree_canon_tests>"
    )
  else()
    add_test(
      NAME unit.tree-canonicalisation
      COMMAND parallelassemblycpp_tree_canon_tests
    )
  endif()
  set_tests_properties(
    unit.tree-canonicalisation
    PROPERTIES
      LABELS "unit"
      TIMEOUT 300
  )

  parallelassemblycpp_add_unit_test(
    cyclic-canonicalisation
    parallelassemblycpp_cyclic_canon_tests
    unitTests/cyclicCanonTester.cpp
  )
  parallelassemblycpp_add_unit_test(
    search-stop-polling
    parallelassemblycpp_search_stop_poll_tests
    unitTests/searchStopPollTester.cpp
  )
  parallelassemblycpp_add_unit_test(
    fragmentation-metadata
    parallelassemblycpp_fragmentation_metadata_tests
    unitTests/fragmentationMetadataTester.cpp
  )
  parallelassemblycpp_add_unit_test(
    parallel-task-transfer
    parallelassemblycpp_parallel_task_transfer_tests
    unitTests/parallelTaskTransferTester.cpp
  )
  set_tests_properties(
    unit.parallel-task-transfer
    PROPERTIES LABELS "unit;parallel"
  )
  target_compile_definitions(
    parallelassemblycpp_parallel_task_transfer_tests
    PRIVATE ASSEMBLY_ENABLE_TELEMETRY
  )
  if(TARGET OpenMP::OpenMP_CXX)
    target_compile_definitions(
      parallelassemblycpp_parallel_task_transfer_tests
      PRIVATE PARALLELASSEMBLYCPP_USE_OPENMP
    )
    target_link_libraries(
      parallelassemblycpp_parallel_task_transfer_tests
      PRIVATE OpenMP::OpenMP_CXX
    )
  endif()
  parallelassemblycpp_add_unit_test(
    molfile-parser
    parallelassemblycpp_molfile_parser_tests
    unitTests/molfileParserTester.cpp
  )
  parallelassemblycpp_add_unit_test(
    native-graph-parser
    parallelassemblycpp_graphio_tests
    unitTests/graphioTester.cpp
  )
  parallelassemblycpp_add_unit_test(
    pathway-generator
    parallelassemblycpp_pathway_generator_tests
    unitTests/pathwayGeneratorTester.cpp
  )
  parallelassemblycpp_add_unit_test(
    string-assembly
    parallelassemblycpp_string_assembly_tests
    unitTests/stringAssemblyTester.cpp
  )
  parallelassemblycpp_add_unit_test(
    library-batch
    parallelassemblycpp_library_batch_tests
    unitTests/libraryBatchTester.cpp
    "${PROJECT_SOURCE_DIR}/unitTests/icosane.mol"
    "${PROJECT_SOURCE_DIR}/unitTests/sucrose.mol"
    "${PROJECT_SOURCE_DIR}/unitTests/butane.mol"
    "${PROJECT_SOURCE_DIR}/unitTests/graphio_test"
  )
  target_link_libraries(
    parallelassemblycpp_library_batch_tests
    PRIVATE ParallelAssemblyCppLibrary
  )

  parallelassemblycpp_add_python_test(
    tooling.test-data-audit
    300
    "${PROJECT_SOURCE_DIR}/unitTests/unitTester.py"
    --audit
  )
  parallelassemblycpp_add_python_test(
    tooling.repository-text-policy
    30
    "${PROJECT_SOURCE_DIR}/tools/check_repository.py"
  )
  parallelassemblycpp_add_python_test(
    tooling.benchmark-runner
    300
    -m unittest discover
    -s "${PROJECT_SOURCE_DIR}/benchmarks"
    -p "test_*.py"
  )
  parallelassemblycpp_add_python_test(
    tooling.benchmark-corpus
    300
    "${PROJECT_SOURCE_DIR}/benchmarks/benchmark.py"
    --list-cases
  )
  parallelassemblycpp_add_python_test(
    tooling.archived-benchmark-help
    30
    "${PROJECT_SOURCE_DIR}/output/benchmarks/amino-acid-cpp-rust-2026-08-26/run_benchmark.py"
    --help
  )

  set(
    PARALLELASSEMBLYCPP_REGRESSION_COMMAND
    "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/unitTests/unitTester.py"
    "$<TARGET_FILE:ParallelAssemblyCpp>"
    --jobs 2
    --timeout 120
  )
  if(NOT PARALLELASSEMBLYCPP_FULL_REGRESSION_TESTS)
    list(APPEND PARALLELASSEMBLYCPP_REGRESSION_COMMAND --limit 20)
  endif()
  add_test(NAME regression.manifest COMMAND ${PARALLELASSEMBLYCPP_REGRESSION_COMMAND})
  set_tests_properties(
    regression.manifest
    PROPERTIES
      LABELS "regression"
      PROCESSORS 2
      TIMEOUT 1800
      WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
  )

  if(TARGET ParallelAssemblyCppTelemetry)
    add_test(
      NAME regression.telemetry
      COMMAND
        "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/unitTests/unitTester.py"
        "$<TARGET_FILE:ParallelAssemblyCppTelemetry>"
        --jobs 1
        --limit 1
        --timeout 120
    )
    set_tests_properties(
      regression.telemetry
      PROPERTIES
        LABELS "regression;telemetry"
        TIMEOUT 300
        WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
    )
  endif()

  if(TARGET ParallelAssemblyCppOMP)
    add_test(
      NAME regression.parallel-solver
      COMMAND
        "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/unitTests/parallelSolverTester.py"
        --serial "$<TARGET_FILE:ParallelAssemblyCpp>"
        --openmp "$<TARGET_FILE:ParallelAssemblyCppOMP>"
        --repetitions 3
        --timeout 30
    )
    set_tests_properties(
      regression.parallel-solver
      PROPERTIES
        LABELS "regression;parallel"
        PROCESSORS 4
        SKIP_RETURN_CODE 77
        TIMEOUT 300
        WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
    )
  endif()

  if(TARGET ParallelAssemblyCppOMPTelemetry)
    add_test(
      NAME regression.parallel-telemetry
      COMMAND
        "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/unitTests/parallelSolverTester.py"
        --serial "$<TARGET_FILE:ParallelAssemblyCpp>"
        --openmp-telemetry "$<TARGET_FILE:ParallelAssemblyCppOMPTelemetry>"
        --telemetry-repetitions 2
        --timeout 30
    )
    set_tests_properties(
      regression.parallel-telemetry
      PROPERTIES
        LABELS "regression;parallel;telemetry"
        PROCESSORS 4
        SKIP_RETURN_CODE 77
        TIMEOUT 300
        WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
    )
  endif()

  set(PARALLELASSEMBLYCPP_DISTRIBUTED_PARITY_OPTIONS)
  if(TARGET ParallelAssemblyCppMPI)
    list(
      APPEND PARALLELASSEMBLYCPP_DISTRIBUTED_PARITY_OPTIONS
      --mpi "$<TARGET_FILE:ParallelAssemblyCppMPI>"
    )
  endif()
  if(TARGET ParallelAssemblyCppHybrid)
    list(
      APPEND PARALLELASSEMBLYCPP_DISTRIBUTED_PARITY_OPTIONS
      --hybrid "$<TARGET_FILE:ParallelAssemblyCppHybrid>"
    )
  endif()
  list(
    LENGTH PARALLELASSEMBLYCPP_DISTRIBUTED_PARITY_OPTIONS
    PARALLELASSEMBLYCPP_DISTRIBUTED_PARITY_OPTION_COUNT
  )
  if(PARALLELASSEMBLYCPP_DISTRIBUTED_PARITY_OPTION_COUNT GREATER 0)
    add_test(
      NAME regression.parallel-distributed
      COMMAND
        "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/unitTests/parallelSolverTester.py"
        --serial "$<TARGET_FILE:ParallelAssemblyCpp>"
        ${PARALLELASSEMBLYCPP_DISTRIBUTED_PARITY_OPTIONS}
        --mpiexec "${MPIEXEC_EXECUTABLE}"
        "--mpiexec-numproc-flag=${MPIEXEC_NUMPROC_FLAG}"
        --timeout 30
    )
    set_tests_properties(
      regression.parallel-distributed
      PROPERTIES
        LABELS "regression;parallel;mpi"
        PROCESSORS 4
        SKIP_RETURN_CODE 77
        TIMEOUT 300
        WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
    )
  endif()

  set(PARALLELASSEMBLYCPP_DISTRIBUTED_TELEMETRY_OPTIONS)
  if(TARGET ParallelAssemblyCppMPITelemetry)
    list(
      APPEND PARALLELASSEMBLYCPP_DISTRIBUTED_TELEMETRY_OPTIONS
      --mpi-telemetry "$<TARGET_FILE:ParallelAssemblyCppMPITelemetry>"
    )
  endif()
  if(TARGET ParallelAssemblyCppHybridTelemetry)
    list(
      APPEND PARALLELASSEMBLYCPP_DISTRIBUTED_TELEMETRY_OPTIONS
      --hybrid-telemetry "$<TARGET_FILE:ParallelAssemblyCppHybridTelemetry>"
    )
  endif()
  list(
    LENGTH PARALLELASSEMBLYCPP_DISTRIBUTED_TELEMETRY_OPTIONS
    PARALLELASSEMBLYCPP_DISTRIBUTED_TELEMETRY_OPTION_COUNT
  )
  if(PARALLELASSEMBLYCPP_DISTRIBUTED_TELEMETRY_OPTION_COUNT GREATER 0)
    add_test(
      NAME regression.parallel-distributed-telemetry
      COMMAND
        "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/unitTests/parallelSolverTester.py"
        --serial "$<TARGET_FILE:ParallelAssemblyCpp>"
        ${PARALLELASSEMBLYCPP_DISTRIBUTED_TELEMETRY_OPTIONS}
        --mpiexec "${MPIEXEC_EXECUTABLE}"
        "--mpiexec-numproc-flag=${MPIEXEC_NUMPROC_FLAG}"
        --timeout 30
    )
    set_tests_properties(
      regression.parallel-distributed-telemetry
      PROPERTIES
        LABELS "regression;parallel;mpi;telemetry"
        PROCESSORS 4
        SKIP_RETURN_CODE 77
        TIMEOUT 300
        WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
    )
  endif()
endif()

set(CPACK_PACKAGE_NAME "parallelassemblycpp")
set(CPACK_PACKAGE_VENDOR "ELIFE-ASU")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "${PROJECT_DESCRIPTION}")
set(CPACK_PACKAGE_HOMEPAGE_URL "${PROJECT_HOMEPAGE_URL}")
set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/License.md")
set(CPACK_RESOURCE_FILE_README "${PROJECT_SOURCE_DIR}/README.md")
set(CPACK_PACKAGE_CHECKSUM SHA256)
set(CPACK_MONOLITHIC_INSTALL ON)
if(CMAKE_VERSION VERSION_GREATER_EQUAL 4.3)
  set(CPACK_ARCHIVE_UID 0)
  set(CPACK_ARCHIVE_GID 0)
endif()
set(
  CPACK_PRE_BUILD_SCRIPTS
  "${PROJECT_SOURCE_DIR}/cmake/NormalizeArchivePermissions.cmake"
)
set(CPACK_VERBATIM_VARIABLES ON)
if(PARALLELASSEMBLYCPP_X86_64_V3)
  set(PARALLELASSEMBLYCPP_PACKAGE_CPU_BASELINE "v3")
else()
  set(PARALLELASSEMBLYCPP_PACKAGE_CPU_BASELINE "portable")
endif()
set(
  CPACK_PACKAGE_FILE_NAME
  "parallelassemblycpp-${PROJECT_VERSION}-${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}-${PARALLELASSEMBLYCPP_PACKAGE_CPU_BASELINE}"
)
if(WIN32)
  set(CPACK_GENERATOR ZIP)
else()
  set(CPACK_GENERATOR TGZ)
endif()
set(CPACK_SOURCE_GENERATOR TGZ)
set(CPACK_SOURCE_PACKAGE_FILE_NAME "parallelassemblycpp-${PROJECT_VERSION}-source")
set(
  CPACK_SOURCE_IGNORE_FILES
  "/\\.git/"
  "/\\.agents/"
  "/\\.codex/"
  "/build/"
  "/dist/"
  "/stage/"
  "/_CPack_Packages/"
  "/\\.ruff_cache/"
  "/__pycache__/"
  "\\.py[cod]$"
  "/\\.vscode/"
  "/CMakeUserPresets\\.json$"
  "/memUsage$"
  "Telemetry\\.json$"
  "/unitTests/.*(Out|Pathway|MAs|Telemetry\\.json)$"
)

include(CPack)
