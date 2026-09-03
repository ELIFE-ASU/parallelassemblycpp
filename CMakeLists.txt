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
  AssemblyCpp
  VERSION 5.0.0
  DESCRIPTION "Compute molecular assembly indices and pathways"
  HOMEPAGE_URL "https://github.com/ELIFE-ASU/assemblycpp-v5"
  LANGUAGES CXX
)

include(GNUInstallDirs)

set(
  ASSEMBLYCPP_EXECUTABLE_PERMISSIONS
  OWNER_READ OWNER_WRITE OWNER_EXECUTE
  GROUP_READ GROUP_EXECUTE
  WORLD_READ WORLD_EXECUTE
)
set(
  ASSEMBLYCPP_FILE_PERMISSIONS
  OWNER_READ OWNER_WRITE
  GROUP_READ
  WORLD_READ
)
set(
  CMAKE_INSTALL_DEFAULT_DIRECTORY_PERMISSIONS
  ${ASSEMBLYCPP_EXECUTABLE_PERMISSIONS}
)

option(
  ASSEMBLYCPP_BUILD_TELEMETRY
  "Build the instrumented AssemblyCppTelemetry executable"
  OFF
)
option(
  ASSEMBLYCPP_X86_64_V3
  "Target x86-64-v3 with POPCNT"
  OFF
)
option(
  ASSEMBLYCPP_ENABLE_LTO
  "Enable link-time optimization"
  OFF
)
option(
  ASSEMBLYCPP_STRICT_WARNINGS
  "Treat high-signal compiler warnings as errors"
  ON
)
set(
  ASSEMBLYCPP_PGO_PHASE
  "OFF"
  CACHE STRING
  "PGO phase: OFF, GENERATE, or USE"
)
set_property(
  CACHE ASSEMBLYCPP_PGO_PHASE
  PROPERTY STRINGS OFF GENERATE USE
)
set(
  ASSEMBLYCPP_PGO_PROFILE_DIR
  "${CMAKE_BINARY_DIR}/pgo-profiles"
  CACHE PATH
  "Shared PGO profile directory"
)
option(
  ASSEMBLYCPP_FULL_REGRESSION_TESTS
  "Run the full regression manifest"
  OFF
)
option(
  ASSEMBLYCPP_BUILD_OPENMP
  "Build the experimental OpenMP executable"
  OFF
)
option(
  ASSEMBLYCPP_BUILD_MPI
  "Build the experimental MPI executable"
  OFF
)
option(
  ASSEMBLYCPP_BUILD_HYBRID
  "Build the experimental MPI/OpenMP executable"
  OFF
)

set(_ASSEMBLYCPP_PGO_PHASES OFF GENERATE USE)
if(NOT ASSEMBLYCPP_PGO_PHASE IN_LIST _ASSEMBLYCPP_PGO_PHASES)
  message(
    FATAL_ERROR
    "ASSEMBLYCPP_PGO_PHASE must be OFF, GENERATE, or USE; got "
    "${ASSEMBLYCPP_PGO_PHASE}"
  )
endif()

if(ASSEMBLYCPP_PGO_PHASE STREQUAL "OFF")
  set(ASSEMBLYCPP_PGO_PROFILE_DIR_ABSOLUTE "")
else()
  if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    message(FATAL_ERROR "PGO requires GCC")
  endif()
  if(ASSEMBLYCPP_PGO_PROFILE_DIR STREQUAL "")
    message(FATAL_ERROR "ASSEMBLYCPP_PGO_PROFILE_DIR must not be empty")
  endif()
  cmake_path(
    ABSOLUTE_PATH ASSEMBLYCPP_PGO_PROFILE_DIR
    BASE_DIRECTORY "${CMAKE_BINARY_DIR}"
    NORMALIZE
    OUTPUT_VARIABLE ASSEMBLYCPP_PGO_PROFILE_DIR_ABSOLUTE
  )
  set(
    ASSEMBLYCPP_PGO_DIRECTORY_MARKER
    "${ASSEMBLYCPP_PGO_PROFILE_DIR_ABSOLUTE}/.assemblycpp-pgo-profile-directory"
  )

  if(ASSEMBLYCPP_PGO_PHASE STREQUAL "GENERATE")
    file(MAKE_DIRECTORY "${ASSEMBLYCPP_PGO_PROFILE_DIR_ABSOLUTE}")
    if(NOT EXISTS "${ASSEMBLYCPP_PGO_DIRECTORY_MARKER}")
      file(
        GLOB ASSEMBLYCPP_PGO_PROFILE_DIRECTORY_ENTRIES
        LIST_DIRECTORIES true
        "${ASSEMBLYCPP_PGO_PROFILE_DIR_ABSOLUTE}/*"
      )
      if(ASSEMBLYCPP_PGO_PROFILE_DIRECTORY_ENTRIES)
        message(
          FATAL_ERROR
          "PGO profile directory is non-empty but unmarked: "
          "${ASSEMBLYCPP_PGO_PROFILE_DIR_ABSOLUTE}"
        )
      endif()
      file(
        WRITE "${ASSEMBLYCPP_PGO_DIRECTORY_MARKER}"
        "Dedicated AssemblyCpp PGO profile directory.\n"
      )
    elseif(
      IS_DIRECTORY "${ASSEMBLYCPP_PGO_DIRECTORY_MARKER}"
      OR IS_SYMLINK "${ASSEMBLYCPP_PGO_DIRECTORY_MARKER}"
    )
      message(FATAL_ERROR "PGO profile marker must be a regular file")
    endif()
  else()
    if(
      NOT EXISTS "${ASSEMBLYCPP_PGO_DIRECTORY_MARKER}"
      OR IS_DIRECTORY "${ASSEMBLYCPP_PGO_DIRECTORY_MARKER}"
      OR IS_SYMLINK "${ASSEMBLYCPP_PGO_DIRECTORY_MARKER}"
    )
      message(
        FATAL_ERROR
        "PGO profile marker is missing or invalid; configure preset "
        "pgo-generate first"
      )
    endif()
    set(
      ASSEMBLYCPP_PGO_COMPLETION_FILE
      "${ASSEMBLYCPP_PGO_PROFILE_DIR_ABSOLUTE}/assemblycpp-pgo-complete.json"
    )
    if(NOT EXISTS "${ASSEMBLYCPP_PGO_COMPLETION_FILE}")
      message(
        FATAL_ERROR
        "PGO training is incomplete; build preset pgo-train first"
      )
    endif()
    file(READ "${ASSEMBLYCPP_PGO_COMPLETION_FILE}" ASSEMBLYCPP_PGO_COMPLETION)
    string(
      JSON ASSEMBLYCPP_PGO_COMPLETION_SCHEMA
      ERROR_VARIABLE ASSEMBLYCPP_PGO_COMPLETION_ERROR
      GET "${ASSEMBLYCPP_PGO_COMPLETION}" schema_version
    )
    if(
      ASSEMBLYCPP_PGO_COMPLETION_ERROR
      OR NOT ASSEMBLYCPP_PGO_COMPLETION_SCHEMA EQUAL 2
    )
      message(
        FATAL_ERROR
        "Invalid PGO training completion record: "
        "${ASSEMBLYCPP_PGO_COMPLETION_FILE}"
      )
    endif()
    string(
      JSON ASSEMBLYCPP_PGO_RECORDED_WEIGHTS_SHA256
      ERROR_VARIABLE ASSEMBLYCPP_PGO_WEIGHTS_ERROR
      GET "${ASSEMBLYCPP_PGO_COMPLETION}" weights sha256
    )
    file(
      SHA256 "${PROJECT_SOURCE_DIR}/benchmarks/pgo-training.tsv"
      ASSEMBLYCPP_PGO_CURRENT_WEIGHTS_SHA256
    )
    if(
      ASSEMBLYCPP_PGO_WEIGHTS_ERROR
      OR NOT ASSEMBLYCPP_PGO_RECORDED_WEIGHTS_SHA256
        STREQUAL ASSEMBLYCPP_PGO_CURRENT_WEIGHTS_SHA256
    )
      message(
        FATAL_ERROR
        "PGO training weights changed; rebuild preset pgo-train"
      )
    endif()
    set(
      ASSEMBLYCPP_PGO_MANIFEST
      "${PROJECT_SOURCE_DIR}/benchmarks/cases.tsv"
    )
    string(
      JSON ASSEMBLYCPP_PGO_RECORDED_MANIFEST_SHA256
      ERROR_VARIABLE ASSEMBLYCPP_PGO_MANIFEST_ERROR
      GET "${ASSEMBLYCPP_PGO_COMPLETION}" corpus manifest sha256
    )
    file(
      SHA256 "${ASSEMBLYCPP_PGO_MANIFEST}"
      ASSEMBLYCPP_PGO_CURRENT_MANIFEST_SHA256
    )
    if(
      ASSEMBLYCPP_PGO_MANIFEST_ERROR
      OR NOT ASSEMBLYCPP_PGO_RECORDED_MANIFEST_SHA256
        STREQUAL ASSEMBLYCPP_PGO_CURRENT_MANIFEST_SHA256
    )
      message(
        FATAL_ERROR
        "Benchmark manifest changed; rebuild preset pgo-train"
      )
    endif()
    string(
      JSON ASSEMBLYCPP_PGO_INPUT_COUNT
      ERROR_VARIABLE ASSEMBLYCPP_PGO_INPUT_COUNT_ERROR
      LENGTH "${ASSEMBLYCPP_PGO_COMPLETION}" corpus inputs
    )
    if(
      ASSEMBLYCPP_PGO_INPUT_COUNT_ERROR
      OR ASSEMBLYCPP_PGO_INPUT_COUNT LESS 1
    )
      message(FATAL_ERROR "Invalid PGO corpus input fingerprints")
    endif()
    math(EXPR ASSEMBLYCPP_PGO_LAST_INPUT "${ASSEMBLYCPP_PGO_INPUT_COUNT} - 1")
    set(ASSEMBLYCPP_PGO_CORPUS_FILES "${ASSEMBLYCPP_PGO_MANIFEST}")
    foreach(ASSEMBLYCPP_PGO_INPUT_INDEX RANGE ${ASSEMBLYCPP_PGO_LAST_INPUT})
      string(
        JSON ASSEMBLYCPP_PGO_RECORDED_INPUT_PATH
        ERROR_VARIABLE ASSEMBLYCPP_PGO_INPUT_PATH_ERROR
        GET "${ASSEMBLYCPP_PGO_COMPLETION}"
        corpus inputs ${ASSEMBLYCPP_PGO_INPUT_INDEX} path
      )
      string(
        JSON ASSEMBLYCPP_PGO_RECORDED_INPUT_SHA256
        ERROR_VARIABLE ASSEMBLYCPP_PGO_INPUT_SHA256_ERROR
        GET "${ASSEMBLYCPP_PGO_COMPLETION}"
        corpus inputs ${ASSEMBLYCPP_PGO_INPUT_INDEX} sha256
      )
      if(
        ASSEMBLYCPP_PGO_INPUT_PATH_ERROR
        OR ASSEMBLYCPP_PGO_INPUT_SHA256_ERROR
      )
        message(FATAL_ERROR "Invalid PGO corpus input fingerprint")
      endif()
      cmake_path(
        ABSOLUTE_PATH ASSEMBLYCPP_PGO_RECORDED_INPUT_PATH
        BASE_DIRECTORY "${PROJECT_SOURCE_DIR}"
        NORMALIZE
        OUTPUT_VARIABLE ASSEMBLYCPP_PGO_CURRENT_INPUT
      )
      if(
        NOT EXISTS "${ASSEMBLYCPP_PGO_CURRENT_INPUT}"
        OR IS_DIRECTORY "${ASSEMBLYCPP_PGO_CURRENT_INPUT}"
      )
        message(
          FATAL_ERROR
          "PGO input is missing: ${ASSEMBLYCPP_PGO_CURRENT_INPUT}; rebuild "
          "preset pgo-train"
        )
      endif()
      file(
        SHA256 "${ASSEMBLYCPP_PGO_CURRENT_INPUT}"
        ASSEMBLYCPP_PGO_CURRENT_INPUT_SHA256
      )
      if(
        NOT ASSEMBLYCPP_PGO_RECORDED_INPUT_SHA256
          STREQUAL ASSEMBLYCPP_PGO_CURRENT_INPUT_SHA256
      )
        message(
          FATAL_ERROR
          "PGO input changed: ${ASSEMBLYCPP_PGO_CURRENT_INPUT}; rebuild "
          "preset pgo-train"
        )
      endif()
      list(APPEND ASSEMBLYCPP_PGO_CORPUS_FILES "${ASSEMBLYCPP_PGO_CURRENT_INPUT}")
    endforeach()
    file(
      GLOB_RECURSE ASSEMBLYCPP_PGO_PROFILE_DATA
      LIST_DIRECTORIES false
      "${ASSEMBLYCPP_PGO_PROFILE_DIR_ABSOLUTE}/*.gcda"
    )
    if(NOT ASSEMBLYCPP_PGO_PROFILE_DATA)
      message(
        FATAL_ERROR
        "PGO data is missing from ${ASSEMBLYCPP_PGO_PROFILE_DIR_ABSOLUTE}; "
        "build preset pgo-train first"
      )
    endif()
    # Reconfigure when training inputs change to reject partial or stale data.
    set_property(
      DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
      "${ASSEMBLYCPP_PGO_COMPLETION_FILE}"
      "${PROJECT_SOURCE_DIR}/benchmarks/pgo-training.tsv"
      ${ASSEMBLYCPP_PGO_CORPUS_FILES}
    )
  endif()
endif()

if(ASSEMBLYCPP_ENABLE_LTO)
  include(CheckIPOSupported)
  check_ipo_supported(
    RESULT ASSEMBLYCPP_LTO_SUPPORTED
    OUTPUT ASSEMBLYCPP_LTO_ERROR
    LANGUAGES CXX
  )
  if(NOT ASSEMBLYCPP_LTO_SUPPORTED)
    message(FATAL_ERROR "LTO is unavailable: ${ASSEMBLYCPP_LTO_ERROR}")
  endif()
endif()

function(assemblycpp_configure_target target)
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

  if(ASSEMBLYCPP_STRICT_WARNINGS)
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
        "ASSEMBLYCPP_STRICT_WARNINGS supports GNU, Clang-compatible, and "
        "MSVC toolchains; disable it only for unsupported compiler diagnosis"
      )
    endif()
  endif()

  if(ASSEMBLYCPP_ENABLE_LTO)
    set_property(
      TARGET ${target}
      PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE
    )

    if(
      CMAKE_CXX_COMPILER_ID STREQUAL "GNU"
      AND ASSEMBLYCPP_PGO_PHASE STREQUAL "OFF"
    )
      # One-partition LTO outperformed GCC's default without PGO.
      set(
        _ASSEMBLYCPP_GNU_LTO_OPTIONS
        -flto-partition=one
        -fipa-pta
        -fno-semantic-interposition
        -fdevirtualize-at-ltrans
        -fno-plt
      )
      target_compile_options(
        ${target}
        PRIVATE ${_ASSEMBLYCPP_GNU_LTO_OPTIONS}
      )
      target_link_options(
        ${target}
        PRIVATE ${_ASSEMBLYCPP_GNU_LTO_OPTIONS}
      )
    endif()
  endif()

  if(WIN32)
    target_compile_definitions(${target} PRIVATE NOMINMAX WIN32_LEAN_AND_MEAN)
  endif()

  if(ASSEMBLYCPP_X86_64_V3)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
      include(CheckCXXCompilerFlag)
      check_cxx_compiler_flag("-march=x86-64-v3" ASSEMBLYCPP_HAS_X86_64_V3)
      check_cxx_compiler_flag("-mpopcnt" ASSEMBLYCPP_HAS_MPOPCNT)
      if(NOT ASSEMBLYCPP_HAS_X86_64_V3 OR NOT ASSEMBLYCPP_HAS_MPOPCNT)
        message(
          FATAL_ERROR
          "ASSEMBLYCPP_X86_64_V3 requires compiler support for "
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
      endif()
    else()
      message(
        FATAL_ERROR
        "ASSEMBLYCPP_X86_64_V3 is supported only by GCC-compatible and "
        "Clang-compatible compilers"
      )
    endif()
  endif()
endfunction()

function(assemblycpp_configure_test_target target)
  assemblycpp_configure_target(${target})

  # Release-style configurations define NDEBUG globally. Keep standard assert
  # checks active in test executables until they use an always-active helper.
  if(MSVC)
    target_compile_options(${target} PRIVATE /UNDEBUG)
  else()
    target_compile_options(${target} PRIVATE -UNDEBUG)
  endif()
endfunction()

function(assemblycpp_configure_pgo_target target)
  if(ASSEMBLYCPP_PGO_PHASE STREQUAL "GENERATE")
    set(
      _ASSEMBLYCPP_PGO_OPTIONS
      "-fprofile-generate=${ASSEMBLYCPP_PGO_PROFILE_DIR_ABSOLUTE}"
      "-fprofile-prefix-path=${CMAKE_BINARY_DIR}"
    )
  elseif(ASSEMBLYCPP_PGO_PHASE STREQUAL "USE")
    set(
      _ASSEMBLYCPP_PGO_OPTIONS
      "-fprofile-use=${ASSEMBLYCPP_PGO_PROFILE_DIR_ABSOLUTE}"
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
  target_compile_options(${target} PRIVATE ${_ASSEMBLYCPP_PGO_OPTIONS})
  target_link_options(${target} PRIVATE ${_ASSEMBLYCPP_PGO_OPTIONS})
endfunction()

function(assemblycpp_add_solver target)
  cmake_parse_arguments(
    PARSE_ARGV 1
    ASSEMBLYCPP_SOLVER
    "MPI;OPENMP;TELEMETRY"
    ""
    ""
  )
  if(ASSEMBLYCPP_SOLVER_UNPARSED_ARGUMENTS)
    message(
      FATAL_ERROR
      "Unknown options for ${target}: ${ASSEMBLYCPP_SOLVER_UNPARSED_ARGUMENTS}"
    )
  endif()

  add_executable(${target} v5/main.cpp)
  set(ASSEMBLYCPP_SOLVER_DEFINITIONS)
  set(ASSEMBLYCPP_SOLVER_LIBRARIES)
  if(ASSEMBLYCPP_SOLVER_MPI)
    list(APPEND ASSEMBLYCPP_SOLVER_DEFINITIONS ASSEMBLYCPP_USE_MPI)
    list(APPEND ASSEMBLYCPP_SOLVER_LIBRARIES MPI::MPI_CXX)
  endif()
  if(ASSEMBLYCPP_SOLVER_OPENMP)
    list(APPEND ASSEMBLYCPP_SOLVER_DEFINITIONS ASSEMBLYCPP_USE_OPENMP)
    list(APPEND ASSEMBLYCPP_SOLVER_LIBRARIES OpenMP::OpenMP_CXX)
  endif()
  if(ASSEMBLYCPP_SOLVER_TELEMETRY)
    list(APPEND ASSEMBLYCPP_SOLVER_DEFINITIONS ASSEMBLY_ENABLE_TELEMETRY)
  endif()
  if(ASSEMBLYCPP_SOLVER_DEFINITIONS)
    target_compile_definitions(
      ${target}
      PRIVATE ${ASSEMBLYCPP_SOLVER_DEFINITIONS}
    )
  endif()
  if(ASSEMBLYCPP_SOLVER_LIBRARIES)
    target_link_libraries(${target} PRIVATE ${ASSEMBLYCPP_SOLVER_LIBRARIES})
  endif()
  assemblycpp_configure_target(${target})
endfunction()

add_library(AssemblyCppLibrary STATIC v5/main.cpp)
add_library(AssemblyCpp::Library ALIAS AssemblyCppLibrary)
target_compile_definitions(
  AssemblyCppLibrary
  PRIVATE ASSEMBLYCPP_LIBRARY_BUILD ASSEMBLYCPP_NO_MAIN
)
target_include_directories(
  AssemblyCppLibrary
  PUBLIC
    "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/v5>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/assemblycpp>"
)
target_compile_features(AssemblyCppLibrary PUBLIC cxx_std_20)
set_target_properties(
  AssemblyCppLibrary
  PROPERTIES
    EXPORT_NAME Library
    OUTPUT_NAME assemblycpp
    POSITION_INDEPENDENT_CODE ON
)
assemblycpp_configure_target(AssemblyCppLibrary)

assemblycpp_add_solver(AssemblyCpp)
assemblycpp_configure_pgo_target(AssemblyCpp)

if(ASSEMBLYCPP_BUILD_OPENMP OR ASSEMBLYCPP_BUILD_HYBRID)
  find_package(OpenMP REQUIRED COMPONENTS CXX)
endif()
if(ASSEMBLYCPP_BUILD_MPI OR ASSEMBLYCPP_BUILD_HYBRID)
  find_package(MPI REQUIRED COMPONENTS CXX)
endif()

if(ASSEMBLYCPP_BUILD_OPENMP)
  assemblycpp_add_solver(AssemblyCppOMP OPENMP)
endif()

if(ASSEMBLYCPP_BUILD_MPI)
  assemblycpp_add_solver(AssemblyCppMPI MPI)
endif()

if(ASSEMBLYCPP_BUILD_HYBRID)
  assemblycpp_add_solver(AssemblyCppHybrid MPI OPENMP)
endif()

if(ASSEMBLYCPP_BUILD_TELEMETRY)
  assemblycpp_add_solver(AssemblyCppTelemetry TELEMETRY)

  if(ASSEMBLYCPP_BUILD_OPENMP)
    assemblycpp_add_solver(AssemblyCppOMPTelemetry OPENMP TELEMETRY)
  endif()

  if(ASSEMBLYCPP_BUILD_MPI)
    assemblycpp_add_solver(AssemblyCppMPITelemetry MPI TELEMETRY)
  endif()

  if(ASSEMBLYCPP_BUILD_HYBRID)
    assemblycpp_add_solver(AssemblyCppHybridTelemetry MPI OPENMP TELEMETRY)
  endif()
endif()

if(ASSEMBLYCPP_X86_64_V3 AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  # Internalize single-file targets while preserving public library symbols.
  foreach(
    ASSEMBLYCPP_WHOLE_PROGRAM_TARGET
    AssemblyCpp
    AssemblyCppLibrary
    AssemblyCppOMP
    AssemblyCppMPI
    AssemblyCppHybrid
    AssemblyCppTelemetry
    AssemblyCppOMPTelemetry
    AssemblyCppMPITelemetry
    AssemblyCppHybridTelemetry
  )
    if(TARGET ${ASSEMBLYCPP_WHOLE_PROGRAM_TARGET})
      target_compile_options(
        ${ASSEMBLYCPP_WHOLE_PROGRAM_TARGET}
        PRIVATE -fwhole-program
      )
    endif()
  endforeach()
endif()

install(
  TARGETS AssemblyCpp AssemblyCppLibrary
  EXPORT AssemblyCppTargets
  RUNTIME
    DESTINATION "${CMAKE_INSTALL_BINDIR}"
    PERMISSIONS ${ASSEMBLYCPP_EXECUTABLE_PERMISSIONS}
  ARCHIVE
    DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    PERMISSIONS ${ASSEMBLYCPP_FILE_PERMISSIONS}
  INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/assemblycpp"
)
install(
  FILES v5/assemblycpp.h
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/assemblycpp"
  PERMISSIONS ${ASSEMBLYCPP_FILE_PERMISSIONS}
)
install(
  FILES README.md License.md
  DESTINATION "${CMAKE_INSTALL_DATAROOTDIR}/doc/assemblycpp"
  PERMISSIONS ${ASSEMBLYCPP_FILE_PERMISSIONS}
)
install(
  EXPORT AssemblyCppTargets
  FILE AssemblyCppTargets.cmake
  NAMESPACE AssemblyCpp::
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/AssemblyCpp"
)
include(CMakePackageConfigHelpers)
configure_package_config_file(
  cmake/AssemblyCppConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/AssemblyCppConfig.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/AssemblyCpp"
)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/AssemblyCppConfigVersion.cmake"
  VERSION "${PROJECT_VERSION}"
  COMPATIBILITY SameMajorVersion
)
install(
  FILES
    "${CMAKE_CURRENT_BINARY_DIR}/AssemblyCppConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/AssemblyCppConfigVersion.cmake"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/AssemblyCpp"
)
if(UNIX AND NOT APPLE)
  # Remove absolute RPATHs injected by Conda compiler wrappers.
  if(IS_ABSOLUTE "${CMAKE_INSTALL_BINDIR}")
    set(
      ASSEMBLYCPP_INSTALLED_EXECUTABLE
      "\$ENV{DESTDIR}${CMAKE_INSTALL_BINDIR}/AssemblyCpp"
    )
  else()
    set(
      ASSEMBLYCPP_INSTALLED_EXECUTABLE
      "\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}/AssemblyCpp"
    )
  endif()
  install(
    CODE
      "file(RPATH_REMOVE FILE \"${ASSEMBLYCPP_INSTALLED_EXECUTABLE}\")"
  )
endif()

include(CTest)
if(ASSEMBLYCPP_FULL_REGRESSION_TESTS AND NOT BUILD_TESTING)
  message(
    FATAL_ERROR
    "ASSEMBLYCPP_FULL_REGRESSION_TESTS requires BUILD_TESTING=ON"
  )
endif()
if(BUILD_TESTING OR ASSEMBLYCPP_PGO_PHASE STREQUAL "GENERATE")
  find_package(Python3 3.10 REQUIRED COMPONENTS Interpreter)
endif()

if(ASSEMBLYCPP_PGO_PHASE STREQUAL "GENERATE")
  add_custom_target(
    pgo-train
    COMMAND
      "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/benchmarks/pgo_training.py"
      --executable "$<TARGET_FILE:AssemblyCpp>"
      --profile-dir "${ASSEMBLYCPP_PGO_PROFILE_DIR_ABSOLUTE}"
      --weights "${PROJECT_SOURCE_DIR}/benchmarks/pgo-training.tsv"
    DEPENDS AssemblyCpp
    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
    USES_TERMINAL
    VERBATIM
  )
endif()

function(assemblycpp_add_unit_test name target source)
  add_executable(${target} ${source})
  assemblycpp_configure_test_target(${target})
  add_test(NAME unit.${name} COMMAND ${target} ${ARGN})
  set_tests_properties(
    unit.${name}
    PROPERTIES
      LABELS "unit"
      TIMEOUT 300
  )
endfunction()

function(assemblycpp_add_python_test name timeout)
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
  assemblycpp_add_unit_test(
    active-word-mask
    assemblycpp_active_word_mask_tests
    unitTests/activeWordMaskTester.cpp
  )
  assemblycpp_add_unit_test(
    assembly-transposition-table
    assemblycpp_transposition_table_tests
    unitTests/assemblyTranspositionTableTester.cpp
  )

  add_executable(
    assemblycpp_tree_canon_tests
    unitTests/treeCanonTester.cpp
  )
  assemblycpp_configure_test_target(assemblycpp_tree_canon_tests)
  if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    add_test(
      NAME unit.tree-canonicalisation
      COMMAND
        sh -c "ulimit -s 64 && exec \"$1\"" _
        "$<TARGET_FILE:assemblycpp_tree_canon_tests>"
    )
  else()
    add_test(
      NAME unit.tree-canonicalisation
      COMMAND assemblycpp_tree_canon_tests
    )
  endif()
  set_tests_properties(
    unit.tree-canonicalisation
    PROPERTIES
      LABELS "unit"
      TIMEOUT 300
  )

  assemblycpp_add_unit_test(
    cyclic-canonicalisation
    assemblycpp_cyclic_canon_tests
    unitTests/cyclicCanonTester.cpp
  )
  assemblycpp_add_unit_test(
    search-stop-polling
    assemblycpp_search_stop_poll_tests
    unitTests/searchStopPollTester.cpp
  )
  assemblycpp_add_unit_test(
    fragmentation-metadata
    assemblycpp_fragmentation_metadata_tests
    unitTests/fragmentationMetadataTester.cpp
  )
  assemblycpp_add_unit_test(
    molfile-parser
    assemblycpp_molfile_parser_tests
    unitTests/molfileParserTester.cpp
  )
  assemblycpp_add_unit_test(
    native-graph-parser
    assemblycpp_graphio_tests
    unitTests/graphioTester.cpp
  )
  assemblycpp_add_unit_test(
    pathway-generator
    assemblycpp_pathway_generator_tests
    unitTests/pathwayGeneratorTester.cpp
  )
  assemblycpp_add_unit_test(
    library-batch
    assemblycpp_library_batch_tests
    unitTests/libraryBatchTester.cpp
    "${PROJECT_SOURCE_DIR}/unitTests/icosane.mol"
    "${PROJECT_SOURCE_DIR}/unitTests/sucrose.mol"
    "${PROJECT_SOURCE_DIR}/unitTests/butane.mol"
    "${PROJECT_SOURCE_DIR}/unitTests/graphio_test"
  )
  target_link_libraries(
    assemblycpp_library_batch_tests
    PRIVATE AssemblyCppLibrary
  )

  assemblycpp_add_python_test(
    tooling.test-data-audit
    300
    "${PROJECT_SOURCE_DIR}/unitTests/unitTester.py"
    --audit
  )
  assemblycpp_add_python_test(
    tooling.repository-text-policy
    30
    "${PROJECT_SOURCE_DIR}/tools/check_repository.py"
  )
  assemblycpp_add_python_test(
    tooling.benchmark-runner
    300
    -m unittest discover
    -s "${PROJECT_SOURCE_DIR}/benchmarks"
    -p "test_*.py"
  )
  assemblycpp_add_python_test(
    tooling.benchmark-corpus
    300
    "${PROJECT_SOURCE_DIR}/benchmarks/benchmark.py"
    --list-cases
  )
  assemblycpp_add_python_test(
    tooling.archived-benchmark-help
    30
    "${PROJECT_SOURCE_DIR}/output/benchmarks/amino-acid-cpp-rust-2026-08-26/run_benchmark.py"
    --help
  )

  set(
    ASSEMBLYCPP_REGRESSION_COMMAND
    "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/unitTests/unitTester.py"
    "$<TARGET_FILE:AssemblyCpp>"
    --jobs 2
    --timeout 120
  )
  if(NOT ASSEMBLYCPP_FULL_REGRESSION_TESTS)
    list(APPEND ASSEMBLYCPP_REGRESSION_COMMAND --limit 20)
  endif()
  add_test(NAME regression.manifest COMMAND ${ASSEMBLYCPP_REGRESSION_COMMAND})
  set_tests_properties(
    regression.manifest
    PROPERTIES
      LABELS "regression"
      PROCESSORS 2
      TIMEOUT 1800
      WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
  )

  if(TARGET AssemblyCppTelemetry)
    add_test(
      NAME regression.telemetry
      COMMAND
        "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/unitTests/unitTester.py"
        "$<TARGET_FILE:AssemblyCppTelemetry>"
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

  if(TARGET AssemblyCppOMP)
    add_test(
      NAME regression.parallel-solver
      COMMAND
        "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/unitTests/parallelSolverTester.py"
        --serial "$<TARGET_FILE:AssemblyCpp>"
        --openmp "$<TARGET_FILE:AssemblyCppOMP>"
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

  if(TARGET AssemblyCppOMPTelemetry)
    add_test(
      NAME regression.parallel-telemetry
      COMMAND
        "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/unitTests/parallelSolverTester.py"
        --serial "$<TARGET_FILE:AssemblyCpp>"
        --openmp-telemetry "$<TARGET_FILE:AssemblyCppOMPTelemetry>"
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

  set(ASSEMBLYCPP_DISTRIBUTED_PARITY_OPTIONS)
  if(TARGET AssemblyCppMPI)
    list(
      APPEND ASSEMBLYCPP_DISTRIBUTED_PARITY_OPTIONS
      --mpi "$<TARGET_FILE:AssemblyCppMPI>"
    )
  endif()
  if(TARGET AssemblyCppHybrid)
    list(
      APPEND ASSEMBLYCPP_DISTRIBUTED_PARITY_OPTIONS
      --hybrid "$<TARGET_FILE:AssemblyCppHybrid>"
    )
  endif()
  list(
    LENGTH ASSEMBLYCPP_DISTRIBUTED_PARITY_OPTIONS
    ASSEMBLYCPP_DISTRIBUTED_PARITY_OPTION_COUNT
  )
  if(ASSEMBLYCPP_DISTRIBUTED_PARITY_OPTION_COUNT GREATER 0)
    add_test(
      NAME regression.parallel-distributed
      COMMAND
        "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/unitTests/parallelSolverTester.py"
        --serial "$<TARGET_FILE:AssemblyCpp>"
        ${ASSEMBLYCPP_DISTRIBUTED_PARITY_OPTIONS}
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

  set(ASSEMBLYCPP_DISTRIBUTED_TELEMETRY_OPTIONS)
  if(TARGET AssemblyCppMPITelemetry)
    list(
      APPEND ASSEMBLYCPP_DISTRIBUTED_TELEMETRY_OPTIONS
      --mpi-telemetry "$<TARGET_FILE:AssemblyCppMPITelemetry>"
    )
  endif()
  if(TARGET AssemblyCppHybridTelemetry)
    list(
      APPEND ASSEMBLYCPP_DISTRIBUTED_TELEMETRY_OPTIONS
      --hybrid-telemetry "$<TARGET_FILE:AssemblyCppHybridTelemetry>"
    )
  endif()
  list(
    LENGTH ASSEMBLYCPP_DISTRIBUTED_TELEMETRY_OPTIONS
    ASSEMBLYCPP_DISTRIBUTED_TELEMETRY_OPTION_COUNT
  )
  if(ASSEMBLYCPP_DISTRIBUTED_TELEMETRY_OPTION_COUNT GREATER 0)
    add_test(
      NAME regression.parallel-distributed-telemetry
      COMMAND
        "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/unitTests/parallelSolverTester.py"
        --serial "$<TARGET_FILE:AssemblyCpp>"
        ${ASSEMBLYCPP_DISTRIBUTED_TELEMETRY_OPTIONS}
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

set(CPACK_PACKAGE_NAME "assemblycpp")
set(CPACK_PACKAGE_VENDOR "University of Glasgow")
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
if(ASSEMBLYCPP_X86_64_V3)
  set(ASSEMBLYCPP_PACKAGE_CPU_BASELINE "v3")
else()
  set(ASSEMBLYCPP_PACKAGE_CPU_BASELINE "portable")
endif()
set(
  CPACK_PACKAGE_FILE_NAME
  "assemblycpp-${PROJECT_VERSION}-${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}-${ASSEMBLYCPP_PACKAGE_CPU_BASELINE}"
)
if(WIN32)
  set(CPACK_GENERATOR ZIP)
else()
  set(CPACK_GENERATOR TGZ)
endif()
set(CPACK_SOURCE_GENERATOR TGZ)
set(CPACK_SOURCE_PACKAGE_FILE_NAME "assemblycpp-${PROJECT_VERSION}-source")
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
