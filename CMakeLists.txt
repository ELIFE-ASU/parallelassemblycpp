cmake_minimum_required(VERSION 3.25)

if(POLICY CMP0206)
  cmake_policy(SET CMP0206 NEW)
endif()

project(
  AssemblyCpp
  VERSION 5.0.0
  DESCRIPTION "Calculate molecular assembly indices and pathways"
  HOMEPAGE_URL "https://github.com/ELIFE-ASU/assemblycpp-v5"
  LANGUAGES CXX
)

include(GNUInstallDirs)

set(
  CMAKE_INSTALL_DEFAULT_DIRECTORY_PERMISSIONS
  OWNER_READ OWNER_WRITE OWNER_EXECUTE
  GROUP_READ GROUP_EXECUTE
  WORLD_READ WORLD_EXECUTE
)

option(
  ASSEMBLYCPP_BUILD_TELEMETRY
  "Build the separately instrumented AssemblyCppTelemetry executable"
  OFF
)
option(
  ASSEMBLYCPP_X86_64_V3
  "Target x86-64-v3 with POPCNT instead of the portable compiler baseline"
  OFF
)
option(
  ASSEMBLYCPP_FULL_REGRESSION_TESTS
  "Run the complete regression manifest instead of its first 20 cases"
  OFF
)

function(assemblycpp_configure_target target)
  target_compile_features(${target} PRIVATE cxx_std_20)
  set_target_properties(
    ${target}
    PROPERTIES
      CXX_EXTENSIONS OFF
      CXX_STANDARD_REQUIRED ON
      INSTALL_REMOVE_ENVIRONMENT_RPATH ON
  )

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
        # The search is intentionally template-heavy. Its filtered and plain
        # matching loops benefit from a slightly larger whole-unit inlining
        # budget in the performance preset.
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

add_executable(AssemblyCpp v5/main.cpp)
assemblycpp_configure_target(AssemblyCpp)

if(ASSEMBLYCPP_BUILD_TELEMETRY)
  add_executable(AssemblyCppTelemetry v5/main.cpp)
  target_compile_definitions(
    AssemblyCppTelemetry
    PRIVATE ASSEMBLY_ENABLE_TELEMETRY
  )
  assemblycpp_configure_target(AssemblyCppTelemetry)
endif()

install(
  TARGETS AssemblyCpp
  RUNTIME
    DESTINATION "${CMAKE_INSTALL_BINDIR}"
    PERMISSIONS
      OWNER_READ OWNER_WRITE OWNER_EXECUTE
      GROUP_READ GROUP_EXECUTE
      WORLD_READ WORLD_EXECUTE
)
install(
  FILES README.md License.md
  DESTINATION "${CMAKE_INSTALL_DATAROOTDIR}/doc/assemblycpp"
  PERMISSIONS
    OWNER_READ OWNER_WRITE
    GROUP_READ
    WORLD_READ
)
if(UNIX AND NOT APPLE)
  # Conda compiler wrappers can inject an absolute environment RPATH outside
  # CMake's normal install-RPATH handling. Never preserve it in an install or
  # archive intended to leave that environment.
  install(
    CODE
      "file(RPATH_REMOVE FILE \"\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}/AssemblyCpp\")"
  )
endif()

include(CTest)
if(BUILD_TESTING)
  find_package(Python3 3.10 REQUIRED COMPONENTS Interpreter)

  add_executable(
    assemblycpp_active_word_mask_tests
    unitTests/activeWordMaskTester.cpp
  )
  assemblycpp_configure_target(assemblycpp_active_word_mask_tests)

  add_executable(
    assemblycpp_tree_canon_tests
    unitTests/treeCanonTester.cpp
  )
  assemblycpp_configure_target(assemblycpp_tree_canon_tests)

  add_executable(
    assemblycpp_cyclic_canon_tests
    unitTests/cyclicCanonTester.cpp
  )
  assemblycpp_configure_target(assemblycpp_cyclic_canon_tests)

  add_executable(
    assemblycpp_search_stop_poll_tests
    unitTests/searchStopPollTester.cpp
  )
  assemblycpp_configure_target(assemblycpp_search_stop_poll_tests)

  add_test(
    NAME unit.active-word-mask
    COMMAND assemblycpp_active_word_mask_tests
  )
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
  add_test(
    NAME unit.cyclic-canonicalisation
    COMMAND assemblycpp_cyclic_canon_tests
  )
  add_test(
    NAME unit.search-stop-polling
    COMMAND assemblycpp_search_stop_poll_tests
  )
  set_tests_properties(
    unit.active-word-mask
    unit.tree-canonicalisation
    unit.cyclic-canonicalisation
    unit.search-stop-polling
    PROPERTIES LABELS "unit"
  )

  add_test(
    NAME tooling.test-data-audit
    COMMAND
      "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/unitTests/unitTester.py"
      --audit
  )
  set_tests_properties(
    tooling.test-data-audit
    PROPERTIES
      LABELS "tooling"
      WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
  )

  add_test(
    NAME tooling.benchmark-runner
    COMMAND
      "${Python3_EXECUTABLE}"
      -m unittest discover
      -s "${PROJECT_SOURCE_DIR}/benchmarks"
      -p "test_*.py"
  )
  set_tests_properties(
    tooling.benchmark-runner
    PROPERTIES
      LABELS "tooling"
      WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
  )

  add_test(
    NAME tooling.benchmark-corpus
    COMMAND
      "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/benchmarks/benchmark.py"
      --list-cases
  )
  set_tests_properties(
    tooling.benchmark-corpus
    PROPERTIES
      LABELS "tooling"
      WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
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
endif()

set(CPACK_PACKAGE_NAME "assemblycpp")
set(CPACK_PACKAGE_VENDOR "University of Glasgow")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "${PROJECT_DESCRIPTION}")
set(CPACK_PACKAGE_HOMEPAGE_URL "${PROJECT_HOMEPAGE_URL}")
set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/License.md")
set(CPACK_RESOURCE_FILE_README "${PROJECT_SOURCE_DIR}/README.md")
set(CPACK_PACKAGE_CHECKSUM SHA256)
set(CPACK_MONOLITHIC_INSTALL ON)
set(CPACK_ARCHIVE_UID 0)
set(CPACK_ARCHIVE_GID 0)
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
