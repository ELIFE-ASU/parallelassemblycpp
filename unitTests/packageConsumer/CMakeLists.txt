cmake_minimum_required(VERSION 3.25)
project(ParallelAssemblyCppPackageConsumer LANGUAGES CXX)

find_package(ParallelAssemblyCpp 0.1.0 CONFIG REQUIRED)

add_executable(package_consumer ../libraryBatchTester.cpp)
target_link_libraries(package_consumer PRIVATE ParallelAssemblyCpp::Library)

enable_testing()
add_test(
  NAME installed-library
  COMMAND
    package_consumer
    "${CMAKE_CURRENT_SOURCE_DIR}/../icosane.mol"
    "${CMAKE_CURRENT_SOURCE_DIR}/../sucrose.mol"
    "${CMAKE_CURRENT_SOURCE_DIR}/../butane.mol"
    "${CMAKE_CURRENT_SOURCE_DIR}/../graphio_test"
)
