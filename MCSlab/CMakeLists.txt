cmake_minimum_required(VERSION 3.14)

project(MCSlab LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/heads/main.zip)
FetchContent_MakeAvailable(googletest)

# Core library
add_library(
  MCLib src/MCSlab.cpp src/Neutron.cpp src/Region.cpp src/RNG.cpp
        src/tinyxml2.cpp # add other source files here
)
target_include_directories(MCLib PUBLIC ${CMAKE_SOURCE_DIR}/include)

# Main executable
add_executable(MCSlab src/main.cpp)
target_link_libraries(MCSlab PRIVATE MCLib)

# Enable testing
enable_testing()

# Tests
add_executable(run_tests tests/tests.cpp tests/neutron_test.cpp
                         tests/region_test.cpp tests/simulation_test.cpp)
target_link_libraries(run_tests PRIVATE MCLib GTest::gtest_main)

include(GoogleTest)
gtest_discover_tests(run_tests)
