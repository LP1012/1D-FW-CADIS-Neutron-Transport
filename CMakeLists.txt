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
  FW-CADIS_lib 
  src/FW-CADIS.cpp 
  src/SN.cpp 
  src/CellBase.cpp 
  src/SNCell.cpp 
  src/Cell.cpp 
  src/MCSlab.cpp 
  src/Neutron.cpp 
  src/Region.cpp 
  src/RNG.cpp
  src/CellQuadrature.cpp
  src/tinyxml2.cpp # add other source files here
)
target_include_directories(FW-CADIS_lib PUBLIC ${CMAKE_SOURCE_DIR}/include)

# Main executable
add_executable(1D-FW-CADIS src/main.cpp)
target_link_libraries(1D-FW-CADIS PRIVATE FW-CADIS_lib)

# Enable testing
enable_testing()

# Tests
add_executable(run_tests tests/neutron_test.cpp tests/region_test.cpp
                         tests/simulation_test.cpp)
target_link_libraries(run_tests PRIVATE FW-CADIS_lib GTest::gtest_main)

include(GoogleTest)
gtest_discover_tests(run_tests)
