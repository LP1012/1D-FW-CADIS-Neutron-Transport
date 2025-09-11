cmake_minimum_required(VERSION 3.10)

# Project name and language
project(MCSlab LANGUAGES CXX)

# Require C++14 (you can stick to 11 if you prefer)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Define the executable and source files
add_executable(MCSlab
    src/MCSlab.cpp
    src/Neutron.cpp
    src/Slab.cpp
    src/main.cpp
    # add other source files here
)

# Tell the compiler where to find headers
target_include_directories(MCSlab
    PUBLIC ${CMAKE_SOURCE_DIR}/include
)

