cmake_minimum_required(VERSION 3.10)
project(xbar_simulator_python)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find Python
find_package(Python COMPONENTS Interpreter Development REQUIRED)

# Add pybind11
add_subdirectory(pybind11)

# Include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include/eigen)

# Define Eigen as a system include to suppress warnings
include_directories(SYSTEM ${CMAKE_CURRENT_SOURCE_DIR}/include/eigen)

# Source files for the crossbar simulator
file(GLOB CORE_SRC "src/core/*.cpp")
file(GLOB CROSSBAR_SRC "src/crossbar_model/*.cpp")
file(GLOB MEMRISTOR_SRC "src/memristor_model/*.cpp")

# Create the pybind11 module
pybind11_add_module(xbar_simulator 
    src/python/xbar_wrapper.cpp
    ${CORE_SRC}
    ${CROSSBAR_SRC}
    ${MEMRISTOR_SRC}
)

# Link against filesystem library if needed (for C++17 std::filesystem)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    target_link_libraries(xbar_simulator PRIVATE stdc++fs)
endif()

# Set output directory
set_target_properties(xbar_simulator PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python
)

# Installation
install(TARGETS xbar_simulator DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/python) 