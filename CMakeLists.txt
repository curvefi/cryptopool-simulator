cmake_minimum_required(VERSION 3.19)
project(CryptoPoolSimulator)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Release build flags matching the original Makefile
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -funroll-loops -fomit-frame-pointer -ffast-math -march=native")

# Debug build flags for simud target
set(CMAKE_CXX_FLAGS_DEBUG "-g -fsanitize=address")

# Find pthreads
find_package(Threads REQUIRED)

# Default target: simu (matches original Makefile default)
add_executable(simu main.cpp)
target_link_libraries(simu Threads::Threads)
set_target_properties(simu PROPERTIES
    CMAKE_BUILD_TYPE Release
)

# Debug version of simu
add_executable(simud main.cpp)
target_link_libraries(simud Threads::Threads)
set_target_properties(simud PROPERTIES
    CMAKE_BUILD_TYPE Debug
)

# Stableswap simulator
add_executable(simus main-stableswap.cpp)
target_link_libraries(simus Threads::Threads)
set_target_properties(simus PROPERTIES
    CMAKE_BUILD_TYPE Release
)

# Dynamic rate stableswap simulator
add_executable(simudyn main-stableswap-dynrate.cpp)
target_link_libraries(simudyn Threads::Threads)
set_target_properties(simudyn PROPERTIES
    CMAKE_BUILD_TYPE Release
)

# Legacy simulator
add_executable(simul main-legacy.cpp)
target_link_libraries(simul Threads::Threads)
set_target_properties(simul PROPERTIES
    CMAKE_BUILD_TYPE Release
)

# Set the default target
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT simu)