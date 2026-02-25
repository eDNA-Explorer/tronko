# CMakeLists.txt for Tronko Tests
cmake_minimum_required(VERSION 3.10)
project(TronkoTests)

# Set C standard
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Find required libraries
find_package(PkgConfig REQUIRED)
pkg_check_modules(ZLIB REQUIRED zlib)

# Include directories
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/../tronko-assign
    ${CMAKE_CURRENT_SOURCE_DIR}/../tronko-build
    ${CMAKE_CURRENT_SOURCE_DIR}/unity/src
)

# Unity framework
add_library(unity STATIC
    unity/src/unity.c
)

# Common source files for testing
set(TRONKO_SOURCES
    ../tronko-assign/readreference.c
    ../tronko-assign/assignment.c
    ../tronko-assign/placement.c
    ../tronko-assign/logger.c
    ../tronko-assign/crash_debug.c
    ../tronko-assign/symbol_resolver.c
    ../tronko-assign/options.c
    # Add other source files as needed
)

# Test executables
set(TEST_SOURCES
    unit/test_readreference.c
    unit/test_assignment.c
    unit/test_crash_debug.c
    integration/test_workflows.c
)

# Create test executables
foreach(test_file ${TEST_SOURCES})
    get_filename_component(test_name ${test_file} NAME_WE)
    add_executable(${test_name}
        ${test_file}
        ${TRONKO_SOURCES}
        test_runner.c
    )
    
    target_link_libraries(${test_name}
        unity
        ${ZLIB_LIBRARIES}
        m
        pthread
        rt
    )
    
    # Add test to CTest
    add_test(NAME ${test_name} COMMAND ${test_name})
endforeach()

# Enable testing
enable_testing()