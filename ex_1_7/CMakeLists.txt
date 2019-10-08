# Copyright Tom Westerhout (c) 2019
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#
#     * Neither the name of Tom Westerhout nor the names of other
#       contributors may be used to endorse or promote products derived
#       from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

cmake_minimum_required(VERSION 3.9)
cmake_policy(SET CMP0074 NEW)

project(ex_1_7 VERSION 0.1
               DESCRIPTION "Solution to Exercise 1.7"
               LANGUAGES C)

# detect if we're being bundled, disable test and install targets in that case
if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
    set(NOT_SUBPROJECT ON)
else()
    set(NOT_SUBPROJECT OFF)
endif()

include(CTest)
include(CheckCCompilerFlag)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_C_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    message(STATUS "Setting build type to RelWithDebInfo, use \
-DCMAKE_BUILD_TYPE= to override."
    )
    set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()
string(TOUPPER "${CMAKE_BUILD_TYPE}" _EX17_BUILD_TYPE)

if(_EX17_BUILD_TYPE STREQUAL "DEBUG")
    set(_EX17_ASSERTS_ARE_OKAY ON)
else()
    set(_EX17_ASSERTS_ARE_OKAY OFF)
endif()

if(BUILD_TESTING AND NOT_SUBPROJECT)
    set(_EX17_TESTS_ARE_OKAY ON)
    set(_EX17_EXAMPLES_ARE_OKAY ON)
else()
    set(_EX17_TESTS_ARE_OKAY OFF)
    set(_EX17_EXAMPLES_ARE_OKAY ON)
endif()

option(EX17_USE_VALGRIND "Perform tests with Valgrind" OFF)
option(EX17_USE_BLAS "Use BLAS for 1D vector operations" OFF)
option(EX17_BUILD_TESTING "Build tests" ${_EX17_TESTS_ARE_OKAY})
option(EX17_BUILD_EXAMPLES "Build examples" ${_EX17_EXAMPLES_ARE_OKAY})
option(EX17_BUILD_BENCHMARKS "Build benchmarks" OFF)
option(EX17_ENABLE_COVERAGE "Generate coverage for codecov.io" OFF)
option(EX17_INSTALL_DOCS "Install documentation alongside library" ON)
option(EX17_DEBUG "Include assertions" ${_EX17_ASSERTS_ARE_OKAY})


# Warnings
################################################################################
if("${CMAKE_C_COMPILER_ID}" STREQUAL "Clang")
    set(EX17_WARNING_FLAGS -Weverything)
elseif("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
    set(EX17_WARNING_FLAGS -pedantic -W -Wall -Wextra)
endif()

find_package(MPI)
if(NOT MPI_C_FOUND)
    message(FATAL_ERROR "MPI is required for compiling and running the code")
endif()


function(setup_ex17_target MY_TARGET)
    target_compile_features(${MY_TARGET} PUBLIC c_std_11)
    target_compile_options(${MY_TARGET} PRIVATE "${EX17_WARNING_FLAGS}")

    if("${CMAKE_C_COMPILER_ID}" STREQUAL "Clang")
        target_compile_options(${MY_TARGET} PUBLIC -fcolor-diagnostics)
    elseif("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
        target_compile_options(${MY_TARGET} PUBLIC -fdiagnostics-color=always)
    endif()

    check_c_compiler_flag("-fvisibility=hidden" EX17_COMPILER_OPT_HIDDEN_SUPPORTED)
    if(EX17_COMPILER_OPT_HIDDEN_SUPPORTED)
        message(STATUS "[ex_1_7] Using -fvisibility=hidden.")
        target_compile_options(${MY_TARGET} PRIVATE "-fvisibility=hidden")
    endif()

    target_include_directories(${MY_TARGET} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
    target_link_libraries(${MY_TARGET} PUBLIC m MPI::MPI_C)
endfunction()

# Main target
set(ex17_HEADERS "include/ex17.h")
set(ex17_SOURCES "src/ex17.c")
# add_library(ex17 SHARED ${ex17_HEADERS} ${ex17_SOURCES})
# target_include_directories(ex17 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
# setup_ex17_target(ex17)
# check_c_compiler_flag("-fvisibility=hidden" EX17_COMPILER_OPT_HIDDEN_SUPPORTED)
# if(EX17_COMPILER_OPT_HIDDEN_SUPPORTED)
#     message(STATUS "[ex_1_7] Using -fvisibility=hidden.")
#     target_compile_options(ex17 PRIVATE "-fvisibility=hidden")
# endif()
# target_link_libraries(ex17 PUBLIC m MPI::MPI_C)

add_executable(generate ${ex17_HEADERS} ${ex17_SOURCES} src/generate.c)
setup_ex17_target(generate)
if(EX17_DEBUG)
    target_compile_definitions(generate PUBLIC EX17_DEBUG=1)
endif()
# target_link_libraries(generate PRIVATE ex17)

add_executable(count include/ex17.h src/count.c src/ex17.c)
target_compile_definitions(count PRIVATE EX17_COUNTERS=1 EX17_DEBUG=0)
setup_ex17_target(count)
# target_include_directories(count PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
# setup_ex17_target(ex17)
# check_c_compiler_flag("-fvisibility=hidden" EX17_COMPILER_OPT_HIDDEN_SUPPORTED)
# if(EX17_COMPILER_OPT_HIDDEN_SUPPORTED)
#     message(STATUS "[ex_1_7] Using -fvisibility=hidden.")
#     target_compile_options(ex17 PRIVATE "-fvisibility=hidden")
# endif()
# target_link_libraries(ex17 PUBLIC m MPI::MPI_C)
# target_link_libraries(PRIVATE ex17)
# target_include_directories(generate PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
# setup_ex17_target(ex17)
# check_c_compiler_flag("-fvisibility=hidden" EX17_COMPILER_OPT_HIDDEN_SUPPORTED)
# if(EX17_COMPILER_OPT_HIDDEN_SUPPORTED)
#     message(STATUS "[ex_1_7] Using -fvisibility=hidden.")
#     target_compile_options(ex17 PRIVATE "-fvisibility=hidden")
# endif()
# target_link_libraries(ex17 PRIVATE m MPI::MPI_C)
# if (LBFGS_USE_BLAS)
#     find_package(CBLAS REQUIRED)
#     target_compile_definitions(lbfgs INTERFACE LBFGS_USE_BLAS=1)
#     target_link_libraries(lbfgs PUBLIC ${CBLAS_LIBRARIES})
#     target_include_directories(lbfgs PUBLIC ${CBLAS_INCLUDE_DIR})
# endif()

# if (LBFGS_BUILD_TESTING)
#     if (NOT TARGET Catch2::Catch2)
#         option(CATCH_BUILD_TESTING OFF)
#         option(CATCH_BUILD_EXAMPLES OFF)
#         option(CATCH_BUILD_EXTRA_TESTS OFF)
#         add_subdirectory(third_party/Catch2)
#     endif()
#     enable_testing()
# endif()

# if (LBFGS_BUILD_EXAMPLES)
#     add_executable(line_search example/line_search.cpp)
#     target_link_libraries(line_search PUBLIC lbfgs)
# 
#     add_executable(lbfgs_example example/lbfgs.cpp)
#     target_link_libraries(lbfgs_example PUBLIC lbfgs)
# endif()
# 
# if (LBFGS_BUILD_TESTING)
#     add_subdirectory(test)
# endif()
