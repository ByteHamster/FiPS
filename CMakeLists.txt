cmake_minimum_required (VERSION 3.25)
cmake_policy(SET CMP0048 NEW)
project(FiPS LANGUAGES CXX)

if(TARGET FiPS)
    return()
endif()

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif ()

if((CMAKE_BUILD_TYPE STREQUAL "Release" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
        AND PROJECT_IS_TOP_LEVEL)
    add_compile_options(-march=native)
endif()

# ---------------------------- Dependencies ----------------------------
if(NOT TARGET tlx)
    add_subdirectory(extlib/tlx)
endif()

if(NOT TARGET ips2ra)
    set(IPS2RA_DISABLE_PARALLEL ON CACHE PATH "ips2ra's FindTBB greps a file that does not exist in recent TBB versions")
    add_subdirectory(extlib/ips2ra SYSTEM)
    find_package(TBB)
    target_compile_options(ips2ra INTERFACE -D_REENTRANT)
    target_link_libraries(ips2ra INTERFACE pthread atomic TBB::tbb)
endif()

if(NOT TARGET ByteHamster::Util)
    add_subdirectory(extlib/util)
endif()

# ---------------------------- Library Setup ----------------------------

add_library(FiPS INTERFACE)
target_include_directories(FiPS INTERFACE include)
target_compile_features(FiPS INTERFACE cxx_std_20)
target_link_libraries(FiPS INTERFACE ips2ra ByteHamster::Util)
add_library(FiPS::fips ALIAS FiPS)

# ---------------------------- Benchmarks ----------------------------
if(PROJECT_IS_TOP_LEVEL)
    add_library(BenchmarkUtils INTERFACE)
    target_include_directories(BenchmarkUtils INTERFACE benchmark)
    target_link_libraries(BenchmarkUtils INTERFACE tlx ByteHamster::Util)

    # Warnings if this is the main project
    target_compile_options(FiPS INTERFACE $<$<COMPILE_LANGUAGE:CXX>:-Wall -Wextra -Wpedantic -Werror -frecord-gcc-switches>)

    add_executable(Benchmark benchmark/benchmark_construction.cpp)
    target_link_libraries(Benchmark PUBLIC BenchmarkUtils FiPS)
endif()
