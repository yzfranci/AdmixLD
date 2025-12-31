cmake_minimum_required(VERSION 3.16)
project(adfinder LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()

# HTSlib via pkg-config
find_package(PkgConfig REQUIRED)
pkg_check_modules(HTS REQUIRED htslib)

# Eigen
find_package(Eigen3 REQUIRED)

add_executable(adfinder
	src/main.cpp
	src/io/vcf_windows.cpp
	src/core/hybrid_index.cpp
	src/core/residualize.cpp
	src/core/scan_blocks.cpp
)

# HTSlib
target_include_directories(adfinder PRIVATE ${HTS_INCLUDE_DIRS})
target_link_directories(adfinder PRIVATE ${HTS_LIBRARY_DIRS})
target_link_libraries(adfinder PRIVATE ${HTS_LIBRARIES})
target_compile_options(adfinder PRIVATE ${HTS_CFLAGS_OTHER})

# Eigen
target_link_libraries(adfinder PRIVATE Eigen3::Eigen)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
	message(STATUS "OpenMP found")
	target_link_libraries(adfinder PRIVATE OpenMP::OpenMP_CXX)
	target_compile_definitions(adfinder PRIVATE ADFINDER_HAS_OPENMP=1)
endif()
