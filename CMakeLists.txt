cmake_minimum_required(VERSION 3.16)
project(admixld LANGUAGES CXX)

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

add_executable(admixld
	src/main.cpp
	src/io/vcf_windows.cpp
	src/io/msp_windows.cpp
	src/io/bed.cpp
	src/io/sample_vector.cpp
	src/core/hybrid_index.cpp
	src/core/residualize.cpp
	src/core/scan_blocks.cpp
	src/core/scan_dynamic_hi.cpp
)

# HTSlib
target_include_directories(admixld PRIVATE ${HTS_INCLUDE_DIRS})
target_link_directories(admixld PRIVATE ${HTS_LIBRARY_DIRS})
target_link_libraries(admixld PRIVATE ${HTS_LIBRARIES})
target_compile_options(admixld PRIVATE ${HTS_CFLAGS_OTHER})

# Eigen
target_link_libraries(admixld PRIVATE Eigen3::Eigen)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
	message(STATUS "OpenMP found")
	target_link_libraries(admixld PRIVATE OpenMP::OpenMP_CXX)
	target_compile_definitions(admixld PRIVATE ADMIXLD_HAS_OPENMP=1)
endif()
