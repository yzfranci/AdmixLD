cmake_minimum_required(VERSION 3.16)
project(admixld LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()

# HTSlib
find_path(HTS_INCLUDE_DIR NAMES htslib/sam.h)
find_library(HTS_LIBRARY NAMES hts)
if(NOT HTS_INCLUDE_DIR OR NOT HTS_LIBRARY)
	message(FATAL_ERROR "htslib not found (need htslib/sam.h and libhts)")
endif()

# Eigen
find_package(Eigen3 REQUIRED)

add_executable(admixld
	src/main.cpp
	src/io/vcf_windows.cpp
	src/io/msp_windows.cpp
	src/io/bed.cpp
	src/io/sample_vector.cpp
	src/io/ref_freq.cpp
	src/core/hybrid_index.cpp
	src/core/residualize.cpp
	src/core/scan_markers.cpp
	src/core/scan_dynamic_hi.cpp
	src/core/empirical_null.cpp
)

# HTSlib
target_include_directories(admixld PRIVATE ${HTS_INCLUDE_DIR})
target_link_libraries(admixld PRIVATE ${HTS_LIBRARY})

# Eigen
target_link_libraries(admixld PRIVATE Eigen3::Eigen)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
	message(STATUS "OpenMP found")
	target_link_libraries(admixld PRIVATE OpenMP::OpenMP_CXX)
	target_compile_definitions(admixld PRIVATE ADMIXLD_HAS_OPENMP=1)
endif()

include(GNUInstallDirs)
install(TARGETS admixld RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
