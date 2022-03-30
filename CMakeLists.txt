cmake_minimum_required(VERSION 2.6)
project ("hwt_pj" CXX)
set(CMAKE_CXX_STANDARD 17)

#find_library(HDF5_LIB hdf5)
set(HDF5_LIB "/usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5.so")

set (CMAKE_CXX_FLAGS "-DSINGLE_CORE -O2 -fopenmp -DNDEBUG -std=c++0x")
find_package(OpenMP REQUIRED)

if(NOT WIN32)
	find_package(PkgConfig QUIET)
endif()

if(PKG_CONFIG_FOUND)
	pkg_check_modules(LIBPMEMOBJ++ REQUIRED libpmemobj++)
else()
	find_package(LIBPMEMOBJ++ REQUIRED)
endif()

include_directories ("${PROJECT_SOURCE_DIR}/include")
include_directories("/usr/include/hdf5/serial/")

link_directories(${LIBPMEMOBJ++_LIBRARY_DIRS})

add_executable (hwt_pj src/node.cpp src/root_node.cpp src/Preprocessing.cpp interface/hwt_interface.cpp)
target_link_libraries (hwt_pj ${HDF5_LIB} "${LIBPMEMKV_LIBRARIES}" ${LIBPMEMOBJ++_LIBRARIES} OpenMP::OpenMP_CXX)

project (linscan)
add_executable (linscan interface/linscan_interface.cpp src/linscan.cpp src/loadVar.cpp)
target_link_libraries (linscan ${HDF5_LIB} "${LIBPMEMKV_LIBRARIES}" ${LIBPMEMOBJ++_LIBRARIES} OpenMP::OpenMP_CXX)


