cmake_minimum_required(VERSION 3.2)
project(discreteGeometry CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(DGtal REQUIRED)
link_directories(${DGTAL_LIBRARY_DIRS})


# require Eigen
#find_package(Eigen3 REQUIRED)
#if(${EIGEN3_FOUND})
#    message(STATUS "lib EIGEN3 found")
#    message(STATUS "  version " ${EIGEN3_VERSION_STRING})
#    message(STATUS "  include " ${EIGEN3_INCLUDE_DIR})
#endif()

set(FILE_NAME TD1)

set(LIBRAIRIES ${DGTAL_LIBRARIES})
set(INCLUDE_DIRECTORIES ${DGTAL_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR}/lib/Eigen3)

add_subdirectory(src)