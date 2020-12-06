cmake_minimum_required(VERSION 3.2)
project(discreteGeometry CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(DGtal REQUIRED)
link_directories(${DGTAL_LIBRARY_DIRS})

set(FILE_NAME TD1)

set(LIBRAIRIES ${DGTAL_LIBRARIES})
set(INCLUDE_DIRECTORIES ${DGTAL_INCLUDE_DIRS})

add_subdirectory(src)