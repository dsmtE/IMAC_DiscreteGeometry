cmake_minimum_required(VERSION 3.12)
project(discreteGeometry C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Setup warnings
if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC") 
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Wall")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall")
endif()

#set output build directory
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/out/${CMAKE_BUILD_TYPE})

# ////////////////////////////// PROJECT FILES ////////////////////////////////
set(MAIN_FILE ${PROJECT_SOURCE_DIR}/src/main.cpp)
file(GLOB_RECURSE SRC_FILES  ${PROJECT_SOURCE_DIR}/src/*.cpp )
# remove main from sources
LIST(REMOVE_ITEM SRC_FILES ${MAIN_FILE})
file(GLOB_RECURSE HEADERS_FILES "src/*.hpp")

# ////////////////////////////// LIBRARIES //////////////////////////////// #
set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR}/install)

# set(Boost_USE_STATIC_LIBS ON)
find_package(Boost REQUIRED)
find_package(zlib REQUIRED)

# ////////////////////////////// EXTERNAL LIBRARIES //////////////////////////////// #
# Use CMake ExternalProject to build dependencies
include (ExternalProject)

set(BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}/external")

set(DGTAL_CONFIG_FLAGS -DBUILD_EXAMPLES=OFF -DBUILD_SHARED_LIBS=OFF -DWITH_GMP=OFF -DWITH_ITK=OFF -DWITH_CAIRO=OFF -DWITH_QGLVIEWER=OFF -DWITH_MAGICK=OFF)
set(DGTAL_TARGET dgtal)
ExternalProject_Add(${DGTAL_TARGET}
    # Download Archive from Web Server
    URL https://github.com/DGtal-team/DGtal/archive/1.1.tar.gz
    URL_HASH MD5=c1f9189ffd42874ccfc613aafb37c527
    # DOWNLOAD_NO_PROGRESS ON
    # Download Source from GIT Repository
    # GIT_REPOSITORY https://github.com/DGtal-team/DGtal
    # GIT_TAG        v1.1   
    # GIT_SHALLOW     ON
    DOWNLOAD_DIR ${BUILD_DIR}/download/dgtal

    PREFIX ${BUILD_DIR}
    BUILD_IN_SOURCE 0
    BUILD_ALWAYS 0
    UPDATE_COMMAND ""
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/dgtal
    BINARY_DIR ${BUILD_DIR}/dgtal_build
    INSTALL_DIR ${CMAKE_INSTALL_PREFIX}

    CONFIGURE_COMMAND ${CMAKE_COMMAND} -DCMAKE_TOOLCHAIN_FILE=C:/vcpkg/scripts/buildsystems/vcpkg.cmake -DVCPKG_TARGET_TRIPLET=x64-windows -A x64 -T host=x64 ${DGTAL_CONFIG_FLAGS} -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR> <SOURCE_DIR>
    BUILD_COMMAND $(MAKE)
    )
set(DGTAL_LIBRARIES ${CMAKE_INSTALL_PREFIX}/lib/DGtald.lib)

set(INCLUDES_DIRS ${Boost_INCLUDE_DIRS} ${ZLIB_INCLUDE_DIRS} ${CMAKE_INSTALL_PREFIX}/include)
set(LINK_DIRS ${Boost_LIBRARY_DIR})
set(LIBRARIES ${Boost_LIBRARIES} ${ZLIB_LIBRARIES} ${DGTAL_LIBRARIES})

set(FILE_NAME discreteGeometry)

add_executable(${FILE_NAME} ${MAIN_FILE} ${HEADERS_FILES} ${SRC_FILES})

add_dependencies(${FILE_NAME} ${DGTAL_TARGET})

target_link_libraries(${FILE_NAME} ${LIBRARIES})
target_include_directories(${FILE_NAME} PUBLIC ${INCLUDES_DIRS})
