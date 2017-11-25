cmake_minimum_required(VERSION 3.1)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS "-D_FILE_OFFSET_BITS=64")

project(Metal)

find_package(ZLIB REQUIRED)
message(STATUS "ZLIB version: " ${ZLIB_VERSION_STRING})
message(STATUS "ZLIB header directory: " ${ZLIB_INCLUDE_DIRS})
message(STATUS "ZLIB library: " ${ZLIB_LIBRARIES})

enable_testing()

add_subdirectory(libsrc)
add_subdirectory(metal)
add_subdirectory(tests)
