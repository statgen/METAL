cmake_minimum_required(VERSION 3.1)

set(CMAKE_C_COMPILER "gcc")
set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS "-O2 -D_FILE_OFFSET_BITS=64")


project(HeavtMetal)

add_subdirectory(libsrc)
add_subdirectory(metal)


