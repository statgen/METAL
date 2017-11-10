include_directories(${CMAKE_SOURCE_DIR}/libsrc/)

find_package(ZLIB REQUIRED)

message(STATUS "ZLIB version: " ${ZLIB_VERSION_STRING})
message(STATUS "ZLIB header directory: " ${ZLIB_INCLUDE_DIRS})
message(STATUS "ZLIB library: " ${ZLIB_LIBRARIES})

include_directories(${ZLIB_INCLUDE_DIRS})

add_executable(metal Main.cpp)
target_link_libraries(metal goncalo ${ZLIB_LIBRARIES})
