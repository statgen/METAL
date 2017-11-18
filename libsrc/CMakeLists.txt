file(GLOB SRCS *.cpp)

include_directories(${ZLIB_INCLUDE_DIRS})
add_definitions(-D__ZLIB_AVAILABLE__)

add_library(goncalo STATIC ${SRCS})
