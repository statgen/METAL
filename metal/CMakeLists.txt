include_directories(${CMAKE_SOURCE_DIR}/libsrc/)

add_definitions(-D__ZLIB_AVAILABLE__)

add_executable(metal Main.cpp)
target_link_libraries(metal goncalo ${ZLIB_LIBRARIES})


install(TARGETS metal DESTINATION ${CMAKE_BINARY_DIR}/bin)