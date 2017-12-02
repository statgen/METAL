include_directories(${CMAKE_SOURCE_DIR}/libsrc/)
add_definitions(-D__ZLIB_AVAILABLE__)

add_executable(metaltest metaltest.cpp)
target_link_libraries(metaltest goncalo ${ZLIB_LIBRARIES})

add_test(NAME NORMP COMMAND metaltest test_normp)
add_test(NAME BINORMP COMMAND metaltest test_binormp)
add_test(NAME BINORMQ COMMAND metaltest test_binormq)

file(COPY test_01 DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
add_test(NAME SAMPLESIZE_NOGC WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test_01
        COMMAND ${CMAKE_COMMAND}
        -DMETAL=${CMAKE_BINARY_DIR}/metal/metal
        -DMETAL_TBL=${CMAKE_CURRENT_SOURCE_DIR}/test_01/METAANALYSIS1.TBL
        -P ${CMAKE_CURRENT_SOURCE_DIR}/compare.cmake)

file(COPY test_02 DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
add_test(NAME SAMPLESIZE_NOGC_OVERLAP WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test_01
        COMMAND ${CMAKE_COMMAND}
        -DMETAL=${CMAKE_BINARY_DIR}/metal/metal
        -DMETAL_TBL=${CMAKE_CURRENT_SOURCE_DIR}/test_02/METAANALYSIS1.TBL
        -P ${CMAKE_CURRENT_SOURCE_DIR}/compare.cmake)
