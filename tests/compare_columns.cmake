
if(EXISTS METAANALYSIS1.TBL)
    file(REMOVE METAANALYSIS1.TBL)
endif()
if(EXISTS METAANALYSIS1.TBL.info)
    file(REMOVE METAANALYSIS1.TBL.info)
endif()

execute_process(COMMAND ${METAL} metal.txt RESULT_VARIABLE metal_exit_code)
if(metal_exit_code)
    message(FATAL_ERROR "METAL failed.")
endif()

execute_process(COMMAND ${METALTEST} test_file_correlation METAANALYSIS1.TBL ${METAL_TBL} ${COLUMN} ${CORR} ${TOL} RESULT_VARIABLE diff_exit_code)
if(diff_exit_code)
    message(FATAL_ERROR "METAL didn't replicate results.")
endif()