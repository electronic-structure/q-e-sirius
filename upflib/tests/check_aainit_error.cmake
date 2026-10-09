# Require the intended QE error, not just a failed launcher or process.
execute_process(COMMAND ${TEST_COMMAND}
    RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(result EQUAL 0 OR NOT "${output}${error}" MATCHES "${EXPECTED_ERROR}")
    message(FATAL_ERROR "Expected aainit error '${EXPECTED_ERROR}', got ${result}:\n${output}${error}")
endif()
