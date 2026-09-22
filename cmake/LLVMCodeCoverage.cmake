include(CMakeParseArguments)

find_program(XCRUN_EXECUTABLE xcrun)
if(NOT XCRUN_EXECUTABLE)
    message(FATAL_ERROR "xcrun not found; LLVM coverage on macOS requires the Xcode command line tools.")
endif()

execute_process(
    COMMAND "${XCRUN_EXECUTABLE}" --find llvm-profdata
    RESULT_VARIABLE LLVM_PROFDATA_RESULT
    OUTPUT_VARIABLE LLVM_PROFDATA_EXECUTABLE
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
)
if(NOT LLVM_PROFDATA_RESULT EQUAL 0 OR NOT LLVM_PROFDATA_EXECUTABLE)
    message(FATAL_ERROR "llvm-profdata not found via xcrun.")
endif()

execute_process(
    COMMAND "${XCRUN_EXECUTABLE}" --find llvm-cov
    RESULT_VARIABLE LLVM_COV_RESULT
    OUTPUT_VARIABLE LLVM_COV_EXECUTABLE
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
)
if(NOT LLVM_COV_RESULT EQUAL 0 OR NOT LLVM_COV_EXECUTABLE)
    message(FATAL_ERROR "llvm-cov not found via xcrun.")
endif()

set(LLVM_COVERAGE_COMPILE_FLAGS "-fprofile-instr-generate -fcoverage-mapping")
set(LLVM_COVERAGE_LINK_FLAGS "-fprofile-instr-generate")

function(append_llvm_coverage_compiler_flags)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${LLVM_COVERAGE_COMPILE_FLAGS}" PARENT_SCOPE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${LLVM_COVERAGE_COMPILE_FLAGS}" PARENT_SCOPE)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LLVM_COVERAGE_LINK_FLAGS}" PARENT_SCOPE)
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${LLVM_COVERAGE_LINK_FLAGS}" PARENT_SCOPE)
    set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} ${LLVM_COVERAGE_LINK_FLAGS}" PARENT_SCOPE)
    message(STATUS "Appending LLVM code coverage compiler flags: ${LLVM_COVERAGE_COMPILE_FLAGS}")
endfunction()

function(setup_target_for_coverage_llvm)
    set(options NONE)
    set(oneValueArgs NAME TARGET EXCLUDE_REGEX)
    set(multiValueArgs DEPENDENCIES)
    cmake_parse_arguments(Coverage "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    if(NOT Coverage_NAME)
        message(FATAL_ERROR "setup_target_for_coverage_llvm requires NAME.")
    endif()
    if(NOT Coverage_TARGET)
        message(FATAL_ERROR "setup_target_for_coverage_llvm requires TARGET.")
    endif()
    if(NOT TARGET ${Coverage_TARGET})
        message(FATAL_ERROR "Coverage target '${Coverage_TARGET}' does not exist.")
    endif()

    set(PROFILE_DIR "${PROJECT_BINARY_DIR}/coverage/profiles")
    set(PROFDATA_FILE "${PROJECT_BINARY_DIR}/coverage/coverage.profdata")
    set(HTML_DIR "${PROJECT_BINARY_DIR}/coverage/html")
    set(MERGE_SCRIPT "${PROJECT_BINARY_DIR}/merge-llvm-coverage.cmake")

    file(WRITE "${MERGE_SCRIPT}"
"file(GLOB LLVM_RAW_PROFILES \"${PROFILE_DIR}/*.profraw\")
if(NOT LLVM_RAW_PROFILES)
    message(FATAL_ERROR \"No LLVM raw coverage profiles were generated.\")
endif()
execute_process(
    COMMAND \"${LLVM_PROFDATA_EXECUTABLE}\" merge -sparse \${LLVM_RAW_PROFILES} -o \"${PROFDATA_FILE}\"
    RESULT_VARIABLE LLVM_PROFILE_MERGE_RESULT
)
if(NOT LLVM_PROFILE_MERGE_RESULT EQUAL 0)
    message(FATAL_ERROR \"llvm-profdata merge failed.\")
endif()
")

    set(LLVM_COV_FILTER_ARGS "")
    if(Coverage_EXCLUDE_REGEX)
        list(APPEND LLVM_COV_FILTER_ARGS "-ignore-filename-regex=${Coverage_EXCLUDE_REGEX}")
    endif()

    add_custom_target(${Coverage_NAME}
        COMMAND "${CMAKE_COMMAND}" -E rm -rf "${PROJECT_BINARY_DIR}/coverage"
        COMMAND "${CMAKE_COMMAND}" -E make_directory "${PROFILE_DIR}"
        COMMAND "${CMAKE_COMMAND}" -E env
                "LLVM_PROFILE_FILE=${PROFILE_DIR}/%p-%m.profraw"
                "${CMAKE_CTEST_COMMAND}" --output-on-failure
        COMMAND "${CMAKE_COMMAND}" -P "${MERGE_SCRIPT}"
        COMMAND "${LLVM_COV_EXECUTABLE}" report
                "$<TARGET_FILE:${Coverage_TARGET}>"
                "-instr-profile=${PROFDATA_FILE}"
                ${LLVM_COV_FILTER_ARGS}
        COMMAND "${CMAKE_COMMAND}" -E make_directory "${HTML_DIR}"
        COMMAND "${LLVM_COV_EXECUTABLE}" show
                "$<TARGET_FILE:${Coverage_TARGET}>"
                "-instr-profile=${PROFDATA_FILE}"
                "-format=html"
                "-output-dir=${HTML_DIR}"
                ${LLVM_COV_FILTER_ARGS}
        WORKING_DIRECTORY "${PROJECT_BINARY_DIR}"
        DEPENDS ${Coverage_DEPENDENCIES} ${Coverage_TARGET}
        VERBATIM
        COMMENT "Running tests and generating LLVM code coverage report."
    )

    add_custom_command(TARGET ${Coverage_NAME} POST_BUILD
        COMMAND "${CMAKE_COMMAND}" -E echo
                "LLVM coverage HTML report: ${HTML_DIR}/index.html"
    )
endfunction()
