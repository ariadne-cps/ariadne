include(CMakeParseArguments)

find_program(LCOV_EXECUTABLE NAMES lcov REQUIRED)
find_program(GENHTML_EXECUTABLE NAMES genhtml REQUIRED)

function(setup_target_for_coverage_gcc)
    set(options NONE)
    set(oneValueArgs NAME)
    set(multiValueArgs DEPENDENCIES)
    cmake_parse_arguments(Coverage "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    if(NOT Coverage_NAME)
        message(FATAL_ERROR "setup_target_for_coverage_gcc requires NAME.")
    endif()

    set(COVERAGE_INFO "${PROJECT_BINARY_DIR}/${Coverage_NAME}.info")
    set(COVERAGE_HTML_DIR "${PROJECT_BINARY_DIR}/${Coverage_NAME}")

    add_custom_target(${Coverage_NAME}
        COMMAND "${CMAKE_CTEST_COMMAND}" --output-on-failure
        COMMAND "${LCOV_EXECUTABLE}"
                --directory "${PROJECT_BINARY_DIR}"
                --capture
                --output-file "${COVERAGE_INFO}"
        COMMAND "${LCOV_EXECUTABLE}"
                --remove "${COVERAGE_INFO}" "/usr/*"
                --output-file "${COVERAGE_INFO}"
        COMMAND "${LCOV_EXECUTABLE}" --list "${COVERAGE_INFO}"
        COMMAND "${GENHTML_EXECUTABLE}"
                --output-directory "${COVERAGE_HTML_DIR}"
                "${COVERAGE_INFO}"
        WORKING_DIRECTORY "${PROJECT_BINARY_DIR}"
        DEPENDS ${Coverage_DEPENDENCIES}
        VERBATIM
        COMMENT "Running tests and generating GCC/lcov code coverage report."
    )

    add_custom_command(TARGET ${Coverage_NAME} POST_BUILD
        COMMAND "${CMAKE_COMMAND}" -E echo
                "Lcov coverage report: ${COVERAGE_HTML_DIR}/index.html"
    )
endfunction()
