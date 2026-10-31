# Prevent accidental dependencies on analysis from evolution sources.
if(NOT DEFINED EVOLUTION_SOURCE_DIR)
    message(FATAL_ERROR "EVOLUTION_SOURCE_DIR is required")
endif()

file(GLOB evolution_sources
    "${EVOLUTION_SOURCE_DIR}/*.hpp"
    "${EVOLUTION_SOURCE_DIR}/*.cpp")

foreach(source_file IN LISTS evolution_sources)
    file(STRINGS "${source_file}" forbidden_includes
        REGEX "^[ \t]*#[ \t]*include[ \t]*[\"<](analysis/|dynamics/)")
    if(forbidden_includes)
        message(FATAL_ERROR "Evolution source includes analysis via ${source_file}: ${forbidden_includes}")
    endif()
endforeach()

file(READ "${EVOLUTION_SOURCE_DIR}/CMakeLists.txt" evolution_cmake)
if(evolution_cmake MATCHES "ariadne-analysis")
    message(FATAL_ERROR "Evolution CMake target must not depend on ariadne-analysis")
endif()
