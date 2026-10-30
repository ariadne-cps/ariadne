# System specifications and dynamics, layered on hybrid set representations.
add_library(ariadne-system OBJECT
    system_function.cpp
    vector_field.cpp
    iterated_map.cpp
    differential_inclusion.cpp
    hybrid_automaton.cpp
    hybrid_automaton-composite.cpp
)
set_target_properties(ariadne-system PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_link_libraries(ariadne-system PUBLIC ariadne-hybrid)

# Every system header must compile as the first and only include of its own
# translation unit. Do not rely on an aggregate Ariadne header's include order.
file(GLOB _ariadne_system_headers CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")
set(_ariadne_system_header_tus)
foreach(_header IN LISTS _ariadne_system_headers)
    get_filename_component(_name "${_header}" NAME)
    set(_tu "${CMAKE_CURRENT_BINARY_DIR}/header-isolation/${_name}.cpp")
    file(GENERATE OUTPUT "${_tu}" CONTENT "#include \"system/${_name}\"\n")
    list(APPEND _ariadne_system_header_tus "${_tu}")
endforeach()
add_library(ariadne-system-header-isolation OBJECT
    ${_ariadne_system_header_tus})
target_link_libraries(ariadne-system-header-isolation PRIVATE ariadne-hybrid)
# Building `ariadne` (including CI Release builds) also compiles these checks.
add_dependencies(ariadne-system ariadne-system-header-isolation)

# An upward include would make the system definitions depend on their
# evolution or analysis. Refuse such dependencies at configure time, even if
# an aggregate build happens to make them compile transitively.
file(GLOB _ariadne_system_sources CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
foreach(_file IN LISTS _ariadne_system_sources)
    file(STRINGS "${_file}" _forbidden_includes
        REGEX "^[ \t]*#[ \t]*include[ \t]*[\"<](evolution|analysis|dynamics)/")
    if(_forbidden_includes)
        message(FATAL_ERROR
            "System source ${_file} includes an upper-layer header:\n"
            "${_forbidden_includes}")
    endif()
endforeach()
