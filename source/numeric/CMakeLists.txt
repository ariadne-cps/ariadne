set(LIBRARY_NAME NUMERIC_SRC)

add_library(${LIBRARY_NAME} OBJECT
    operators.cpp
    int.cpp
    double.cpp
    integer.cpp
    dyadic.cpp
    decimal.cpp
    rational.cpp
    real.cpp
    complex.cpp
    number.cpp
    floatdp.cpp
    floatmp.cpp
    float_approximation.cpp
    float_lower_bound.cpp
    float_upper_bound.cpp
    float_bounds.cpp
    float_ball.cpp
    float_error.cpp
    float_literals.cpp

    concepts.cpp
    logical_sequence.cpp
)

if(NOT TARGET numeric)
    add_library(numeric INTERFACE)
    target_include_directories(numeric INTERFACE
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/source>
        $<INSTALL_INTERFACE:include/ariadne-numeric>
    )
    target_link_libraries(numeric INTERFACE
        foundation
        utility
        GMP::GMP
        MPFR::MPFR
    )
endif()

target_link_libraries(${LIBRARY_NAME} PRIVATE numeric)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU" OR NOT APPLE)
    target_compile_options(${LIBRARY_NAME} PRIVATE -frounding-math)
endif()

if(NOT TARGET ariadne-numeric)
    setup_project_library(
        TARGET ariadne-numeric
        OBJECTS
            NUMERIC_SRC
            FOUNDATION_SRC
            UTILITY_SRC
        LINK_LIBRARIES
            numeric
    )
endif()
