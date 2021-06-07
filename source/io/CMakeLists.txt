set(LIBRARY_NAME ariadne-io)

add_library(${LIBRARY_NAME} OBJECT
    logging.cpp
    figure.cpp
    cairo.cpp
    gnuplot.cpp
    graphics_manager.cpp
    drawer.cpp
    command_line_interface.cpp
    configuration_file.cpp
)

foreach(WARN ${LIBRARY_EXCLUSIVE_WARN})
    target_compile_options(${LIBRARY_NAME} PUBLIC -W${WARN})
endforeach()

if(COVERAGE)
    target_link_libraries(${LIBRARY_NAME} PUBLIC coverage_config)
endif()

