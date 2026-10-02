set(LIBRARY_NAME ariadne-io)

add_library(${LIBRARY_NAME} OBJECT
    figure.cpp
    cairo.cpp
    graphics_manager.cpp
    drawer.cpp
    command_line_interface.cpp
)



