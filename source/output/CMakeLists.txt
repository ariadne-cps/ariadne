add_library(ariadne-output OBJECT
    logging.cpp
    graphics.cpp
    drawer.cpp
)

# GTK2 headers on gcc mess with cotire
if(NOT((DEFINED GTK2_FOUND) AND (CMAKE_CXX_COMPILER_ID MATCHES "GNU")))
    cotire(ariadne-output CONFIGURATIONS None)
endif()
