add_library(ariadne-output OBJECT
    logging.cpp
    graphics.cpp
    drawer.cpp
)

if(NOT((DEFINED GTK2_FOUND)))
    cotire(ariadne-output CONFIGURATIONS None)
else()
    message(STATUS "cotire disabled for the output module due to conflicts with the GTK2 library.")
endif()
