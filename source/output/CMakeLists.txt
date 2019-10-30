add_library(ariadne-output OBJECT
    logging.cpp
    graphics.cpp
    drawer.cpp
)

if(NOT GTK2_FOUND)
    cotire(ariadne-output CONFIGURATIONS None)
else()
    message(STATUS "cotire disabled for the output module due to conflicts with the GTK2 library.")
endif()

if(COVERAGE)
    target_link_libraries(ariadne-output PUBLIC coverage_config)
endif()
