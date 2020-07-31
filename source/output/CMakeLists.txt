set(LIBRARY_NAME ariadne-output)

add_library(${LIBRARY_NAME} OBJECT
    logging.cpp
    graphics.cpp
    drawer.cpp
)

if(NOT GTK2_FOUND)
    cotire(${LIBRARY_NAME} CONFIGURATIONS None)
else()
    message(STATUS "cotire disabled for the output module due to conflicts with the GTK2 library.")
endif()

foreach(WARN ${LIBRARY_EXCLUSIVE_WARN})
    target_compile_options(${LIBRARY_NAME} PUBLIC -W${WARN})
endforeach()

if(COVERAGE)
    target_link_libraries(${LIBRARY_NAME} PUBLIC coverage_config)
endif()

