set(LIBRARY_NAME ariadne-foundations)

add_library(${LIBRARY_NAME} OBJECT
    logical.cpp
)

target_link_libraries(${LIBRARY_NAME} PUBLIC utility)
