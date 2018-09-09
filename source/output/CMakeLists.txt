add_library(ariadne-output OBJECT
    logging.cpp
    graphics.cpp
    drawer.cpp
)

cotire(ariadne-output CONFIGURATIONS None)
