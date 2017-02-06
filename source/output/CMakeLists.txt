add_library(ariadne-output OBJECT
    logging.cpp
    graphics.cpp
    drawer.cpp
    textplot.cpp
)

cotire(ariadne-output)
