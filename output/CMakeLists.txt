add_library(ariadne-output OBJECT
    logging.cc
    graphics.cc
    drawer.cc
    textplot.cc
)

cotire(ariadne-output)
