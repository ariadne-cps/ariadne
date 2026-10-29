# Staging interface: numerical evolution is still housed in dynamics/hybrid.
add_library(ariadne-evolution INTERFACE)
target_link_libraries(ariadne-evolution INTERFACE
    ariadne-system
    integration
    labelling
    equation_solving
)
