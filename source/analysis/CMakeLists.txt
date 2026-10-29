# Staging interface: reachability remains in dynamics/hybrid.
add_library(ariadne-analysis INTERFACE)
target_link_libraries(ariadne-analysis INTERFACE
    ariadne-evolution
    smt
    sat
)
