# Staging interface: system definitions will be migrated from dynamics/hybrid.
# System must not link to labelling, evolution, or analysis.
add_library(ariadne-system INTERFACE)
target_link_libraries(ariadne-system INTERFACE symbolic)
