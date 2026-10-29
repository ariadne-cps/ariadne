# Staging target: the evolvers still compile in the legacy source directories.
add_library(ariadne-evolution INTERFACE)
target_link_libraries(ariadne-evolution INTERFACE ariadne-system)
