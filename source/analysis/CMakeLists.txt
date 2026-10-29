# Staging target: reachability code still compiles in the legacy directories.
add_library(ariadne-analysis INTERFACE)
target_link_libraries(ariadne-analysis INTERFACE ariadne-evolution)
