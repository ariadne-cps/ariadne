# The external hybrid submodule is initialized before local modules by ss/.
# Internal dependency order: system -> evolution -> analysis.
add_subdirectory(system)
add_subdirectory(evolution)
add_subdirectory(analysis)
