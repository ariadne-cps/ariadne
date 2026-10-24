# Exclusions for the Ariadne documentation-only Doxygen target.
# All paths are relative to the Ariadne repository root. Add C++ source
# exclusions here, not to Configuration or the generated Doxyfile.
set(ARIADNE_DOXYGEN_EXCLUSIONS
  doc/site
  doc/shared-configuration

  # Legacy implementation exclusions (absent files are harmless).
  source/prototyping
  source/experimental
  source/algebra/vector-crtp.hpp
  source/algebra/vector-sfinae.hpp
  source/algebra/matrix-crtp.hpp
  source/algebra/matrix-sfinae.hpp

  # Pages embedding C++/Python example code (snippet/dontinclude).
  doc/hybrid_evolution_tutorial.dox
  doc/rigorous_numerics_tutorial.dox
  doc/python_module_tutorials.dox
  doc/python_module_demonstrations.dox
  doc/python_example_tutorials.dox
  doc/neural_network.dox
)
