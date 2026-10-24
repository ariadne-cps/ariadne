cmake_minimum_required(VERSION 3.15)
# Documentation-only project: intentionally independent of Ariadne's
# older build-time Configuration pin and its CMake compatibility contracts.
project(Ariadne VERSION 2.5.3 LANGUAGES NONE)

get_filename_component(ARIADNE_ROOT "${CMAKE_CURRENT_LIST_DIR}/../.." ABSOLUTE)
set(_doc_support "${ARIADNE_ROOT}/doc/shared-configuration/cmake/UseDoxygen.cmake")
if(NOT EXISTS "${_doc_support}")
  message(FATAL_ERROR "Initialize documentation submodule: git submodule update --init doc/shared-configuration")
endif()
include("${_doc_support}")
ariadne_setup_documentation(
  INPUT
    "${ARIADNE_ROOT}/source"
    "${ARIADNE_ROOT}/python"
    "${ARIADNE_ROOT}/tutorials"
    "${ARIADNE_ROOT}/doc"
  IMAGE_PATH "${ARIADNE_ROOT}/doc"
)

# Retain Ariadne's project-specific Doxygen semantics while Configuration
# continues to own the base Doxyfile and XML-to-Markdown conversion.
# This extends the generated config instead of copying Doxyfile.in.
if(TARGET doxygen)
  set(_doxyfile "${CMAKE_CURRENT_BINARY_DIR}/Doxyfile")
  file(APPEND "${_doxyfile}" "
PROJECT_NAME = Ariadne
PROJECT_BRIEF = \"Ariadne library and mathematical documentation\"
STRIP_FROM_PATH = \"${ARIADNE_ROOT}\"
LAYOUT_FILE = \"${ARIADNE_ROOT}/doc/DoxygenLayout.xml\"
MATHJAX_CODEFILE = \"${ARIADNE_ROOT}/doc/macros.js\"
EXAMPLE_PATH = \"${ARIADNE_ROOT}/tutorials\" \"${ARIADNE_ROOT}/python\"
EXCLUDE = \"${ARIADNE_ROOT}/source/prototyping\" \"${ARIADNE_ROOT}/source/experimental\" \"${ARIADNE_ROOT}/python/pybind11\" \"${ARIADNE_ROOT}/source/algebra/vector-crtp.hpp\" \"${ARIADNE_ROOT}/source/algebra/vector-sfinae.hpp\" \"${ARIADNE_ROOT}/source/algebra/matrix-crtp.hpp\" \"${ARIADNE_ROOT}/source/algebra/matrix-sfinae.hpp\"
EXTRA_PACKAGES = amsmath amsfonts amssymb \"${ARIADNE_ROOT}/doc/macros\"
")
endif()
