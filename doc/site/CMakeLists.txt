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
    "${ARIADNE_ROOT}/doc"
  IMAGE_PATH "${ARIADNE_ROOT}/doc"
)

# Retain Ariadne's project-specific Doxygen semantics while Configuration
# continues to own the base Doxyfile and XML-to-Markdown conversion.
# This extends the generated config instead of copying Doxyfile.in.
if(TARGET doxygen)
  set(_doxyfile "${CMAKE_CURRENT_BINARY_DIR}/Doxyfile")
  # Centralize all Ariadne-only exclusions under doc/site.
  include("${ARIADNE_ROOT}/doc/site/exclusions.cmake")
  set(_doxygen_excludes "")
  foreach(_excluded IN LISTS ARIADNE_DOXYGEN_EXCLUSIONS)
    if(IS_ABSOLUTE "${_excluded}")
      message(FATAL_ERROR "Doxygen exclusion must be relative: ${_excluded}")
    endif()
    if(EXISTS "${ARIADNE_ROOT}/${_excluded}")
      string(APPEND _doxygen_excludes " \"${ARIADNE_ROOT}/${_excluded}\"")
    endif()
  endforeach()
  file(APPEND "${_doxyfile}" "
PROJECT_NAME = Ariadne
PROJECT_BRIEF = \"Ariadne library and mathematical documentation\"
STRIP_FROM_PATH = \"${ARIADNE_ROOT}\"
LAYOUT_FILE =
MATHJAX_CODEFILE = \"${ARIADNE_ROOT}/doc/macros.js\"
EXAMPLE_PATH =
EXCLUDE =${_doxygen_excludes}
EXTRA_PACKAGES = amsmath amsfonts amssymb \"${ARIADNE_ROOT}/doc/macros\"
")
endif()
