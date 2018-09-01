# Try to find the Ariadne library for rigorous numerics.
#
# Once done this will define:
#
#  ARIADNE_FOUND - system has Ariadne
#  ARIADNE_INCLUDE_DIRS - the Ariadne include directories
#  ARIADNE_LIBRARIES - Link these to use Ariadne
#
# Copyright (C) 2018, Luca Geretti <luca.geretti@univr.it>
#
# Redistribution and use is allowed according to the terms of the GPLv2 license.

if(ARIADNE_INCLUDE_DIRS AND ARIADNE_LIBRARIES)
  message("Ariadne includes and libraries already identified.")

  set(ARIADNE_FOUND TRUE)

else()

    find_library(ARIADNE_LIBRARY ariadne)
    set(ARIADNE_LIBRARIES "${ARIADNE_LIBRARY}")

    find_path(ARIADNE_INCLUDE_DIR ariadne.hpp PATH_SUFFIXES ariadne)
    set(ARIADNE_INCLUDE_DIRS "${ARIADNE_INCLUDE_DIR}")

    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(Ariadne DEFAULT_MSG ARIADNE_LIBRARIES ARIADNE_INCLUDE_DIRS)

endif()

mark_as_advanced(
  ARIADNE_INCLUDE_DIRS
  ARIADNE_LIBRARIES
)
