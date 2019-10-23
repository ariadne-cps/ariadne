# ##############################################################################
# Copyright (C) 2014 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH    # #
# This software is distributed under the terms of the             # GNU Lesser
# General Public Licence (LGPL) version 3,             # copied verbatim in the
# file "LICENSE"                       #
# ##############################################################################
# * Try to find GLPK instalation Once done this will define
#

message(STATUS "Looking for GLPK ...")

find_path(
  GLPK_INCLUDE_DIR
  NAMES glpk.h
  PATHS
  HINTS ${GLPK_INCLUDE_DIR} /usr/local/include /usr/local/include/glpk/include
  NO_DEFAULT_PATH)

find_path(
  GLPK_LIBRARIES_DIR
  NAMES libglpk.so
  PATHS ${SIMPATH}/basics/glpk/lib ${SIMPATH}/local/lib ${SIMPATH}/local/lib
        /usr/local/lib /usr/lib/x86_64-linux-gnu/ /usr/local/include/glpk/lib
  NO_DEFAULT_PATH)

find_library(
  GLPK_LIBRARIES
  NAMES libglpk.so
  HINTS ${GLPK_LIBRARIES_DIR} /usr/lib/libglpk.so)

if(GLPK_INCLUDE_DIR AND GLPK_LIBRARIES_DIR)
  set(GLPK_FOUND TRUE)
endif(GLPK_INCLUDE_DIR AND GLPK_LIBRARIES_DIR)

if(GLPK_FOUND)
  if(NOT GLPK_FOUND_QUIETLY)
    message(STATUS "Looking for GLPK... - found ${GLPK_LIBRARIES_DIR}")
    # message(STATUS "Found PLUTO: ${PLUTO_LIBRARY_DIR}")
    set(LD_LIBRARY_PATH ${LD_LIBRARY_PATH} ${GLPK_LIBRARIES_DIR})
  endif(NOT GLPK_FOUND_QUIETLY)
else(GLPK_FOUND)
  if(GLPK_FOUND_REQUIRED)
    message(FATAL_ERROR "Looking for GLPK... - Not found")
  endif(GLPK_FOUND_REQUIRED)
endif(GLPK_FOUND)
