VERSION = 0.0.0
NAME = ariadne

# set the following variables 
DEBUG = no
USE_OCTAVE = no
FASTER_BUT_DIRTIER = yes

CC = gcc
CXX = g++
AR = ar

MAKE = make

SRCDIR = src
INCLUDEDIR = include

OCTAVE_HEADERS = /usr/include/octave-2.1.57

OBJECTS = approx_type.o arc.o automaton.o basic_set_list.o \
	linear_algebra.o location.o maintain.o map.o set.o \
	variable.o vectorfield.o polyhedron.o

TARGET = lib${NAME}.a

CXXFLAGS = -O3 -Wall -fPIC -I../${INCLUDEDIR}

ifeq ($(DEBUG),yes)
	CXXFLAGS += -g
endif

ifeq ($(USE_OCTAVE),yes)
	CPPFLAGS += -DUSE_OCTAVE
	CXXFLAGS += -I${OCTAVE_HEADERS}
endif

ifeq ($(FASTER_BUT_DIRTIER),yes)
	CPPFLAGS += -DFASTER_BUT_DIRTIER
endif

