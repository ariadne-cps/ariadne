VERSION = 0.0.0
NAME = ariadne

# set the following variables 
DEBUG = no
VERBATIM = yes

PREFIX=${HOME}

CC = gcc
CXX = g++

MAKE = make

INCLUDEDIR = include
SRCDIR = src
TESTDIR = test
EXAMPLEDIR = examples
PYTHONDIR=python

LIBGMPXX = -lgmpxx -lgmp

SUBDIRS=${INCLUDEDIR} ${SRCDIR} ${TESTDIR} ${EXAMPLEDIR} ${PYTHONDIR}

PYTHONVERSION=2.4
PYTHONINCLUDEDIR=/usr/include/python$(PYTHONVERSION)

PYTHONINSTALLDIR=${PREFIX}/python/

CXXFLAGS =  -I../../${INCLUDEDIR}
LDXXFLAGS = 
RPATH = ${PREFIX}/lib
RPATHFLAGS = -Wl,-rpath -Wl,${RPATH}

ifeq ($(DEBUG),yes)
	CXXFLAGS += -g -DDEBUG -Wall
	LIBS = libariadne.a
else
	CXXFLAGS += -O2 -fPIC -Wall
	LIBS = libariadne.so
endif

ifeq ($(VERBATIM),yes)
	CXXFLAGS += -DVERBATIM
endif
