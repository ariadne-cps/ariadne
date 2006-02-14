VERSION = 0.0.0
NAME = ariadne

# set the following variables 
DEBUG = no
VERBATIM = yes

PREFIX=${HOME}

LOCALDIR = ../..

CC = gcc
CXX = g++

MAKE = make

INCLUDEDIR = include
SRCDIR = src
TESTDIR = test
EXAMPLEDIR = examples
PYTHONDIR=python

LIBGMPXX = -lgmpxx -lgmp

SUBDIRS=${INCLUDEDIR} ${SRCDIR} ${TESTDIR} ${EXAMPLEDIR}

PYTHONVERSION=2.4
PYTHONINCLUDEDIR=/usr/include/python$(PYTHONVERSION)

CXXFLAGS =  -I${LOCALDIR}/${INCLUDEDIR}
LDXXFLAGS = -L${LOCALDIR}/${SRCDIR} 

PYTHONINSTALLDIR=${PREFIX}/python/

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
