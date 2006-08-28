VERSION = 0.0.0
NAME = ariadne

# set the following variables 
DEBUG = no
VERBATIM = no

PREFIX=${HOME}
LIBPREFIX=${PREFIX}/lib/
PYTHONPREFIX=${PREFIX}/python/

LOCALDIR = ../..

CC = gcc
CXX = g++

MAKE = make

INCLUDEDIR = include
SRCDIR = src
TESTDIR = test
EXAMPLEDIR = examples
PYTHONDIR=python

LIBPPL = -lppl
LIBGMPXX = -lgmpxx -lgmp

REAL_TYPE=MPFLOAT

CXXFLAGS = -D${REAL_TYPE}_REAL 

SUBDIRS=${INCLUDEDIR} ${SRCDIR} ${TESTDIR} ${EXAMPLEDIR}

PYTHONVERSION=2.4
PYTHONINCLUDEDIR=/usr/include/python$(PYTHONVERSION)

RPATH = ${PREFIX}/lib
RPATHFLAGS = -Wl,-rpath -Wl,${RPATH}

ifeq ($(DEBUG),yes)
	CXXFLAGS += -g -DDEBUG -Wall -Wextra
	LIBS = libariadne.so
else
	CXXFLAGS += -O2 -fPIC -Wall  
	#CXXFLAGS += -O2 -fPIC -Wall -march=athlon64 -pipe 
	CXXFLAGS += -fomit-frame-pointer -funswitch-loops -fgcse-after-reload 
	LIBS = libariadne.so
endif

ifeq ($(VERBATIM),yes)
	CXXFLAGS = -DVERBATIM
endif

CXXFLAGS +=  -I${LOCALDIR}/${INCLUDEDIR} 
LDXXFLAGS = -L${LOCALDIR}/${SRCDIR}
