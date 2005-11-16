VERSION = 0.0.0
NAME = ariadne

# set the following variables 
DEBUG = no
VERBATIM = yes

CC = gcc
CXX = g++

MAKE = make

INCLUDEDIR = include
SRCDIR = src
TESTDIR = test
EXAMPLEDIR = examples
PYTHONDIR=python

SUBDIRS=${INCLUDEDIR} ${SRCDIR} ${TESTDIR} ${EXAMPLEDIR} ${PYTHONDIR}

PYTHONVERSION=2.4
PYTHONINCLUDEDIR=/usr/include/python$(PYTHONVERSION)

CXXFLAGS =  -I../${INCLUDEDIR}
LDXXFLAGS = -lppl -lgmpxx -lgmp

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
