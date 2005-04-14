VERSION = 0.0.0
NAME = ariadne

# set the following variables 
DEBUG = no
VERBATIM = yes

CC = gcc
CXX = g++

MAKE = make

EXAMPLEDIR = examples
INCLUDEDIR = include

CXXFLAGS = -I../${INCLUDEDIR}

LDXXFLAGS = -lppl -lgmpxx

ifeq ($(DEBUG),yes)
	CXXFLAGS += -g -DDEBUG -Wall
else
#	CXXFLAGS += -O3 -fPIC
endif

ifeq ($(VERBATIM),yes)
	CXXFLAGS += -DVERBATIM
endif

