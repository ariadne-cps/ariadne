INCLUDES = -I$(top_srcdir)/include
LIBS = $(GTK_LIBS)
LDADD = ../src/libariadne.la $(LIBS)

bin_PROGRAMS=tutorial watertank bouncingball ecoli
tutorial_SOURCES=tutorial.cc
watertank_SOURCES=watertank.cc
bouncingball_SOURCES = bouncingball.cc
ecoli_SOURCES=ecoli.cc

examples: $(bin_PROGRAMS)

