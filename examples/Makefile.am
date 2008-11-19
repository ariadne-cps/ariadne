INCLUDES = -I$(top_srcdir)/include
LIBS = $(GTK_LIBS)
LDADD = ../src/libariadne.la $(LIBS)

bin_PROGRAMS=tutorial watertank
tutorial_SOURCES=tutorial.cc
watertank_SOURCES=watertank.cc

examples: $(bin_PROGRAMS)

