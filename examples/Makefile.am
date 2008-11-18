INCLUDES = -I$(top_srcdir)/include
LIBS = $(GTK_LIBS)
LDADD = ../src/libariadne.la $(LIBS)

bin_PROGRAMS=tutorial simple
tutorial_SOURCES=tutorial.cc
simple_SOURCES=simple.cc

examples: $(bin_PROGRAMS)

