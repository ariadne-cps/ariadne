AM_CXXFLAGS = -Wall

INCLUDES = -I$(top_srcdir)/include
FLEX_LIBS = -lfl
LIBS +=  $(BOOST_LIBS) $(CAIRO_LIBS) $(GTK_LIBS) $(FLEX_LIBS)
LDADD = ../src/libariadne.la $(LIBS)

BUILT_SOURCES = cif_bison.h 
CLEANFILES = $(BUILT_SOURCES) cif_bison.cc cif_flex.cc
AM_YFLAGS = -d -v
bin_PROGRAMS=cifparser 

cifparser_SOURCES=cifparser.cc ciftypes.cc cif_bison.yy cif_flex.ll

cif2ariadne: $(bin_PROGRAMS)

