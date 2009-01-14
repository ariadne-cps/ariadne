#AM_CXXFLAGS = -O2 -Wall

phony: lib

lib: libariadne.la

INCLUDES = -I$(top_srcdir)/include $(CAIRO_CFLAGS) $(GTK_CFLAGS)
LDADD = $(GTK_LIBS)

SRC = \
	numeric.cc \
	vector.cc \
	matrix.cc \
	linear_programming.cc \
	dense_differential.cc \
	sparse_differential.cc \
	taylor_variable.cc  \
	approximate_taylor_model.cc  \
	differential_calculus.cc \
	point.cc \
	box.cc \
	zonotope.cc \
	polytope.cc \
	polyhedron.cc \
	curve.cc \
	function_set.cc \
	taylor_set.cc \
	grid_set.cc \
	approximate_set.cc \
	hybrid_set.cc \
	geometry.cc \
	hybrid_automaton.cc \
	orbit.cc \
	hybrid_evolver.cc \
	discrete_evolver.cc \
	reachability_analyser.cc \
	logging.cc \
	graphics.cc 



lib_LTLIBRARIES = libariadne.la 
libariadne_la_SOURCES = $(SRC)
nodist_EXTRA_libariadne_la_SOURCES = dummy.cc

EXTRA_DIST =  
bin_PROGRAMS = 
include_HEADERS =

rounding_CXXFLAGS = -O0
