#AM_CXXFLAGS = -O2 -Wall

phony: lib

lib: libariadne.la libcalculus.la

AM_CXXFLAGS = -Wall

INCLUDES = -I$(top_srcdir)/include $(CAIRO_CFLAGS) $(GTK_CFLAGS)
LIBS =  $(CAIRO_LIBS) $(GTK_LIBS)
LDADD =  $(LIBS)

SRC = \
	numeric.cc \
	vector.cc \
	matrix.cc \
	linear_programming.cc \
	polynomial.cc \
	differential.cc \
	taylor_model.cc  \
	taylor_expression.cc  \
	formula.cc \
	function.cc \
	taylor_function.cc  \
	simulation_toolbox.cc \
	taylor_calculus.cc \
	newton_solver.cc \
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
	map_evolver.cc \
	vector_field_evolver.cc \
	hybrid_evolver.cc \
	stable_hybrid_evolver.cc \
	discretiser.cc \
	reachability_analyser.cc \
	logging.cc \
	graphics.cc \
	textplot.cc

CALCULUS_SRC = \
	numeric.cc \
	vector.cc \
	matrix.cc \
	function.cc \
	polynomial.cc \
	differential.cc \
	taylor_model.cc  \
	taylor_expression.cc  \
	taylor_set.cc  \
	taylor_function.cc  \
	taylor_calculus.cc

#	hybrid_simulator.cc

#	taylor_model.cc
#	taylor_set.cc

lib_LTLIBRARIES = libariadne.la libcalculus.la
libariadne_la_SOURCES = $(SRC)
libcalculus_la_SOURCES = $(CALCULUS_SRC)
nodist_EXTRA_libariadne_la_SOURCES = dummy.cc

EXTRA_DIST =
bin_PROGRAMS =
include_HEADERS =

rounding_CXXFLAGS = -O0
