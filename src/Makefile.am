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
	differential.cc \
	propagator.cc \
	taylor_model.cc  \
	expression.cc \
	affine.cc \
	polynomial.cc \
	function.cc \
	nonlinear_programming.cc \
	taylor_function.cc  \
	simulation_toolbox.cc \
	taylor_calculus.cc \
	solver.cc \
	integrator.cc \
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
	affine_set.cc \
	hybrid_set.cc \
	hybrid_enclosure.cc \
	geometry.cc \
	hybrid_automaton-monolithic.cc \
	hybrid_automaton-composite.cc \
	orbit.cc \
	constraint_solver.cc \
	map_evolver.cc \
	vector_field_evolver.cc \
	hybrid_evolver-stable.cc \
	hybrid_evolver-image.cc \
	hybrid_evolver-constrained.cc \
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
	taylor_set.cc  \
	taylor_function.cc  \
	taylor_calculus.cc

#	hybrid_simulator.cc
#	hybrid_system.cc


lib_LTLIBRARIES = libariadne.la libcalculus.la
libariadne_la_SOURCES = $(SRC)
libcalculus_la_SOURCES = $(CALCULUS_SRC)
nodist_EXTRA_libariadne_la_SOURCES = dummy.cc

EXTRA_DIST =
bin_PROGRAMS =
include_HEADERS =

rounding_CXXFLAGS = -O0
