phony: lib extra

lib: libariadne.la libcalculus.la
extra: libariadne-extra.la

AM_CXXFLAGS = $(ARIADNE_CXXFLAGS)

INCLUDES = -I$(top_srcdir)/include $(CAIRO_CFLAGS) $(GTK_CFLAGS)
LIBS =  $(CAIRO_LIBS) $(GTK_LIBS)
LDADD =  $(LIBS)

ARIADNE_SRC = \
	numeric.cc \
	vector.cc \
	matrix.cc \
	differential.cc \
	affine.cc \
	polynomial.cc \
	propagator.cc \
	expression.cc \
	function.cc \
	taylor_model.cc  \
	taylor_function.cc  \
	solver.cc \
	integrator.cc \
	constraint_solver.cc \
	linear_programming.cc \
	nonlinear_programming.cc \
	point.cc \
	box.cc \
	curve.cc \
	orbit.cc \
	function_set.cc \
	taylor_set.cc \
	grid_set.cc \
	affine_set.cc \
	hybrid_set.cc \
	hybrid_enclosure.cc \
	hybrid_automaton-monolithic.cc \
	hybrid_automaton-composite.cc \
	map_evolver.cc \
	vector_field_evolver.cc \
	hybrid_evolver.cc \
	discretiser.cc \
	reachability_analyser.cc \
	logging.cc \
	graphics.cc

EXTRA_SRC = \
	geometry.cc \
	zonotope.cc \
	polytope.cc \
	polyhedron.cc \
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
	taylor_function.cc

lib_LTLIBRARIES = libariadne.la libariadne-extra.la libcalculus.la
libariadne_la_SOURCES = $(ARIADNE_SRC)
libariadne_extra_la_SOURCES = $(EXTRA_SRC)
libcalculus_la_SOURCES = $(CALCULUS_SRC)
nodist_EXTRA_libariadne_la_SOURCES = dummy.cc

EXTRA_DIST =
bin_PROGRAMS =
include_HEADERS =

rounding_CXXFLAGS = -O0
