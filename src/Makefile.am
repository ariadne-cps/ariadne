all: lib lib-extras

.PHONY: all lib lib-kernel lib-geometry lib-dynamics lib-hybrid lib-extras show-env

show-env:
	@echo 'INCLUDES =' $(INCLUDES)
	@echo 'LIBS =' $(LIBS)


lib: libariadne.la libcalculus.la
lib-kernel: libariadne-kernel.la
lib-geometry: libariadne-geometry.la
lib-dynamics: libariadne-dynamics.la
lib-hybrid: libariadne-hybrid.la
lib-extras: libariadne-extras.la

AM_CXXFLAGS = $(ARIADNE_CXXFLAGS)
AM_LDFLAGS = $(EXEC_LDFLAGS)

INCLUDES = -iquote$(top_srcdir)/include $(CAIRO_CFLAGS) $(GTK_CFLAGS)
LIBS =  $(BOOST_SERIALIZATION_LIBS) $(CAIRO_LIBS) $(GTK_LIBS)
LDADD =  $(LIBS)

ARIADNE_KERNEL_SRC = \
	logging.cc \
	integer.cc \
	rational.cc \
	float.cc \
	interval.cc \
	real.cc \
	vector.cc \
	matrix.cc \
	differential.cc \
	affine.cc \
	polynomial.cc \
	operators.cc \
	function.cc \
	algebra.cc \
	taylor_model.cc  \
	affine_model.cc  \
	taylor_function.cc  \
	solver.cc \
	integrator.cc \
	constraint_solver.cc \
	simplex_algorithm.cc \
	linear_programming.cc \
	nonlinear_programming.cc

ARIADNE_GEOMETRY_SRC = \
	point.cc \
	box.cc \
	curve.cc \
	orbit.cc \
	function_set.cc \
	grid_set.cc \
	paver.cc \
	affine_set.cc \
	enclosure.cc \
	drawer.cc \
	graphics.cc

ARIADNE_DYNAMICS_SRC = \
	map_evolver.cc \
	vector_field_evolver.cc

ARIADNE_HYBRID_SRC = \
	expression.cc \
	expression_set.cc \
	hybrid_set.cc \
	hybrid_enclosure.cc \
	hybrid_automaton-monolithic.cc \
	hybrid_automaton-composite.cc \
	hybrid_automaton-restrictive.cc \
	hybrid_simulator.cc \
	hybrid_evolver.cc \
	hybrid_reachability_analyser.cc \
	hybrid_graphics.cc

ARIADNE_SRC = $(ARIADNE_KERNEL_SRC) $(ARIADNE_GEOMETRY_SRC) $(ARIADNE_DYNAMICS_SRC) $(ARIADNE_HYBRID_SRC)

EXTRAS_SRC = \
	runge_kutta_integrator.cc \
	geometry.cc \
	zonotope.cc \
	polytope.cc \
	polyhedron.cc \
	textplot.cc

CALCULUS_SRC = \
	float.cc \
	interval.cc \
	real.cc \
	vector.cc \
	matrix.cc \
	function.cc \
	polynomial.cc \
	differential.cc \
	taylor_model.cc  \
	taylor_function.cc

lib_LTLIBRARIES = libariadne.la libariadne-kernel.la libariadne-geometry.la libariadne-dynamics.la libariadne-hybrid.la libariadne-extras.la libcalculus.la
libariadne_la_SOURCES = $(ARIADNE_SRC)
libariadne_kernel_la_SOURCES = $(ARIADNE_KERNEL_SRC)
libariadne_geometry_la_SOURCES = $(ARIADNE_GEOMETRY_SRC)
libariadne_dynamics_la_SOURCES = $(ARIADNE_DYNAMICS_SRC)
libariadne_hybrid_la_SOURCES = $(ARIADNE_HYBRID_SRC)
libariadne_extras_la_SOURCES = $(EXTRAS_SRC)
libcalculus_la_SOURCES = $(CALCULUS_SRC)
nodist_EXTRA_libariadne_la_SOURCES = dummy.cc

EXTRA_DIST =
bin_PROGRAMS =
include_HEADERS =

rounding_CXXFLAGS = -O0
