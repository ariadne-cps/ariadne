#AM_CXXFLAGS = -O2 -Wall

phony: lib

lib: libgeometric.la libariadne.la

INCLUDES = -I$(top_srcdir)/include 
LDADD = -lmpfr -lppl -lgmpxx

#SUBDIRS = base numeric linear_algebra linear_programming combinatoric function geometry system evaluation output input python

SUBDIRS = . python

BASE_SRC = 

NUMERIC_SRC = \
	numeric/integer.cc numeric/rational.cc \
	numeric/float64.cc numeric/floatmp.cc \
	numeric/interval.cc

LINEAR_ALGEBRA_SRC = \
	linear_algebra/vector.cc linear_algebra/matrix.cc \
	linear_algebra/tensor.cc \
	linear_algebra/matrix_function.cc 

LINEAR_PROGRAMMING_SRC = \
	linear_programming/linear_program.cc 

COMBINATORIC_SRC = \
	combinatoric/array_operations.cc \
	combinatoric/binary_word.cc combinatoric/binary_tree.cc \
	combinatoric/lattice_set.cc combinatoric/lattice_map.cc \
	combinatoric/subdivision_tree_set.cc 

FUNCTION_SRC = \
	function/interpreted_function.cc \
	function/affine_model.cc \
	function/polynomial_model.cc \
	function/polynomial.cc

GEOMETRY_SRC = \
	geometry/point.cc geometry/point_list.cc \
	geometry/rectangle.cc \
	geometry/polytope.cc geometry/polyhedron.cc \
	geometry/zonotope.cc geometry/parallelotope.cc geometry/simplex.cc \
	geometry/ellipsoid.cc geometry/sphere.cc \
	geometry/taylor_set.cc \
	geometry/list_set.cc \
	geometry/grid.cc \
	geometry/grid_approximation.cc \
	geometry/grid_cell.cc geometry/grid_block.cc \
	geometry/grid_cell_list_set.cc geometry/grid_mask_set.cc \
	geometry/irregular_grid.cc geometry/irregular_grid_set.cc \
	geometry/partition_tree_set.cc \
	geometry/curve.cc \
	geometry/constraint.cc geometry/poincare_section.cc \
	geometry/constraint_set.cc geometry/level_set.cc \
	geometry/hybrid_set.cc \
	geometry/ppl_polyhedron.cc

SYSTEM_SRC = \
	system/map.cc system/function_map.cc system/affine_map.cc system/polynomial_map.cc \
	system/vector_field.cc system/affine_vector_field.cc \
	system/multimap.cc system/affine_multimap.cc \
	system/set_based_hybrid_automaton.cc system/constraint_based_hybrid_automaton.cc \
	system/grid_multimap.cc 


EVALUATION_SRC = \
	evaluation/evolution_parameters.cc \
	evaluation/solver.cc \
	evaluation/newton_solver.cc \
	evaluation/euler_integrator.cc \
	evaluation/affine_integrator.cc \
	evaluation/lohner_integrator.cc \
	evaluation/bounder.cc \
	evaluation/detector.cc \
	evaluation/flow.cc \
	evaluation/applicator.cc \
	evaluation/integrator.cc \
	evaluation/set_based_hybrid_evolver.cc \
	evaluation/constraint_based_hybrid_evolver.cc \
	evaluation/constraint_based_hybrid_evolver_plugin.cc 

OUTPUT_SRC = \
	output/logging.cc \
	output/colour.cc \
	output/geometry2d.cc \
	output/epsstream.cc \
	output/textstream.cc \
	output/chompstream.cc \
	output/dotstream.cc \
	output/texstream.cc 

INPUT_SRC = \
	input/modelica.cc 


noinst_LTLIBRARIES = libbase.la libnumeric.la liblinearalgebra.la liblinearprogramming.la \
	libcombinatoric.la libfunction.la libgeometry.la libsystem.la libevaluation.la \
	liboutput.la libinput.la

lib_LTLIBRARIES = libariadne.la libgeometric.la 
libariadne_la_LIBADD = $(noinst_LTLIBRARIES) -lmpfr -lppl -lgmpxx
libariadne_la_SOURCES = 
nodist_EXTRA_libariadne_la_SOURCES = dummy.cc

libgeometric_la_LIBADD = libbase.la libnumeric.la liblinearalgebra.la liblinearprogramming.la \
	libcombinatoric.la libfunction.la libgeometry.la -lmpfr -lppl -lgmpxx
libgeometric_la_SOURCES = 
nodist_EXTRA_libgeometric_la_SOURCES = dummy.cc


libbase_la_SOURCES = $(BASE_SRC)
libnumeric_la_SOURCES = $(NUMERIC_SRC)
liblinearalgebra_la_SOURCES = $(LINEAR_ALGEBRA_SRC)
liblinearprogramming_la_SOURCES = $(LINEAR_PROGRAMMING_SRC)
libcombinatoric_la_SOURCES = $(COMBINATORIC_SRC)
libfunction_la_SOURCES = $(FUNCTION_SRC)
libgeometry_la_SOURCES = $(GEOMETRY_SRC)
libsystem_la_SOURCES = $(SYSTEM_SRC)
libevaluation_la_SOURCES = $(EVALUATION_SRC)
liboutput_la_SOURCES = $(OUTPUT_SRC)
libinput_la_SOURCES = $(INPUT_SRC)


EXTRA_DIST =  
bin_PROGRAMS = 
include_HEADERS =
