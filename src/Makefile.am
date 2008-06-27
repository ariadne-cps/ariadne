#AM_CXXFLAGS = -O2 -Wall

phony: lib

lib: libariadne.la

INCLUDES = -I$(top_srcdir)/include 
LDADD = -lmpfr -lgmpxx

#SUBDIRS = base numeric linear_algebra linear_programming combinatoric function geometry system evaluation output input python

SUBDIRS = . python

BASE_SRC = 

NUMERIC_SRC = \
	numeric/double.cc \
	numeric/integer.cc \
	numeric/rational.cc \
	numeric/float.cc \
	numeric/interval.cc

LINEAR_ALGEBRA_SRC = \
	linear_algebra/vector.cc \
	linear_algebra/covector.cc \
	linear_algebra/matrix.cc \
	linear_algebra/matrix_function.cc 

LINEAR_PROGRAMMING_SRC = \
	linear_programming/linear_program.cc 

COMBINATORIC_SRC = \
	combinatoric/array_operations.cc \
	combinatoric/binary_word.cc \
	combinatoric/binary_tree.cc \
	combinatoric/lattice_set.cc \
	combinatoric/lattice_map.cc \
	combinatoric/subdivision_tree_set.cc 

DIFFERENTIATION_SRC = \
	differentiation/affine_variable.cc \
	differentiation/taylor_series.cc \
	differentiation/taylor_variable.cc \
	differentiation/taylor_derivative.cc \
	differentiation/sparse_differential.cc

FUNCTION_SRC = \
	function/affine_model.cc \
	function/taylor_model.cc  \
	function/approximate_taylor_model.cc  \
	function/constant_function.cc \
	function/identity_function.cc \
	function/affine_function.cc \
	function/polynomial_function.cc


GEOMETRY_SRC = \
	geometry/point.cc \
	geometry/interval_set.cc \
	geometry/box.cc \
	geometry/box_list_set.cc \
	geometry/point_list.cc \
	geometry/rectangle.cc \
	geometry/polyhedron.cc \
	geometry/polytope.cc \
	geometry/simplex.cc \
	geometry/zonotope.cc \
	geometry/ellipsoid.cc \
	geometry/sphere.cc \
	geometry/taylor_set.cc \
	geometry/list_set.cc \
	geometry/curve.cc \
	geometry/constraint.cc \
	geometry/poincare_section.cc \
	geometry/constraint_set.cc \
	geometry/level_set.cc \
	geometry/hybrid_set.cc 

STORAGE_SRC = \
	geometry/grid.cc \
	geometry/grid_approximation.cc \
	geometry/grid_cell.cc \
	geometry/grid_block.cc \
	geometry/grid_cell_list_set.cc \
	geometry/grid_mask_set.cc \
	geometry/partition_tree_set.cc 

SYSTEM_SRC = \
	system/numerical_system.cc \
	system/transition_system.cc \
	system/map.cc \
	system/vector_field.cc \
	system/hybrid_automaton.cc

#	system/constraint_based_hybrid_automaton.cc \
#	system/grid_multimap.cc 


EVALUATION_SRC = \
	evaluation/evolution_parameters.cc \
	\
	evaluation/solver.cc \
	\
	evaluation/approximator.cc \
	evaluation/satisfier.cc \
	evaluation/reducer.cc \
	evaluation/subdivider.cc \
	\
	evaluation/applicator.cc \
	evaluation/affine_applicator.cc \
	\
	evaluation/bounder.cc \
	evaluation/flower.cc \
	evaluation/integrator.cc \
	evaluation/euler_integrator.cc \
	evaluation/affine_integrator.cc \
	\
	evaluation/evolver.cc \
	evaluation/discrete_evolver.cc \
	evaluation/discretiser.cc \
	evaluation/reachability_analyser.cc \
	evaluation/model_checker.cc


#	evaluation/map_evolver.cc \
#	evaluation/vector_field_evolver.cc \
#	evaluation/set_based_hybrid_evolver.cc \
#	evaluation/constraint_based_hybrid_evolver.cc  

OUTPUT_SRC = \
	output/logging.cc \
	output/colour.cc \
	output/geometry2d.cc \
	output/textstream.cc \
	output/latexstream.cc \
	output/epsstream.cc \
	output/chompstream.cc \
	output/dotstream.cc 

INPUT_SRC = 



noinst_LTLIBRARIES = libbase.la libnumeric.la liblinearalgebra.la liblinearprogramming.la \
	libcombinatoric.la libdifferentiation.la \
	libfunction.la libgeometry.la libstorage.la libsystem.la libevaluation.la \
	liboutput.la libinput.la

lib_LTLIBRARIES = libariadne.la libgeometric.la 
libariadne_la_LIBADD = $(noinst_LTLIBRARIES) -lmpfr -lgmpxx
libariadne_la_SOURCES = 
nodist_EXTRA_libariadne_la_SOURCES = dummy.cc

libgeometric_la_LIBADD = libbase.la libnumeric.la liblinearalgebra.la liblinearprogramming.la \
	libcombinatoric.la libdifferentiation.la libfunction.la libgeometry.la -lmpfr -lgmpxx
libgeometric_la_SOURCES = 
nodist_EXTRA_libgeometric_la_SOURCES = dummy.cc


libbase_la_SOURCES = $(BASE_SRC)
libnumeric_la_SOURCES = $(NUMERIC_SRC)
liblinearalgebra_la_SOURCES = $(LINEAR_ALGEBRA_SRC)
liblinearprogramming_la_SOURCES = $(LINEAR_PROGRAMMING_SRC)
libcombinatoric_la_SOURCES = $(COMBINATORIC_SRC)
libdifferentiation_la_SOURCES = $(DIFFERENTIATION_SRC)
libfunction_la_SOURCES = $(FUNCTION_SRC)
libgeometry_la_SOURCES = $(GEOMETRY_SRC)
libstorage_la_SOURCES = $(STORAGE_SRC)
libsystem_la_SOURCES = $(SYSTEM_SRC)
libevaluation_la_SOURCES = $(EVALUATION_SRC)
liboutput_la_SOURCES = $(OUTPUT_SRC)
libinput_la_SOURCES = $(INPUT_SRC)


EXTRA_DIST =  
bin_PROGRAMS = 
include_HEADERS =

rounding_CXXFLAGS = -O0
