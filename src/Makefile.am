#AM_CXXFLAGS = -O2 -Wall

phony: lib

lib: libariadne.la

INCLUDES = -I$(top_srcdir)/include 
LDADD = -lmpfr -lppl -lgmpxx

#SUBDIRS = base numeric linear_algebra combinatoric geometry system evaluation output python

SUBDIRS = . python

NUMERIC_SRC = \
	numeric/approximation.cc numeric/arithmetic.cc numeric/function.cc 

LINEAR_ALGEBRA_SRC = \
	linear_algebra/vector.cc linear_algebra/matrix.cc \
	linear_algebra/tensor.cc \
	linear_algebra/linear_program.cc \
	linear_algebra/matrix_function.cc 

COMBINATORIC_SRC = \
	combinatoric/array_operations.cc \
	combinatoric/binary_tree.cc \
	combinatoric/lattice_set.cc combinatoric/lattice_map.cc \
	combinatoric/subdivision_tree_set.cc 

GEOMETRY_SRC = \
	geometry/ppl_polyhedron.cc \
	geometry/point.cc geometry/point_list.cc \
	geometry/rectangle.cc \
	geometry/polytope.cc geometry/polyhedron.cc \
	geometry/zonotope.cc geometry/parallelotope.cc geometry/simplex.cc \
	geometry/ellipsoid.cc geometry/sphere.cc \
	geometry/list_set.cc \
	geometry/grid.cc geometry/grid_set.cc \
	geometry/partition_tree_set.cc 

SYSTEM_SRC = \
	system/map.cc system/affine_map.cc system/polynomial_map.cc \
	system/vector_field.cc system/affine_vector_field.cc \
	system/multimap.cc system/affine_multimap.cc 

EVALUATION_SRC = \
	evaluation/newton.cc \
	evaluation/apply.cc 	evaluation/integrator.cc evaluation/lohner_integrator.cc 

OUTPUT_SRC = 

noinst_LTLIBRARIES = libnumeric.la liblinearalgebra.la libcombinatoric.la libgeometry.la libsystem.la libevaluation.la liboutput.la

lib_LTLIBRARIES = libariadne.la
libariadne_la_LIBADD = $(noinst_LTLIBRARIES) -lmpfr -lppl -lgmpxx
libariadne_la_SOURCES = 
nodist_EXTRA_libariadne_la_SOURCES = dummy.cc


libnumeric_la_SOURCES = $(NUMERIC_SRC)
liblinearalgebra_la_SOURCES = $(LINEAR_ALGEBRA_SRC)
libcombinatoric_la_SOURCES = $(COMBINATORIC_SRC)
libgeometry_la_SOURCES = $(GEOMETRY_SRC)
libsystem_la_SOURCES = $(SYSTEM_SRC)
libevaluation_la_SOURCES = $(EVALUATION_SRC)
liboutput_la_SOURCES = $(OUTPUT_SRC)


EXTRA_DIST =  
bin_PROGRAMS = 
include_HEADERS = 
