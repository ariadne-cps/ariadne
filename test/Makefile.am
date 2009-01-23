all: tests

#TESTS_ENVIRONMENT = (function run_tst { $1 > $1.log; }; run_tst )
TESTS_ENVIRONMENT = ./run_test.sh


INCLUDES = -I$(top_srcdir)/include
LIBS = $(GTK_LIBS)
LDADD = ../src/libariadne.la $(LIBS)

check_PROGRAMS = \
	test_float \
	test_interval \
	test_vector \
	test_matrix \
	test_differential \
	test_differential_vector \
	test_approximate_taylor_model \
	test_taylor_function \
	test_taylor_set \
	test_grid_set \
	test_hybrid_set \
	test_discrete_evolution \
	test_continuous_evolution \
	test_hybrid_evolution \
	test_discretised_evolution \
	test_reachability_analysis \
	test_serialization \
	test_graphics 

tests: $(check_PROGRAMS)

#test_integer_LDADD = -lgmp $(LIBS)
#test_rational_LDADD = -lgmp $(LIBS)

# Explicitly list test sources
test_float_SOURCES=test_float.cc
test_interval_SOURCES=test_interval.cc

test_vector_SOURCES=test_vector.cc
test_matrix_SOURCES=test_matrix.cc

test_differential_SOURCES=test_differential.cc
test_differential_vector_SOURCES=test_differential_vector.cc

test_approximate_taylor_model_SOURCES=test_approximate_taylor_model.cc
test_taylor_function_SOURCES=test_taylor_function.cc

test_discrete_evolution_SOURCES=test_discrete_evolution.cc
test_continuous_evolution_SOURCES=test_continuous_evolution.cc
test_hybrid_evolution_SOURCES=test_hybrid_evolution.cc

test_grid_set_SOURCES=test_grid_set.cc
test_hybrid_set_SOURCES=test_hybrid_set.cc

test_discretised_evolution_SOURCES=test_discretised_evolution.cc
test_reachability_analysis_SOURCES=test_reachability_analysis.cc

test_serialization_SOURCES=test_serialization.cc
test_graphics_SOURCES=test_graphics.cc

EXTRA_DIST =
SUBDIRS =

TESTS = $(check_PROGRAMS)
XFAIL_TESTS = test_interval test_taylor_set test_taylor_function test_hybrid_evolution test_reachability_analysis
testsubdir = ./
