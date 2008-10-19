all: tests

#TESTS_ENVIRONMENT = (function run_tst { $1 > $1.log; }; run_tst )
TESTS_ENVIRONMENT = ./run_test.sh


INCLUDES = -I$(top_srcdir)/include
LIBS = -lmpfr -lgmpxx -lgmp $(GTK_LIBS) -lboost_serialization 
LDADD = ../src/libariadne.la $(LIBS)

check_PROGRAMS = \
	test_differential \
	test_approximate_taylor_model \
	test_discrete_evolution \
	test_grid_set \
	test_serialization \
	test_graphics

tests: $(check_PROGRAMS)

#test_integer_LDADD = -lgmp $(LIBS)
#test_rational_LDADD = -lgmp $(LIBS)

# Explicitly list test sources
#test_tribool_SOURCES=test_tribool.cc
#test_array_SOURCES=test_array.cc
#test_boost_SOURCES=test_boost.cc
#test_mpfr_SOURCES=test_mpfr.cc

#test_integer_SOURCES=test_integer.cc
#test_rational_SOURCES=test_rational.cc
#test_float_SOURCES=test_float.cc
#test_interval_SOURCES=test_interval.cc

#test_vector_SOURCES=test_vector.cc
#test_matrix_SOURCES=test_matrix.cc

#test_multi_index_SOURCES=test_multi_index.cc
test_differential_SOURCES=test_differential.cc

test_approximate_taylor_model_SOURCES=test_approximate_taylor_model.cc

test_serialization_SOURCES=test_serialization.cc
test_graphics_SOURCES=test_graphics.cc

EXTRA_DIST =
SUBDIRS =

TESTS = $(check_PROGRAMS)
XFAIL_TESTS = 
testsubdir = ./
