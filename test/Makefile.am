#TESTS_ENVIRONMENT = (function run_tst { $1 > $1.log; }; run_tst )
TESTS_ENVIRONMENT = ./run_test.sh


INCLUDES = -I$(top_srcdir)/include
LDADD   = ../src/libariadne.la -lppl -lmpfr -lgmpxx 

check_PROGRAMS = \
	test_array \
	test_ddconv test_mpfr \
	test_integer \
	test_arithmetic test_function test_interval \
	test_vector test_matrix test_tensor \
	test_linear_algebra test_linear_program \
	test_binary_word test_lattice test_subdivision_tree \
	test_point \
	test_rectangle test_zonotope test_parallelotope test_simplex \
	test_polytope test_polyhedron \
	test_list_set test_grid_set test_partition_tree_set \
	test_hybrid_automaton \
	test_newton test_apply test_integration_step test_integrate test_chainreach \
	test_epsfstream test_chompfstream 

XFAIL_TESTS = 

# Explicitly list test sources
# This should not be necessary!!
test_array_SOURCES=test_array.cc
test_ddconv_SOURCES=test_ddconv.cc
test_mpfr_SOURCES=test_mpfr.cc
test_integer_SOURCES=test_integer.cc
test_arithmetic_SOURCES=test_arithmetic.cc
test_function_SOURCES=test_function.cc
test_interval_SOURCES=test_interval.cc
test_vector_SOURCES=test_vector.cc
test_matrix_SOURCES=test_matrix.cc
test_tensor_SOURCES=test_tensor.cc
test_linear_algebra_SOURCES=test_linear_algebra.cc
test_linear_program_SOURCES=test_linear_program.cc
test_binary_word_SOURCES=test_binary_word.cc
test_lattice_SOURCES=test_lattice.cc
test_subdivision_tree_SOURCES=test_subdivision_tree.cc
test_point_SOURCES=test_point.cc
test_list_set_SOURCES=test_list_set.cc
test_grid_set_SOURCES=test_grid_set.cc
test_partition_tree_set_SOURCES=test_partition_tree_set.cc
test_rectangle_SOURCES=test_rectangle.cc
test_zonotope_SOURCES=test_zonotope.cc
test_parallelotope_SOURCES=test_parallelotope.cc
test_simplex_SOURCES=test_simplex.cc
test_polytope_SOURCES=test_polytope.cc
test_polyhedron_SOURCES=test_polyhedron.cc
test_hybrid_automaton_SOURCES=test_hybrid_automaton.cc
test_newton_SOURCES=test_newton.cc
test_apply_SOURCES=test_apply.cc
test_integration_step_SOURCES=test_integration_step.cc
test_integrate_SOURCES=test_integrate.cc
test_chainreach_SOURCES=test_chainreach.cc
test_epsfstream_SOURCES=test_epsfstream.cc
test_chompfstream_SOURCES=test_chompfstream.cc

EXTRA_DIST =
SUBDIRS =

TESTS = $(check_PROGRAMS)
testsubdir = ./
