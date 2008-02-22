all: tests

#TESTS_ENVIRONMENT = (function run_tst { $1 > $1.log; }; run_tst )
TESTS_ENVIRONMENT = ./run_test.sh


INCLUDES = -I$(top_srcdir)/include
LDADD   = ../src/libariadne.la -lmpfr -lgmpxx 

BASE_TESTS = \
	test_tribool \
	test_array \
	test_ddconv \
	test_boost \
	test_mpfr 

NUMERIC_TESTS = \
	test_integer \
	test_rational \
	test_float \
	test_interval 

LINEAR_ALGEBRA_TESTS = \
	test_vector \
	test_matrix \
	test_linear_algebra

LINEAR_PROGRAMMING_TESTS = \
	test_linear_program \
	test_lpslv \
	test_lpstp \
	test_lpfsp \
	test_lptst 

COMBINATORIC_TESTS = \
	test_binary_word \
	test_lattice \
	test_subdivision_tree 

FUNCTION_TESTS = \
	test_index \
	test_interpreted_function \
	test_polynomial_function \
	test_affine_variable \
	test_affine_model \
	test_taylor_series \
	test_taylor_variable \
	test_taylor_derivative \
	test_taylor_model 

GEOMETRY_TESTS = \
	test_point \
	test_box \
	test_list_set \
	test_grid_set \
	test_partition_tree_set \
	test_rectangle \
	test_zonotope \
	test_simplex \
	test_polytope \
	test_polyhedron \
	test_constraint_set \
	test_level_set 

SYSTEM_TESTS = \
	test_hybrid_automaton 

EVALUATION_TESTS = \
	test_newton_solver \
	test_lohner_integrator \
	test_affine_integrator \
	test_map_evolver \
	test_vector_field_evolver \
	test_set_based_hybrid_evolver 

OUTPUT_TESTS = \
	test_textstream \
        test_latexstream \
        test_epsstream \
	test_chompstream 

check_PROGRAMS = \
	$(BASE_TESTS) \
	$(NUMERIC_TESTS) \
	$(LINEAR_ALGEBRA_TESTS) \
	$(LINEAR_PROGRAMMING_TESTS) \
	$(COMBINATORIC_TESTS) \
	$(FUNCTION_TESTS) \
	$(GEOMETRY_TESTS) \
	$(SYSTEM_TESTS) \
	$(EVALUATION_TESTS) \
	$(OUTPUT_TESTS)

tests: $(check_PROGRAMS)

test_integer_LDADD = ../src/libnumeric.la -lmpfr -lgmpxx -lgmp
test_rational_LDADD = ../src/libnumeric.la -lmpfr -lgmpxx -lgmp
test_float_LDADD = ../src/libnumeric.la -lmpfr -lgmpxx -lgmp
test_interval_LDADD = ../src/libnumeric.la -lmpfr -lgmpxx -lgmp
check_numericdir = ./
check_numeric_PROGRAMS = $(NUMERIC_TESTS) 
check_numeric: $(check_numeric_PROGRAMS)

#test_vector_LDADD = ../src/libnumeric.la ../src/liblinearalgebra.la -lmpfr -lgmpxx -lgmp
#test_matrix_LDADD = ../src/libnumeric.la ../src/liblinearalgebra.la -lmpfr -lgmpxx -lgmp
#test_tensor_LDADD = ../src/libnumeric.la ../src/liblinearalgebra.la -lmpfr -lgmpxx -lgmp
#test_linear_algebra_LDADD = ../src/libnumeric.la ../src/liblinearalgebra.la -lmpfr -lgmpxx -lgmp
check_linear_algebradir = ./
check_linear_algebra_PROGRAMS = $(LINEAR_ALGEBRA_TESTS) 
check_linear_algebra: $(check_linear_algebra_PROGRAMS)




# Explicitly list test sources
# This should not be necessary!!
test_tribool_SOURCES=test_tribool.cc
test_array_SOURCES=test_array.cc
test_ddconv_SOURCES=test_ddconv.cc
test_boost_SOURCES=test_boost.cc
test_mpfr_SOURCES=test_mpfr.cc

test_integer_SOURCES=test_integer.cc
test_rational_SOURCES=test_rational.cc
test_float_SOURCES=test_float.cc
test_interval_SOURCES=test_interval.cc

test_vector_SOURCES=test_vector.cc
test_matrix_SOURCES=test_matrix.cc
test_tensor_SOURCES=test_tensor.cc
test_linear_algebra_SOURCES=test_linear_algebra.cc

test_linear_program_SOURCES=test_linear_program.cc
test_lpslv_SOURCES=test_lpslv.cc
test_lpstp_SOURCES=test_lpstp.cc
test_lpfsp_SOURCES=test_lpfsp.cc
test_lptst_SOURCES=test_lptst.cc

test_binary_word_SOURCES=test_binary_word.cc
test_lattice_SOURCES=test_lattice.cc
test_subdivision_tree_SOURCES=test_subdivision_tree.cc

test_interpreted_function_SOURCES=test_interpreted_function.cc
test_polynomial_function_SOURCES=test_polynomial_function.cc
test_index_SOURCES=test_index.cc
test_taylor_series_SOURCES=test_taylor_series.cc
test_taylor_variable_SOURCES=test_taylor_variable.cc
test_taylor_derivative_SOURCES=test_taylor_derivative.cc
test_affine_model_SOURCES=test_affine_model.cc
test_taylor_model_SOURCES=test_taylor_model.cc

test_point_SOURCES=test_point.cc
test_list_set_SOURCES=test_list_set.cc
test_grid_set_SOURCES=test_grid_set.cc
test_partition_tree_set_SOURCES=test_partition_tree_set.cc
test_rectangle_SOURCES=test_rectangle.cc
test_zonotope_SOURCES=test_zonotope.cc
test_simplex_SOURCES=test_simplex.cc
test_polytope_SOURCES=test_polytope.cc
test_polyhedron_SOURCES=test_polyhedron.cc
test_constraint_set_SOURCES=test_constraint_set.cc
test_level_set_SOURCES=test_level_set.cc

test_hybrid_automaton_SOURCES=test_hybrid_automaton.cc

test_newton_solver_SOURCES=test_newton_solver.cc
test_lohner_integrator_SOURCES=test_lohner_integrator.cc
test_affine_integrator_SOURCES=test_affine_integrator.cc
test_map_evolver_SOURCES=test_map_evolver.cc
test_vector_field_evolver_SOURCES=test_vector_field_evolver.cc
test_set_based_hybrid_evolver_SOURCES=test_set_based_hybrid_evolver.cc

test_epsstream_SOURCES=test_epsstream.cc 
test_latexstream_SOURCES=test_latexstream.cc
test_textstream_SOURCES=test_textstream.cc
test_chompstream_SOURCES=test_chompstream.cc


EXTRA_DIST =
SUBDIRS =

TESTS = $(check_PROGRAMS)
XFAIL_TESTS = test_boost test_linear_program test_float
testsubdir = ./
