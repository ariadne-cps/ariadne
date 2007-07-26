#TESTS_ENVIRONMENT = (function run_tst { $1 > $1.log; }; run_tst )
TESTS_ENVIRONMENT = ./run_test.sh


INCLUDES = -I$(top_srcdir)/include
LDADD   = ../src/libariadne.la -lppl -lmpfr -lgmpxx 

check_PROGRAMS = \
	test_array \
	test_ddconv \
	test_boost \
	test_mpfr \
	test_integer \
	test_rational \
	test_float \
	test_interval \
	test_differential \
	test_vector \
	test_matrix \
	test_tensor \
	test_linear_algebra \
	test_linear_program \
	test_lpslv \
	test_lpstp \
	test_lpfsp \
	test_lptst \
	test_binary_word \
	test_lattice \
	test_subdivision_tree \
	test_point \
	test_list_set \
	test_grid_set \
	test_partition_tree_set \
	test_rectangle \
	test_zonotope \
	test_parallelotope \
	test_simplex \
	test_polytope \
	test_polyhedron \
	test_constraint \
	test_constraint_set test_level_set \
	test_hybrid_automaton \
	test_constraint_hybrid_automaton \
	test_affine_model \
	test_function \
	test_newton \
	test_apply \
	test_lohner_integrator \
	test_affine_integrator \
	test_integrate \
	test_chainreach \
	test_hybrid_evolution \
	test_constraint_hybrid_evolution \
        test_epsfstream \
	test_txtfstream \
	test_chompfstream 

XFAIL_TESTS = 

tests: $(check_PROGRAMS)

# Explicitly list test sources
# This should not be necessary!!
test_array_SOURCES=test_array.cc
test_ddconv_SOURCES=test_ddconv.cc
test_boost_SOURCES=test_boost.cc
test_mpfr_SOURCES=test_mpfr.cc
test_integer_SOURCES=test_integer.cc
test_rational_SOURCES=test_rational.cc
test_float_SOURCES=test_float.cc
test_interval_SOURCES=test_interval.cc
test_differential_SOURCES=test_differential.cc
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
test_constraint_SOURCES=test_constraint.cc
test_constraint_set_SOURCES=test_constraint_set.cc
test_level_set_SOURCES=test_level_set.cc
test_affine_model_SOURCES=test_affine_model.cc
test_hybrid_automaton_SOURCES=test_hybrid_automaton.cc
test_constraint_hybrid_automaton_SOURCES=test_constraint_hybrid_automaton.cc
test_function_SOURCES=test_function.cc
test_newton_SOURCES=test_newton.cc
test_apply_SOURCES=test_apply.cc
test_lohner_integrator_SOURCES=test_lohner_integrator.cc
test_affine_integrator_SOURCES=test_affine_integrator.cc
test_integrate_SOURCES=test_integrate.cc
test_chainreach_SOURCES=test_chainreach.cc
test_hybrid_evolution_SOURCES=test_hybrid_evolution.cc
test_constraint_hybrid_evolution_SOURCES=test_constraint_hybrid_evolution.cc
test_epsfstream_SOURCES=test_epsfstream.cc 
test_txtfstream_SOURCES=test_txtfstream.cc
test_chompfstream_SOURCES=test_chompfstream.cc

EXTRA_DIST =
SUBDIRS =

TESTS = $(check_PROGRAMS)
testsubdir = ./
