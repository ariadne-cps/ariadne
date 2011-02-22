all: tests profiles

profile-lib: $(bin_PROGRAMS)
	for file in $(bin_PROGRAMS); do ./run_profile.sh $$file; done

TESTS_ENVIRONMENT = ./run_test.sh
PROFILES_ENVIRONMENT = ./run_profile.sh

AM_CXXFLAGS = $(ARIADNE_CXXFLAGS)
AM_LDFLAGS = $(EXEC_LDFLAGS)

INCLUDES = -iquote$(top_srcdir)/include $(CAIRO_CFLAGS) $(GTK_CFLAGS)
LIBS += $(CAIRO_LIBS) $(GTK_LIBS)
LDADD = ../src/libariadne.la $(LIBS)

bin_PROGRAMS = \
	profile_arithmetic \
	profile_taylor_model

#test_misc_LDADD = $(LDADD) -lboost_python -lpython2.5

check_PROGRAMS = \
	test_float \
	test_interval \
	test_vector \
	test_matrix \
	test_multi_index \
	test_expansion \
	test_polynomial \
	test_differential \
	test_procedure \
	test_function \
	test_expression \
	test_taylor_model \
	test_taylor_function \
	test_function_sets \
	test_affine_sets \
	test_taylor_sets \
	test_grid_set \
	test_linear_programming \
	test_nonlinear_programming \
	test_integrator \
	test_constraint_solver \
	test_solvers \
	test_discrete_evolution \
	test_continuous_evolution \
	test_discretiser \
	test_graphics \
	test_serialization \
	test_hybrid_set \
	test_hybrid_automaton \
	test_hybrid_evolver \
	test_hybrid_evolution \
	test_hybrid_discretiser \
	test_hybrid_reachability_analyser

EXTRA_TESTS = \
	test_textplot \
	test_misc

tests: $(check_PROGRAMS)
profiles: $(bin_PROGRAMS)

#test_integer_LDADD = -lgmp $(LIBS)
#test_rational_LDADD = -lgmp $(LIBS)

# Explicitly list profile sources
profile_arithmetic_SOURCES=profile_arithmetic.cc
profile_taylor_model_SOURCES=profile_taylor_model.cc

# Explicitly list test sources
test_float_SOURCES=test_float.cc
test_interval_SOURCES=test_interval.cc
test_vector_SOURCES=test_vector.cc
test_matrix_SOURCES=test_matrix.cc
test_multi_index_SOURCES=test_multi_index.cc
test_expansion_SOURCES=test_expansion.cc
test_polynomial_SOURCES=test_polynomial.cc
test_differential_SOURCES=test_differential.cc
test_taylor_model_SOURCES=test_taylor_model.cc
test_taylor_function_SOURCES=test_taylor_function.cc
test_function_sets_SOURCES=test_function_sets.cc
test_affine_sets_SOURCES=test_affine_sets.cc
test_taylor_sets_SOURCES=test_taylor_sets.cc
test_expression_SOURCES=test_expression.cc
test_function_SOURCES=test_function.cc
test_grid_set_SOURCES=test_grid_set.cc
test_nonlinear_programming_SOURCES=test_nonlinear_programming.cc
test_integrator_SOURCES=test_integrator.cc
test_constraint_solver_SOURCES=test_constraint_solver.cc
test_solvers_SOURCES=test_solvers.cc
test_discrete_evolution_SOURCES=test_discrete_evolution.cc
test_continuous_evolution_SOURCES=test_continuous_evolution.cc
test_discretiser_SOURCES=test_discretiser.cc
test_hybrid_set_SOURCES=test_hybrid_set.cc
test_hybrid_automaton_SOURCES=test_hybrid_automaton.cc
test_hybrid_evolver_SOURCES=test_hybrid_evolver.cc
test_hybrid_evolution_SOURCES=test_hybrid_evolution.cc
test_hybrid_discretiser_SOURCES=test_hybrid_discretiser.cc
test_hybrid_reachability_analyser_SOURCES=test_hybrid_reachability_analyser.cc
test_serialization_SOURCES=test_serialization.cc
test_graphics_SOURCES=test_graphics.cc
test_textplot_SOURCES=test_textplot.cc
test_misc_SOURCES=test_misc.cc

EXTRA_DIST =
SUBDIRS =

TESTS = $(check_PROGRAMS)
XFAIL_TESTS = test_integrator test_discretiser
testsubdir = ./
