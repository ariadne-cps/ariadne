LIBS = -lboost_python -lmpfr -lgmpxx
AM_LDFLAGS = -module
CXX_FLAGS = -Wall -O0 -fpic
INCLUDES = -I$(top_srcdir)/include 

pythondir=${pkgpyexecdir}
#pyexecdir=${prefix}/lib/python2.4
libdir=${pyexecdir}

python_DATA = __init__.py
python_LIBS = ariadne.la
python_MODULES = .libs/ariadne.so

lib_LTLIBRARIES = $(python_LIBS)

python: $(python_LIBS)

install-exec-local:
	echo "Using the installer " $(INSTALL)
	echo $(pyexecdir)
	echo $(pkgpyexecdir)
	echo $(python_MODULES)
	echo $(python_LIBS)
	mkdir -p $(pkgpyexecdir)
	$(INSTALL) $(python_MODULES) $(pkgpyexecdir)/
	$(INSTALL) $(python_DATA) $(pkgpyexecdir)/


ariadne_la_SOURCES=$(ARIADNE_SRC)
ariadne_la_LIBADD = ../libariadne.la

ARIADNE_SRC = ariadne_module.cc \
	export_debug.cc \
	export_logging.cc \
	export_exceptions.cc \
	\
	export_tribool.cc \
	export_array.cc \
	 \
	export_binary_tree.cc \
	export_lattice_set.cc \
	export_lattice_map.cc \
	 \
	export_integer.cc \
	export_rational.cc \
	export_float.cc \
	export_interval.cc \
	 \
	export_vector.cc \
	export_covector.cc \
	export_matrix.cc \
	export_tensor.cc \
	export_linear_program.cc \
	 \
	export_function.cc \
	export_taylor_series.cc \
	export_taylor_variable.cc \
	export_taylor_derivative.cc \
	export_polynomial_variable.cc \
	 \
	export_point.cc \
	export_box.cc \
	export_set.cc \
	export_interval_set.cc \
	export_rectangle.cc \
	export_orthotope.cc \
	export_zonotope.cc \
	export_parallelotope.cc \
	export_polytope.cc \
	export_simplex.cc \
	export_polyhedron.cc \
	export_empty_set.cc \
	export_rectangular_set.cc \
	export_polyhedral_set.cc \
	export_list_set.cc \
	export_grid_set.cc \
	export_partition_tree_set.cc \
	export_constraint_set.cc \
	export_level_set.cc \
	export_hybrid_set.cc \
	 \
	export_constraint.cc \
	 \
	export_map.cc \
	export_affine_map.cc \
	export_polynomial_map.cc \
	export_affine_multimap.cc \
	export_vector_field.cc \
	export_affine_vector_field.cc \
	export_set_based_hybrid_automaton.cc \
	export_constraint_based_hybrid_automaton.cc \
	 \
	export_evolution_parameters.cc \
	export_solver.cc \
	export_applicator.cc \
	export_integrator.cc \
	export_detector.cc \
	export_map_evolver.cc \
	export_vector_field_evolver.cc \
	export_set_based_hybrid_evolver.cc \
	export_constraint_based_hybrid_evolver.cc \
	 \
	export_text_output.cc \
	export_latex_output.cc \
	export_postscript_output.cc \
	 \
	export_models.cc
