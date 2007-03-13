LIBS = -lboost_python -lmpfr -lppl -lgmpxx
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
	\
	export_tribool.cc \
	export_array.cc \
	 \
	export_binary_tree.cc \
	export_lattice_set.cc \
	export_lattice_map.cc \
	 \
	export_numeric.cc \
	export_interval.cc \
	export_function.cc \
	 \
	export_vector.cc \
	export_matrix.cc \
	export_tensor.cc \
	export_linear_program.cc \
	 \
	export_point.cc \
	export_set.cc \
	export_rectangle.cc \
	export_zonotope.cc \
	export_parallelotope.cc \
	export_polytope.cc \
	export_simplex.cc \
	export_polyhedron.cc \
	export_polyhedral_set.cc \
	export_list_set.cc \
	export_grid_set.cc \
	export_partition_tree_set.cc \
	export_hybrid_set.cc \
	 \
	export_map.cc \
	export_affine_map.cc \
	export_polynomial_map.cc \
	export_affine_multimap.cc \
	export_vector_field.cc \
	export_affine_vector_field.cc \
	export_hybrid_automaton.cc \
	 \
	export_solve.cc \
	export_apply.cc \
	export_integrate.cc \
	export_hybrid_evolver.cc \
	 \
	export_postscript_output.cc \
	export_models.cc


