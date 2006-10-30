LIBS = -lboost_python -lmpfr -lppl -lgmpxx
AM_LDFLAGS = -module

pythondir=${prefix}/python
pyexecdir=${prefix}/python
libdir=${pyexecdir}

python_LIBS = base.la numeric.la linear_algebra.la \
	combinatoric.la \
	geometry.la system.la evaluation.la \
	output.la \
	models.la

python_MODULES = .libs/base.so .libs/numeric.so .libs/linear_algebra.so \
	.libs/combinatoric.so \
	.libs/geometry.so .libs/system.so .libs/evaluation.so \
	.libs/output.so \
	.libs/models.so

pkglib_LTLIBRARIES = $(python_LIBS)
pkgpython_DATA = __init__.py

#install-exec-local:
#	echo "Using the installer " $(INSTALL)
#	echo $(pyexecdir)
#	echo $(pkgpyexecdir)
#	echo $(python_MODULES)
#	echo $(python_LIBS)
#	mkdir -p $(pkgpyexecdir)
#	$(INSTALL) $(python_MODULES) $(pkgpyexecdir)/
#	$(INSTALL) __init__.py $(pkgpyexecdir)/


INCLUDES = -I$(top_srcdir)/include 
CXX_FLAGS = -fpic

EXTRA_DIST = __init__.py


base_la_SOURCES=$(BASE_MODULE_SRC)
numeric_la_SOURCES=$(NUMERIC_MODULE_SRC)
linear_algebra_la_SOURCES=$(LINEAR_ALGEBRA_MODULE_SRC)
combinatoric_la_SOURCES=$(COMBINATORIC_MODULE_SRC)
geometry_la_SOURCES=$(GEOMETRY_MODULE_SRC)
system_la_SOURCES=$(SYSTEM_MODULE_SRC)
evaluation_la_SOURCES=$(EVALUATION_MODULE_SRC)
output_la_SOURCES=$(OUTPUT_MODULE_SRC)
models_la_SOURCES=$(MODELS_SRC)

# name
base_la_LIBADD = ../libariadne.la
numeric_la_LIBADD = ../libariadne.la
linear_algebra_la_LIBADD = ../libariadne.la
combinatoric_la_LIBADD = ../libariadne.la
geometry_la_LIBADD = ../libariadne.la
system_la_LIBADD = ../libariadne.la
evaluation_la_LIBADD = ../libariadne.la
output_la_LIBADD = ../libariadne.la
models_la_LIBADD = ../libariadne.la

BASE_MODULE_SRC=base_module.cc \
	export_array.cc

COMBINATORIC_MODULE_SRC=combinatoric_module.cc \
	export_binary_tree.cc export_lattice_set.cc	export_lattice_map.cc


NUMERIC_MODULE_SRC=numeric_module.cc \
	export_numeric.cc export_interval.cc export_function.cc 

LINEAR_ALGEBRA_MODULE_SRC = linear_algebra_module.cc \
	export_vector.cc export_matrix.cc export_tensor.cc export_linear_program.cc

GEOMETRY_MODULE_SRC = geometry_module.cc \
	export_point.cc \
	export_rectangle.cc \
	export_zonotope.cc export_parallelotope.cc \
	export_polytope.cc export_simplex.cc \
	export_polyhedron.cc \
	export_list_set.cc \
	export_grid_set.cc \
	export_partition_tree_set.cc 

SYSTEM_MODULE_SRC = system_module.cc \
	export_map.cc export_affine_map.cc export_polynomial_map.cc \
	export_affine_multimap.cc \
	export_vector_field.cc export_affine_vector_field.cc 

EVALUATION_MODULE_SRC = evaluation_module.cc \
	export_solve.cc \
	export_apply.cc \
	export_integrate.cc

OUTPUT_MODULE_SRC = output_module.cc \
	export_postscript_output.cc

MODELS_SRC = models.cc \
	export_models.cc
