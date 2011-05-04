AM_CXXFLAGS = -Wall -fpic
AM_LDFLAGS = -module

INCLUDES = -I$(top_srcdir)/include $(PYTHON_INCLUDES)
LIBS += $(BOOST_PYTHON_LIB) $(PYTHON_LIB) $(CAIRO_LIBS) $(GTK_LIBS)

pythondir=${prefix}/lib/python$(PYTHON_VERSION)
libdir=${pyexecdir}

pythoninstalldir=${prefix}/lib/python$(PYTHON_VERSION)/dist-packages/ariadne

python_DATA = __init__.py
python_LIBS = ariadne.la interval.la differentiation.la
python_MODULES = .libs/ariadne.so .libs/interval.so .libs/differentiation.so

python_LTLIBRARIES =
if WANT_PYTHON_INTERFACE
  python_LTLIBRARIES += $(python_LIBS)
endif

wrap: python

python: $(python_LTLIBRARIES)

install-exec-local:
	echo "Using the installer " $(INSTALL)
	echo $(pyexecdir)
	echo $(pythoninstalldir)
	echo $(python_MODULES)
	echo $(python_LIBS)
	mkdir -p $(pythoninstalldir)
	$(INSTALL) $(python_MODULES) $(pythoninstalldir)/
	$(INSTALL) $(python_DATA) $(pythoninstalldir)/


interval_la_SOURCES=$(INTERVAL_SRC)
interval_la_LIBADD = ../src/numeric.lo

differentiation_la_SOURCES=$(DIFFERENTIATION_SRC)
differentiation_la_LIBADD = ../src/numeric.lo ../src/vector.lo ../src/matrix.lo

ariadne_la_SOURCES=$(ARIADNE_SRC)
ariadne_la_LIBADD = ../src/libariadne.la

INTERVAL_SRC = interval_module.cc numeric_submodule.cc linear_algebra_submodule.cc

DIFFERENTIATION_SRC = differentiation_module.cc differentiation_submodule.cc

ARIADNE_SRC = utilities.cc \
	ariadne_module.cc \
	numeric_submodule.cc \
	linear_algebra_submodule.cc \
	optimization_submodule.cc \
	differentiation_submodule.cc \
	function_submodule.cc \
	calculus_submodule.cc \
	geometry_submodule.cc \
	solver_submodule.cc \
	storage_submodule.cc \
	evolution_submodule.cc \
	graphics_submodule.cc

