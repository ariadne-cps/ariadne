#TESTS_ENVIRONMENT = "export PYTHONPATH=$PYTHONPATH:$(top_srcdir)/python/.libs/"

pythondir=${prefix}/python
pkgpythondir=${pythondir}/ariadne

pkgpython_DATA = src/hybrid_system.py  src/hybrid_system_evolution.py

check_SCRIPTS = test/test_import.py

TESTS = $(check_SCRIPTS)

