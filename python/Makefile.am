#TESTS_ENVIRONMENT = "export PYTHONPATH=$PYTHONPATH:$(top_srcdir)/python/.libs/"
TESTS_ENVIRONMENT = export PYTHONPATH=$PYTHONPATH:$(top_srcdir)/src/python/.libs/ ; $(PYTHON)
#pythondir=${prefix}/lib/python2.4
#pyexecdir=${prefix}/lib/python2.4

#pkgpythondir=${pythondir}/ariadne

pkgpython_DATA = 

check_SCRIPTS = test/test_import.py test/test_numeric.py test/test_linear_algebra.py

TESTS = $(check_SCRIPTS)
