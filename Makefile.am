SUBDIRS = include src test examples python

EXTRA_DIST = Doxyfile

CXXFLAGS += -Wall 

.PHONY: shiny lib python examples check-lib 
shiny: maintainer-clean
	-rm -rf *~ */*~
	-rm -f *.log */*.log */*.eps
	-rm -rf html latexb
	-rm -rf */.deps */.libs
	-rm -f aclocal.m4 config.* configure depcomp install-sh libtool ltmain.sh missing mkinstalldirs stamp-h* .tm_project.cache
	-rm -rf autom4te.cache
	-rm -f Makefile.in */Makefile.in
	-rm -f */Makefile
	-rm -f Makefile

lib: 
	(cd src; make lib)

lib-numeric:
	(cd src; make libnumeric.la)

python: lib
	(cd src/python; make python)

tests: lib
	(cd test; make tests)

examples: lib
	(cd examples; make examples)

check-lib: lib
	(cd test; make check)

check-numeric: lib-numeric
	(cd test; make check_numeric)

check-python: python
	(cd python; make check)

install-lib: lib python
	(cd src; make install)