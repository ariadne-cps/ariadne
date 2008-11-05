SUBDIRS = include src test wrap examples 

EXTRA_DIST = Doxyfile

CXXFLAGS += -Wall 

.PHONY: shiny doc all lib wrap examples check-lib check-python 

all: lib wrap tests examples

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

doc:
	doxygen

lib: 
	(cd src; make lib)

wrap: lib
	(cd wrap; make wrap)

python: lib
	(cd wrap; make python)

tests: lib
	(cd test; make tests)

examples: lib
	(cd examples; make examples)

check-lib: lib
	(cd test; make check)

check-python: python
	(cd python; make check)

install-lib: lib
	(cd src; make install)