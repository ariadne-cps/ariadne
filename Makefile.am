SUBDIRS = include src test python

EXTRA_DIST = Doxyfile

CXXFLAGS += -Wall 

.PHONY: shiny lib check-lib
shiny: maintainer-clean
	-rm -rf *~ */*~
	-rm -f *.log */*.log */*.eps
	-rm -rf html latex
	-rm -rf */.deps */.libs
	-rm -f aclocal.m4 config.* configure depcomp install-sh libtool ltmain.sh missing mkinstalldirs stamp-h* .tm_project.cache
	-rm -rf autom4te.cache
	-rm -f Makefile.in */Makefile.in
	-rm -f */Makefile
	-rm -f Makefile

lib: 
	(cd src; make lib)

check-lib: lib
	(cd test; make check)

check-python:
	(cd python; make check)
