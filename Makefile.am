SUBDIRS = include src test wrap examples cifparser

EXTRA_DIST = Doxyfile

AM_CXXFLAGS = -Wall -fno-strict-aliasing

ACLOCAL_AMFLAGS = -I m4

.PHONY: shiny doc all lib wrap examples profiles check-lib check-python

all: lib wrap tests examples profiles doc 

# Show the values of the environment variables
evn:
	@echo 'INCLUDES =' $(INCLUDES)
	@echo 'LIBS =' $(LIBS)

shiny: maintainer-clean
	-rm -rf *~ */*~
	-rm -f *.log */*.log */*.eps
	-rm -rf docs 
	-rm -rf */.deps */.libs
	-rm -f aclocal.m4 config.* configure depcomp install-sh libtool ltmain.sh missing mkinstalldirs stamp-h* .tm_project.cache
	-rm -rf autom4te.cache
	-rm -f Makefile.in */Makefile.in
	-rm -f */Makefile
	-rm -f Makefile

doc: docs docs/refman.pdf docs/html

docs: 
	doxygen

docs/html: docs

docs/refman.pdf: docs
	(cd docs/latex; make)
	cp docs/latex/refman.pdf docs/refman.pdf 
lib:
	(cd src; make lib)

cifparser: lib
	(cd cifparser; make)

wrap: lib
	(cd wrap; make wrap)

python: lib
	(cd wrap; make python)

tests: lib
	(cd test; make tests)

profiles: lib
	(cd test; make profiles)

check-lib: lib
	(cd test; make check)

check-python: python
	(cd python; make check)

profile-lib: lib profiles
	(cd test; make profile-lib)

install-lib: lib
	(cd src; make install)

install-doc: doc
	install -d $(DESTDIR)/usr/share/doc/libariadne-doc/
	cp -rp docs/refman.pdf $(DESTDIR)/usr/share/doc/libariadne-doc/ariadne-manual.pdf
	cp -rp docs/html $(DESTDIR)/usr/share/doc/libariadne-doc/

install-examples: install-cxx-examples install-python-examples

install-cxx-examples:
	install -d $(DESTDIR)/usr/share/doc/libariadne-doc/examples/cxx
	cp -rp examples/*.h $(DESTDIR)/usr/share/doc/libariadne-doc/cxx/examples/
	cp -rp examples/*.cc $(DESTDIR)/usr/share/doc/libariadne-doc/cxx/examples/
	cp -rp examples/*.plt $(DESTDIR)/usr/share/doc/libariadne-doc/cxx/examples/

install-python-examples:
	install -d $(DESTDIR)/usr/share/doc/libariadne-doc/examples/python
	cp -rp python/*.py $(DESTDIR)/usr/share/doc/libariadne-doc/python/examples/


install-am: install-doc 

