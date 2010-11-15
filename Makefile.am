SUBDIRS = include src test wrap examples

EXTRA_DIST = Doxyfile

.PHONY: shiny doc all lib lib-extras wrap examples profiles check-lib check-python

all: lib lib-extras wrap tests examples profiles

# No point in setting AM_CXXFLAGS here as these do not propagate to subdirectories
# Global CXXFLAGS are defined in configure.ac in the ARIADNE_CXXFLAGS variable.

# Show the values of the environment variables
evn:
	@echo 'INCLUDES =' $(INCLUDES)
	@echo 'LIBS =' $(LIBS)

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
	(cd latex; make)

lib:
	(cd src; make lib)

lib-extras:
	(cd src; make lib-extras)

wrap: lib lib-extras
	(cd wrap; make wrap)

python: lib lib-extras
	(cd wrap; make python)

tests: lib
	(cd test; make tests)

examples: lib lib-extras
	(cd examples; make examples)

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
