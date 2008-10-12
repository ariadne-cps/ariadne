.PHONY: doc lib tests wrappers clean install

all: lib tests wrappers
	(cd src/; make plotting)

doc:
	doxygen;

lib:
	(cd src/; make lib)

main:
	(cd src/; make main)

tests:
	(cd test/; make tests)

wrappers:
	(cd wrap/; make wrappers)

clean:
	(rm -rf html/ latex/ *~)
	(cd src/; make clean)
	(cd test/; make clean)
	(cd wrap/; make clean)

dep: 
	(touch src/Makefile.dep test/Makefile.dep wrap/Makefile.dep)
	(cd src/; make dep)
	(cd test/; make dep)
	(cd wrap/; make dep)

depclean: 
	(cd ${SRCDIR}; ${MAKE} depclean);
	(cd ${TESTDIR}; ${MAKE} depclean);
