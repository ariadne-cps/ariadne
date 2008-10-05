.PHONY: doc lib tests wrappers clean install

all: lib tests wrappers

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
	(cd src/; make dep);
	(cd test/; make dep);
	(cd wrap/; make dep);

depclean: 
	(cd ${SRCDIR}; ${MAKE} depclean);
	(cd ${TESTDIR}; ${MAKE} depclean);
