include ./config.mk

.PHONY: doc lib python clean install install-lib 

all: lib python 

doc:
	doxygen;

lib: 
	(cd ${SRCDIR}; ${MAKE} lib);

python: lib
	(cd ${SRCDIR}/${PYTHONDIR}; ${MAKE} all);

clean:                                                                          
	rm -f *~ *% core 
	rm -f *.eps test_*.log
	rm -rf latex html
	(cd ${SRCDIR}; ${MAKE} clean; cd ..; )
	(cd ${TESTDIR}; ${MAKE} clean; cd ..; )
	(cd ${PYTHONDIR}; ${MAKE} clean; cd ..; )

install: lib python
	(cd ${SRCDIR}; ${MAKE} install);
	(cd ${PYTHONDIR}; ${MAKE} install);

install-lib: lib
	(cd ${SRCDIR}; ${MAKE} install-lib);

check: install
	(cd ${TESTDIR}; ${MAKE} check);
	(cd ${PYTHONDIR}; ${MAKE} check);

check-lib: lib
	(cd ${TESTDIR}; ${MAKE} check);

check-python: install
	(cd ${PYTHONDIR}; ${MAKE} check);

dep: 
	(cd ${SRCDIR}; ${MAKE} dep);
	(cd ${TESTDIR}; ${MAKE} dep);

depclean: 
	(cd ${SRCDIR}; ${MAKE} depclean);
	(cd ${TESTDIR}; ${MAKE} depclean);

