include ./config.mk

all:  python 

lib: 
	(cd ${SRCDIR}; ${MAKE} lib);

python: lib
	(cd ${SRCDIR}/${PYTHONDIR}; ${MAKE} all);

clean:                                                                          
	rm -f *.bak *.o *~ *% core 
	rm -f *.eps
	rm -rf latex html
	(cd ${SRCDIR}; ${MAKE} clean; cd ..; )
	(cd ${EXAMPLEDIR}; ${MAKE} clean; cd ..; )
	(cd ${TESTDIR}; ${MAKE} clean; cd ..; )

install:
	(cd ${SRCDIR}; ${MAKE} install);

check:
	(cd ${TESTDIR}; ${MAKE});

dep: 
	(cd ${SRCDIR}; ${MAKE} dep);
