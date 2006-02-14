include ./config.mk

all:  python 

lib: 
	(cd ${SRCDIR}; ${MAKE} lib);

python: lib
	(cd ${SRCDIR}/${PYTHONDIR}; ${MAKE} all);

clean:                                                                          
	rm -f *.bak *.o *~ *% core
	(cd ${SRCDIR}; ${MAKE} clean; cd ..; )
	(cd ${EXAMPLEDIR}; ${MAKE} clean; cd ..; )
	(cd ${TESTDIR}; ${MAKE} clean; cd ..; )

install:
#	cp src/libariadne.so ${PREFIX}/lib/libariadne.so
	ln -s ${PWD}/src/libariadne.so ${PREFIX}/lib/libariadne.so 

check:
	(cd ${TESTDIR}; ${MAKE});

dep: 
	(cd ${SRCDIR}; ${MAKE} dep);
