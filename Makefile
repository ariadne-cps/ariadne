include ./config.mk

.PHONY: lib python check clean

all:  ${SRCDIR}/Makefile 
	(cd ${SRCDIR}; ${MAKE} all);

lib: 
	(cd ${SRCDIR}; ${MAKE} lib);

python: 
	(cd ${SRCDIR}; ${MAKE} lib);
	(cd ${SRCDIR}; ${MAKE} python);

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
