include ./config.mk

.PHONY: lib python check clean

all:  ${SRCDIR}/Makefile  ${PYTHONDIR}/Makefile
	(cd ${SRCDIR}; ${MAKE} all);
	(cd ${PYTHONDIR}; ${MAKE} all);

lib: 
	(cd ${SRCDIR}; ${MAKE} all);

python: 
	(cd ${SRCDIR}; ${MAKE} all);
	(cd ${PYTHONDIR}; ${MAKE} all);

clean:                                                                          
	rm -f *.bak *.o *~ *% core
	(cd ${SRCDIR}; ${MAKE} clean; cd ..; )
	(cd ${EXAMPLEDIR}; ${MAKE} clean; cd ..; )
	(cd ${TESTDIR}; ${MAKE} clean; cd ..; )
	(cd ${PYTHONDIR}; ${MAKE} clean; cd ..; )

install:
#	cp src/libariadne.so ${PREFIX}/lib/libariadne.so
	ln -s ${PWD}/src/libariadne.so ${PREFIX}/lib/libariadne.so 
   
check:
	(cd ${TESTDIR}; ${MAKE});

dep: 
	(cd ${SRCDIR}; ${MAKE} dep);
	(cd ${PYTHONDIR}; ${MAKE} dep);
 
