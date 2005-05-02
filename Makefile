include ./config.mk

all:  ${EXAMPLEDIR}/Makefile
	(cd ${EXAMPLEDIR}; ${MAKE} all);

clean:                                                                          
	rm -f *.bak *.o *~ *% core
	(cd ${EXAMPLEDIR}; ${MAKE} clean; cd ..; cd ${TESTDIR}; ${MAKE} clean);

check:
	(cd ${TESTDIR}; ${MAKE});
