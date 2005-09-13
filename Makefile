include ./config.mk

.PHONY: lib python check clean

all:  ${SRCDIR}/Makefile  ${EXAMPLEDIR}/Makefile ${PYTHONDIR}/Makefile
	(cd ${SRCDIR}; ${MAKE} all);
	(cd ${EXAMPLEDIR}; ${MAKE} all);
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

check:
	(cd ${TESTDIR}; ${MAKE});
