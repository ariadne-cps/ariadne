include ./config.mk

all:  ${SRCDIR}/Makefile
	(cd ${SRCDIR}; ${MAKE} all);

clean:                                                                          
	rm -f *.bak *.o *~ *% core
	(cd ${SRCDIR}; ${MAKE} clean);

dep: ${SRCDIR}/Makefile
	(cd ${SRCDIR}; ${MAKE} dep);
