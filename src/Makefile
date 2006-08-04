include ../config.mk

BASEDIR=base/
NUMERICDIR=numeric/
LINEARALGEBRADIR=linear_algebra/
GEOMETRYDIR=geometry/
SYSTEMDIR=system/
EVALUATIONDIR=evaluation/

.PHONY: lib base numeric linear_algebra geometry system evaluation python check clean

all: lib python

lib: libariadne.so

base: 
	(cd base; ${MAKE} all)

numeric: 
	(cd numeric; ${MAKE} all)

linear_algebra:
	(cd linear_algebra; ${MAKE} all)

geometry:
	(cd geometry; ${MAKE} all)

system:
	(cd system; ${MAKE} all)

evaluation:
	(cd evaluation; ${MAKE} all)

libariadne.so: base numeric linear_algebra geometry system evaluation
	${CXX} ${LDXXFLAGS} -shared -o libariadne.so base/*.o numeric/*.o linear_algebra/*.o geometry/*.o system/*.o evaluation/*.o ${LIBPPL} ${LIBGMPXX} ${RPATHFLAGS}

libariadne.a: ${OBJ}
	${CXX} ${LDXXFLAGS} -o libariadne.a ${OBJ} -lppl -lgmpxx -lgmp

python: lib
	(cd ${PYTHONDIR}; ${MAKE} all)

install: lib python
	mkdir -p ${LIBPREFIX}
	cp -f libariadne.so ${LIBPREFIX}/libariadne.so
	(cd ${PYTHONDIR}; ${MAKE} install)

install-lib: lib
	mkdir -p ${LIBPREFIX}
	cp -f libariadne.so ${LIBPREFIX}/libariadne.so


clean:                                                                          
	rm -f *.bak *.o *~ *% core libariadne.*
	(cd ${BASEDIR}; ${MAKE} clean )
	(cd ${NUMERICDIR}; ${MAKE} clean )
	(cd ${LINEARALGEBRADIR}; ${MAKE} clean )
	(cd ${GEOMETRYDIR}; ${MAKE} clean )
	(cd ${SYSTEMDIR}; ${MAKE} clean )
	(cd ${EVALUATIONDIR}; ${MAKE} clean )
	(cd ${PYTHONDIR}; ${MAKE} clean )


dep:  
	(cd ${BASEDIR}; ${MAKE} dep )
	(cd ${NUMERICDIR}; ${MAKE} dep )
	(cd ${LINEARALGEBRADIR}; ${MAKE} dep )
	(cd ${GEOMETRYDIR}; ${MAKE} dep )
	(cd ${SYSTEMDIR}; ${MAKE} dep )
	(cd ${EVALUATIONDIR}; ${MAKE} dep )
	(cd ${PYTHONDIR}; ${MAKE} dep )
