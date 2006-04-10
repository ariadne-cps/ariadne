include ../config.mk

BASEDIR=base/
NUMERICDIR=numeric/
LINEARALGEBRADIR=linear_algebra/
GEOMETRYDIR=geometry/
EVALUATIONDIR=evaluation/

.PHONY: lib base linear_algebra geometry evaluation python check clean

all: lib python

lib: libariadne.so

base: 
	(cd base; ${MAKE} all)

numeric: 
	(cd numeric; ${MAKE} all)

linear_algebra: numeric
	(cd linear_algebra; ${MAKE} all)

geometry: base numeric linear_algebra
	(cd geometry; ${MAKE} all)

evaluation: base numeric linear_algebra geometry 
	(cd evaluation; ${MAKE} all)

libariadne.so: base numeric linear_algebra geometry evaluation
	${CXX} ${LDXXFLAGS} -shared -o libariadne.so base/*.o linear_algebra/*.o geometry/*.o evaluation/*.o ${LIBPPL} ${LIBGMPXX} ${RPATHFLAGS}

libariadne.a: ${OBJ}
	${CXX} ${LDXXFLAGS} -o libariadne.a ${OBJ} -lppl -lgmpxx -lgmp

python: lib
	(cd ${PYTHONDIR}; ${MAKE} all)

install: lib python
	cp -f libariadne.so ${PREFIX}/lib/libariadne.so
	(cd ${PYTHONDIR}; ${MAKE} install)


clean:                                                                          
	rm -f *.bak *.o *~ *% core libariadne.*
	(cd ${BASEDIR}; ${MAKE} clean )
	(cd ${LINEARALGEBRADIR}; ${MAKE} clean )
	(cd ${GEOMETRYDIR}; ${MAKE} clean )
	(cd ${EVALUATIONDIR}; ${MAKE} clean )
	(cd ${PYTHONDIR}; ${MAKE} clean )


dep:  
	(cd ${BASEDIR}; ${MAKE} dep);
	(cd ${LINEARALGEBRADIR}; ${MAKE} dep);
	(cd ${GEOMETRYDIR}; ${MAKE} dep);
	(cd ${EVALUATIONDIR}; ${MAKE} dep; )
	(cd ${PYTHONDIR}; ${MAKE} dep);
