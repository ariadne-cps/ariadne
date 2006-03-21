include ../config.mk

GEOMETRYDIR=geometry/
BASEDIR = base/
LINEARALGEBRADIR=linear_algebra/

.PHONY: lib base geometry linear_algebra python check clean

all: lib python

lib: libariadne.so

base: 
	(cd base; ${MAKE} all)

linear_algebra: base 
	(cd linear_algebra; ${MAKE} all)

geometry: base 
	(cd geometry; ${MAKE} all)

libariadne.so: base geometry linear_algebra
	${CXX} ${LDXXFLAGS} -shared -o libariadne.so base/*.o linear_algebra/*.o geometry/*.o ${LIBPPL} ${LIBGMPXX} ${RPATHFLAGS}

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
	(cd ${PYTHONDIR}; ${MAKE} clean )


dep:  
	(cd ${BASEDIR}; ${MAKE} dep);
	(cd ${LINEARALGEBRADIR}; ${MAKE} dep);
	(cd ${GEOMETRYDIR}; ${MAKE} dep);
	(cd ${PYTHONDIR}; ${MAKE} dep);
