include ../config.mk

GEOMETRYDIR=geometry/
BASEDIR = base/

.PHONY: lib base geometry python check clean

all: lib python

lib: libariadne.so

base: 
	(cd base/; ${MAKE} all)

geometry: 
	(cd geometry/; ${MAKE} all)

libariadne.so: base geometry
	${CXX} ${LDXXFLAGS} -shared -o libariadne.so base/*.o geometry/*.o ${LIBPPL} ${LIBGMPXX} ${RPATHFLAGS}

libariadne.a: ${OBJ}
	${CXX} ${LDXXFLAGS} -o libariadne.a ${OBJ} -lppl -lgmpxx -lgmp

python: lib
	(cd python/; ${MAKE} all)

clean:                                                                          
	rm -f *.bak *.o *~ *% core libariadne.*
	(cd ${BASEDIR}; ${MAKE} clean; cd ..; )
	(cd ${GEOMETRYDIR}; ${MAKE} clean; cd ..; )
	(cd ${PYTHONDIR}; ${MAKE} clean; cd ..; )


dep:  
	(cd ${BASEDIR}; ${MAKE} dep);
	(cd ${GEOMETRYDIR}; ${MAKE} dep);
	(cd ${PYTHONDIR}; ${MAKE} dep);
