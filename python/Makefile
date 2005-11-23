include ../config.mk

CXXFLAGS = -I../${INCLUDEDIR} -O2 -fpic -Wall
LDXXFLAGS += -lariadne -lboost_python -lsuperlu -lblas -lmpfr -lgmpxx -lgmp

ARIADNE_OBJ=ariadne_module.o export_numeric.o export_function.o 

GEOMETRY_OBJ = geometry_module.o \
	export_state.o export_interval.o export_rectangle.o export_polyhedron.o \
	export_list_set.o export_grid.o export_partition_tree.o

LINEAR_ALGEBRA_OBJ = linear_algebra_module.o \
	export_vector.o export_matrix.o
  
all: python

python: ariadne.so geometry.so linear_algebra.so

ariadne.so: ${ARIADNE_OBJ}
	${CXX} -shared ${ARIADNE_OBJ} ${LDXXFLAGS} -o ariadne.so

geometry.so: ${GEOMETRY_OBJ}
	${CXX} -shared ${GEOMETRY_OBJ} ${LDXXFLAGS} -o geometry.so

linear_algebra.so: ${LINEAR_ALGEBRA_OBJ}
	${CXX} -shared ${LINEAR_ALGEBRA_OBJ} ${LDXXFLAGS} -o linear_algebra.so

.cc.o: ariadne.cc function.cc geometry.cc grid.cc interval.cc list_set.cc partition_tree.cc rectangle.cc state.cc
	${CXX} ${CXXFLAGS} -I${PYTHONINCLUDEDIR} -c $<

clean:
	rm -f *.o *.so
