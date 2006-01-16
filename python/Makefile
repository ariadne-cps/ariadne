include ../config.mk

CXXFLAGS = -I../${INCLUDEDIR} -O2 -fpic -Wall
LDXXFLAGS += -L../src
LDLIBS = -lariadne -lboost_python ${LIBSUPERLU} -lblas -lmpfr ${LIBGMPXX}

ARIADNE_OBJ=ariadne_module.o export_numeric.o export_function.o 

GEOMETRY_OBJ = geometry_module.o \
	export_state.o export_interval.o \
	export_rectangle.o export_parallelopiped.o export_simplex.o export_polyhedron.o \
	export_list_set.o export_grid.o export_partition_tree.o

LINEAR_ALGEBRA_OBJ = linear_algebra_module.o \
	export_vector.o export_matrix.o
  
all: python

python: ariadne.so geometry.so linear_algebra.so

ariadne.so: ${ARIADNE_OBJ}
	${CXX} -shared -o ariadne.so ${LDXXFLAGS} ${ARIADNE_OBJ} ${LDLIBS} ${RPATHFLAGS}

geometry.so: ${GEOMETRY_OBJ}
	${CXX} -shared -o geometry.so ${LDXXFLAGS} ${GEOMETRY_OBJ} ${LDLIBS} ${RPATHFLAGS}

linear_algebra.so: ${LINEAR_ALGEBRA_OBJ}
	${CXX} -shared -o linear_algebra.so ${LDXXFLAGS} ${LINEAR_ALGEBRA_OBJ} ${LDLIBS} ${RPATHFLAGS}

.cc.o: ariadne.cc function.cc geometry.cc grid.cc interval.cc list_set.cc partition_tree.cc rectangle.cc state.cc
	${CXX} ${CXXFLAGS} -I${PYTHONINCLUDEDIR} -c $<

clean:
	rm -f *.o *.so
