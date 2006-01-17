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

dep:  
	@sed -e '/^### Do Not edit this line$$/,$$d' Makefile \
		> MMM.$$$$ && mv MMM.$$$$ Makefile
	@echo "### Do Not edit this line" >> Makefile
	@echo "### Everything below is auto-generated" >> Makefile
	@for ff in ${ARIADNE_OBJ} ${GEOMETRY_OBJ} ${LINEAR_ALGEBRA_OBJ}; \
	do f=`basename $$ff .o`.cc ; \
	g++ -MM ${CXXFLAGS} -I${PYTHONINCLUDEDIR} $$f; echo ""; \
	done >> Makefile

### Do Not edit this line
### Everything below is auto-generated
ariadne_module.o: ariadne_module.cc ../include/numerical_type.h \
  ../include/interval.h ../include/state.h ../include/utility.h \
  ../include/linear_algebra.h ../include/numerical_type.h \
  ../include/rectangle.h ../include/ariadne.h ../include/utility.h \
  ../include/interval.h ../include/state.h ../include/binary_word.h \
  ../include/utility.h ../include/array.h ../include/utility.h \
  ../include/list_set.h

export_numeric.o: export_numeric.cc ../include/numerical_type.h

export_function.o: export_function.cc ../include/function.h \
  ../include/numerical_type.h real_typedef.h

geometry_module.o: geometry_module.cc ../include/numerical_type.h \
  ../include/interval.h ../include/state.h ../include/utility.h \
  ../include/linear_algebra.h ../include/numerical_type.h \
  ../include/rectangle.h ../include/ariadne.h ../include/utility.h \
  ../include/interval.h ../include/state.h ../include/binary_word.h \
  ../include/utility.h ../include/array.h ../include/utility.h \
  ../include/parallelopiped.h ../include/ariadne.h ../include/utility.h \
  ../include/linear_algebra.h ../include/interval.h ../include/state.h \
  ../include/rectangle.h ../include/list_set.h ../include/constraint.h \
  ../include/ariadne.h ../include/linear_algebra.h \
  ../include/polyhedron.h ../include/constraint.h ../include/generator.h \
  ../include/state.h ../include/rectangle.h ../include/linear_algebra.h \
  ../include/grid_set.h ../include/array.h ../include/interval.h \
  ../include/binary_word.h ../include/rectangle.h ../include/state.h \
  ../include/grid_operations.h ../include/array.h ../include/rectangle.h \
  ../include/list_set.h ../include/binary_tree.h ../include/array.h \
  ../include/sequence.h ../include/utility.h ../include/binary_word.h \
  ../include/rectangle.h ../include/partition_tree_set.h \
  ../include/array.h ../include/sequence.h ../include/utility.h \
  ../include/binary_word.h ../include/binary_tree.h \
  ../include/rectangle.h ../include/list_set.h real_typedef.h \
  python_utilities.h

export_state.o: export_state.cc ../include/numerical_type.h \
  ../include/state.h ../include/utility.h ../include/linear_algebra.h \
  ../include/numerical_type.h real_typedef.h python_utilities.h

export_interval.o: export_interval.cc ../include/numerical_type.h \
  ../include/interval.h real_typedef.h

export_rectangle.o: export_rectangle.cc ../include/numerical_type.h \
  ../include/rectangle.h ../include/ariadne.h ../include/utility.h \
  ../include/interval.h ../include/state.h ../include/utility.h \
  ../include/linear_algebra.h ../include/numerical_type.h \
  ../include/binary_word.h ../include/utility.h ../include/array.h \
  ../include/utility.h ../include/list_set.h real_typedef.h \
  python_utilities.h

export_parallelopiped.o: export_parallelopiped.cc \
  ../include/numerical_type.h ../include/rectangle.h ../include/ariadne.h \
  ../include/utility.h ../include/interval.h ../include/state.h \
  ../include/utility.h ../include/linear_algebra.h \
  ../include/numerical_type.h ../include/binary_word.h \
  ../include/utility.h ../include/array.h ../include/utility.h \
  ../include/parallelopiped.h ../include/ariadne.h ../include/utility.h \
  ../include/linear_algebra.h ../include/interval.h ../include/state.h \
  ../include/rectangle.h ../include/list_set.h ../include/constraint.h \
  ../include/ariadne.h ../include/linear_algebra.h \
  ../include/polyhedron.h ../include/constraint.h ../include/generator.h \
  ../include/state.h ../include/rectangle.h ../include/linear_algebra.h \
  real_typedef.h python_utilities.h

export_simplex.o: export_simplex.cc ../include/numerical_type.h \
  ../include/rectangle.h ../include/ariadne.h ../include/utility.h \
  ../include/interval.h ../include/state.h ../include/utility.h \
  ../include/linear_algebra.h ../include/numerical_type.h \
  ../include/binary_word.h ../include/utility.h ../include/array.h \
  ../include/utility.h ../include/simplex.h ../include/ariadne.h \
  ../include/array.h ../include/utility.h ../include/linear_algebra.h \
  ../include/interval.h ../include/state.h ../include/rectangle.h \
  ../include/polyhedron.h ../include/constraint.h ../include/ariadne.h \
  ../include/linear_algebra.h ../include/generator.h ../include/state.h \
  ../include/rectangle.h ../include/linear_algebra.h \
  ../include/list_set.h real_typedef.h python_utilities.h

export_polyhedron.o: export_polyhedron.cc ../include/numerical_type.h \
  ../include/rectangle.h ../include/ariadne.h ../include/utility.h \
  ../include/interval.h ../include/state.h ../include/utility.h \
  ../include/linear_algebra.h ../include/numerical_type.h \
  ../include/binary_word.h ../include/utility.h ../include/array.h \
  ../include/utility.h ../include/parallelopiped.h ../include/ariadne.h \
  ../include/utility.h ../include/linear_algebra.h ../include/interval.h \
  ../include/state.h ../include/rectangle.h ../include/list_set.h \
  ../include/constraint.h ../include/ariadne.h \
  ../include/linear_algebra.h ../include/polyhedron.h \
  ../include/constraint.h ../include/generator.h ../include/state.h \
  ../include/rectangle.h ../include/linear_algebra.h real_typedef.h \
  python_utilities.h

export_list_set.o: export_list_set.cc ../include/numerical_type.h \
  ../include/rectangle.h ../include/ariadne.h ../include/utility.h \
  ../include/interval.h ../include/state.h ../include/utility.h \
  ../include/linear_algebra.h ../include/numerical_type.h \
  ../include/binary_word.h ../include/utility.h ../include/array.h \
  ../include/utility.h ../include/list_set.h real_typedef.h \
  python_utilities.h

export_grid.o: export_grid.cc ../include/numerical_type.h \
  ../include/grid_set.h ../include/array.h ../include/utility.h \
  ../include/interval.h ../include/binary_word.h ../include/utility.h \
  ../include/array.h ../include/rectangle.h ../include/ariadne.h \
  ../include/utility.h ../include/interval.h ../include/state.h \
  ../include/utility.h ../include/linear_algebra.h \
  ../include/numerical_type.h ../include/binary_word.h ../include/state.h \
  ../include/grid_operations.h ../include/array.h ../include/rectangle.h \
  ../include/list_set.h real_typedef.h

export_partition_tree.o: export_partition_tree.cc \
  ../include/numerical_type.h ../include/rectangle.h ../include/ariadne.h \
  ../include/utility.h ../include/interval.h ../include/state.h \
  ../include/utility.h ../include/linear_algebra.h \
  ../include/numerical_type.h ../include/binary_word.h \
  ../include/utility.h ../include/array.h ../include/utility.h \
  ../include/binary_tree.h ../include/array.h ../include/sequence.h \
  ../include/utility.h ../include/binary_word.h ../include/rectangle.h \
  ../include/partition_tree_set.h ../include/array.h \
  ../include/sequence.h ../include/utility.h ../include/binary_word.h \
  ../include/binary_tree.h ../include/rectangle.h ../include/list_set.h \
  real_typedef.h

linear_algebra_module.o: linear_algebra_module.cc \
  ../include/linear_algebra.h ../include/numerical_type.h real_typedef.h

export_vector.o: export_vector.cc ../include/numerical_type.h \
  ../include/linear_algebra.h ../include/numerical_type.h real_typedef.h \
  python_utilities.h

export_matrix.o: export_matrix.cc ../include/numerical_type.h \
  ../include/linear_algebra.h ../include/numerical_type.h real_typedef.h \
  python_utilities.h

