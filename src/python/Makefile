include ../../config.mk

CXXFLAGS += -O2 -fpic -Wall -I/usr/include/python2.4/
LDXXFLAGS += -L../src
LDLIBS = -lariadne -lboost_python ${LIBGMPXX} -lppl

BASE_OBJ=base_module.o \
	export_numeric.o export_interval.o export_function.o 

GEOMETRY_OBJ = geometry_module.o \
	export_point.o \
	export_rectangle.o export_parallelopiped.o export_simplex.o export_polyhedron.o \
	export_list_set.o export_grid.o export_partition_tree.o \
	export_postscript_output.o \

EVALUATION_OBJ = evaluation_module.o \
	export_map.o export_affine_map.o export_polynomial_map.o \
	export_apply.o \
	export_henon_map.o

LINEAR_ALGEBRA_OBJ = linear_algebra_module.o \
	export_vector.o export_matrix.o

all: python

python: base.so geometry.so evaluation.so linear_algebra.so

base.so: ${BASE_OBJ}
	${CXX} -shared -o base.so ${LDXXFLAGS} ${BASE_OBJ} ${LDLIBS} ${RPATHFLAGS}

geometry.so: ${GEOMETRY_OBJ}
	${CXX} -shared -o geometry.so ${LDXXFLAGS} ${GEOMETRY_OBJ} ${LDLIBS} ${RPATHFLAGS}

evaluation.so: ${EVALUATION_OBJ}
	${CXX} -shared -o evaluation.so ${LDXXFLAGS} ${EVALUATION_OBJ} ${LDLIBS} ${RPATHFLAGS}

linear_algebra.so: ${LINEAR_ALGEBRA_OBJ}
	${CXX} -shared -o linear_algebra.so ${LDXXFLAGS} ${LINEAR_ALGEBRA_OBJ} ${LDLIBS} ${RPATHFLAGS}

install: python
	mkdir ${PYTHONINSTALLDIR}/ariadne
	cp -f *.so __init__.py ${PYTHONINSTALLDIR}/ariadne/

clean:
	rm -f *.o *.so ariadne/*.so

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
geometry_module.o: geometry_module.cc ../../include/base/binary_word.h \
  ../../include/base/utility.h ../../include/base/array.h \
  ../../include/base/basic_type.h ../../include/base/sequence.h \
  ../../include/base/binary_tree.h ../../include/geometry/rectangle.h \
  ../../include/ariadne.h ../../include/base/interval.h \
  ../../include/geometry/point.h ../../include/linear_algebra/vector.h \
  ../../include/base/numerical_type.h ../../include/base/dyadic.h \
  ../../include/geometry/geometry_declarations.h \
  ../../include/geometry/partition_tree_set.h \
  ../../include/geometry/list_set.h ../../include/python/real_typedef.h \
  ../../include/python/python_utilities.h

export_point.o: export_point.cc ../../include/base/numerical_type.h \
  ../../include/base/dyadic.h ../../include/geometry/point.h \
  ../../include/base/utility.h ../../include/base/array.h \
  ../../include/linear_algebra/vector.h ../../include/base/basic_type.h \
  ../../include/base/sequence.h ../../include/base/interval.h \
  ../../include/geometry/geometry_declarations.h \
  ../../include/python/real_typedef.h

export_rectangle.o: export_rectangle.cc \
  ../../include/base/numerical_type.h ../../include/base/dyadic.h \
  ../../include/geometry/rectangle.h ../../include/ariadne.h \
  ../../include/base/utility.h ../../include/base/array.h \
  ../../include/base/interval.h ../../include/base/binary_word.h \
  ../../include/base/basic_type.h ../../include/base/sequence.h \
  ../../include/geometry/point.h ../../include/linear_algebra/vector.h \
  ../../include/geometry/geometry_declarations.h \
  ../../include/geometry/list_set.h ../../include/python/real_typedef.h

export_parallelopiped.o: export_parallelopiped.cc \
  ../../include/base/numerical_type.h ../../include/base/dyadic.h \
  ../../include/geometry/rectangle.h ../../include/ariadne.h \
  ../../include/base/utility.h ../../include/base/array.h \
  ../../include/base/interval.h ../../include/base/binary_word.h \
  ../../include/base/basic_type.h ../../include/base/sequence.h \
  ../../include/geometry/point.h ../../include/linear_algebra/vector.h \
  ../../include/geometry/geometry_declarations.h \
  ../../include/geometry/parallelopiped.h \
  ../../include/linear_algebra/linear_algebra.h \
  ../../include/linear_algebra/matrix.h \
  ../../include/linear_algebra/constraint.h \
  ../../include/geometry/list_set.h \
  ../../include/geometry/grid_operations.h \
  ../../include/geometry/polyhedron.h ../../include/geometry/generator.h \
  ../../include/python/real_typedef.h \
  ../../include/python/python_utilities.h

export_simplex.o: export_simplex.cc ../../include/base/numerical_type.h \
  ../../include/base/dyadic.h ../../include/geometry/rectangle.h \
  ../../include/ariadne.h ../../include/base/utility.h \
  ../../include/base/array.h ../../include/base/interval.h \
  ../../include/base/binary_word.h ../../include/base/basic_type.h \
  ../../include/base/sequence.h ../../include/geometry/point.h \
  ../../include/linear_algebra/vector.h \
  ../../include/geometry/geometry_declarations.h \
  ../../include/geometry/simplex.h ../../include/geometry/polyhedron.h \
  ../../include/linear_algebra/constraint.h \
  ../../include/linear_algebra/matrix.h \
  ../../include/geometry/generator.h ../../include/geometry/list_set.h \
  ../../include/python/real_typedef.h

export_polyhedron.o: export_polyhedron.cc \
  ../../include/base/numerical_type.h ../../include/base/dyadic.h \
  ../../include/geometry/rectangle.h ../../include/ariadne.h \
  ../../include/base/utility.h ../../include/base/array.h \
  ../../include/base/interval.h ../../include/base/binary_word.h \
  ../../include/base/basic_type.h ../../include/base/sequence.h \
  ../../include/geometry/point.h ../../include/linear_algebra/vector.h \
  ../../include/geometry/geometry_declarations.h \
  ../../include/geometry/parallelopiped.h \
  ../../include/linear_algebra/linear_algebra.h \
  ../../include/linear_algebra/matrix.h \
  ../../include/linear_algebra/constraint.h \
  ../../include/geometry/list_set.h \
  ../../include/geometry/grid_operations.h \
  ../../include/geometry/polyhedron.h ../../include/geometry/generator.h \
  ../../include/python/real_typedef.h

export_list_set.o: export_list_set.cc ../../include/base/numerical_type.h \
  ../../include/base/dyadic.h ../../include/geometry/rectangle.h \
  ../../include/ariadne.h ../../include/base/utility.h \
  ../../include/base/array.h ../../include/base/interval.h \
  ../../include/base/binary_word.h ../../include/base/basic_type.h \
  ../../include/base/sequence.h ../../include/geometry/point.h \
  ../../include/linear_algebra/vector.h \
  ../../include/geometry/geometry_declarations.h \
  ../../include/geometry/parallelopiped.h \
  ../../include/linear_algebra/linear_algebra.h \
  ../../include/linear_algebra/matrix.h \
  ../../include/linear_algebra/constraint.h \
  ../../include/geometry/list_set.h \
  ../../include/geometry/grid_operations.h \
  ../../include/geometry/polyhedron.h ../../include/geometry/generator.h \
  ../../include/geometry/grid_set.h ../../include/python/real_typedef.h \
  ../../include/python/python_utilities.h

export_grid.o: export_grid.cc ../../include/base/numerical_type.h \
  ../../include/base/dyadic.h \
  ../../include/linear_algebra/linear_algebra.h \
  ../../include/linear_algebra/vector.h ../../include/base/basic_type.h \
  ../../include/base/array.h ../../include/base/sequence.h \
  ../../include/base/interval.h ../../include/linear_algebra/matrix.h \
  ../../include/linear_algebra/constraint.h \
  ../../include/geometry/grid_set.h ../../include/base/binary_word.h \
  ../../include/base/utility.h \
  ../../include/geometry/geometry_declarations.h \
  ../../include/geometry/rectangle.h ../../include/ariadne.h \
  ../../include/geometry/point.h ../../include/geometry/grid_operations.h \
  ../../include/geometry/list_set.h \
  ../../include/geometry/parallelopiped.h \
  ../../include/geometry/polyhedron.h ../../include/geometry/generator.h \
  ../../include/python/real_typedef.h

export_partition_tree.o: export_partition_tree.cc \
  ../../include/base/numerical_type.h ../../include/base/dyadic.h \
  ../../include/base/binary_word.h ../../include/base/utility.h \
  ../../include/base/array.h ../../include/base/basic_type.h \
  ../../include/base/sequence.h ../../include/base/binary_tree.h \
  ../../include/geometry/rectangle.h ../../include/ariadne.h \
  ../../include/base/interval.h ../../include/geometry/point.h \
  ../../include/linear_algebra/vector.h \
  ../../include/geometry/geometry_declarations.h \
  ../../include/geometry/partition_tree_set.h \
  ../../include/geometry/list_set.h ../../include/python/real_typedef.h

export_postscript_output.o: export_postscript_output.cc \
  ../../include/geometry/rectangle.h ../../include/ariadne.h \
  ../../include/base/utility.h ../../include/base/array.h \
  ../../include/base/interval.h ../../include/base/binary_word.h \
  ../../include/base/basic_type.h ../../include/base/sequence.h \
  ../../include/geometry/point.h ../../include/linear_algebra/vector.h \
  ../../include/base/numerical_type.h ../../include/base/dyadic.h \
  ../../include/geometry/geometry_declarations.h \
  ../../include/geometry/parallelopiped.h \
  ../../include/linear_algebra/linear_algebra.h \
  ../../include/linear_algebra/matrix.h \
  ../../include/linear_algebra/constraint.h \
  ../../include/geometry/list_set.h \
  ../../include/geometry/grid_operations.h \
  ../../include/geometry/polyhedron.h ../../include/geometry/generator.h \
  ../../include/geometry/grid_set.h \
  ../../include/geometry/partition_tree_set.h \
  ../../include/base/binary_tree.h ../../include/utility/epsfstream.h \
  ../../include/python/real_typedef.h

linear_algebra_module.o: linear_algebra_module.cc

export_vector.o: export_vector.cc ../../include/base/numerical_type.h \
  ../../include/base/dyadic.h ../../include/base/interval.h \
  ../../include/linear_algebra/linear_algebra.h \
  ../../include/linear_algebra/vector.h ../../include/base/basic_type.h \
  ../../include/base/array.h ../../include/base/sequence.h \
  ../../include/linear_algebra/matrix.h \
  ../../include/linear_algebra/constraint.h \
  ../../include/python/real_typedef.h

export_matrix.o: export_matrix.cc ../../include/base/numerical_type.h \
  ../../include/base/dyadic.h ../../include/base/interval.h \
  ../../include/linear_algebra/matrix.h ../../include/base/basic_type.h \
  ../../include/base/array.h ../../include/base/sequence.h \
  ../../include/python/real_typedef.h \
  ../../include/python/python_utilities.h

