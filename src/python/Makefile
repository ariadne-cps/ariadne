include ../../config.mk

CXXFLAGS += -O2 -fpic -Wall -I/usr/include/python2.4/
LDXXFLAGS += -L../src
LDLIBS = -lariadne -lboost_python ${LIBGMPXX} -lppl

BASE_OBJ=base_module.o \
	export_numeric.o export_interval.o export_function.o \
	export_array.o export_binary_tree.o 

LINEAR_ALGEBRA_OBJ = linear_algebra_module.o \
	export_vector.o export_matrix.o export_linear_program.o

GEOMETRY_OBJ = geometry_module.o \
	export_point.o export_rectangle.o export_parallelotope.o \
	export_zonotope.o export_simplex.o export_polyhedron.o \
	export_list_set.o export_grid_set.o export_partition_tree_set.o \
	export_postscript_output.o \

EVALUATION_OBJ = evaluation_module.o \
	export_map.o export_affine_map.o export_polynomial_map.o \
	export_vector_field.o export_affine_vector_field.o \
	export_apply.o \
	export_integrate.o \
	export_henon_map.o export_lorenz_system.o

all: python

python: base.so linear_algebra.so geometry.so evaluation.so 

base.so: ${BASE_OBJ}
	${CXX} -shared -o base.so ${LDXXFLAGS} ${BASE_OBJ} ${LDLIBS} ${RPATHFLAGS}

linear_algebra.so: ${LINEAR_ALGEBRA_OBJ}
	${CXX} -shared -o linear_algebra.so ${LDXXFLAGS} ${LINEAR_ALGEBRA_OBJ} ${LDLIBS} ${RPATHFLAGS}

geometry.so: ${GEOMETRY_OBJ}
	${CXX} -shared -o geometry.so ${LDXXFLAGS} ${GEOMETRY_OBJ} ${LDLIBS} ${RPATHFLAGS}

evaluation.so: ${EVALUATION_OBJ}
	${CXX} -shared -o evaluation.so ${LDXXFLAGS} ${EVALUATION_OBJ} ${LDLIBS} ${RPATHFLAGS}


install: python
	mkdir -p ${PYTHONINSTALLDIR}/ariadne
	cp -f *.so __init__.py ${PYTHONINSTALLDIR}/ariadne/

clean:
	rm -f *.o *.so ariadne/*.so

dep:  
	@sed -e '/^### Do Not edit this line$$/,$$d' Makefile \
		> MMM.$$$$ && mv MMM.$$$$ Makefile
	@echo "### Do Not edit this line" >> Makefile
	@echo "### Everything below is auto-generated" >> Makefile
	@for ff in ${BASE_OBJ} ${GEOMETRY_OBJ} ${EVALUATION_OBJ} ${LINEAR_ALGEBRA_OBJ}; \
	do f=`basename $$ff .o`.cc ; \
	g++ -MM ${CXXFLAGS} -I${PYTHONINCLUDEDIR} $$f; echo ""; \
	done >> Makefile

### Do Not edit this line
### Everything below is auto-generated
base_module.o: base_module.cc

export_numeric.o: export_numeric.cc ../../include/base/numerical_type.h \
  ../../include/base/../base/dyadic.h ../../include/base/arithmetic.h \
  ../../include/base/../base/numerical_type.h \
  ../../include/base/../base/approximation.h \
  ../../include/base/../base/../base/numerical_type.h

export_interval.o: export_interval.cc ../../include/base/numerical_type.h \
  ../../include/base/../base/dyadic.h ../../include/base/interval.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../base/numerical_type.h \
  ../../include/python/../declarations.h \
  ../../include/python/../base/basic_type.h

export_function.o: export_function.cc ../../include/base/function.h \
  /ufs/collins/include/mpfr.h ../../include/base/../base/numerical_type.h \
  ../../include/base/../base/../base/dyadic.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../base/numerical_type.h \
  ../../include/python/../declarations.h \
  ../../include/python/../base/basic_type.h

export_array.o: export_array.cc ../../include/base/array.h \
  ../../include/base/basic_type.h ../../include/base/array_operations.h \
  ../../include/base/../base/array.h \
  ../../include/base/../base/basic_type.h ../../include/base/utility.h \
  ../../include/base/basic_type.h ../../include/python/typedefs.h \
  ../../include/python/../real_typedef.h \
  ../../include/python/../base/numerical_type.h \
  ../../include/python/../base/../base/dyadic.h \
  ../../include/python/../declarations.h \
  ../../include/python/../base/basic_type.h \
  ../../include/python/python_utilities.h

export_binary_tree.o: export_binary_tree.cc \
  ../../include/base/binary_word.h \
  ../../include/base/../base/basic_type.h \
  ../../include/base/../base/array.h \
  ../../include/base/../base/basic_type.h \
  ../../include/base/binary_tree.h \
  ../../include/base/../base/binary_word.h \
  ../../include/base/basic_type.h ../../include/python/typedefs.h \
  ../../include/python/../real_typedef.h \
  ../../include/python/../base/numerical_type.h \
  ../../include/python/../base/../base/dyadic.h \
  ../../include/python/../declarations.h \
  ../../include/python/../base/basic_type.h \
  ../../include/python/python_utilities.h

geometry_module.o: geometry_module.cc

export_point.o: export_point.cc ../../include/base/numerical_type.h \
  ../../include/base/../base/dyadic.h ../../include/geometry/point.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../base/utility.h \
  ../../include/geometry/../base/../base/array.h \
  ../../include/geometry/../base/../base/basic_type.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/../base/basic_type.h \
  ../../include/geometry/../linear_algebra/../base/numerical_type.h \
  ../../include/linear_algebra/vector.h ../../include/python/typedefs.h \
  ../../include/python/../real_typedef.h \
  ../../include/python/../base/numerical_type.h \
  ../../include/python/../declarations.h \
  ../../include/python/../base/basic_type.h

export_rectangle.o: export_rectangle.cc \
  ../../include/base/numerical_type.h ../../include/base/../base/dyadic.h \
  ../../include/geometry/rectangle.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../base/utility.h \
  ../../include/geometry/../base/../base/array.h \
  ../../include/geometry/../base/../base/basic_type.h \
  ../../include/geometry/../base/interval.h \
  ../../include/geometry/../base/binary_word.h \
  ../../include/geometry/../base/../base/basic_type.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../base/utility.h \
  ../../include/geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../base/basic_type.h \
  ../../include/geometry/../geometry/../linear_algebra/../base/numerical_type.h \
  ../../include/geometry/list_set.h ../../include/python/typedefs.h \
  ../../include/python/../real_typedef.h \
  ../../include/python/../base/numerical_type.h \
  ../../include/python/../declarations.h \
  ../../include/python/../base/basic_type.h

export_parallelotope.o: export_parallelotope.cc \
  ../../include/base/numerical_type.h ../../include/base/../base/dyadic.h \
  ../../include/geometry/rectangle.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../base/utility.h \
  ../../include/geometry/../base/../base/array.h \
  ../../include/geometry/../base/../base/basic_type.h \
  ../../include/geometry/../base/interval.h \
  ../../include/geometry/../base/binary_word.h \
  ../../include/geometry/../base/../base/basic_type.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../base/utility.h \
  ../../include/geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../base/basic_type.h \
  ../../include/geometry/../geometry/../linear_algebra/../base/numerical_type.h \
  ../../include/geometry/list_set.h \
  ../../include/geometry/parallelotope.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../linear_algebra/../base/basic_type.h \
  ../../include/geometry/../linear_algebra/../base/numerical_type.h \
  ../../include/geometry/../linear_algebra/../base/arithmetic.h \
  ../../include/geometry/../linear_algebra/../base/../base/numerical_type.h \
  ../../include/geometry/../linear_algebra/../base/../base/approximation.h \
  ../../include/geometry/../linear_algebra/../base/../base/../base/numerical_type.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/list_set.h \
  ../../include/geometry/../geometry/polyhedron.h \
  /ufs/collins/include/ppl.hh \
  ../../include/geometry/../geometry/../linear_algebra/constraint.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../geometry/generator.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/zonotope.h \
  ../../include/geometry/../geometry/../linear_algebra/zonotopic_vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../base/utility.h \
  ../../include/geometry/../geometry/../linear_algebra/../base/interval.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/interval_vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/../base/basic_type.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/../base/numerical_type.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/../base/interval.h \
  ../../include/geometry/../geometry/../base/interval.h \
  ../../include/geometry/../geometry/../geometry/list_set.h \
  ../../include/geometry/../geometry/../geometry/parallelotope.h \
  ../../include/geometry/../geometry/../geometry/polyhedron.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../base/numerical_type.h \
  ../../include/python/../declarations.h \
  ../../include/python/../base/basic_type.h \
  ../../include/python/python_utilities.h

export_zonotope.o: export_zonotope.cc ../../include/base/numerical_type.h \
  ../../include/base/../base/dyadic.h ../../include/geometry/rectangle.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../base/utility.h \
  ../../include/geometry/../base/../base/array.h \
  ../../include/geometry/../base/../base/basic_type.h \
  ../../include/geometry/../base/interval.h \
  ../../include/geometry/../base/binary_word.h \
  ../../include/geometry/../base/../base/basic_type.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../base/utility.h \
  ../../include/geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../base/basic_type.h \
  ../../include/geometry/../geometry/../linear_algebra/../base/numerical_type.h \
  ../../include/geometry/list_set.h \
  ../../include/geometry/parallelotope.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../linear_algebra/../base/basic_type.h \
  ../../include/geometry/../linear_algebra/../base/numerical_type.h \
  ../../include/geometry/../linear_algebra/../base/arithmetic.h \
  ../../include/geometry/../linear_algebra/../base/../base/numerical_type.h \
  ../../include/geometry/../linear_algebra/../base/../base/approximation.h \
  ../../include/geometry/../linear_algebra/../base/../base/../base/numerical_type.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/list_set.h \
  ../../include/geometry/../geometry/polyhedron.h \
  /ufs/collins/include/ppl.hh \
  ../../include/geometry/../geometry/../linear_algebra/constraint.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../geometry/generator.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/zonotope.h \
  ../../include/geometry/../geometry/../linear_algebra/zonotopic_vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../base/utility.h \
  ../../include/geometry/../geometry/../linear_algebra/../base/interval.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/interval_vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/../base/basic_type.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/../base/numerical_type.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/../base/interval.h \
  ../../include/geometry/../geometry/../base/interval.h \
  ../../include/geometry/../geometry/../geometry/list_set.h \
  ../../include/geometry/../geometry/../geometry/parallelotope.h \
  ../../include/geometry/../geometry/../geometry/polyhedron.h \
  ../../include/geometry/zonotope.h ../../include/python/typedefs.h \
  ../../include/python/../real_typedef.h \
  ../../include/python/../base/numerical_type.h \
  ../../include/python/../declarations.h \
  ../../include/python/../base/basic_type.h \
  ../../include/python/python_utilities.h

export_simplex.o: export_simplex.cc ../../include/base/numerical_type.h \
  ../../include/base/../base/dyadic.h ../../include/geometry/rectangle.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../base/utility.h \
  ../../include/geometry/../base/../base/array.h \
  ../../include/geometry/../base/../base/basic_type.h \
  ../../include/geometry/../base/interval.h \
  ../../include/geometry/../base/binary_word.h \
  ../../include/geometry/../base/../base/basic_type.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../base/utility.h \
  ../../include/geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../base/basic_type.h \
  ../../include/geometry/../geometry/../linear_algebra/../base/numerical_type.h \
  ../../include/geometry/simplex.h ../../include/geometry/../base/array.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/polyhedron.h \
  /ufs/collins/include/ppl.hh \
  ../../include/geometry/../geometry/../linear_algebra/constraint.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/../base/basic_type.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/../base/numerical_type.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/../base/arithmetic.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/../base/../base/numerical_type.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/../base/../base/approximation.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/../base/../base/../base/numerical_type.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../geometry/generator.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/list_set.h \
  ../../include/geometry/list_set.h ../../include/python/typedefs.h \
  ../../include/python/../real_typedef.h \
  ../../include/python/../base/numerical_type.h \
  ../../include/python/../declarations.h \
  ../../include/python/../base/basic_type.h

export_polyhedron.o: export_polyhedron.cc \
  ../../include/base/numerical_type.h ../../include/base/../base/dyadic.h \
  ../../include/geometry/rectangle.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../base/utility.h \
  ../../include/geometry/../base/../base/array.h \
  ../../include/geometry/../base/../base/basic_type.h \
  ../../include/geometry/../base/interval.h \
  ../../include/geometry/../base/binary_word.h \
  ../../include/geometry/../base/../base/basic_type.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../base/utility.h \
  ../../include/geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../base/basic_type.h \
  ../../include/geometry/../geometry/../linear_algebra/../base/numerical_type.h \
  ../../include/geometry/parallelotope.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../linear_algebra/../base/basic_type.h \
  ../../include/geometry/../linear_algebra/../base/numerical_type.h \
  ../../include/geometry/../linear_algebra/../base/arithmetic.h \
  ../../include/geometry/../linear_algebra/../base/../base/numerical_type.h \
  ../../include/geometry/../linear_algebra/../base/../base/approximation.h \
  ../../include/geometry/../linear_algebra/../base/../base/../base/numerical_type.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/list_set.h \
  ../../include/geometry/../geometry/polyhedron.h \
  /ufs/collins/include/ppl.hh \
  ../../include/geometry/../geometry/../linear_algebra/constraint.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../geometry/generator.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/zonotope.h \
  ../../include/geometry/../geometry/../linear_algebra/zonotopic_vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../base/utility.h \
  ../../include/geometry/../geometry/../linear_algebra/../base/interval.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/interval_vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/../base/basic_type.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/../base/numerical_type.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/../base/interval.h \
  ../../include/geometry/../geometry/../base/interval.h \
  ../../include/geometry/../geometry/../geometry/list_set.h \
  ../../include/geometry/../geometry/../geometry/parallelotope.h \
  ../../include/geometry/../geometry/../geometry/polyhedron.h \
  ../../include/geometry/polyhedron.h ../../include/python/typedefs.h \
  ../../include/python/../real_typedef.h \
  ../../include/python/../base/numerical_type.h \
  ../../include/python/../declarations.h \
  ../../include/python/../base/basic_type.h

export_list_set.o: export_list_set.cc ../../include/base/numerical_type.h \
  ../../include/base/../base/dyadic.h ../../include/geometry/rectangle.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../base/utility.h \
  ../../include/geometry/../base/../base/array.h \
  ../../include/geometry/../base/../base/basic_type.h \
  ../../include/geometry/../base/interval.h \
  ../../include/geometry/../base/binary_word.h \
  ../../include/geometry/../base/../base/basic_type.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../base/utility.h \
  ../../include/geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../base/basic_type.h \
  ../../include/geometry/../geometry/../linear_algebra/../base/numerical_type.h \
  ../../include/geometry/parallelotope.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../linear_algebra/../base/basic_type.h \
  ../../include/geometry/../linear_algebra/../base/numerical_type.h \
  ../../include/geometry/../linear_algebra/../base/arithmetic.h \
  ../../include/geometry/../linear_algebra/../base/../base/numerical_type.h \
  ../../include/geometry/../linear_algebra/../base/../base/approximation.h \
  ../../include/geometry/../linear_algebra/../base/../base/../base/numerical_type.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/list_set.h \
  ../../include/geometry/../geometry/polyhedron.h \
  /ufs/collins/include/ppl.hh \
  ../../include/geometry/../geometry/../linear_algebra/constraint.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../geometry/generator.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/zonotope.h \
  ../../include/geometry/../geometry/../linear_algebra/zonotopic_vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../base/utility.h \
  ../../include/geometry/../geometry/../linear_algebra/../base/interval.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/interval_vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/../base/basic_type.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/../base/numerical_type.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/../base/interval.h \
  ../../include/geometry/../geometry/../base/interval.h \
  ../../include/geometry/../geometry/../geometry/list_set.h \
  ../../include/geometry/../geometry/../geometry/parallelotope.h \
  ../../include/geometry/../geometry/../geometry/polyhedron.h \
  ../../include/geometry/list_set.h ../../include/geometry/grid_set.h \
  ../../include/geometry/../base/array.h \
  ../../include/geometry/../geometry/grid.h \
  ../../include/geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../base/binary_word.h \
  ../../include/geometry/../geometry/../geometry/lattice_set.h \
  ../../include/geometry/../geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../geometry/../base/iterator.h \
  ../../include/geometry/../geometry/../geometry/../base/../base/array.h \
  ../../include/geometry/../geometry/../geometry/../base/basic_type.h \
  ../../include/geometry/../geometry/../geometry/../base/interval.h \
  ../../include/geometry/../geometry/../geometry/../base/array_operations.h \
  ../../include/geometry/../geometry/../geometry/../base/../base/basic_type.h \
  ../../include/geometry/../geometry/lattice_set.h \
  ../../include/geometry/partition_tree_set.h \
  ../../include/geometry/../base/basic_type.h \
  ../../include/geometry/../base/binary_tree.h \
  ../../include/geometry/../base/../base/binary_word.h \
  ../../include/geometry/../base/iterator.h \
  ../../include/geometry/../geometry/subdivision_tree_set.h \
  ../../include/geometry/../geometry/../base/basic_type.h \
  ../../include/geometry/../geometry/../base/sequence.h \
  ../../include/geometry/../geometry/../base/binary_tree.h \
  ../../include/geometry/../geometry/../base/iterator.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../base/numerical_type.h \
  ../../include/python/../declarations.h \
  ../../include/python/../base/basic_type.h \
  ../../include/python/python_utilities.h

export_grid_set.o: export_grid_set.cc ../../include/base/numerical_type.h \
  ../../include/base/../base/dyadic.h \
  ../../include/linear_algebra/linear_algebra.h \
  ../../include/linear_algebra/../linear_algebra/vector.h \
  ../../include/linear_algebra/../linear_algebra/../base/basic_type.h \
  ../../include/linear_algebra/../linear_algebra/../base/numerical_type.h \
  ../../include/linear_algebra/../linear_algebra/matrix.h \
  ../../include/linear_algebra/../linear_algebra/../base/arithmetic.h \
  ../../include/linear_algebra/../linear_algebra/../base/../base/numerical_type.h \
  ../../include/linear_algebra/../linear_algebra/../base/../base/approximation.h \
  ../../include/linear_algebra/../linear_algebra/../base/../base/../base/numerical_type.h \
  ../../include/linear_algebra/../linear_algebra/../linear_algebra/vector.h \
  ../../include/linear_algebra/../linear_algebra/constraint.h \
  /ufs/collins/include/ppl.hh \
  ../../include/linear_algebra/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/parallelotope.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../base/utility.h \
  ../../include/geometry/../base/../base/array.h \
  ../../include/geometry/../base/../base/basic_type.h \
  ../../include/geometry/../base/interval.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../base/utility.h \
  ../../include/geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/../base/interval.h \
  ../../include/geometry/../geometry/../base/binary_word.h \
  ../../include/geometry/../geometry/../base/../base/basic_type.h \
  ../../include/geometry/../geometry/../base/../base/array.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/list_set.h \
  ../../include/geometry/../geometry/polyhedron.h \
  ../../include/geometry/../geometry/../linear_algebra/constraint.h \
  ../../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../geometry/generator.h \
  ../../include/geometry/../geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/zonotope.h \
  ../../include/geometry/../geometry/../linear_algebra/zonotopic_vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../base/utility.h \
  ../../include/geometry/../geometry/../linear_algebra/../base/interval.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/interval_vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/../base/basic_type.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/../base/numerical_type.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/../base/interval.h \
  ../../include/geometry/../geometry/../geometry/list_set.h \
  ../../include/geometry/../geometry/../geometry/parallelotope.h \
  ../../include/geometry/../geometry/../geometry/polyhedron.h \
  ../../include/geometry/list_set.h ../../include/geometry/grid_set.h \
  ../../include/geometry/../base/array.h \
  ../../include/geometry/../base/binary_word.h \
  ../../include/geometry/../geometry/grid.h \
  ../../include/geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../geometry/lattice_set.h \
  ../../include/geometry/../geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../geometry/../base/iterator.h \
  ../../include/geometry/../geometry/../geometry/../base/../base/array.h \
  ../../include/geometry/../geometry/../geometry/../base/basic_type.h \
  ../../include/geometry/../geometry/../geometry/../base/interval.h \
  ../../include/geometry/../geometry/../geometry/../base/array_operations.h \
  ../../include/geometry/../geometry/../geometry/../base/../base/basic_type.h \
  ../../include/geometry/../geometry/lattice_set.h \
  ../../include/geometry/partition_tree_set.h \
  ../../include/geometry/../base/basic_type.h \
  ../../include/geometry/../base/binary_tree.h \
  ../../include/geometry/../base/../base/basic_type.h \
  ../../include/geometry/../base/../base/binary_word.h \
  ../../include/geometry/../base/iterator.h \
  ../../include/geometry/../geometry/subdivision_tree_set.h \
  ../../include/geometry/../geometry/../base/basic_type.h \
  ../../include/geometry/../geometry/../base/sequence.h \
  ../../include/geometry/../geometry/../base/binary_tree.h \
  ../../include/geometry/../geometry/../base/iterator.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../base/numerical_type.h \
  ../../include/python/../declarations.h \
  ../../include/python/../base/basic_type.h

export_partition_tree_set.o: export_partition_tree_set.cc \
  ../../include/base/numerical_type.h ../../include/base/../base/dyadic.h \
  ../../include/base/binary_word.h \
  ../../include/base/../base/basic_type.h \
  ../../include/base/../base/array.h \
  ../../include/base/../base/basic_type.h \
  ../../include/base/binary_tree.h \
  ../../include/base/../base/binary_word.h \
  ../../include/geometry/rectangle.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../base/utility.h \
  ../../include/geometry/../base/../base/array.h \
  ../../include/geometry/../base/interval.h \
  ../../include/geometry/../base/binary_word.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../base/utility.h \
  ../../include/geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../base/basic_type.h \
  ../../include/geometry/../geometry/../linear_algebra/../base/numerical_type.h \
  ../../include/geometry/grid_set.h \
  ../../include/geometry/../base/array.h \
  ../../include/geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/grid.h \
  ../../include/geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../base/interval.h \
  ../../include/geometry/../geometry/../base/binary_word.h \
  ../../include/geometry/../geometry/../geometry/lattice_set.h \
  ../../include/geometry/../geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../geometry/../base/iterator.h \
  ../../include/geometry/../geometry/../geometry/../base/../base/array.h \
  ../../include/geometry/../geometry/../geometry/../base/basic_type.h \
  ../../include/geometry/../geometry/../geometry/../base/interval.h \
  ../../include/geometry/../geometry/../geometry/../base/array_operations.h \
  ../../include/geometry/../geometry/../geometry/../base/../base/basic_type.h \
  ../../include/geometry/../geometry/lattice_set.h \
  ../../include/geometry/../geometry/list_set.h \
  ../../include/geometry/partition_tree_set.h \
  ../../include/geometry/../base/basic_type.h \
  ../../include/geometry/../base/binary_tree.h \
  ../../include/geometry/../base/iterator.h \
  ../../include/geometry/../geometry/subdivision_tree_set.h \
  ../../include/geometry/../geometry/../base/basic_type.h \
  ../../include/geometry/../geometry/../base/sequence.h \
  ../../include/geometry/../geometry/../base/binary_tree.h \
  ../../include/geometry/../geometry/../base/iterator.h \
  ../../include/geometry/../geometry/../geometry/rectangle.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../base/numerical_type.h \
  ../../include/python/../declarations.h \
  ../../include/python/../base/basic_type.h

export_postscript_output.o: export_postscript_output.cc \
  ../../include/geometry/rectangle.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../base/utility.h \
  ../../include/geometry/../base/../base/array.h \
  ../../include/geometry/../base/../base/basic_type.h \
  ../../include/geometry/../base/interval.h \
  ../../include/geometry/../base/binary_word.h \
  ../../include/geometry/../base/../base/basic_type.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../base/utility.h \
  ../../include/geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../base/basic_type.h \
  ../../include/geometry/../geometry/../linear_algebra/../base/numerical_type.h \
  ../../include/geometry/../geometry/../linear_algebra/../base/../base/dyadic.h \
  ../../include/geometry/parallelotope.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../linear_algebra/../base/basic_type.h \
  ../../include/geometry/../linear_algebra/../base/numerical_type.h \
  ../../include/geometry/../linear_algebra/../base/arithmetic.h \
  ../../include/geometry/../linear_algebra/../base/../base/numerical_type.h \
  ../../include/geometry/../linear_algebra/../base/../base/approximation.h \
  ../../include/geometry/../linear_algebra/../base/../base/../base/numerical_type.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/list_set.h \
  ../../include/geometry/../geometry/polyhedron.h \
  /ufs/collins/include/ppl.hh \
  ../../include/geometry/../geometry/../linear_algebra/constraint.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../geometry/generator.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/zonotope.h \
  ../../include/geometry/../geometry/../linear_algebra/zonotopic_vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../base/utility.h \
  ../../include/geometry/../geometry/../linear_algebra/../base/interval.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/interval_vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/../base/basic_type.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/../base/numerical_type.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/../base/interval.h \
  ../../include/geometry/../geometry/../base/interval.h \
  ../../include/geometry/../geometry/../geometry/list_set.h \
  ../../include/geometry/../geometry/../geometry/parallelotope.h \
  ../../include/geometry/../geometry/../geometry/polyhedron.h \
  ../../include/geometry/list_set.h ../../include/geometry/grid_set.h \
  ../../include/geometry/../base/array.h \
  ../../include/geometry/../geometry/grid.h \
  ../../include/geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../base/binary_word.h \
  ../../include/geometry/../geometry/../geometry/lattice_set.h \
  ../../include/geometry/../geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../geometry/../base/iterator.h \
  ../../include/geometry/../geometry/../geometry/../base/../base/array.h \
  ../../include/geometry/../geometry/../geometry/../base/basic_type.h \
  ../../include/geometry/../geometry/../geometry/../base/interval.h \
  ../../include/geometry/../geometry/../geometry/../base/array_operations.h \
  ../../include/geometry/../geometry/../geometry/../base/../base/basic_type.h \
  ../../include/geometry/../geometry/lattice_set.h \
  ../../include/geometry/partition_tree_set.h \
  ../../include/geometry/../base/basic_type.h \
  ../../include/geometry/../base/binary_tree.h \
  ../../include/geometry/../base/../base/binary_word.h \
  ../../include/geometry/../base/iterator.h \
  ../../include/geometry/../geometry/subdivision_tree_set.h \
  ../../include/geometry/../geometry/../base/basic_type.h \
  ../../include/geometry/../geometry/../base/sequence.h \
  ../../include/geometry/../geometry/../base/binary_tree.h \
  ../../include/geometry/../geometry/../base/iterator.h \
  ../../include/utility/epsfstream.h \
  ../../include/utility/../geometry/rectangle.h \
  ../../include/utility/../geometry/list_set.h \
  ../../include/utility/../geometry/grid_set.h \
  ../../include/utility/../geometry/parallelotope.h \
  ../../include/utility/../geometry/partition_tree_set.h \
  ../../include/utility/../linear_algebra/matrix.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../base/numerical_type.h \
  ../../include/python/../declarations.h \
  ../../include/python/../base/basic_type.h

evaluation_module.o: evaluation_module.cc

export_map.o: export_map.cc ../../include/evaluation/map.h \
  ../../include/evaluation/../declarations.h \
  ../../include/evaluation/../linear_algebra/vector.h \
  ../../include/evaluation/../linear_algebra/../base/basic_type.h \
  ../../include/evaluation/../linear_algebra/../base/numerical_type.h \
  ../../include/evaluation/../linear_algebra/../base/../base/dyadic.h \
  ../../include/evaluation/../linear_algebra/matrix.h \
  ../../include/evaluation/../linear_algebra/../base/arithmetic.h \
  ../../include/evaluation/../linear_algebra/../base/../base/numerical_type.h \
  ../../include/evaluation/../linear_algebra/../base/../base/approximation.h \
  ../../include/evaluation/../linear_algebra/../base/../base/../base/numerical_type.h \
  ../../include/evaluation/../linear_algebra/../linear_algebra/vector.h \
  ../../include/evaluation/../linear_algebra/interval_matrix.h \
  ../../include/evaluation/../linear_algebra/../declarations.h \
  ../../include/evaluation/../linear_algebra/../base/interval.h \
  ../../include/evaluation/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/evaluation/../linear_algebra/../linear_algebra/interval_vector.h \
  ../../include/evaluation/../linear_algebra/../linear_algebra/../base/basic_type.h \
  ../../include/evaluation/../linear_algebra/../linear_algebra/../base/numerical_type.h \
  ../../include/evaluation/../linear_algebra/../linear_algebra/../base/interval.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../base/numerical_type.h \
  ../../include/python/../declarations.h \
  ../../include/python/../base/basic_type.h

export_affine_map.o: export_affine_map.cc \
  ../../include/base/numerical_type.h ../../include/base/../base/dyadic.h \
  ../../include/geometry/point.h ../../include/geometry/../declarations.h \
  ../../include/geometry/../base/utility.h \
  ../../include/geometry/../base/../base/array.h \
  ../../include/geometry/../base/../base/basic_type.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/../base/basic_type.h \
  ../../include/geometry/../linear_algebra/../base/numerical_type.h \
  ../../include/geometry/rectangle.h \
  ../../include/geometry/../base/interval.h \
  ../../include/geometry/../base/binary_word.h \
  ../../include/geometry/../base/../base/basic_type.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/parallelotope.h \
  ../../include/geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../linear_algebra/../base/arithmetic.h \
  ../../include/geometry/../linear_algebra/../base/../base/numerical_type.h \
  ../../include/geometry/../linear_algebra/../base/../base/approximation.h \
  ../../include/geometry/../linear_algebra/../base/../base/../base/numerical_type.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/list_set.h \
  ../../include/geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/polyhedron.h \
  /ufs/collins/include/ppl.hh \
  ../../include/geometry/../geometry/../linear_algebra/constraint.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../geometry/generator.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/zonotope.h \
  ../../include/geometry/../geometry/../linear_algebra/zonotopic_vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../base/utility.h \
  ../../include/geometry/../geometry/../linear_algebra/../base/interval.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/interval_vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/../base/basic_type.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/../base/numerical_type.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/../base/interval.h \
  ../../include/geometry/../geometry/../base/utility.h \
  ../../include/geometry/../geometry/../base/interval.h \
  ../../include/geometry/../geometry/../geometry/list_set.h \
  ../../include/geometry/../geometry/../geometry/parallelotope.h \
  ../../include/geometry/../geometry/../geometry/polyhedron.h \
  ../../include/geometry/simplex.h ../../include/geometry/../base/array.h \
  ../../include/geometry/zonotope.h ../../include/geometry/polyhedron.h \
  ../../include/evaluation/affine_map.h \
  ../../include/evaluation/../linear_algebra/vector.h \
  ../../include/evaluation/../linear_algebra/matrix.h \
  ../../include/evaluation/../linear_algebra/interval_matrix.h \
  ../../include/evaluation/../linear_algebra/../declarations.h \
  ../../include/evaluation/../linear_algebra/../base/numerical_type.h \
  ../../include/evaluation/../linear_algebra/../base/interval.h \
  ../../include/evaluation/../linear_algebra/../linear_algebra/vector.h \
  ../../include/evaluation/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/evaluation/../linear_algebra/../linear_algebra/interval_vector.h \
  ../../include/evaluation/../geometry/point.h \
  ../../include/evaluation/../geometry/simplex.h \
  ../../include/evaluation/../geometry/rectangle.h \
  ../../include/evaluation/../geometry/parallelotope.h \
  ../../include/evaluation/../geometry/zonotope.h \
  ../../include/evaluation/../geometry/polyhedron.h \
  ../../include/evaluation/../evaluation/map.h \
  ../../include/evaluation/../evaluation/../declarations.h \
  ../../include/evaluation/../evaluation/../linear_algebra/vector.h \
  ../../include/evaluation/../evaluation/../linear_algebra/matrix.h \
  ../../include/evaluation/../evaluation/../linear_algebra/interval_matrix.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../base/numerical_type.h \
  ../../include/python/../declarations.h \
  ../../include/python/../base/basic_type.h \
  ../../include/python/python_utilities.h

export_polynomial_map.o: export_polynomial_map.cc \
  ../../include/base/numerical_type.h ../../include/base/../base/dyadic.h \
  ../../include/geometry/point.h ../../include/geometry/../declarations.h \
  ../../include/geometry/../base/utility.h \
  ../../include/geometry/../base/../base/array.h \
  ../../include/geometry/../base/../base/basic_type.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/../base/basic_type.h \
  ../../include/geometry/../linear_algebra/../base/numerical_type.h \
  ../../include/evaluation/polynomial_map.h \
  ../../include/evaluation/../linear_algebra/vector.h \
  ../../include/evaluation/../linear_algebra/matrix.h \
  ../../include/evaluation/../linear_algebra/../base/basic_type.h \
  ../../include/evaluation/../linear_algebra/../base/numerical_type.h \
  ../../include/evaluation/../linear_algebra/../base/arithmetic.h \
  ../../include/evaluation/../linear_algebra/../base/../base/numerical_type.h \
  ../../include/evaluation/../linear_algebra/../base/../base/approximation.h \
  ../../include/evaluation/../linear_algebra/../base/../base/../base/numerical_type.h \
  ../../include/evaluation/../linear_algebra/../linear_algebra/vector.h \
  ../../include/evaluation/../base/utility.h \
  ../../include/evaluation/../base/array.h \
  ../../include/evaluation/../base/interval.h \
  ../../include/evaluation/../geometry/point.h \
  ../../include/evaluation/../geometry/rectangle.h \
  ../../include/evaluation/../geometry/../declarations.h \
  ../../include/evaluation/../geometry/../base/utility.h \
  ../../include/evaluation/../geometry/../base/interval.h \
  ../../include/evaluation/../geometry/../base/binary_word.h \
  ../../include/evaluation/../geometry/../base/../base/basic_type.h \
  ../../include/evaluation/../geometry/../base/../base/array.h \
  ../../include/evaluation/../geometry/../geometry/point.h \
  ../../include/evaluation/../geometry/polyhedron.h \
  /ufs/collins/include/ppl.hh \
  ../../include/evaluation/../geometry/../linear_algebra/constraint.h \
  ../../include/evaluation/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/evaluation/../geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/evaluation/../geometry/../linear_algebra/vector.h \
  ../../include/evaluation/../geometry/../linear_algebra/matrix.h \
  ../../include/evaluation/../geometry/../geometry/generator.h \
  ../../include/evaluation/../geometry/../geometry/rectangle.h \
  ../../include/evaluation/../geometry/list_set.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../base/numerical_type.h \
  ../../include/python/../declarations.h \
  ../../include/python/../base/basic_type.h

export_vector_field.o: export_vector_field.cc \
  ../../include/evaluation/vector_field.h \
  ../../include/evaluation/../declarations.h \
  ../../include/evaluation/../base/interval.h \
  ../../include/evaluation/../linear_algebra/vector.h \
  ../../include/evaluation/../linear_algebra/../base/basic_type.h \
  ../../include/evaluation/../linear_algebra/../base/numerical_type.h \
  ../../include/evaluation/../linear_algebra/../base/../base/dyadic.h \
  ../../include/evaluation/../linear_algebra/matrix.h \
  ../../include/evaluation/../linear_algebra/../base/arithmetic.h \
  ../../include/evaluation/../linear_algebra/../base/../base/numerical_type.h \
  ../../include/evaluation/../linear_algebra/../base/../base/approximation.h \
  ../../include/evaluation/../linear_algebra/../base/../base/../base/numerical_type.h \
  ../../include/evaluation/../linear_algebra/../linear_algebra/vector.h \
  ../../include/evaluation/../linear_algebra/interval_vector.h \
  ../../include/evaluation/../linear_algebra/../base/interval.h \
  ../../include/evaluation/../linear_algebra/interval_matrix.h \
  ../../include/evaluation/../linear_algebra/../declarations.h \
  ../../include/evaluation/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/evaluation/../linear_algebra/../linear_algebra/interval_vector.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../base/numerical_type.h \
  ../../include/python/../declarations.h \
  ../../include/python/../base/basic_type.h

export_affine_vector_field.o: export_affine_vector_field.cc \
  ../../include/base/numerical_type.h ../../include/base/../base/dyadic.h \
  ../../include/base/interval.h ../../include/geometry/point.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../base/utility.h \
  ../../include/geometry/../base/../base/array.h \
  ../../include/geometry/../base/../base/basic_type.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/../base/basic_type.h \
  ../../include/geometry/../linear_algebra/../base/numerical_type.h \
  ../../include/geometry/rectangle.h \
  ../../include/geometry/../base/interval.h \
  ../../include/geometry/../base/binary_word.h \
  ../../include/geometry/../base/../base/basic_type.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/parallelotope.h \
  ../../include/geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../linear_algebra/../base/arithmetic.h \
  ../../include/geometry/../linear_algebra/../base/../base/numerical_type.h \
  ../../include/geometry/../linear_algebra/../base/../base/approximation.h \
  ../../include/geometry/../linear_algebra/../base/../base/../base/numerical_type.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/list_set.h \
  ../../include/geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/polyhedron.h \
  /ufs/collins/include/ppl.hh \
  ../../include/geometry/../geometry/../linear_algebra/constraint.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../geometry/generator.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/zonotope.h \
  ../../include/geometry/../geometry/../linear_algebra/zonotopic_vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../base/utility.h \
  ../../include/geometry/../geometry/../linear_algebra/../base/interval.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/interval_vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/../base/basic_type.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/../base/numerical_type.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/../base/interval.h \
  ../../include/geometry/../geometry/../base/utility.h \
  ../../include/geometry/../geometry/../base/interval.h \
  ../../include/geometry/../geometry/../geometry/list_set.h \
  ../../include/geometry/../geometry/../geometry/parallelotope.h \
  ../../include/geometry/../geometry/../geometry/polyhedron.h \
  ../../include/geometry/simplex.h ../../include/geometry/../base/array.h \
  ../../include/geometry/zonotope.h ../../include/geometry/polyhedron.h \
  ../../include/evaluation/affine_vector_field.h \
  ../../include/evaluation/../linear_algebra/vector.h \
  ../../include/evaluation/../linear_algebra/matrix.h \
  ../../include/evaluation/../linear_algebra/interval_vector.h \
  ../../include/evaluation/../linear_algebra/interval_matrix.h \
  ../../include/evaluation/../linear_algebra/../declarations.h \
  ../../include/evaluation/../linear_algebra/../base/numerical_type.h \
  ../../include/evaluation/../linear_algebra/../base/interval.h \
  ../../include/evaluation/../linear_algebra/../linear_algebra/vector.h \
  ../../include/evaluation/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/evaluation/../linear_algebra/../linear_algebra/interval_vector.h \
  ../../include/evaluation/../evaluation/vector_field.h \
  ../../include/evaluation/../evaluation/../declarations.h \
  ../../include/evaluation/../evaluation/../base/interval.h \
  ../../include/evaluation/../evaluation/../linear_algebra/vector.h \
  ../../include/evaluation/../evaluation/../linear_algebra/matrix.h \
  ../../include/evaluation/../evaluation/../linear_algebra/interval_vector.h \
  ../../include/evaluation/../evaluation/../linear_algebra/interval_matrix.h \
  ../../include/evaluation/../evaluation/affine_map.h \
  ../../include/evaluation/../evaluation/../geometry/point.h \
  ../../include/evaluation/../evaluation/../geometry/simplex.h \
  ../../include/evaluation/../evaluation/../geometry/rectangle.h \
  ../../include/evaluation/../evaluation/../geometry/parallelotope.h \
  ../../include/evaluation/../evaluation/../geometry/zonotope.h \
  ../../include/evaluation/../evaluation/../geometry/polyhedron.h \
  ../../include/evaluation/../evaluation/../evaluation/map.h \
  ../../include/evaluation/../evaluation/../evaluation/../declarations.h \
  ../../include/evaluation/../evaluation/../evaluation/../linear_algebra/vector.h \
  ../../include/evaluation/../evaluation/../evaluation/../linear_algebra/matrix.h \
  ../../include/evaluation/../evaluation/../evaluation/../linear_algebra/interval_matrix.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../base/numerical_type.h \
  ../../include/python/../declarations.h \
  ../../include/python/../base/basic_type.h \
  ../../include/python/python_utilities.h

export_apply.o: export_apply.cc ../../include/geometry/rectangle.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../base/utility.h \
  ../../include/geometry/../base/../base/array.h \
  ../../include/geometry/../base/../base/basic_type.h \
  ../../include/geometry/../base/interval.h \
  ../../include/geometry/../base/binary_word.h \
  ../../include/geometry/../base/../base/basic_type.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../base/utility.h \
  ../../include/geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../base/basic_type.h \
  ../../include/geometry/../geometry/../linear_algebra/../base/numerical_type.h \
  ../../include/geometry/../geometry/../linear_algebra/../base/../base/dyadic.h \
  ../../include/geometry/parallelotope.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../linear_algebra/../base/basic_type.h \
  ../../include/geometry/../linear_algebra/../base/numerical_type.h \
  ../../include/geometry/../linear_algebra/../base/arithmetic.h \
  ../../include/geometry/../linear_algebra/../base/../base/numerical_type.h \
  ../../include/geometry/../linear_algebra/../base/../base/approximation.h \
  ../../include/geometry/../linear_algebra/../base/../base/../base/numerical_type.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/list_set.h \
  ../../include/geometry/../geometry/polyhedron.h \
  /ufs/collins/include/ppl.hh \
  ../../include/geometry/../geometry/../linear_algebra/constraint.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../geometry/generator.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/zonotope.h \
  ../../include/geometry/../geometry/../linear_algebra/zonotopic_vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../base/utility.h \
  ../../include/geometry/../geometry/../linear_algebra/../base/interval.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/interval_vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/../base/basic_type.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/../base/numerical_type.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/../base/interval.h \
  ../../include/geometry/../geometry/../base/interval.h \
  ../../include/geometry/../geometry/../geometry/list_set.h \
  ../../include/geometry/../geometry/../geometry/parallelotope.h \
  ../../include/geometry/../geometry/../geometry/polyhedron.h \
  ../../include/geometry/list_set.h ../../include/geometry/grid_set.h \
  ../../include/geometry/../base/array.h \
  ../../include/geometry/../geometry/grid.h \
  ../../include/geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../base/binary_word.h \
  ../../include/geometry/../geometry/../geometry/lattice_set.h \
  ../../include/geometry/../geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../geometry/../base/iterator.h \
  ../../include/geometry/../geometry/../geometry/../base/../base/array.h \
  ../../include/geometry/../geometry/../geometry/../base/basic_type.h \
  ../../include/geometry/../geometry/../geometry/../base/interval.h \
  ../../include/geometry/../geometry/../geometry/../base/array_operations.h \
  ../../include/geometry/../geometry/../geometry/../base/../base/basic_type.h \
  ../../include/geometry/../geometry/lattice_set.h \
  ../../include/evaluation/apply.h \
  ../../include/evaluation/../declarations.h \
  ../../include/evaluation/../evaluation/map.h \
  ../../include/evaluation/../evaluation/../declarations.h \
  ../../include/evaluation/../evaluation/../linear_algebra/vector.h \
  ../../include/evaluation/../evaluation/../linear_algebra/matrix.h \
  ../../include/evaluation/../evaluation/../linear_algebra/interval_matrix.h \
  ../../include/evaluation/../evaluation/../linear_algebra/../declarations.h \
  ../../include/evaluation/../evaluation/../linear_algebra/../base/numerical_type.h \
  ../../include/evaluation/../evaluation/../linear_algebra/../base/interval.h \
  ../../include/evaluation/../evaluation/../linear_algebra/../linear_algebra/vector.h \
  ../../include/evaluation/../evaluation/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/evaluation/../evaluation/../linear_algebra/../linear_algebra/interval_vector.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../base/numerical_type.h \
  ../../include/python/../declarations.h \
  ../../include/python/../base/basic_type.h

export_integrate.o: export_integrate.cc \
  ../../include/geometry/rectangle.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../base/utility.h \
  ../../include/geometry/../base/../base/array.h \
  ../../include/geometry/../base/../base/basic_type.h \
  ../../include/geometry/../base/interval.h \
  ../../include/geometry/../base/binary_word.h \
  ../../include/geometry/../base/../base/basic_type.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../base/utility.h \
  ../../include/geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../base/basic_type.h \
  ../../include/geometry/../geometry/../linear_algebra/../base/numerical_type.h \
  ../../include/geometry/../geometry/../linear_algebra/../base/../base/dyadic.h \
  ../../include/geometry/parallelotope.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../linear_algebra/../base/basic_type.h \
  ../../include/geometry/../linear_algebra/../base/numerical_type.h \
  ../../include/geometry/../linear_algebra/../base/arithmetic.h \
  ../../include/geometry/../linear_algebra/../base/../base/numerical_type.h \
  ../../include/geometry/../linear_algebra/../base/../base/approximation.h \
  ../../include/geometry/../linear_algebra/../base/../base/../base/numerical_type.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/list_set.h \
  ../../include/geometry/../geometry/polyhedron.h \
  /ufs/collins/include/ppl.hh \
  ../../include/geometry/../geometry/../linear_algebra/constraint.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../geometry/generator.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/zonotope.h \
  ../../include/geometry/../geometry/../linear_algebra/zonotopic_vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../base/utility.h \
  ../../include/geometry/../geometry/../linear_algebra/../base/interval.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/interval_vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/../base/basic_type.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/../base/numerical_type.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/../base/interval.h \
  ../../include/geometry/../geometry/../base/interval.h \
  ../../include/geometry/../geometry/../geometry/list_set.h \
  ../../include/geometry/../geometry/../geometry/parallelotope.h \
  ../../include/geometry/../geometry/../geometry/polyhedron.h \
  ../../include/geometry/zonotope.h ../../include/geometry/grid_set.h \
  ../../include/geometry/../base/array.h \
  ../../include/geometry/../geometry/grid.h \
  ../../include/geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../base/binary_word.h \
  ../../include/geometry/../geometry/../geometry/lattice_set.h \
  ../../include/geometry/../geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../geometry/../base/iterator.h \
  ../../include/geometry/../geometry/../geometry/../base/../base/array.h \
  ../../include/geometry/../geometry/../geometry/../base/basic_type.h \
  ../../include/geometry/../geometry/../geometry/../base/interval.h \
  ../../include/geometry/../geometry/../geometry/../base/array_operations.h \
  ../../include/geometry/../geometry/../geometry/../base/../base/basic_type.h \
  ../../include/geometry/../geometry/lattice_set.h \
  ../../include/evaluation/integrate.h \
  ../../include/evaluation/../declarations.h \
  ../../include/evaluation/../evaluation/vector_field.h \
  ../../include/evaluation/../evaluation/../declarations.h \
  ../../include/evaluation/../evaluation/../base/interval.h \
  ../../include/evaluation/../evaluation/../linear_algebra/vector.h \
  ../../include/evaluation/../evaluation/../linear_algebra/matrix.h \
  ../../include/evaluation/../evaluation/../linear_algebra/interval_vector.h \
  ../../include/evaluation/../evaluation/../linear_algebra/interval_matrix.h \
  ../../include/evaluation/../evaluation/../linear_algebra/../declarations.h \
  ../../include/evaluation/../evaluation/../linear_algebra/../base/numerical_type.h \
  ../../include/evaluation/../evaluation/../linear_algebra/../base/interval.h \
  ../../include/evaluation/../evaluation/../linear_algebra/../linear_algebra/vector.h \
  ../../include/evaluation/../evaluation/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/evaluation/../evaluation/../linear_algebra/../linear_algebra/interval_vector.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../base/numerical_type.h \
  ../../include/python/../declarations.h \
  ../../include/python/../base/basic_type.h

export_henon_map.o: export_henon_map.cc \
  ../../include/evaluation/henon_map.h \
  ../../include/evaluation/../linear_algebra/matrix.h \
  ../../include/evaluation/../linear_algebra/../base/basic_type.h \
  ../../include/evaluation/../linear_algebra/../base/numerical_type.h \
  ../../include/evaluation/../linear_algebra/../base/../base/dyadic.h \
  ../../include/evaluation/../linear_algebra/../base/arithmetic.h \
  ../../include/evaluation/../linear_algebra/../base/../base/numerical_type.h \
  ../../include/evaluation/../linear_algebra/../base/../base/approximation.h \
  ../../include/evaluation/../linear_algebra/../base/../base/../base/numerical_type.h \
  ../../include/evaluation/../linear_algebra/../linear_algebra/vector.h \
  ../../include/evaluation/../linear_algebra/../linear_algebra/../base/basic_type.h \
  ../../include/evaluation/../linear_algebra/../linear_algebra/../base/numerical_type.h \
  ../../include/evaluation/../geometry/point.h \
  ../../include/evaluation/../geometry/../declarations.h \
  ../../include/evaluation/../geometry/../base/utility.h \
  ../../include/evaluation/../geometry/../base/../base/array.h \
  ../../include/evaluation/../geometry/../base/../base/basic_type.h \
  ../../include/evaluation/../geometry/../linear_algebra/vector.h \
  ../../include/evaluation/../geometry/rectangle.h \
  ../../include/evaluation/../geometry/../base/interval.h \
  ../../include/evaluation/../geometry/../base/binary_word.h \
  ../../include/evaluation/../geometry/../base/../base/basic_type.h \
  ../../include/evaluation/../geometry/../geometry/point.h \
  ../../include/evaluation/../geometry/parallelotope.h \
  ../../include/evaluation/../geometry/../linear_algebra/matrix.h \
  ../../include/evaluation/../geometry/../geometry/rectangle.h \
  ../../include/evaluation/../geometry/../geometry/list_set.h \
  ../../include/evaluation/../geometry/../geometry/../declarations.h \
  ../../include/evaluation/../geometry/../geometry/polyhedron.h \
  /ufs/collins/include/ppl.hh \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/constraint.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/vector.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/generator.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/point.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/rectangle.h \
  ../../include/evaluation/../geometry/../geometry/zonotope.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/zonotopic_vector.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/../base/utility.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/../base/interval.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/../linear_algebra/interval_vector.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/../linear_algebra/../base/basic_type.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/../linear_algebra/../base/numerical_type.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/../linear_algebra/../base/interval.h \
  ../../include/evaluation/../geometry/../geometry/../base/utility.h \
  ../../include/evaluation/../geometry/../geometry/../base/interval.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/list_set.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/parallelotope.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/polyhedron.h \
  ../../include/evaluation/../evaluation/map.h \
  ../../include/evaluation/../evaluation/../declarations.h \
  ../../include/evaluation/../evaluation/../linear_algebra/vector.h \
  ../../include/evaluation/../evaluation/../linear_algebra/matrix.h \
  ../../include/evaluation/../evaluation/../linear_algebra/interval_matrix.h \
  ../../include/evaluation/../evaluation/../linear_algebra/../declarations.h \
  ../../include/evaluation/../evaluation/../linear_algebra/../base/numerical_type.h \
  ../../include/evaluation/../evaluation/../linear_algebra/../base/interval.h \
  ../../include/evaluation/../evaluation/../linear_algebra/../linear_algebra/vector.h \
  ../../include/evaluation/../evaluation/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/evaluation/../evaluation/../linear_algebra/../linear_algebra/interval_vector.h \
  ../../include/evaluation/../evaluation/apply.h \
  ../../include/evaluation/../evaluation/../evaluation/map.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../base/numerical_type.h \
  ../../include/python/../declarations.h \
  ../../include/python/../base/basic_type.h

export_lorenz_system.o: export_lorenz_system.cc \
  ../../include/evaluation/lorenz_system.h \
  ../../include/evaluation/../linear_algebra/vector.h \
  ../../include/evaluation/../linear_algebra/../base/basic_type.h \
  ../../include/evaluation/../linear_algebra/../base/numerical_type.h \
  ../../include/evaluation/../linear_algebra/../base/../base/dyadic.h \
  ../../include/evaluation/../linear_algebra/interval_vector.h \
  ../../include/evaluation/../linear_algebra/../base/interval.h \
  ../../include/evaluation/../linear_algebra/matrix.h \
  ../../include/evaluation/../linear_algebra/../base/arithmetic.h \
  ../../include/evaluation/../linear_algebra/../base/../base/numerical_type.h \
  ../../include/evaluation/../linear_algebra/../base/../base/approximation.h \
  ../../include/evaluation/../linear_algebra/../base/../base/../base/numerical_type.h \
  ../../include/evaluation/../linear_algebra/../linear_algebra/vector.h \
  ../../include/evaluation/../linear_algebra/interval_matrix.h \
  ../../include/evaluation/../linear_algebra/../declarations.h \
  ../../include/evaluation/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/evaluation/../linear_algebra/../linear_algebra/interval_vector.h \
  ../../include/evaluation/../geometry/point.h \
  ../../include/evaluation/../geometry/../declarations.h \
  ../../include/evaluation/../geometry/../base/utility.h \
  ../../include/evaluation/../geometry/../base/../base/array.h \
  ../../include/evaluation/../geometry/../base/../base/basic_type.h \
  ../../include/evaluation/../geometry/../linear_algebra/vector.h \
  ../../include/evaluation/../geometry/rectangle.h \
  ../../include/evaluation/../geometry/../base/interval.h \
  ../../include/evaluation/../geometry/../base/binary_word.h \
  ../../include/evaluation/../geometry/../base/../base/basic_type.h \
  ../../include/evaluation/../geometry/../geometry/point.h \
  ../../include/evaluation/../geometry/parallelotope.h \
  ../../include/evaluation/../geometry/../linear_algebra/matrix.h \
  ../../include/evaluation/../geometry/../geometry/rectangle.h \
  ../../include/evaluation/../geometry/../geometry/list_set.h \
  ../../include/evaluation/../geometry/../geometry/../declarations.h \
  ../../include/evaluation/../geometry/../geometry/polyhedron.h \
  /ufs/collins/include/ppl.hh \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/constraint.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/vector.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/generator.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/point.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/rectangle.h \
  ../../include/evaluation/../geometry/../geometry/zonotope.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/zonotopic_vector.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/../base/utility.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/../base/interval.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/../linear_algebra/interval_vector.h \
  ../../include/evaluation/../geometry/../geometry/../base/utility.h \
  ../../include/evaluation/../geometry/../geometry/../base/interval.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/list_set.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/parallelotope.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/polyhedron.h \
  ../../include/evaluation/../evaluation/vector_field.h \
  ../../include/evaluation/../evaluation/../declarations.h \
  ../../include/evaluation/../evaluation/../base/interval.h \
  ../../include/evaluation/../evaluation/../linear_algebra/vector.h \
  ../../include/evaluation/../evaluation/../linear_algebra/matrix.h \
  ../../include/evaluation/../evaluation/../linear_algebra/interval_vector.h \
  ../../include/evaluation/../evaluation/../linear_algebra/interval_matrix.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../base/numerical_type.h \
  ../../include/python/../declarations.h \
  ../../include/python/../base/basic_type.h

linear_algebra_module.o: linear_algebra_module.cc

export_vector.o: export_vector.cc ../../include/base/numerical_type.h \
  ../../include/base/../base/dyadic.h ../../include/base/interval.h \
  ../../include/linear_algebra/vector.h \
  ../../include/linear_algebra/../base/basic_type.h \
  ../../include/linear_algebra/../base/numerical_type.h \
  ../../include/linear_algebra/interval_vector.h \
  ../../include/linear_algebra/../base/interval.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../base/numerical_type.h \
  ../../include/python/../declarations.h \
  ../../include/python/../base/basic_type.h

export_matrix.o: export_matrix.cc ../../include/base/numerical_type.h \
  ../../include/base/../base/dyadic.h ../../include/base/interval.h \
  ../../include/linear_algebra/matrix.h \
  ../../include/linear_algebra/../base/basic_type.h \
  ../../include/linear_algebra/../base/numerical_type.h \
  ../../include/linear_algebra/../base/arithmetic.h \
  ../../include/linear_algebra/../base/../base/numerical_type.h \
  ../../include/linear_algebra/../base/../base/approximation.h \
  ../../include/linear_algebra/../base/../base/../base/numerical_type.h \
  ../../include/linear_algebra/../linear_algebra/vector.h \
  ../../include/linear_algebra/../linear_algebra/../base/basic_type.h \
  ../../include/linear_algebra/../linear_algebra/../base/numerical_type.h \
  ../../include/linear_algebra/interval_matrix.h \
  ../../include/linear_algebra/../declarations.h \
  ../../include/linear_algebra/../base/interval.h \
  ../../include/linear_algebra/../linear_algebra/matrix.h \
  ../../include/linear_algebra/../linear_algebra/interval_vector.h \
  ../../include/linear_algebra/../linear_algebra/../base/interval.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../base/numerical_type.h \
  ../../include/python/../declarations.h \
  ../../include/python/../base/basic_type.h \
  ../../include/python/python_utilities.h

export_linear_program.o: export_linear_program.cc \
  ../../include/base/numerical_type.h ../../include/base/../base/dyadic.h \
  ../../include/linear_algebra/linear_program.h \
  ../../include/linear_algebra/../base/basic_type.h \
  ../../include/linear_algebra/../linear_algebra/matrix.h \
  ../../include/linear_algebra/../linear_algebra/../base/basic_type.h \
  ../../include/linear_algebra/../linear_algebra/../base/numerical_type.h \
  ../../include/linear_algebra/../linear_algebra/../base/arithmetic.h \
  ../../include/linear_algebra/../linear_algebra/../base/../base/numerical_type.h \
  ../../include/linear_algebra/../linear_algebra/../base/../base/approximation.h \
  ../../include/linear_algebra/../linear_algebra/../base/../base/../base/numerical_type.h \
  ../../include/linear_algebra/../linear_algebra/../linear_algebra/vector.h \
  ../../include/linear_algebra/../linear_algebra/../linear_algebra/../base/basic_type.h \
  ../../include/linear_algebra/../linear_algebra/../linear_algebra/../base/numerical_type.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../base/numerical_type.h \
  ../../include/python/../declarations.h \
  ../../include/python/../base/basic_type.h \
  ../../include/python/python_utilities.h

