include ../../config.mk

CXXFLAGS += -I/usr/include/python2.4/
LDXXFLAGS += -L../src
LDLIBS = -lariadne -lboost_python ${LIBGMPXX} -lppl

BASE_OBJ=base_module.o \
	export_array.o export_binary_tree.o 

NUMERIC_OBJ=numeric_module.o \
	export_numeric.o export_interval.o export_function.o \

LINEAR_ALGEBRA_OBJ = linear_algebra_module.o \
	export_vector.o export_matrix.o export_linear_program.o

GEOMETRY_OBJ = geometry_module.o \
	export_point.o export_rectangle.o export_parallelotope.o \
	export_zonotope.o export_simplex.o export_polyhedron.o \
	export_lattice_set.o export_list_set.o export_grid_set.o \
	export_partition_tree_set.o export_postscript_output.o

SYSTEM_OBJ = system_module.o \
	export_lattice_map.o \
	export_map.o export_affine_map.o export_affine_multimap.o \
	export_polynomial_map.o \
	export_vector_field.o export_affine_vector_field.o \
	export_henon_map.o export_lorenz_system.o 

EVALUATION_OBJ = evaluation_module.o \
	export_apply.o \
	export_integrate.o

all: python 

python: base.so numeric.so linear_algebra.so geometry.so system.so evaluation.so 

init_file:
	@echo -n '__all__=["base","numeric", '> __init__.py
	@echo '"linear_algebra","geometry","system","evaluation"]' >> __init__.py
	@echo >> __init__.py
	@echo 'from base import *' >> __init__.py
	@echo 'from numeric import *' >> __init__.py
	@echo -n 'Real=' >> __init__.py
	if [ "${REAL_TYPE}" == "RATIONAL" ]; then \
	  echo 'Rational' >> __init__.py; \
	elif [ "${REAL_TYPE}" == "MPFLOAT" ]; then \
	  echo 'MPFloat' >> __init__.py; \
	else \
	  echo "Unknown REAL type"; \
	  false; \
	fi 
	@echo "numeric.Real=Real" >> __init__.py

base.so: ${BASE_OBJ}
	${CXX} -shared -o base.so ${LDXXFLAGS} ${BASE_OBJ} ${LDLIBS} ${RPATHFLAGS}

numeric.so: ${NUMERIC_OBJ}
	${CXX} -shared -o numeric.so ${LDXXFLAGS} ${NUMERIC_OBJ} ${LDLIBS} ${RPATHFLAGS}

linear_algebra.so: ${LINEAR_ALGEBRA_OBJ}
	${CXX} -shared -o linear_algebra.so ${LDXXFLAGS} ${LINEAR_ALGEBRA_OBJ} ${LDLIBS} ${RPATHFLAGS}

geometry.so: ${GEOMETRY_OBJ}
	${CXX} -shared -o geometry.so ${LDXXFLAGS} ${GEOMETRY_OBJ} ${LDLIBS} ${RPATHFLAGS}

system.so: ${SYSTEM_OBJ}
	${CXX} -shared -o system.so ${LDXXFLAGS} ${SYSTEM_OBJ} ${LDLIBS} ${RPATHFLAGS}

evaluation.so: ${EVALUATION_OBJ}
	${CXX} -shared -o evaluation.so ${LDXXFLAGS} ${EVALUATION_OBJ} ${LDLIBS} ${RPATHFLAGS}


install: python init_file  
	mkdir -p ${PYTHONPREFIX}/ariadne
	rm -rf ${PYTHONPREFIX}/ariadne/*.py
	rm -rf ${PYTHONPREFIX}/ariadne/*.pyc
	rm -rf ${PYTHONPREFIX}/ariadne/*.so
	cp -f *.so *.py ${PYTHONPREFIX}/ariadne/

clean:
	rm -f *.o *.so ariadne/*.so

dep:  
	@sed -e '/^### Do Not edit this line$$/,$$d' Makefile \
		> MMM.$$$$ && mv MMM.$$$$ Makefile
	@echo "### Do Not edit this line" >> Makefile
	@echo "### Everything below is auto-generated" >> Makefile
	@for ff in ${BASE_OBJ} ${GEOMETRY_OBJ} ${EVALUATION_OBJ} ${SYSTEM_OBJ} ${LINEAR_ALGEBRA_OBJ}; \
	do f=`basename $$ff .o`.cc ; \
	g++ -MM ${CXXFLAGS} -I${PYTHONINCLUDEDIR} $$f; echo ""; \
	done >> Makefile

### Do Not edit this line
### Everything below is auto-generated
base_module.o: base_module.cc

export_array.o: export_array.cc ../../include/base/array.h \
  ../../include/base/../declarations.h \
  ../../include/base/array_operations.h ../../include/utility/stlio.h \
  ../../include/utility/../base/array.h ../../include/python/typedefs.h \
  ../../include/python/../real_typedef.h \
  ../../include/python/../numeric/numerical_types.h \
  ../../include/python/../numeric/../declarations.h \
  ../../include/python/../numeric/../numeric/integer.h \
  ../../include/python/../numeric/../numeric/../declarations.h \
  ../../include/python/../numeric/../numeric/../utility/stlio.h \
  ../../include/python/../numeric/../numeric/dyadic.h \
  ../../include/python/../numeric/../numeric/numerical_traits.h \
  ../../include/python/../numeric/../utility/stlio.h \
  ../../include/python/../declarations.h \
  ../../include/python/python_utilities.h

export_binary_tree.o: export_binary_tree.cc \
  ../../include/base/binary_word.h ../../include/base/../declarations.h \
  ../../include/base/binary_tree.h ../../include/base/../base/array.h \
  ../../include/base/../base/../declarations.h \
  ../../include/base/../base/binary_word.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../numeric/numerical_types.h \
  ../../include/python/../numeric/../declarations.h \
  ../../include/python/../numeric/../numeric/integer.h \
  ../../include/python/../numeric/../numeric/../declarations.h \
  ../../include/python/../numeric/../numeric/../utility/stlio.h \
  ../../include/python/../numeric/../numeric/../utility/../base/array.h \
  ../../include/python/../numeric/../numeric/dyadic.h \
  ../../include/python/../numeric/../numeric/numerical_traits.h \
  ../../include/python/../numeric/../utility/stlio.h \
  ../../include/python/../declarations.h \
  ../../include/python/python_utilities.h

geometry_module.o: geometry_module.cc

export_point.o: export_point.cc ../../include/geometry/point.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../base/array.h \
  ../../include/geometry/../base/../declarations.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../linear_algebra/../numeric/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/../utility/stlio.h \
  ../../include/geometry/../linear_algebra/../numeric/../utility/../base/array.h \
  ../../include/geometry/rectangle.h \
  ../../include/geometry/../numeric/arithmetic.h \
  ../../include/geometry/../numeric/../numeric/numerical_types.h \
  ../../include/geometry/../numeric/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/../numeric/integer.h \
  ../../include/geometry/../numeric/../numeric/../numeric/dyadic.h \
  ../../include/geometry/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../numeric/../numeric/../utility/stlio.h \
  ../../include/geometry/../numeric/../numeric/approximation.h \
  ../../include/geometry/../numeric/../numeric/../numeric/numerical_types.h \
  ../../include/geometry/../numeric/interval.h \
  ../../include/geometry/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../linear_algebra/interval_vector.h \
  ../../include/geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/linear_algebra/vector.h ../../include/python/typedefs.h \
  ../../include/python/../real_typedef.h \
  ../../include/python/../numeric/numerical_types.h \
  ../../include/python/../declarations.h

export_rectangle.o: export_rectangle.cc \
  ../../include/geometry/rectangle.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../numeric/arithmetic.h \
  ../../include/geometry/../numeric/../numeric/numerical_types.h \
  ../../include/geometry/../numeric/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/../numeric/integer.h \
  ../../include/geometry/../numeric/../numeric/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/../numeric/../utility/stlio.h \
  ../../include/geometry/../numeric/../numeric/../numeric/../utility/../base/array.h \
  ../../include/geometry/../numeric/../numeric/../numeric/../utility/../base/../declarations.h \
  ../../include/geometry/../numeric/../numeric/../numeric/dyadic.h \
  ../../include/geometry/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../numeric/../numeric/../utility/stlio.h \
  ../../include/geometry/../numeric/../numeric/approximation.h \
  ../../include/geometry/../numeric/../numeric/../numeric/numerical_types.h \
  ../../include/geometry/../numeric/interval.h \
  ../../include/geometry/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../linear_algebra/interval_vector.h \
  ../../include/geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/zonotope.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../linear_algebra/transformation_system.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/interval_vector.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/interval_matrix.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/list_set.h \
  ../../include/geometry/../base/utility.h \
  ../../include/geometry/../base/../declarations.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../numeric/numerical_types.h \
  ../../include/python/../declarations.h

export_parallelotope.o: export_parallelotope.cc \
  ../../include/geometry/parallelotope.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../numeric/interval.h \
  ../../include/geometry/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../linear_algebra/../numeric/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/../utility/stlio.h \
  ../../include/geometry/../linear_algebra/../numeric/../utility/../base/array.h \
  ../../include/geometry/../linear_algebra/../numeric/../utility/../base/../declarations.h \
  ../../include/geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/polyhedron.h \
  /ufs/collins/include/ppl.hh \
  ../../include/geometry/../geometry/../linear_algebra/constraint_system.h \
  ../../include/geometry/../geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/numerical_types.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/../declarations.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/../numeric/integer.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/../numeric/dyadic.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/../utility/stlio.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../linear_algebra/interval_vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../geometry/../linear_algebra/generator_system.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../geometry/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../geometry/../numeric/../numeric/numerical_types.h \
  ../../include/geometry/../geometry/../geometry/../numeric/../numeric/approximation.h \
  ../../include/geometry/../geometry/../geometry/../numeric/../numeric/../declarations.h \
  ../../include/geometry/../geometry/../geometry/../numeric/../numeric/../numeric/numerical_types.h \
  ../../include/geometry/../geometry/../geometry/../numeric/interval.h \
  ../../include/geometry/../geometry/../geometry/../linear_algebra/interval_vector.h \
  ../../include/geometry/../geometry/../geometry/../geometry/point.h \
  ../../include/geometry/rectangle.h ../../include/geometry/zonotope.h \
  ../../include/geometry/../linear_algebra/transformation_system.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/interval_vector.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/interval_matrix.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/polyhedron.h ../../include/geometry/list_set.h \
  ../../include/geometry/../base/utility.h \
  ../../include/geometry/../base/../declarations.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../numeric/numerical_types.h \
  ../../include/python/../declarations.h \
  ../../include/python/python_utilities.h

export_zonotope.o: export_zonotope.cc ../../include/geometry/zonotope.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../linear_algebra/../numeric/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/../utility/stlio.h \
  ../../include/geometry/../linear_algebra/../numeric/../utility/../base/array.h \
  ../../include/geometry/../linear_algebra/../numeric/../utility/../base/../declarations.h \
  ../../include/geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../linear_algebra/transformation_system.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/interval_vector.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../numeric/../declarations.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/interval_matrix.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../numeric/interval.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/rectangle.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../numeric/arithmetic.h \
  ../../include/geometry/../numeric/../numeric/numerical_types.h \
  ../../include/geometry/../numeric/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/../numeric/integer.h \
  ../../include/geometry/../numeric/../numeric/../numeric/dyadic.h \
  ../../include/geometry/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../numeric/../numeric/../utility/stlio.h \
  ../../include/geometry/../numeric/../numeric/approximation.h \
  ../../include/geometry/../numeric/../numeric/../numeric/numerical_types.h \
  ../../include/geometry/../linear_algebra/interval_vector.h \
  ../../include/geometry/parallelotope.h \
  ../../include/geometry/../geometry/polyhedron.h \
  /ufs/collins/include/ppl.hh \
  ../../include/geometry/../geometry/../linear_algebra/constraint_system.h \
  ../../include/geometry/../geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/numerical_types.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../linear_algebra/interval_vector.h \
  ../../include/geometry/../geometry/../linear_algebra/generator_system.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../geometry/rectangle.h \
  ../../include/geometry/polyhedron.h ../../include/geometry/list_set.h \
  ../../include/geometry/../base/utility.h \
  ../../include/geometry/../base/../declarations.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../numeric/numerical_types.h \
  ../../include/python/../declarations.h \
  ../../include/python/python_utilities.h

export_simplex.o: export_simplex.cc ../../include/geometry/simplex.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../base/array.h \
  ../../include/geometry/../base/../declarations.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/../declarations.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/../utility/stlio.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/../utility/../base/array.h \
  ../../include/geometry/rectangle.h \
  ../../include/geometry/../numeric/arithmetic.h \
  ../../include/geometry/../numeric/../numeric/numerical_types.h \
  ../../include/geometry/../numeric/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/../numeric/integer.h \
  ../../include/geometry/../numeric/../numeric/../numeric/dyadic.h \
  ../../include/geometry/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../numeric/../numeric/../utility/stlio.h \
  ../../include/geometry/../numeric/../numeric/approximation.h \
  ../../include/geometry/../numeric/../numeric/../numeric/numerical_types.h \
  ../../include/geometry/../numeric/interval.h \
  ../../include/geometry/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../linear_algebra/interval_vector.h \
  ../../include/geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/polyhedron.h /ufs/collins/include/ppl.hh \
  ../../include/geometry/../linear_algebra/constraint_system.h \
  ../../include/geometry/../linear_algebra/../numeric/numerical_types.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../linear_algebra/generator_system.h \
  ../../include/geometry/../geometry/rectangle.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../numeric/numerical_types.h \
  ../../include/python/../declarations.h

export_polyhedron.o: export_polyhedron.cc ../../include/utility/stlio.h \
  ../../include/utility/../base/array.h \
  ../../include/utility/../base/../declarations.h \
  ../../include/geometry/rectangle.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../numeric/arithmetic.h \
  ../../include/geometry/../numeric/../numeric/numerical_types.h \
  ../../include/geometry/../numeric/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/../numeric/integer.h \
  ../../include/geometry/../numeric/../numeric/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/../numeric/../utility/stlio.h \
  ../../include/geometry/../numeric/../numeric/../numeric/dyadic.h \
  ../../include/geometry/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../numeric/../numeric/../utility/stlio.h \
  ../../include/geometry/../numeric/../numeric/approximation.h \
  ../../include/geometry/../numeric/../numeric/../numeric/numerical_types.h \
  ../../include/geometry/../numeric/interval.h \
  ../../include/geometry/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../linear_algebra/interval_vector.h \
  ../../include/geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/parallelotope.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../geometry/polyhedron.h \
  /ufs/collins/include/ppl.hh \
  ../../include/geometry/../geometry/../linear_algebra/constraint_system.h \
  ../../include/geometry/../geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/numerical_types.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../linear_algebra/interval_vector.h \
  ../../include/geometry/../geometry/../linear_algebra/generator_system.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../geometry/rectangle.h \
  ../../include/geometry/polyhedron.h ../../include/python/typedefs.h \
  ../../include/python/../real_typedef.h \
  ../../include/python/../numeric/numerical_types.h \
  ../../include/python/../declarations.h

export_lattice_set.o: export_lattice_set.cc \
  ../../include/geometry/lattice_set.h \
  ../../include/geometry/../base/array.h \
  ../../include/geometry/../base/../declarations.h \
  ../../include/geometry/../base/iterator.h \
  ../../include/geometry/../base/../base/array.h \
  ../../include/geometry/../numeric/interval.h \
  ../../include/geometry/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../base/array_operations.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../numeric/numerical_types.h \
  ../../include/python/../numeric/../declarations.h \
  ../../include/python/../numeric/../numeric/integer.h \
  ../../include/python/../numeric/../numeric/../declarations.h \
  ../../include/python/../numeric/../numeric/../utility/stlio.h \
  ../../include/python/../numeric/../numeric/../utility/../base/array.h \
  ../../include/python/../numeric/../numeric/dyadic.h \
  ../../include/python/../numeric/../numeric/numerical_traits.h \
  ../../include/python/../numeric/../utility/stlio.h \
  ../../include/python/../declarations.h

export_list_set.o: export_list_set.cc ../../include/geometry/rectangle.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../numeric/arithmetic.h \
  ../../include/geometry/../numeric/../numeric/numerical_types.h \
  ../../include/geometry/../numeric/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/../numeric/integer.h \
  ../../include/geometry/../numeric/../numeric/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/../numeric/../utility/stlio.h \
  ../../include/geometry/../numeric/../numeric/../numeric/../utility/../base/array.h \
  ../../include/geometry/../numeric/../numeric/../numeric/../utility/../base/../declarations.h \
  ../../include/geometry/../numeric/../numeric/../numeric/dyadic.h \
  ../../include/geometry/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../numeric/../numeric/../utility/stlio.h \
  ../../include/geometry/../numeric/../numeric/approximation.h \
  ../../include/geometry/../numeric/../numeric/../numeric/numerical_types.h \
  ../../include/geometry/../numeric/interval.h \
  ../../include/geometry/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../linear_algebra/interval_vector.h \
  ../../include/geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/parallelotope.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../geometry/polyhedron.h \
  /ufs/collins/include/ppl.hh \
  ../../include/geometry/../geometry/../linear_algebra/constraint_system.h \
  ../../include/geometry/../geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/numerical_types.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../linear_algebra/interval_vector.h \
  ../../include/geometry/../geometry/../linear_algebra/generator_system.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../geometry/rectangle.h \
  ../../include/geometry/zonotope.h \
  ../../include/geometry/../linear_algebra/transformation_system.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/interval_vector.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/interval_matrix.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/polyhedron.h ../../include/geometry/list_set.h \
  ../../include/geometry/../base/utility.h \
  ../../include/geometry/../base/../declarations.h \
  ../../include/geometry/grid_set.h \
  ../../include/geometry/../base/array.h \
  ../../include/geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/grid.h \
  ../../include/geometry/../geometry/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../geometry/lattice_set.h \
  ../../include/geometry/../geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../geometry/../base/iterator.h \
  ../../include/geometry/../geometry/../geometry/../base/../base/array.h \
  ../../include/geometry/../geometry/../geometry/../numeric/interval.h \
  ../../include/geometry/../geometry/../geometry/../base/array_operations.h \
  ../../include/geometry/../geometry/../geometry/../base/../declarations.h \
  ../../include/geometry/../geometry/lattice_set.h \
  ../../include/geometry/partition_tree_set.h \
  ../../include/geometry/../base/binary_word.h \
  ../../include/geometry/../base/binary_tree.h \
  ../../include/geometry/../base/../base/array.h \
  ../../include/geometry/../base/../base/binary_word.h \
  ../../include/geometry/../base/iterator.h \
  ../../include/geometry/../geometry/subdivision_tree_set.h \
  ../../include/geometry/../geometry/../base/sequence.h \
  ../../include/geometry/../geometry/../base/binary_word.h \
  ../../include/geometry/../geometry/../base/binary_tree.h \
  ../../include/geometry/../geometry/../base/iterator.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../numeric/numerical_types.h \
  ../../include/python/../declarations.h \
  ../../include/python/python_utilities.h

export_grid_set.o: export_grid_set.cc \
  ../../include/linear_algebra/linear_algebra.h \
  ../../include/linear_algebra/../linear_algebra/vector.h \
  ../../include/linear_algebra/../linear_algebra/../declarations.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/integer.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/../declarations.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/../utility/stlio.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/../utility/../base/array.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/../utility/../base/../declarations.h \
  ../../include/linear_algebra/../linear_algebra/matrix.h \
  ../../include/linear_algebra/../linear_algebra/constraint_system.h \
  /ufs/collins/include/ppl.hh \
  ../../include/linear_algebra/../linear_algebra/../numeric/numerical_types.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/../numeric/integer.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/../numeric/dyadic.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/linear_algebra/../linear_algebra/../linear_algebra/vector.h \
  ../../include/linear_algebra/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/polyhedron.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../linear_algebra/constraint_system.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../linear_algebra/interval_vector.h \
  ../../include/geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../linear_algebra/../numeric/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/generator_system.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../numeric/../numeric/numerical_types.h \
  ../../include/geometry/../geometry/../numeric/../numeric/approximation.h \
  ../../include/geometry/../geometry/../numeric/../numeric/../declarations.h \
  ../../include/geometry/../geometry/../numeric/../numeric/../numeric/numerical_types.h \
  ../../include/geometry/../geometry/../numeric/interval.h \
  ../../include/geometry/../geometry/../linear_algebra/interval_vector.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/parallelotope.h \
  ../../include/geometry/../numeric/interval.h \
  ../../include/geometry/../geometry/polyhedron.h \
  ../../include/geometry/zonotope.h \
  ../../include/geometry/../linear_algebra/transformation_system.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/interval_vector.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/interval_matrix.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/rectangle.h ../../include/geometry/list_set.h \
  ../../include/geometry/../base/utility.h \
  ../../include/geometry/../base/../declarations.h \
  ../../include/geometry/grid_set.h \
  ../../include/geometry/../base/array.h \
  ../../include/geometry/../geometry/grid.h \
  ../../include/geometry/../geometry/../geometry/lattice_set.h \
  ../../include/geometry/../geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../geometry/../base/iterator.h \
  ../../include/geometry/../geometry/../geometry/../base/../base/array.h \
  ../../include/geometry/../geometry/../geometry/../numeric/interval.h \
  ../../include/geometry/../geometry/../geometry/../base/array_operations.h \
  ../../include/geometry/../geometry/../geometry/../base/../declarations.h \
  ../../include/geometry/../geometry/lattice_set.h \
  ../../include/geometry/partition_tree_set.h \
  ../../include/geometry/../base/binary_word.h \
  ../../include/geometry/../base/binary_tree.h \
  ../../include/geometry/../base/../base/array.h \
  ../../include/geometry/../base/../base/binary_word.h \
  ../../include/geometry/../base/iterator.h \
  ../../include/geometry/../geometry/subdivision_tree_set.h \
  ../../include/geometry/../geometry/../base/sequence.h \
  ../../include/geometry/../geometry/../base/binary_word.h \
  ../../include/geometry/../geometry/../base/binary_tree.h \
  ../../include/geometry/../geometry/../base/iterator.h \
  ../../include/geometry/../geometry/../geometry/rectangle.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../numeric/numerical_types.h \
  ../../include/python/../declarations.h \
  ../../include/python/python_utilities.h

export_partition_tree_set.o: export_partition_tree_set.cc \
  ../../include/base/binary_word.h ../../include/base/../declarations.h \
  ../../include/base/binary_tree.h ../../include/base/../base/array.h \
  ../../include/base/../base/../declarations.h \
  ../../include/base/../base/binary_word.h \
  ../../include/geometry/rectangle.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../numeric/arithmetic.h \
  ../../include/geometry/../numeric/../numeric/numerical_types.h \
  ../../include/geometry/../numeric/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/../numeric/integer.h \
  ../../include/geometry/../numeric/../numeric/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/../numeric/../utility/stlio.h \
  ../../include/geometry/../numeric/../numeric/../numeric/../utility/../base/array.h \
  ../../include/geometry/../numeric/../numeric/../numeric/dyadic.h \
  ../../include/geometry/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../numeric/../numeric/../utility/stlio.h \
  ../../include/geometry/../numeric/../numeric/approximation.h \
  ../../include/geometry/../numeric/../numeric/../numeric/numerical_types.h \
  ../../include/geometry/../numeric/interval.h \
  ../../include/geometry/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../linear_algebra/interval_vector.h \
  ../../include/geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/grid_set.h \
  ../../include/geometry/../base/array.h \
  ../../include/geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/grid.h \
  ../../include/geometry/../geometry/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../geometry/lattice_set.h \
  ../../include/geometry/../geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../geometry/../base/iterator.h \
  ../../include/geometry/../geometry/../geometry/../base/../base/array.h \
  ../../include/geometry/../geometry/../geometry/../numeric/interval.h \
  ../../include/geometry/../geometry/../geometry/../base/array_operations.h \
  ../../include/geometry/../geometry/../geometry/../base/../declarations.h \
  ../../include/geometry/../geometry/lattice_set.h \
  ../../include/geometry/partition_tree_set.h \
  ../../include/geometry/../base/binary_word.h \
  ../../include/geometry/../base/binary_tree.h \
  ../../include/geometry/../base/iterator.h \
  ../../include/geometry/../geometry/subdivision_tree_set.h \
  ../../include/geometry/../geometry/../base/sequence.h \
  ../../include/geometry/../geometry/../base/binary_word.h \
  ../../include/geometry/../geometry/../base/binary_tree.h \
  ../../include/geometry/../geometry/../base/iterator.h \
  ../../include/geometry/../geometry/../geometry/rectangle.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../numeric/numerical_types.h \
  ../../include/python/../declarations.h

export_postscript_output.o: export_postscript_output.cc \
  ../../include/geometry/rectangle.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../numeric/arithmetic.h \
  ../../include/geometry/../numeric/../numeric/numerical_types.h \
  ../../include/geometry/../numeric/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/../numeric/integer.h \
  ../../include/geometry/../numeric/../numeric/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/../numeric/../utility/stlio.h \
  ../../include/geometry/../numeric/../numeric/../numeric/../utility/../base/array.h \
  ../../include/geometry/../numeric/../numeric/../numeric/../utility/../base/../declarations.h \
  ../../include/geometry/../numeric/../numeric/../numeric/dyadic.h \
  ../../include/geometry/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../numeric/../numeric/../utility/stlio.h \
  ../../include/geometry/../numeric/../numeric/approximation.h \
  ../../include/geometry/../numeric/../numeric/../numeric/numerical_types.h \
  ../../include/geometry/../numeric/interval.h \
  ../../include/geometry/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../linear_algebra/interval_vector.h \
  ../../include/geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/parallelotope.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../geometry/polyhedron.h \
  /ufs/collins/include/ppl.hh \
  ../../include/geometry/../geometry/../linear_algebra/constraint_system.h \
  ../../include/geometry/../geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/numerical_types.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../linear_algebra/interval_vector.h \
  ../../include/geometry/../geometry/../linear_algebra/generator_system.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../geometry/rectangle.h \
  ../../include/geometry/zonotope.h \
  ../../include/geometry/../linear_algebra/transformation_system.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/interval_vector.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/interval_matrix.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/polyhedron.h ../../include/geometry/list_set.h \
  ../../include/geometry/../base/utility.h \
  ../../include/geometry/../base/../declarations.h \
  ../../include/geometry/grid_set.h \
  ../../include/geometry/../base/array.h \
  ../../include/geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/grid.h \
  ../../include/geometry/../geometry/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../geometry/lattice_set.h \
  ../../include/geometry/../geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../geometry/../base/iterator.h \
  ../../include/geometry/../geometry/../geometry/../base/../base/array.h \
  ../../include/geometry/../geometry/../geometry/../numeric/interval.h \
  ../../include/geometry/../geometry/../geometry/../base/array_operations.h \
  ../../include/geometry/../geometry/../geometry/../base/../declarations.h \
  ../../include/geometry/../geometry/lattice_set.h \
  ../../include/geometry/partition_tree_set.h \
  ../../include/geometry/../base/binary_word.h \
  ../../include/geometry/../base/binary_tree.h \
  ../../include/geometry/../base/../base/array.h \
  ../../include/geometry/../base/../base/binary_word.h \
  ../../include/geometry/../base/iterator.h \
  ../../include/geometry/../geometry/subdivision_tree_set.h \
  ../../include/geometry/../geometry/../base/sequence.h \
  ../../include/geometry/../geometry/../base/binary_word.h \
  ../../include/geometry/../geometry/../base/binary_tree.h \
  ../../include/geometry/../geometry/../base/iterator.h \
  ../../include/utility/epsfstream.h \
  ../../include/utility/../geometry/point.h \
  ../../include/utility/../geometry/rectangle.h \
  ../../include/utility/../geometry/list_set.h \
  ../../include/utility/../geometry/grid_set.h \
  ../../include/utility/../geometry/parallelotope.h \
  ../../include/utility/../geometry/zonotope.h \
  ../../include/utility/../geometry/polyhedron.h \
  ../../include/utility/../geometry/partition_tree_set.h \
  ../../include/utility/../linear_algebra/matrix.h \
  ../../include/utility/../system/affine_map.h \
  ../../include/utility/../system/../numeric/numerical_types.h \
  ../../include/utility/../system/../linear_algebra/vector.h \
  ../../include/utility/../system/../linear_algebra/matrix.h \
  ../../include/utility/../system/../linear_algebra/interval_matrix.h \
  ../../include/utility/../system/../geometry/point.h \
  ../../include/utility/../system/../geometry/simplex.h \
  ../../include/utility/../system/../geometry/../declarations.h \
  ../../include/utility/../system/../geometry/../base/array.h \
  ../../include/utility/../system/../geometry/../geometry/point.h \
  ../../include/utility/../system/../geometry/rectangle.h \
  ../../include/utility/../system/../geometry/parallelotope.h \
  ../../include/utility/../system/../geometry/zonotope.h \
  ../../include/utility/../system/../geometry/polyhedron.h \
  ../../include/utility/../system/../geometry/list_set.h \
  ../../include/utility/../system/../geometry/grid_set.h \
  ../../include/utility/../system/../system/map.h \
  ../../include/utility/../system/../system/../declarations.h \
  ../../include/utility/../system/../system/../linear_algebra/vector.h \
  ../../include/utility/../system/../system/../linear_algebra/matrix.h \
  ../../include/utility/../system/../system/../linear_algebra/interval_matrix.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../numeric/numerical_types.h \
  ../../include/python/../declarations.h

evaluation_module.o: evaluation_module.cc

export_apply.o: export_apply.cc ../../include/geometry/rectangle.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../numeric/arithmetic.h \
  ../../include/geometry/../numeric/../numeric/numerical_types.h \
  ../../include/geometry/../numeric/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/../numeric/integer.h \
  ../../include/geometry/../numeric/../numeric/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/../numeric/../utility/stlio.h \
  ../../include/geometry/../numeric/../numeric/../numeric/../utility/../base/array.h \
  ../../include/geometry/../numeric/../numeric/../numeric/../utility/../base/../declarations.h \
  ../../include/geometry/../numeric/../numeric/../numeric/dyadic.h \
  ../../include/geometry/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../numeric/../numeric/../utility/stlio.h \
  ../../include/geometry/../numeric/../numeric/approximation.h \
  ../../include/geometry/../numeric/../numeric/../numeric/numerical_types.h \
  ../../include/geometry/../numeric/interval.h \
  ../../include/geometry/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../linear_algebra/interval_vector.h \
  ../../include/geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/parallelotope.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../geometry/polyhedron.h \
  /ufs/collins/include/ppl.hh \
  ../../include/geometry/../geometry/../linear_algebra/constraint_system.h \
  ../../include/geometry/../geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/numerical_types.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../linear_algebra/interval_vector.h \
  ../../include/geometry/../geometry/../linear_algebra/generator_system.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../geometry/rectangle.h \
  ../../include/geometry/list_set.h \
  ../../include/geometry/../base/utility.h \
  ../../include/geometry/../base/../declarations.h \
  ../../include/geometry/grid_set.h \
  ../../include/geometry/../base/array.h \
  ../../include/geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/grid.h \
  ../../include/geometry/../geometry/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../geometry/lattice_set.h \
  ../../include/geometry/../geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../geometry/../base/iterator.h \
  ../../include/geometry/../geometry/../geometry/../base/../base/array.h \
  ../../include/geometry/../geometry/../geometry/../numeric/interval.h \
  ../../include/geometry/../geometry/../geometry/../base/array_operations.h \
  ../../include/geometry/../geometry/../geometry/../base/../declarations.h \
  ../../include/geometry/../geometry/lattice_set.h \
  ../../include/evaluation/apply.h \
  ../../include/evaluation/../declarations.h \
  ../../include/evaluation/../system/map.h \
  ../../include/evaluation/../system/../declarations.h \
  ../../include/evaluation/../system/../linear_algebra/vector.h \
  ../../include/evaluation/../system/../linear_algebra/matrix.h \
  ../../include/evaluation/../system/../linear_algebra/interval_matrix.h \
  ../../include/evaluation/../system/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/evaluation/../system/../linear_algebra/../declarations.h \
  ../../include/evaluation/../system/../linear_algebra/../numeric/interval.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../numeric/numerical_types.h \
  ../../include/python/../declarations.h

export_integrate.o: export_integrate.cc \
  ../../include/geometry/rectangle.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../numeric/arithmetic.h \
  ../../include/geometry/../numeric/../numeric/numerical_types.h \
  ../../include/geometry/../numeric/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/../numeric/integer.h \
  ../../include/geometry/../numeric/../numeric/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/../numeric/../utility/stlio.h \
  ../../include/geometry/../numeric/../numeric/../numeric/../utility/../base/array.h \
  ../../include/geometry/../numeric/../numeric/../numeric/../utility/../base/../declarations.h \
  ../../include/geometry/../numeric/../numeric/../numeric/dyadic.h \
  ../../include/geometry/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../numeric/../numeric/../utility/stlio.h \
  ../../include/geometry/../numeric/../numeric/approximation.h \
  ../../include/geometry/../numeric/../numeric/../numeric/numerical_types.h \
  ../../include/geometry/../numeric/interval.h \
  ../../include/geometry/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../linear_algebra/interval_vector.h \
  ../../include/geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/parallelotope.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../geometry/polyhedron.h \
  /ufs/collins/include/ppl.hh \
  ../../include/geometry/../geometry/../linear_algebra/constraint_system.h \
  ../../include/geometry/../geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/numerical_types.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../linear_algebra/interval_vector.h \
  ../../include/geometry/../geometry/../linear_algebra/generator_system.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../geometry/rectangle.h \
  ../../include/geometry/zonotope.h \
  ../../include/geometry/../linear_algebra/transformation_system.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/interval_vector.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/interval_matrix.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/grid_set.h \
  ../../include/geometry/../base/array.h \
  ../../include/geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/grid.h \
  ../../include/geometry/../geometry/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../geometry/lattice_set.h \
  ../../include/geometry/../geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../geometry/../base/iterator.h \
  ../../include/geometry/../geometry/../geometry/../base/../base/array.h \
  ../../include/geometry/../geometry/../geometry/../numeric/interval.h \
  ../../include/geometry/../geometry/../geometry/../base/array_operations.h \
  ../../include/geometry/../geometry/../geometry/../base/../declarations.h \
  ../../include/geometry/../geometry/lattice_set.h \
  ../../include/geometry/list_set.h \
  ../../include/geometry/../base/utility.h \
  ../../include/geometry/../base/../declarations.h \
  ../../include/system/vector_field.h \
  ../../include/system/../declarations.h \
  ../../include/system/../numeric/interval.h \
  ../../include/system/../linear_algebra/vector.h \
  ../../include/system/../linear_algebra/matrix.h \
  ../../include/system/../linear_algebra/interval_vector.h \
  ../../include/system/../linear_algebra/interval_matrix.h \
  ../../include/evaluation/integrator.h \
  ../../include/evaluation/../declarations.h \
  ../../include/evaluation/lohner_integrator.h \
  ../../include/evaluation/../evaluation/integrator.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../numeric/numerical_types.h \
  ../../include/python/../declarations.h

system_module.o: system_module.cc

export_lattice_map.o: export_lattice_map.cc \
  ../../include/system/lattice_map.h \
  ../../include/system/../geometry/lattice_set.h \
  ../../include/system/../geometry/../base/array.h \
  ../../include/system/../geometry/../base/../declarations.h \
  ../../include/system/../geometry/../base/iterator.h \
  ../../include/system/../geometry/../base/../base/array.h \
  ../../include/system/../geometry/../numeric/interval.h \
  ../../include/system/../geometry/../numeric/../declarations.h \
  ../../include/system/../geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/system/../geometry/../base/array_operations.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../numeric/numerical_types.h \
  ../../include/python/../numeric/../declarations.h \
  ../../include/python/../numeric/../numeric/integer.h \
  ../../include/python/../numeric/../numeric/../declarations.h \
  ../../include/python/../numeric/../numeric/../utility/stlio.h \
  ../../include/python/../numeric/../numeric/../utility/../base/array.h \
  ../../include/python/../numeric/../numeric/dyadic.h \
  ../../include/python/../numeric/../numeric/numerical_traits.h \
  ../../include/python/../numeric/../utility/stlio.h \
  ../../include/python/../declarations.h

export_map.o: export_map.cc ../../include/system/map.h \
  ../../include/system/../declarations.h \
  ../../include/system/../linear_algebra/vector.h \
  ../../include/system/../linear_algebra/../declarations.h \
  ../../include/system/../linear_algebra/../numeric/integer.h \
  ../../include/system/../linear_algebra/../numeric/../declarations.h \
  ../../include/system/../linear_algebra/../numeric/../utility/stlio.h \
  ../../include/system/../linear_algebra/../numeric/../utility/../base/array.h \
  ../../include/system/../linear_algebra/../numeric/../utility/../base/../declarations.h \
  ../../include/system/../linear_algebra/matrix.h \
  ../../include/system/../linear_algebra/interval_matrix.h \
  ../../include/system/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/system/../linear_algebra/../numeric/interval.h \
  ../../include/system/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../numeric/numerical_types.h \
  ../../include/python/../numeric/../declarations.h \
  ../../include/python/../numeric/../numeric/integer.h \
  ../../include/python/../numeric/../numeric/dyadic.h \
  ../../include/python/../numeric/../numeric/numerical_traits.h \
  ../../include/python/../numeric/../utility/stlio.h \
  ../../include/python/../declarations.h

export_affine_map.o: export_affine_map.cc ../../include/geometry/point.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../base/array.h \
  ../../include/geometry/../base/../declarations.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../linear_algebra/../numeric/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/../utility/stlio.h \
  ../../include/geometry/../linear_algebra/../numeric/../utility/../base/array.h \
  ../../include/geometry/rectangle.h \
  ../../include/geometry/../numeric/arithmetic.h \
  ../../include/geometry/../numeric/../numeric/numerical_types.h \
  ../../include/geometry/../numeric/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/../numeric/integer.h \
  ../../include/geometry/../numeric/../numeric/../numeric/dyadic.h \
  ../../include/geometry/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../numeric/../numeric/../utility/stlio.h \
  ../../include/geometry/../numeric/../numeric/approximation.h \
  ../../include/geometry/../numeric/../numeric/../numeric/numerical_types.h \
  ../../include/geometry/../numeric/interval.h \
  ../../include/geometry/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../linear_algebra/interval_vector.h \
  ../../include/geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/parallelotope.h \
  ../../include/geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/polyhedron.h \
  /ufs/collins/include/ppl.hh \
  ../../include/geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../linear_algebra/constraint_system.h \
  ../../include/geometry/../geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/numerical_types.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../linear_algebra/interval_vector.h \
  ../../include/geometry/../geometry/../linear_algebra/generator_system.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../geometry/rectangle.h \
  ../../include/geometry/simplex.h ../../include/geometry/zonotope.h \
  ../../include/geometry/../linear_algebra/transformation_system.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/interval_vector.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/interval_matrix.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/polyhedron.h ../../include/system/affine_map.h \
  ../../include/system/../numeric/numerical_types.h \
  ../../include/system/../linear_algebra/vector.h \
  ../../include/system/../linear_algebra/matrix.h \
  ../../include/system/../linear_algebra/interval_matrix.h \
  ../../include/system/../geometry/point.h \
  ../../include/system/../geometry/simplex.h \
  ../../include/system/../geometry/rectangle.h \
  ../../include/system/../geometry/parallelotope.h \
  ../../include/system/../geometry/zonotope.h \
  ../../include/system/../geometry/polyhedron.h \
  ../../include/system/../geometry/list_set.h \
  ../../include/system/../geometry/../declarations.h \
  ../../include/system/../geometry/../base/utility.h \
  ../../include/system/../geometry/../base/../declarations.h \
  ../../include/system/../geometry/grid_set.h \
  ../../include/system/../geometry/../base/array.h \
  ../../include/system/../geometry/../geometry/rectangle.h \
  ../../include/system/../geometry/../geometry/grid.h \
  ../../include/system/../geometry/../geometry/../declarations.h \
  ../../include/system/../geometry/../geometry/../base/array.h \
  ../../include/system/../geometry/../geometry/../numeric/arithmetic.h \
  ../../include/system/../geometry/../geometry/../geometry/lattice_set.h \
  ../../include/system/../geometry/../geometry/../geometry/../base/array.h \
  ../../include/system/../geometry/../geometry/../geometry/../base/iterator.h \
  ../../include/system/../geometry/../geometry/../geometry/../base/../base/array.h \
  ../../include/system/../geometry/../geometry/../geometry/../numeric/interval.h \
  ../../include/system/../geometry/../geometry/../geometry/../base/array_operations.h \
  ../../include/system/../geometry/../geometry/../geometry/../base/../declarations.h \
  ../../include/system/../geometry/../geometry/lattice_set.h \
  ../../include/system/../system/map.h \
  ../../include/system/../system/../declarations.h \
  ../../include/system/../system/../linear_algebra/vector.h \
  ../../include/system/../system/../linear_algebra/matrix.h \
  ../../include/system/../system/../linear_algebra/interval_matrix.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../numeric/numerical_types.h \
  ../../include/python/../declarations.h \
  ../../include/python/python_utilities.h

export_affine_multimap.o: export_affine_multimap.cc \
  ../../include/geometry/point.h ../../include/geometry/../declarations.h \
  ../../include/geometry/../base/array.h \
  ../../include/geometry/../base/../declarations.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../linear_algebra/../numeric/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/../utility/stlio.h \
  ../../include/geometry/../linear_algebra/../numeric/../utility/../base/array.h \
  ../../include/geometry/rectangle.h \
  ../../include/geometry/../numeric/arithmetic.h \
  ../../include/geometry/../numeric/../numeric/numerical_types.h \
  ../../include/geometry/../numeric/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/../numeric/integer.h \
  ../../include/geometry/../numeric/../numeric/../numeric/dyadic.h \
  ../../include/geometry/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../numeric/../numeric/../utility/stlio.h \
  ../../include/geometry/../numeric/../numeric/approximation.h \
  ../../include/geometry/../numeric/../numeric/../numeric/numerical_types.h \
  ../../include/geometry/../numeric/interval.h \
  ../../include/geometry/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../linear_algebra/interval_vector.h \
  ../../include/geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/parallelotope.h \
  ../../include/geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/polyhedron.h \
  /ufs/collins/include/ppl.hh \
  ../../include/geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../linear_algebra/constraint_system.h \
  ../../include/geometry/../geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/numerical_types.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../linear_algebra/interval_vector.h \
  ../../include/geometry/../geometry/../linear_algebra/generator_system.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../geometry/rectangle.h \
  ../../include/geometry/simplex.h ../../include/geometry/zonotope.h \
  ../../include/geometry/../linear_algebra/transformation_system.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/interval_vector.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/interval_matrix.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/polyhedron.h \
  ../../include/system/affine_multimap.h \
  ../../include/system/../linear_algebra/vector.h \
  ../../include/system/../linear_algebra/matrix.h \
  ../../include/system/../linear_algebra/interval_matrix.h \
  ../../include/system/../geometry/point.h \
  ../../include/system/../geometry/simplex.h \
  ../../include/system/../geometry/rectangle.h \
  ../../include/system/../geometry/parallelotope.h \
  ../../include/system/../geometry/zonotope.h \
  ../../include/system/../geometry/polyhedron.h \
  ../../include/system/../geometry/list_set.h \
  ../../include/system/../geometry/../declarations.h \
  ../../include/system/../geometry/../base/utility.h \
  ../../include/system/../geometry/../base/../declarations.h \
  ../../include/system/../geometry/grid_set.h \
  ../../include/system/../geometry/../base/array.h \
  ../../include/system/../geometry/../geometry/rectangle.h \
  ../../include/system/../geometry/../geometry/grid.h \
  ../../include/system/../geometry/../geometry/../declarations.h \
  ../../include/system/../geometry/../geometry/../base/array.h \
  ../../include/system/../geometry/../geometry/../numeric/arithmetic.h \
  ../../include/system/../geometry/../geometry/../geometry/lattice_set.h \
  ../../include/system/../geometry/../geometry/../geometry/../base/array.h \
  ../../include/system/../geometry/../geometry/../geometry/../base/iterator.h \
  ../../include/system/../geometry/../geometry/../geometry/../base/../base/array.h \
  ../../include/system/../geometry/../geometry/../geometry/../numeric/interval.h \
  ../../include/system/../geometry/../geometry/../geometry/../base/array_operations.h \
  ../../include/system/../geometry/../geometry/../geometry/../base/../declarations.h \
  ../../include/system/../geometry/../geometry/lattice_set.h \
  ../../include/system/../system/affine_map.h \
  ../../include/system/../system/../numeric/numerical_types.h \
  ../../include/system/../system/../linear_algebra/vector.h \
  ../../include/system/../system/../linear_algebra/matrix.h \
  ../../include/system/../system/../linear_algebra/interval_matrix.h \
  ../../include/system/../system/../geometry/point.h \
  ../../include/system/../system/../geometry/simplex.h \
  ../../include/system/../system/../geometry/rectangle.h \
  ../../include/system/../system/../geometry/parallelotope.h \
  ../../include/system/../system/../geometry/zonotope.h \
  ../../include/system/../system/../geometry/polyhedron.h \
  ../../include/system/../system/../geometry/list_set.h \
  ../../include/system/../system/../geometry/grid_set.h \
  ../../include/system/../system/../system/map.h \
  ../../include/system/../system/../system/../declarations.h \
  ../../include/system/../system/../system/../linear_algebra/vector.h \
  ../../include/system/../system/../system/../linear_algebra/matrix.h \
  ../../include/system/../system/../system/../linear_algebra/interval_matrix.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../numeric/numerical_types.h \
  ../../include/python/../declarations.h \
  ../../include/python/python_utilities.h

export_polynomial_map.o: export_polynomial_map.cc \
  ../../include/geometry/point.h ../../include/geometry/../declarations.h \
  ../../include/geometry/../base/array.h \
  ../../include/geometry/../base/../declarations.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../linear_algebra/../numeric/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/../utility/stlio.h \
  ../../include/geometry/../linear_algebra/../numeric/../utility/../base/array.h \
  ../../include/geometry/rectangle.h \
  ../../include/geometry/../numeric/arithmetic.h \
  ../../include/geometry/../numeric/../numeric/numerical_types.h \
  ../../include/geometry/../numeric/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/../numeric/integer.h \
  ../../include/geometry/../numeric/../numeric/../numeric/dyadic.h \
  ../../include/geometry/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../numeric/../numeric/../utility/stlio.h \
  ../../include/geometry/../numeric/../numeric/approximation.h \
  ../../include/geometry/../numeric/../numeric/../numeric/numerical_types.h \
  ../../include/geometry/../numeric/interval.h \
  ../../include/geometry/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../linear_algebra/interval_vector.h \
  ../../include/geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/system/polynomial_map.h \
  ../../include/system/../system/map.h \
  ../../include/system/../system/../declarations.h \
  ../../include/system/../system/../linear_algebra/vector.h \
  ../../include/system/../system/../linear_algebra/matrix.h \
  ../../include/system/../system/../linear_algebra/../declarations.h \
  ../../include/system/../system/../linear_algebra/../numeric/integer.h \
  ../../include/system/../system/../linear_algebra/interval_matrix.h \
  ../../include/system/../system/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/system/../system/../linear_algebra/../numeric/interval.h \
  ../../include/system/../linear_algebra/matrix.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../numeric/numerical_types.h \
  ../../include/python/../declarations.h

export_vector_field.o: export_vector_field.cc \
  ../../include/system/vector_field.h \
  ../../include/system/../declarations.h \
  ../../include/system/../numeric/interval.h \
  ../../include/system/../numeric/../declarations.h \
  ../../include/system/../numeric/../numeric/numerical_traits.h \
  ../../include/system/../linear_algebra/vector.h \
  ../../include/system/../linear_algebra/../declarations.h \
  ../../include/system/../linear_algebra/../numeric/integer.h \
  ../../include/system/../linear_algebra/../numeric/../declarations.h \
  ../../include/system/../linear_algebra/../numeric/../utility/stlio.h \
  ../../include/system/../linear_algebra/../numeric/../utility/../base/array.h \
  ../../include/system/../linear_algebra/../numeric/../utility/../base/../declarations.h \
  ../../include/system/../linear_algebra/matrix.h \
  ../../include/system/../linear_algebra/interval_vector.h \
  ../../include/system/../linear_algebra/../numeric/interval.h \
  ../../include/system/../linear_algebra/../linear_algebra/vector.h \
  ../../include/system/../linear_algebra/interval_matrix.h \
  ../../include/system/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../numeric/numerical_types.h \
  ../../include/python/../numeric/../declarations.h \
  ../../include/python/../numeric/../numeric/integer.h \
  ../../include/python/../numeric/../numeric/dyadic.h \
  ../../include/python/../numeric/../numeric/numerical_traits.h \
  ../../include/python/../numeric/../utility/stlio.h \
  ../../include/python/../declarations.h

export_affine_vector_field.o: export_affine_vector_field.cc \
  ../../include/numeric/interval.h \
  ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/numerical_traits.h \
  ../../include/geometry/point.h ../../include/geometry/../declarations.h \
  ../../include/geometry/../base/array.h \
  ../../include/geometry/../base/../declarations.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../linear_algebra/../numeric/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/../utility/stlio.h \
  ../../include/geometry/../linear_algebra/../numeric/../utility/../base/array.h \
  ../../include/geometry/rectangle.h \
  ../../include/geometry/../numeric/arithmetic.h \
  ../../include/geometry/../numeric/../numeric/numerical_types.h \
  ../../include/geometry/../numeric/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/../numeric/integer.h \
  ../../include/geometry/../numeric/../numeric/../numeric/dyadic.h \
  ../../include/geometry/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../numeric/../numeric/../utility/stlio.h \
  ../../include/geometry/../numeric/../numeric/approximation.h \
  ../../include/geometry/../numeric/../numeric/../numeric/numerical_types.h \
  ../../include/geometry/../numeric/interval.h \
  ../../include/geometry/../linear_algebra/interval_vector.h \
  ../../include/geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/parallelotope.h \
  ../../include/geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/polyhedron.h \
  /ufs/collins/include/ppl.hh \
  ../../include/geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../linear_algebra/constraint_system.h \
  ../../include/geometry/../geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/numerical_types.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../linear_algebra/interval_vector.h \
  ../../include/geometry/../geometry/../linear_algebra/generator_system.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../geometry/rectangle.h \
  ../../include/geometry/simplex.h ../../include/geometry/zonotope.h \
  ../../include/geometry/../linear_algebra/transformation_system.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/interval_vector.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/interval_matrix.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/polyhedron.h \
  ../../include/system/affine_vector_field.h \
  ../../include/system/../linear_algebra/vector.h \
  ../../include/system/../linear_algebra/matrix.h \
  ../../include/system/../linear_algebra/interval_vector.h \
  ../../include/system/../linear_algebra/interval_matrix.h \
  ../../include/system/../system/vector_field.h \
  ../../include/system/../system/../declarations.h \
  ../../include/system/../system/../numeric/interval.h \
  ../../include/system/../system/../linear_algebra/vector.h \
  ../../include/system/../system/../linear_algebra/matrix.h \
  ../../include/system/../system/../linear_algebra/interval_vector.h \
  ../../include/system/../system/../linear_algebra/interval_matrix.h \
  ../../include/system/../system/affine_map.h \
  ../../include/system/../system/../numeric/numerical_types.h \
  ../../include/system/../system/../geometry/point.h \
  ../../include/system/../system/../geometry/simplex.h \
  ../../include/system/../system/../geometry/rectangle.h \
  ../../include/system/../system/../geometry/parallelotope.h \
  ../../include/system/../system/../geometry/zonotope.h \
  ../../include/system/../system/../geometry/polyhedron.h \
  ../../include/system/../system/../geometry/list_set.h \
  ../../include/system/../system/../geometry/../declarations.h \
  ../../include/system/../system/../geometry/../base/utility.h \
  ../../include/system/../system/../geometry/../base/../declarations.h \
  ../../include/system/../system/../geometry/grid_set.h \
  ../../include/system/../system/../geometry/../base/array.h \
  ../../include/system/../system/../geometry/../geometry/rectangle.h \
  ../../include/system/../system/../geometry/../geometry/grid.h \
  ../../include/system/../system/../geometry/../geometry/../declarations.h \
  ../../include/system/../system/../geometry/../geometry/../base/array.h \
  ../../include/system/../system/../geometry/../geometry/../numeric/arithmetic.h \
  ../../include/system/../system/../geometry/../geometry/../geometry/lattice_set.h \
  ../../include/system/../system/../geometry/../geometry/../geometry/../base/array.h \
  ../../include/system/../system/../geometry/../geometry/../geometry/../base/iterator.h \
  ../../include/system/../system/../geometry/../geometry/../geometry/../base/../base/array.h \
  ../../include/system/../system/../geometry/../geometry/../geometry/../numeric/interval.h \
  ../../include/system/../system/../geometry/../geometry/../geometry/../base/array_operations.h \
  ../../include/system/../system/../geometry/../geometry/../geometry/../base/../declarations.h \
  ../../include/system/../system/../geometry/../geometry/lattice_set.h \
  ../../include/system/../system/../system/map.h \
  ../../include/system/../system/../system/../declarations.h \
  ../../include/system/../system/../system/../linear_algebra/vector.h \
  ../../include/system/../system/../system/../linear_algebra/matrix.h \
  ../../include/system/../system/../system/../linear_algebra/interval_matrix.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../numeric/numerical_types.h \
  ../../include/python/../declarations.h \
  ../../include/python/python_utilities.h

export_henon_map.o: export_henon_map.cc ../../include/system/henon_map.h \
  ../../include/system/../linear_algebra/matrix.h \
  ../../include/system/../linear_algebra/../declarations.h \
  ../../include/system/../linear_algebra/../numeric/integer.h \
  ../../include/system/../linear_algebra/../numeric/../declarations.h \
  ../../include/system/../linear_algebra/../numeric/../utility/stlio.h \
  ../../include/system/../linear_algebra/../numeric/../utility/../base/array.h \
  ../../include/system/../linear_algebra/../numeric/../utility/../base/../declarations.h \
  ../../include/system/../geometry/point.h \
  ../../include/system/../geometry/../declarations.h \
  ../../include/system/../geometry/../base/array.h \
  ../../include/system/../geometry/../linear_algebra/vector.h \
  ../../include/system/../geometry/../linear_algebra/../declarations.h \
  ../../include/system/../geometry/../linear_algebra/../numeric/integer.h \
  ../../include/system/../geometry/rectangle.h \
  ../../include/system/../geometry/../numeric/arithmetic.h \
  ../../include/system/../geometry/../numeric/../numeric/numerical_types.h \
  ../../include/system/../geometry/../numeric/../numeric/../declarations.h \
  ../../include/system/../geometry/../numeric/../numeric/../numeric/integer.h \
  ../../include/system/../geometry/../numeric/../numeric/../numeric/dyadic.h \
  ../../include/system/../geometry/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/system/../geometry/../numeric/../numeric/../utility/stlio.h \
  ../../include/system/../geometry/../numeric/../numeric/approximation.h \
  ../../include/system/../geometry/../numeric/../numeric/../numeric/numerical_types.h \
  ../../include/system/../geometry/../numeric/interval.h \
  ../../include/system/../geometry/../numeric/../declarations.h \
  ../../include/system/../geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/system/../geometry/../linear_algebra/interval_vector.h \
  ../../include/system/../geometry/../linear_algebra/../numeric/interval.h \
  ../../include/system/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/system/../geometry/../geometry/point.h \
  ../../include/system/../geometry/parallelotope.h \
  ../../include/system/../geometry/../linear_algebra/matrix.h \
  ../../include/system/../geometry/../geometry/polyhedron.h \
  /ufs/collins/include/ppl.hh \
  ../../include/system/../geometry/../geometry/../declarations.h \
  ../../include/system/../geometry/../geometry/../linear_algebra/constraint_system.h \
  ../../include/system/../geometry/../geometry/../linear_algebra/../declarations.h \
  ../../include/system/../geometry/../geometry/../linear_algebra/../numeric/numerical_types.h \
  ../../include/system/../geometry/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/system/../geometry/../geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/system/../geometry/../geometry/../linear_algebra/vector.h \
  ../../include/system/../geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/system/../geometry/../geometry/../linear_algebra/interval_vector.h \
  ../../include/system/../geometry/../geometry/../linear_algebra/generator_system.h \
  ../../include/system/../geometry/../geometry/../geometry/point.h \
  ../../include/system/../geometry/../geometry/../geometry/rectangle.h \
  ../../include/system/../system/map.h \
  ../../include/system/../system/../declarations.h \
  ../../include/system/../system/../linear_algebra/vector.h \
  ../../include/system/../system/../linear_algebra/matrix.h \
  ../../include/system/../system/../linear_algebra/interval_matrix.h \
  ../../include/system/../system/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/system/../system/../linear_algebra/../declarations.h \
  ../../include/system/../system/../linear_algebra/../numeric/interval.h \
  ../../include/system/../evaluation/apply.h \
  ../../include/system/../evaluation/../declarations.h \
  ../../include/system/../evaluation/../system/map.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../numeric/numerical_types.h \
  ../../include/python/../declarations.h

export_lorenz_system.o: export_lorenz_system.cc \
  ../../include/system/lorenz_system.h \
  ../../include/system/../linear_algebra/vector.h \
  ../../include/system/../linear_algebra/../declarations.h \
  ../../include/system/../linear_algebra/../numeric/integer.h \
  ../../include/system/../linear_algebra/../numeric/../declarations.h \
  ../../include/system/../linear_algebra/../numeric/../utility/stlio.h \
  ../../include/system/../linear_algebra/../numeric/../utility/../base/array.h \
  ../../include/system/../linear_algebra/../numeric/../utility/../base/../declarations.h \
  ../../include/system/../linear_algebra/interval_vector.h \
  ../../include/system/../linear_algebra/../numeric/interval.h \
  ../../include/system/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/system/../linear_algebra/../linear_algebra/vector.h \
  ../../include/system/../linear_algebra/matrix.h \
  ../../include/system/../linear_algebra/interval_matrix.h \
  ../../include/system/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/system/../geometry/point.h \
  ../../include/system/../geometry/../declarations.h \
  ../../include/system/../geometry/../base/array.h \
  ../../include/system/../geometry/../linear_algebra/vector.h \
  ../../include/system/../geometry/rectangle.h \
  ../../include/system/../geometry/../numeric/arithmetic.h \
  ../../include/system/../geometry/../numeric/../numeric/numerical_types.h \
  ../../include/system/../geometry/../numeric/../numeric/../declarations.h \
  ../../include/system/../geometry/../numeric/../numeric/../numeric/integer.h \
  ../../include/system/../geometry/../numeric/../numeric/../numeric/dyadic.h \
  ../../include/system/../geometry/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/system/../geometry/../numeric/../numeric/../utility/stlio.h \
  ../../include/system/../geometry/../numeric/../numeric/approximation.h \
  ../../include/system/../geometry/../numeric/../numeric/../numeric/numerical_types.h \
  ../../include/system/../geometry/../numeric/interval.h \
  ../../include/system/../geometry/../linear_algebra/interval_vector.h \
  ../../include/system/../geometry/../geometry/point.h \
  ../../include/system/../geometry/parallelotope.h \
  ../../include/system/../geometry/../linear_algebra/matrix.h \
  ../../include/system/../geometry/../geometry/polyhedron.h \
  /ufs/collins/include/ppl.hh \
  ../../include/system/../geometry/../geometry/../declarations.h \
  ../../include/system/../geometry/../geometry/../linear_algebra/constraint_system.h \
  ../../include/system/../geometry/../geometry/../linear_algebra/../declarations.h \
  ../../include/system/../geometry/../geometry/../linear_algebra/../numeric/numerical_types.h \
  ../../include/system/../geometry/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/system/../geometry/../geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/system/../geometry/../geometry/../linear_algebra/vector.h \
  ../../include/system/../geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/system/../geometry/../geometry/../linear_algebra/interval_vector.h \
  ../../include/system/../geometry/../geometry/../linear_algebra/generator_system.h \
  ../../include/system/../geometry/../geometry/../geometry/point.h \
  ../../include/system/../geometry/../geometry/../geometry/rectangle.h \
  ../../include/system/../system/vector_field.h \
  ../../include/system/../system/../declarations.h \
  ../../include/system/../system/../numeric/interval.h \
  ../../include/system/../system/../linear_algebra/vector.h \
  ../../include/system/../system/../linear_algebra/matrix.h \
  ../../include/system/../system/../linear_algebra/interval_vector.h \
  ../../include/system/../system/../linear_algebra/interval_matrix.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../numeric/numerical_types.h \
  ../../include/python/../declarations.h

linear_algebra_module.o: linear_algebra_module.cc

export_vector.o: export_vector.cc ../../include/numeric/numerical_types.h \
  ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/integer.h \
  ../../include/numeric/../numeric/../declarations.h \
  ../../include/numeric/../numeric/../utility/stlio.h \
  ../../include/numeric/../numeric/../utility/../base/array.h \
  ../../include/numeric/../numeric/../utility/../base/../declarations.h \
  ../../include/numeric/../numeric/dyadic.h \
  ../../include/numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../utility/stlio.h \
  ../../include/numeric/interval.h ../../include/linear_algebra/vector.h \
  ../../include/linear_algebra/../declarations.h \
  ../../include/linear_algebra/../numeric/integer.h \
  ../../include/linear_algebra/interval_vector.h \
  ../../include/linear_algebra/../numeric/interval.h \
  ../../include/linear_algebra/../linear_algebra/vector.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../numeric/numerical_types.h \
  ../../include/python/../declarations.h

export_matrix.o: export_matrix.cc ../../include/numeric/interval.h \
  ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/numerical_traits.h \
  ../../include/numeric/integer.h \
  ../../include/numeric/../utility/stlio.h \
  ../../include/numeric/../utility/../base/array.h \
  ../../include/numeric/../utility/../base/../declarations.h \
  ../../include/linear_algebra/vector.h \
  ../../include/linear_algebra/../declarations.h \
  ../../include/linear_algebra/../numeric/integer.h \
  ../../include/linear_algebra/matrix.h \
  ../../include/linear_algebra/interval_matrix.h \
  ../../include/linear_algebra/../linear_algebra/matrix.h \
  ../../include/linear_algebra/../numeric/interval.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../numeric/numerical_types.h \
  ../../include/python/../numeric/../declarations.h \
  ../../include/python/../numeric/../numeric/integer.h \
  ../../include/python/../numeric/../numeric/dyadic.h \
  ../../include/python/../numeric/../numeric/numerical_traits.h \
  ../../include/python/../numeric/../utility/stlio.h \
  ../../include/python/../declarations.h \
  ../../include/python/python_utilities.h

export_linear_program.o: export_linear_program.cc \
  ../../include/linear_algebra/vector.h \
  ../../include/linear_algebra/../declarations.h \
  ../../include/linear_algebra/../numeric/integer.h \
  ../../include/linear_algebra/../numeric/../declarations.h \
  ../../include/linear_algebra/../numeric/../utility/stlio.h \
  ../../include/linear_algebra/../numeric/../utility/../base/array.h \
  ../../include/linear_algebra/../numeric/../utility/../base/../declarations.h \
  ../../include/linear_algebra/linear_program.h \
  ../../include/linear_algebra/../linear_algebra/matrix.h \
  ../../include/linear_algebra/../linear_algebra/../declarations.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/integer.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../numeric/numerical_types.h \
  ../../include/python/../numeric/../declarations.h \
  ../../include/python/../numeric/../numeric/integer.h \
  ../../include/python/../numeric/../numeric/dyadic.h \
  ../../include/python/../numeric/../numeric/numerical_traits.h \
  ../../include/python/../numeric/../utility/stlio.h \
  ../../include/python/../declarations.h

