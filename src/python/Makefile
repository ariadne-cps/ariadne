include ../../config.mk

# Must seach ../ for libariadne.so first since otherwise it may find
# old version. To see this, try putting 
# PYLDXXFLAGS = -L$HOME/lib -L../
# rm -f $HOME/lib/libariadne.so; touch $HOME/lib/libariadne.so
# and then remake. The correct ariande is not used for linking
CXXFLAGS += -I/usr/include/python2.4/
PYLDXXFLAGS = -L../ ${LDXXFLAGS}
LDLIBS = -lariadne ${LIBBOOSTPYTHON} ${LIBMPFR} ${LIBPPL} ${LIBGMPXX} 

BASE_OBJ=base_module.o \
	export_array.o export_binary_tree.o 

NUMERIC_OBJ=numeric_module.o \
	export_numeric.o export_interval.o export_function.o \

LINEAR_ALGEBRA_OBJ = linear_algebra_module.o \
	export_vector.o export_matrix.o export_linear_program.o

GEOMETRY_OBJ = geometry_module.o \
	export_point.o export_rectangle.o export_parallelotope.o \
	export_zonotope.o export_polytope.o export_simplex.o export_polyhedron.o \
	export_lattice_set.o export_list_set.o export_grid_set.o \
	export_partition_tree_set.o 

SYSTEM_OBJ = system_module.o \
	export_lattice_map.o \
	export_map.o export_affine_map.o export_affine_multimap.o \
	export_polynomial_map.o \
	export_vector_field.o export_affine_vector_field.o \
	export_henon_map.o export_lorenz_system.o 

EVALUATION_OBJ = evaluation_module.o \
	export_apply.o \
	export_integrate.o

OUTPUT_OBJ = output_module.o \
	export_postscript_output.o

all: python 

python: base.so numeric.so linear_algebra.so \
	geometry.so system.so evaluation.so \
	output.so

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
	${CXX} -shared -o base.so ${PYLDXXFLAGS} ${BASE_OBJ} ${LDLIBS} ${RPATHFLAGS}

numeric.so: ${NUMERIC_OBJ}
	${CXX} -shared -o numeric.so ${PYLDXXFLAGS} ${NUMERIC_OBJ} ${LDLIBS} ${RPATHFLAGS}

linear_algebra.so: ${LINEAR_ALGEBRA_OBJ}
	${CXX} -shared -o linear_algebra.so ${PYLDXXFLAGS} ${LINEAR_ALGEBRA_OBJ} ${LDLIBS} ${RPATHFLAGS}

geometry.so: ${GEOMETRY_OBJ}
	${CXX} -shared -o geometry.so ${PYLDXXFLAGS} ${GEOMETRY_OBJ} ${LDLIBS} ${RPATHFLAGS}

system.so: ${SYSTEM_OBJ}
	${CXX} -shared -o system.so ${PYLDXXFLAGS} ${SYSTEM_OBJ} ${LDLIBS} ${RPATHFLAGS}

evaluation.so: ${EVALUATION_OBJ}
	${CXX} -shared -o evaluation.so ${PYLDXXFLAGS} ${EVALUATION_OBJ} ${LDLIBS} ${RPATHFLAGS}

output.so: ${OUTPUT_OBJ}
	${CXX} -shared -o output.so ${PYLDXXFLAGS} ${OUTPUT_OBJ} ${LDLIBS} ${RPATHFLAGS}


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

depclean:
	@sed -e '/^### Do Not edit this line$$/,$$d' Makefile \
		> MMM.$$$$ && mv MMM.$$$$ Makefile
	@echo "### Do Not edit this line" >> Makefile

### Do Not edit this line
### Everything below is auto-generated
base_module.o: base_module.cc

export_array.o: export_array.cc ../../include/base/array.h \
  ../../include/base/../declarations.h \
  ../../include/base/array_operations.h ../../include/base/stlio.h \
  ../../include/base/../base/array.h ../../include/python/typedefs.h \
  ../../include/python/../real_typedef.h \
  ../../include/python/../numeric/numerical_types.h \
  ../../include/python/../numeric/mpfloat.h /home/collins/include/mpfr.h \
  ../../include/python/../numeric/../declarations.h \
  ../../include/python/../numeric/../numeric/numerical_traits.h \
  ../../include/python/../numeric/../numeric/approximation.h \
  ../../include/python/../numeric/../numeric/function.h \
  ../../include/python/../numeric/../numeric/rational.h \
  ../../include/python/../numeric/../numeric/../declarations.h \
  ../../include/python/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/python/../numeric/../numeric/../numeric/arithmetic.h \
  ../../include/python/../numeric/../numeric/../numeric/function.h \
  ../../include/python/../numeric/../numeric/../numeric/integer.h \
  ../../include/python/../numeric/../numeric/../numeric/../declarations.h \
  ../../include/python/../numeric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/python/../numeric/../numeric/../numeric/../numeric/arithmetic.h \
  ../../include/python/../numeric/float64.h \
  ../../include/python/../numeric/../numeric/integer.h \
  ../../include/python/../numeric/mpfloat.h \
  ../../include/python/../numeric/dyadic.h \
  ../../include/python/../numeric/../base/dyadic.h \
  ../../include/python/../numeric/rational.h \
  ../../include/python/../declarations.h \
  ../../include/python/python_utilities.h

export_binary_tree.o: export_binary_tree.cc \
  ../../include/combinatoric/binary_word.h \
  ../../include/combinatoric/../declarations.h \
  ../../include/combinatoric/binary_tree.h \
  ../../include/combinatoric/../base/array.h \
  ../../include/combinatoric/../base/../declarations.h \
  ../../include/combinatoric/../combinatoric/binary_word.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../numeric/numerical_types.h \
  ../../include/python/../numeric/mpfloat.h /home/collins/include/mpfr.h \
  ../../include/python/../numeric/../declarations.h \
  ../../include/python/../numeric/../numeric/numerical_traits.h \
  ../../include/python/../numeric/../numeric/approximation.h \
  ../../include/python/../numeric/../numeric/function.h \
  ../../include/python/../numeric/../numeric/rational.h \
  ../../include/python/../numeric/../numeric/../declarations.h \
  ../../include/python/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/python/../numeric/../numeric/../numeric/arithmetic.h \
  ../../include/python/../numeric/../numeric/../numeric/function.h \
  ../../include/python/../numeric/../numeric/../numeric/integer.h \
  ../../include/python/../numeric/../numeric/../numeric/../declarations.h \
  ../../include/python/../numeric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/python/../numeric/../numeric/../numeric/../numeric/arithmetic.h \
  ../../include/python/../numeric/float64.h \
  ../../include/python/../numeric/../numeric/integer.h \
  ../../include/python/../numeric/mpfloat.h \
  ../../include/python/../numeric/dyadic.h \
  ../../include/python/../numeric/../base/dyadic.h \
  ../../include/python/../numeric/rational.h \
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
  ../../include/geometry/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/approximation.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../../include/geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/rounding.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/function.h \
  ../../include/geometry/../linear_algebra/matrix.h \
  ../../include/geometry/point_list.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/rectangle.h \
  ../../include/geometry/../numeric/arithmetic.h \
  ../../include/geometry/../numeric/function.h \
  ../../include/geometry/../numeric/interval.h \
  ../../include/geometry/../geometry/ppl_polyhedron.h \
  ../../include/geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../numeric/rational.h \
  ../../include/geometry/../geometry/../numeric/../declarations.h \
  ../../include/geometry/../geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../geometry/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../numeric/../numeric/function.h \
  ../../include/geometry/../geometry/../numeric/../numeric/integer.h \
  ../../include/geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/rectangle_expression.h \
  ../../include/linear_algebra/vector.h ../../include/python/typedefs.h \
  ../../include/python/../real_typedef.h \
  ../../include/python/../numeric/numerical_types.h \
  ../../include/python/../numeric/mpfloat.h /home/collins/include/mpfr.h \
  ../../include/python/../numeric/../declarations.h \
  ../../include/python/../numeric/../numeric/numerical_traits.h \
  ../../include/python/../numeric/../numeric/function.h \
  ../../include/python/../numeric/../numeric/rational.h \
  ../../include/python/../numeric/float64.h \
  ../../include/python/../numeric/../numeric/integer.h \
  ../../include/python/../numeric/mpfloat.h \
  ../../include/python/../numeric/dyadic.h \
  ../../include/python/../numeric/../base/dyadic.h \
  ../../include/python/../numeric/rational.h \
  ../../include/python/../declarations.h

export_rectangle.o: export_rectangle.cc \
  ../../include/geometry/rectangle.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../numeric/arithmetic.h \
  ../../include/geometry/../numeric/function.h \
  ../../include/geometry/../numeric/interval.h \
  ../../include/geometry/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../numeric/../numeric/approximation.h \
  ../../include/geometry/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../numeric/../numeric/rounding.h \
  ../../include/geometry/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../numeric/../numeric/../numeric/function.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../linear_algebra/../numeric/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../../include/geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../geometry/ppl_polyhedron.h \
  ../../include/geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../numeric/rational.h \
  ../../include/geometry/../geometry/../numeric/../declarations.h \
  ../../include/geometry/../geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../geometry/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../numeric/../numeric/function.h \
  ../../include/geometry/../geometry/../numeric/../numeric/integer.h \
  ../../include/geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../base/../declarations.h \
  ../../include/geometry/../geometry/rectangle_expression.h \
  ../../include/geometry/zonotope.h \
  ../../include/geometry/../linear_algebra/matrix.h \
  ../../include/geometry/list_set.h \
  ../../include/geometry/../base/utility.h \
  ../../include/geometry/../base/../declarations.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../numeric/numerical_types.h \
  ../../include/python/../numeric/mpfloat.h /home/collins/include/mpfr.h \
  ../../include/python/../numeric/../declarations.h \
  ../../include/python/../numeric/../numeric/numerical_traits.h \
  ../../include/python/../numeric/../numeric/function.h \
  ../../include/python/../numeric/../numeric/rational.h \
  ../../include/python/../numeric/float64.h \
  ../../include/python/../numeric/../numeric/integer.h \
  ../../include/python/../numeric/mpfloat.h \
  ../../include/python/../numeric/dyadic.h \
  ../../include/python/../numeric/../base/dyadic.h \
  ../../include/python/../numeric/rational.h \
  ../../include/python/../declarations.h

export_parallelotope.o: export_parallelotope.cc \
  ../../include/geometry/parallelotope.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../numeric/interval.h \
  ../../include/geometry/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../numeric/../numeric/approximation.h \
  ../../include/geometry/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../numeric/../numeric/rounding.h \
  ../../include/geometry/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../numeric/../numeric/../numeric/function.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../linear_algebra/../numeric/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../../include/geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../base/../declarations.h \
  ../../include/geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/zonotope.h \
  ../../include/geometry/../geometry/../numeric/interval.h \
  ../../include/geometry/../geometry/../geometry/ppl_polyhedron.h \
  ../../include/geometry/../geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../geometry/../numeric/rational.h \
  ../../include/geometry/../geometry/../geometry/../numeric/../declarations.h \
  ../../include/geometry/../geometry/../geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../geometry/../geometry/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../geometry/../numeric/../numeric/function.h \
  ../../include/geometry/../geometry/../geometry/../numeric/../numeric/integer.h \
  ../../include/geometry/../geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/rectangle.h \
  ../../include/geometry/../numeric/arithmetic.h \
  ../../include/geometry/../numeric/function.h \
  ../../include/geometry/../geometry/ppl_polyhedron.h \
  ../../include/geometry/../geometry/rectangle_expression.h \
  ../../include/geometry/zonotope.h ../../include/geometry/polyhedron.h \
  ../../include/geometry/../geometry/point_list.h \
  ../../include/geometry/../geometry/rectangle.h \
  ../../include/geometry/list_set.h \
  ../../include/geometry/../base/utility.h \
  ../../include/geometry/../base/../declarations.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../numeric/numerical_types.h \
  ../../include/python/../numeric/mpfloat.h /home/collins/include/mpfr.h \
  ../../include/python/../numeric/../declarations.h \
  ../../include/python/../numeric/../numeric/numerical_traits.h \
  ../../include/python/../numeric/../numeric/function.h \
  ../../include/python/../numeric/../numeric/rational.h \
  ../../include/python/../numeric/float64.h \
  ../../include/python/../numeric/../numeric/integer.h \
  ../../include/python/../numeric/mpfloat.h \
  ../../include/python/../numeric/dyadic.h \
  ../../include/python/../numeric/../base/dyadic.h \
  ../../include/python/../numeric/rational.h \
  ../../include/python/../declarations.h \
  ../../include/python/python_utilities.h

export_zonotope.o: export_zonotope.cc ../../include/geometry/zonotope.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../linear_algebra/../numeric/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/approximation.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../../include/geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/rounding.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/function.h \
  ../../include/geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../numeric/interval.h \
  ../../include/geometry/../geometry/ppl_polyhedron.h \
  ../../include/geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../numeric/rational.h \
  ../../include/geometry/../geometry/../numeric/../declarations.h \
  ../../include/geometry/../geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../geometry/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../numeric/../numeric/function.h \
  ../../include/geometry/../geometry/../numeric/../numeric/integer.h \
  ../../include/geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../base/../declarations.h \
  ../../include/geometry/rectangle.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../numeric/arithmetic.h \
  ../../include/geometry/../numeric/function.h \
  ../../include/geometry/../geometry/rectangle_expression.h \
  ../../include/geometry/parallelotope.h \
  ../../include/geometry/../geometry/zonotope.h \
  ../../include/geometry/polyhedron.h \
  ../../include/geometry/../geometry/point_list.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/rectangle.h \
  ../../include/geometry/list_set.h \
  ../../include/geometry/../base/utility.h \
  ../../include/geometry/../base/../declarations.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../numeric/numerical_types.h \
  ../../include/python/../numeric/mpfloat.h /home/collins/include/mpfr.h \
  ../../include/python/../numeric/../declarations.h \
  ../../include/python/../numeric/../numeric/numerical_traits.h \
  ../../include/python/../numeric/../numeric/function.h \
  ../../include/python/../numeric/../numeric/rational.h \
  ../../include/python/../numeric/float64.h \
  ../../include/python/../numeric/../numeric/integer.h \
  ../../include/python/../numeric/mpfloat.h \
  ../../include/python/../numeric/dyadic.h \
  ../../include/python/../numeric/../base/dyadic.h \
  ../../include/python/../numeric/rational.h \
  ../../include/python/../declarations.h \
  ../../include/python/python_utilities.h

export_polytope.o: export_polytope.cc ../../include/geometry/rectangle.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../numeric/arithmetic.h \
  ../../include/geometry/../numeric/function.h \
  ../../include/geometry/../numeric/interval.h \
  ../../include/geometry/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../numeric/../numeric/approximation.h \
  ../../include/geometry/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../numeric/../numeric/rounding.h \
  ../../include/geometry/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../numeric/../numeric/../numeric/function.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../linear_algebra/../numeric/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../../include/geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../geometry/ppl_polyhedron.h \
  ../../include/geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../numeric/rational.h \
  ../../include/geometry/../geometry/../numeric/../declarations.h \
  ../../include/geometry/../geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../geometry/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../numeric/../numeric/function.h \
  ../../include/geometry/../geometry/../numeric/../numeric/integer.h \
  ../../include/geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../base/../declarations.h \
  ../../include/geometry/../geometry/rectangle_expression.h \
  ../../include/geometry/polytope.h \
  ../../include/geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/point_list.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/rectangle.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../numeric/numerical_types.h \
  ../../include/python/../numeric/mpfloat.h /home/collins/include/mpfr.h \
  ../../include/python/../numeric/../declarations.h \
  ../../include/python/../numeric/../numeric/numerical_traits.h \
  ../../include/python/../numeric/../numeric/function.h \
  ../../include/python/../numeric/../numeric/rational.h \
  ../../include/python/../numeric/float64.h \
  ../../include/python/../numeric/../numeric/integer.h \
  ../../include/python/../numeric/mpfloat.h \
  ../../include/python/../numeric/dyadic.h \
  ../../include/python/../numeric/../base/dyadic.h \
  ../../include/python/../numeric/rational.h \
  ../../include/python/../declarations.h

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
  ../../include/geometry/../geometry/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/../numeric/approximation.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/../numeric/rounding.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/../numeric/../numeric/function.h \
  ../../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/point_list.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/polytope.h \
  ../../include/geometry/../geometry/../geometry/ppl_polyhedron.h \
  ../../include/geometry/../geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../geometry/../numeric/rational.h \
  ../../include/geometry/../geometry/../geometry/../numeric/../declarations.h \
  ../../include/geometry/../geometry/../geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../geometry/../geometry/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../geometry/../numeric/../numeric/function.h \
  ../../include/geometry/../geometry/../geometry/../numeric/../numeric/integer.h \
  ../../include/geometry/../geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../geometry/point_list.h \
  ../../include/geometry/../geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/../geometry/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../geometry/../numeric/function.h \
  ../../include/geometry/../geometry/../geometry/../numeric/interval.h \
  ../../include/geometry/../geometry/../geometry/../geometry/ppl_polyhedron.h \
  ../../include/geometry/../geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../geometry/../geometry/rectangle_expression.h \
  ../../include/geometry/../geometry/../geometry/../geometry/../declarations.h \
  ../../include/geometry/rectangle.h ../../include/geometry/polyhedron.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/ppl_polyhedron.h \
  ../../include/geometry/../geometry/rectangle.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../numeric/numerical_types.h \
  ../../include/python/../numeric/mpfloat.h /home/collins/include/mpfr.h \
  ../../include/python/../numeric/../declarations.h \
  ../../include/python/../numeric/../numeric/numerical_traits.h \
  ../../include/python/../numeric/../numeric/function.h \
  ../../include/python/../numeric/../numeric/rational.h \
  ../../include/python/../numeric/float64.h \
  ../../include/python/../numeric/../numeric/integer.h \
  ../../include/python/../numeric/mpfloat.h \
  ../../include/python/../numeric/dyadic.h \
  ../../include/python/../numeric/../base/dyadic.h \
  ../../include/python/../numeric/rational.h \
  ../../include/python/../declarations.h

export_polyhedron.o: export_polyhedron.cc \
  ../../include/geometry/rectangle.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../numeric/arithmetic.h \
  ../../include/geometry/../numeric/function.h \
  ../../include/geometry/../numeric/interval.h \
  ../../include/geometry/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../numeric/../numeric/approximation.h \
  ../../include/geometry/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../numeric/../numeric/rounding.h \
  ../../include/geometry/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../numeric/../numeric/../numeric/function.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../linear_algebra/../numeric/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../../include/geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../geometry/ppl_polyhedron.h \
  ../../include/geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../numeric/rational.h \
  ../../include/geometry/../geometry/../numeric/../declarations.h \
  ../../include/geometry/../geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../geometry/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../numeric/../numeric/function.h \
  ../../include/geometry/../geometry/../numeric/../numeric/integer.h \
  ../../include/geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../base/../declarations.h \
  ../../include/geometry/../geometry/rectangle_expression.h \
  ../../include/geometry/polyhedron.h \
  ../../include/geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/point_list.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/rectangle.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../numeric/numerical_types.h \
  ../../include/python/../numeric/mpfloat.h /home/collins/include/mpfr.h \
  ../../include/python/../numeric/../declarations.h \
  ../../include/python/../numeric/../numeric/numerical_traits.h \
  ../../include/python/../numeric/../numeric/function.h \
  ../../include/python/../numeric/../numeric/rational.h \
  ../../include/python/../numeric/float64.h \
  ../../include/python/../numeric/../numeric/integer.h \
  ../../include/python/../numeric/mpfloat.h \
  ../../include/python/../numeric/dyadic.h \
  ../../include/python/../numeric/../base/dyadic.h \
  ../../include/python/../numeric/rational.h \
  ../../include/python/../declarations.h

export_lattice_set.o: export_lattice_set.cc \
  ../../include/combinatoric/lattice_set.h \
  ../../include/combinatoric/../base/array.h \
  ../../include/combinatoric/../base/../declarations.h \
  ../../include/combinatoric/../base/iterator.h \
  ../../include/combinatoric/../base/../base/array.h \
  ../../include/combinatoric/../base/array_operations.h \
  ../../include/combinatoric/../numeric/interval.h \
  ../../include/combinatoric/../numeric/../declarations.h \
  ../../include/combinatoric/../numeric/../numeric/numerical_traits.h \
  ../../include/combinatoric/../numeric/../numeric/approximation.h \
  ../../include/combinatoric/../numeric/../numeric/arithmetic.h \
  ../../include/combinatoric/../numeric/../numeric/rounding.h \
  ../../include/combinatoric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/combinatoric/../numeric/../numeric/../numeric/function.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../numeric/numerical_types.h \
  ../../include/python/../numeric/mpfloat.h /home/collins/include/mpfr.h \
  ../../include/python/../numeric/../declarations.h \
  ../../include/python/../numeric/../numeric/numerical_traits.h \
  ../../include/python/../numeric/../numeric/function.h \
  ../../include/python/../numeric/../numeric/rational.h \
  ../../include/python/../numeric/../numeric/../declarations.h \
  ../../include/python/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/python/../numeric/../numeric/../numeric/arithmetic.h \
  ../../include/python/../numeric/../numeric/../numeric/function.h \
  ../../include/python/../numeric/../numeric/../numeric/integer.h \
  ../../include/python/../numeric/../numeric/../numeric/../declarations.h \
  ../../include/python/../numeric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/python/../numeric/../numeric/../numeric/../numeric/arithmetic.h \
  ../../include/python/../numeric/float64.h \
  ../../include/python/../numeric/../numeric/integer.h \
  ../../include/python/../numeric/mpfloat.h \
  ../../include/python/../numeric/dyadic.h \
  ../../include/python/../numeric/../base/dyadic.h \
  ../../include/python/../numeric/rational.h \
  ../../include/python/../declarations.h

export_list_set.o: export_list_set.cc ../../include/geometry/rectangle.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../numeric/arithmetic.h \
  ../../include/geometry/../numeric/function.h \
  ../../include/geometry/../numeric/interval.h \
  ../../include/geometry/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../numeric/../numeric/approximation.h \
  ../../include/geometry/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../numeric/../numeric/rounding.h \
  ../../include/geometry/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../numeric/../numeric/../numeric/function.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../linear_algebra/../numeric/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../../include/geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../geometry/ppl_polyhedron.h \
  ../../include/geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../numeric/rational.h \
  ../../include/geometry/../geometry/../numeric/../declarations.h \
  ../../include/geometry/../geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../geometry/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../numeric/../numeric/function.h \
  ../../include/geometry/../geometry/../numeric/../numeric/integer.h \
  ../../include/geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../base/../declarations.h \
  ../../include/geometry/../geometry/rectangle_expression.h \
  ../../include/geometry/parallelotope.h \
  ../../include/geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/zonotope.h \
  ../../include/geometry/../geometry/../numeric/interval.h \
  ../../include/geometry/../geometry/../geometry/ppl_polyhedron.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/zonotope.h ../../include/geometry/polyhedron.h \
  ../../include/geometry/../geometry/point_list.h \
  ../../include/geometry/../geometry/rectangle.h \
  ../../include/geometry/list_set.h \
  ../../include/geometry/../base/utility.h \
  ../../include/geometry/../base/../declarations.h \
  ../../include/geometry/grid_set.h \
  ../../include/geometry/../base/array.h \
  ../../include/geometry/../combinatoric/lattice_set.h \
  ../../include/geometry/../combinatoric/../base/array.h \
  ../../include/geometry/../combinatoric/../base/iterator.h \
  ../../include/geometry/../combinatoric/../base/../base/array.h \
  ../../include/geometry/../combinatoric/../base/array_operations.h \
  ../../include/geometry/../combinatoric/../base/../declarations.h \
  ../../include/geometry/../combinatoric/../numeric/interval.h \
  ../../include/geometry/../geometry/grid.h \
  ../../include/geometry/../geometry/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../combinatoric/lattice_set.h \
  ../../include/geometry/partition_tree_set.h \
  ../../include/geometry/../base/iterator.h \
  ../../include/geometry/../combinatoric/binary_word.h \
  ../../include/geometry/../combinatoric/../declarations.h \
  ../../include/geometry/../combinatoric/binary_tree.h \
  ../../include/geometry/../combinatoric/../combinatoric/binary_word.h \
  ../../include/geometry/../combinatoric/subdivision_tree_set.h \
  ../../include/geometry/../combinatoric/../base/sequence.h \
  ../../include/geometry/../combinatoric/../combinatoric/binary_tree.h \
  ../../include/geometry/../combinatoric/../geometry/rectangle.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../numeric/numerical_types.h \
  ../../include/python/../numeric/mpfloat.h /home/collins/include/mpfr.h \
  ../../include/python/../numeric/../declarations.h \
  ../../include/python/../numeric/../numeric/numerical_traits.h \
  ../../include/python/../numeric/../numeric/function.h \
  ../../include/python/../numeric/../numeric/rational.h \
  ../../include/python/../numeric/float64.h \
  ../../include/python/../numeric/../numeric/integer.h \
  ../../include/python/../numeric/mpfloat.h \
  ../../include/python/../numeric/dyadic.h \
  ../../include/python/../numeric/../base/dyadic.h \
  ../../include/python/../numeric/rational.h \
  ../../include/python/../declarations.h \
  ../../include/python/python_utilities.h

export_grid_set.o: export_grid_set.cc \
  ../../include/linear_algebra/linear_algebra.h \
  ../../include/linear_algebra/../linear_algebra/vector.h \
  ../../include/linear_algebra/../linear_algebra/../declarations.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/integer.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/../declarations.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/../numeric/approximation.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/../numeric/arithmetic.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/arithmetic.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/interval.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/../numeric/rounding.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/../numeric/../numeric/function.h \
  ../../include/linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/polytope.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/ppl_polyhedron.h \
  ../../include/geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../numeric/rational.h \
  ../../include/geometry/../geometry/../numeric/../declarations.h \
  ../../include/geometry/../geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../geometry/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../numeric/../numeric/function.h \
  ../../include/geometry/../geometry/../numeric/../numeric/integer.h \
  ../../include/geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../base/../declarations.h \
  ../../include/geometry/../geometry/point_list.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../numeric/function.h \
  ../../include/geometry/../geometry/../numeric/interval.h \
  ../../include/geometry/../geometry/../geometry/ppl_polyhedron.h \
  ../../include/geometry/../geometry/../geometry/rectangle_expression.h \
  ../../include/geometry/../geometry/../geometry/../declarations.h \
  ../../include/geometry/parallelotope.h \
  ../../include/geometry/../numeric/interval.h \
  ../../include/geometry/../geometry/zonotope.h \
  ../../include/geometry/zonotope.h ../../include/geometry/rectangle.h \
  ../../include/geometry/list_set.h \
  ../../include/geometry/../base/utility.h \
  ../../include/geometry/../base/../declarations.h \
  ../../include/geometry/grid_set.h \
  ../../include/geometry/../base/array.h \
  ../../include/geometry/../combinatoric/lattice_set.h \
  ../../include/geometry/../combinatoric/../base/array.h \
  ../../include/geometry/../combinatoric/../base/iterator.h \
  ../../include/geometry/../combinatoric/../base/../base/array.h \
  ../../include/geometry/../combinatoric/../base/array_operations.h \
  ../../include/geometry/../combinatoric/../base/../declarations.h \
  ../../include/geometry/../combinatoric/../numeric/interval.h \
  ../../include/geometry/../geometry/grid.h \
  ../../include/geometry/../geometry/../combinatoric/lattice_set.h \
  ../../include/geometry/partition_tree_set.h \
  ../../include/geometry/../base/iterator.h \
  ../../include/geometry/../combinatoric/binary_word.h \
  ../../include/geometry/../combinatoric/../declarations.h \
  ../../include/geometry/../combinatoric/binary_tree.h \
  ../../include/geometry/../combinatoric/../combinatoric/binary_word.h \
  ../../include/geometry/../combinatoric/subdivision_tree_set.h \
  ../../include/geometry/../combinatoric/../base/sequence.h \
  ../../include/geometry/../combinatoric/../combinatoric/binary_tree.h \
  ../../include/geometry/../combinatoric/../geometry/rectangle.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../numeric/numerical_types.h \
  ../../include/python/../numeric/mpfloat.h /home/collins/include/mpfr.h \
  ../../include/python/../numeric/../declarations.h \
  ../../include/python/../numeric/../numeric/numerical_traits.h \
  ../../include/python/../numeric/../numeric/function.h \
  ../../include/python/../numeric/../numeric/rational.h \
  ../../include/python/../numeric/float64.h \
  ../../include/python/../numeric/../numeric/integer.h \
  ../../include/python/../numeric/mpfloat.h \
  ../../include/python/../numeric/dyadic.h \
  ../../include/python/../numeric/../base/dyadic.h \
  ../../include/python/../numeric/rational.h \
  ../../include/python/../declarations.h \
  ../../include/python/python_utilities.h

export_partition_tree_set.o: export_partition_tree_set.cc \
  ../../include/combinatoric/binary_word.h \
  ../../include/combinatoric/../declarations.h \
  ../../include/combinatoric/binary_tree.h \
  ../../include/combinatoric/../base/array.h \
  ../../include/combinatoric/../base/../declarations.h \
  ../../include/combinatoric/../combinatoric/binary_word.h \
  ../../include/geometry/rectangle.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../numeric/arithmetic.h \
  ../../include/geometry/../numeric/function.h \
  ../../include/geometry/../numeric/interval.h \
  ../../include/geometry/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../numeric/../numeric/approximation.h \
  ../../include/geometry/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../numeric/../numeric/rounding.h \
  ../../include/geometry/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../numeric/../numeric/../numeric/function.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../linear_algebra/../numeric/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../../include/geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../geometry/ppl_polyhedron.h \
  ../../include/geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../numeric/rational.h \
  ../../include/geometry/../geometry/../numeric/../declarations.h \
  ../../include/geometry/../geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../geometry/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../numeric/../numeric/function.h \
  ../../include/geometry/../geometry/../numeric/../numeric/integer.h \
  ../../include/geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/rectangle_expression.h \
  ../../include/geometry/grid_set.h \
  ../../include/geometry/../base/array.h \
  ../../include/geometry/../combinatoric/lattice_set.h \
  ../../include/geometry/../combinatoric/../base/array.h \
  ../../include/geometry/../combinatoric/../base/iterator.h \
  ../../include/geometry/../combinatoric/../base/../base/array.h \
  ../../include/geometry/../combinatoric/../base/array_operations.h \
  ../../include/geometry/../combinatoric/../base/../declarations.h \
  ../../include/geometry/../combinatoric/../numeric/interval.h \
  ../../include/geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/grid.h \
  ../../include/geometry/../geometry/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../combinatoric/lattice_set.h \
  ../../include/geometry/partition_tree_set.h \
  ../../include/geometry/../base/iterator.h \
  ../../include/geometry/../combinatoric/binary_word.h \
  ../../include/geometry/../combinatoric/binary_tree.h \
  ../../include/geometry/../combinatoric/subdivision_tree_set.h \
  ../../include/geometry/../combinatoric/../declarations.h \
  ../../include/geometry/../combinatoric/../base/sequence.h \
  ../../include/geometry/../combinatoric/../combinatoric/binary_word.h \
  ../../include/geometry/../combinatoric/../combinatoric/binary_tree.h \
  ../../include/geometry/../combinatoric/../geometry/rectangle.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../numeric/numerical_types.h \
  ../../include/python/../numeric/mpfloat.h /home/collins/include/mpfr.h \
  ../../include/python/../numeric/../declarations.h \
  ../../include/python/../numeric/../numeric/numerical_traits.h \
  ../../include/python/../numeric/../numeric/function.h \
  ../../include/python/../numeric/../numeric/rational.h \
  ../../include/python/../numeric/float64.h \
  ../../include/python/../numeric/../numeric/integer.h \
  ../../include/python/../numeric/mpfloat.h \
  ../../include/python/../numeric/dyadic.h \
  ../../include/python/../numeric/../base/dyadic.h \
  ../../include/python/../numeric/rational.h \
  ../../include/python/../declarations.h

evaluation_module.o: evaluation_module.cc

export_apply.o: export_apply.cc ../../include/geometry/rectangle.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../numeric/arithmetic.h \
  ../../include/geometry/../numeric/function.h \
  ../../include/geometry/../numeric/interval.h \
  ../../include/geometry/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../numeric/../numeric/approximation.h \
  ../../include/geometry/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../numeric/../numeric/rounding.h \
  ../../include/geometry/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../numeric/../numeric/../numeric/function.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../linear_algebra/../numeric/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../../include/geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../geometry/ppl_polyhedron.h \
  ../../include/geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../numeric/rational.h \
  ../../include/geometry/../geometry/../numeric/../declarations.h \
  ../../include/geometry/../geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../geometry/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../numeric/../numeric/function.h \
  ../../include/geometry/../geometry/../numeric/../numeric/integer.h \
  ../../include/geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../base/../declarations.h \
  ../../include/geometry/../geometry/rectangle_expression.h \
  ../../include/geometry/parallelotope.h \
  ../../include/geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/zonotope.h \
  ../../include/geometry/../geometry/../numeric/interval.h \
  ../../include/geometry/../geometry/../geometry/ppl_polyhedron.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/list_set.h \
  ../../include/geometry/../base/utility.h \
  ../../include/geometry/../base/../declarations.h \
  ../../include/geometry/grid_set.h \
  ../../include/geometry/../base/array.h \
  ../../include/geometry/../combinatoric/lattice_set.h \
  ../../include/geometry/../combinatoric/../base/array.h \
  ../../include/geometry/../combinatoric/../base/iterator.h \
  ../../include/geometry/../combinatoric/../base/../base/array.h \
  ../../include/geometry/../combinatoric/../base/array_operations.h \
  ../../include/geometry/../combinatoric/../base/../declarations.h \
  ../../include/geometry/../combinatoric/../numeric/interval.h \
  ../../include/geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/grid.h \
  ../../include/geometry/../geometry/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../combinatoric/lattice_set.h \
  ../../include/system/map.h ../../include/system/../declarations.h \
  ../../include/evaluation/apply.h \
  ../../include/evaluation/../declarations.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../numeric/numerical_types.h \
  ../../include/python/../numeric/mpfloat.h /home/collins/include/mpfr.h \
  ../../include/python/../numeric/../declarations.h \
  ../../include/python/../numeric/../numeric/numerical_traits.h \
  ../../include/python/../numeric/../numeric/function.h \
  ../../include/python/../numeric/../numeric/rational.h \
  ../../include/python/../numeric/float64.h \
  ../../include/python/../numeric/../numeric/integer.h \
  ../../include/python/../numeric/mpfloat.h \
  ../../include/python/../numeric/dyadic.h \
  ../../include/python/../numeric/../base/dyadic.h \
  ../../include/python/../numeric/rational.h \
  ../../include/python/../declarations.h

export_integrate.o: export_integrate.cc \
  ../../include/geometry/rectangle.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../numeric/arithmetic.h \
  ../../include/geometry/../numeric/function.h \
  ../../include/geometry/../numeric/interval.h \
  ../../include/geometry/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../numeric/../numeric/approximation.h \
  ../../include/geometry/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../numeric/../numeric/rounding.h \
  ../../include/geometry/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../numeric/../numeric/../numeric/function.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../linear_algebra/../numeric/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../../include/geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../geometry/ppl_polyhedron.h \
  ../../include/geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../numeric/rational.h \
  ../../include/geometry/../geometry/../numeric/../declarations.h \
  ../../include/geometry/../geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../geometry/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../numeric/../numeric/function.h \
  ../../include/geometry/../geometry/../numeric/../numeric/integer.h \
  ../../include/geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../base/../declarations.h \
  ../../include/geometry/../geometry/rectangle_expression.h \
  ../../include/geometry/parallelotope.h \
  ../../include/geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/zonotope.h \
  ../../include/geometry/../geometry/../numeric/interval.h \
  ../../include/geometry/../geometry/../geometry/ppl_polyhedron.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/zonotope.h ../../include/geometry/grid_set.h \
  ../../include/geometry/../base/array.h \
  ../../include/geometry/../combinatoric/lattice_set.h \
  ../../include/geometry/../combinatoric/../base/array.h \
  ../../include/geometry/../combinatoric/../base/iterator.h \
  ../../include/geometry/../combinatoric/../base/../base/array.h \
  ../../include/geometry/../combinatoric/../base/array_operations.h \
  ../../include/geometry/../combinatoric/../base/../declarations.h \
  ../../include/geometry/../combinatoric/../numeric/interval.h \
  ../../include/geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/grid.h \
  ../../include/geometry/../geometry/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../combinatoric/lattice_set.h \
  ../../include/geometry/list_set.h \
  ../../include/geometry/../base/utility.h \
  ../../include/geometry/../base/../declarations.h \
  ../../include/system/vector_field.h \
  ../../include/system/../declarations.h \
  ../../include/evaluation/integrator.h \
  ../../include/evaluation/../declarations.h \
  ../../include/evaluation/lohner_integrator.h \
  ../../include/evaluation/../evaluation/integrator.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../numeric/numerical_types.h \
  ../../include/python/../numeric/mpfloat.h /home/collins/include/mpfr.h \
  ../../include/python/../numeric/../declarations.h \
  ../../include/python/../numeric/../numeric/numerical_traits.h \
  ../../include/python/../numeric/../numeric/function.h \
  ../../include/python/../numeric/../numeric/rational.h \
  ../../include/python/../numeric/float64.h \
  ../../include/python/../numeric/../numeric/integer.h \
  ../../include/python/../numeric/mpfloat.h \
  ../../include/python/../numeric/dyadic.h \
  ../../include/python/../numeric/../base/dyadic.h \
  ../../include/python/../numeric/rational.h \
  ../../include/python/../declarations.h

system_module.o: system_module.cc

export_lattice_map.o: export_lattice_map.cc \
  ../../include/combinatoric/lattice_map.h \
  ../../include/combinatoric/../combinatoric/lattice_set.h \
  ../../include/combinatoric/../combinatoric/../base/array.h \
  ../../include/combinatoric/../combinatoric/../base/../declarations.h \
  ../../include/combinatoric/../combinatoric/../base/iterator.h \
  ../../include/combinatoric/../combinatoric/../base/../base/array.h \
  ../../include/combinatoric/../combinatoric/../base/array_operations.h \
  ../../include/combinatoric/../combinatoric/../numeric/interval.h \
  ../../include/combinatoric/../combinatoric/../numeric/../declarations.h \
  ../../include/combinatoric/../combinatoric/../numeric/../numeric/numerical_traits.h \
  ../../include/combinatoric/../combinatoric/../numeric/../numeric/approximation.h \
  ../../include/combinatoric/../combinatoric/../numeric/../numeric/arithmetic.h \
  ../../include/combinatoric/../combinatoric/../numeric/../numeric/rounding.h \
  ../../include/combinatoric/../combinatoric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/combinatoric/../combinatoric/../numeric/../numeric/../numeric/function.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../numeric/numerical_types.h \
  ../../include/python/../numeric/mpfloat.h /home/collins/include/mpfr.h \
  ../../include/python/../numeric/../declarations.h \
  ../../include/python/../numeric/../numeric/numerical_traits.h \
  ../../include/python/../numeric/../numeric/function.h \
  ../../include/python/../numeric/../numeric/rational.h \
  ../../include/python/../numeric/../numeric/../declarations.h \
  ../../include/python/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/python/../numeric/../numeric/../numeric/arithmetic.h \
  ../../include/python/../numeric/../numeric/../numeric/function.h \
  ../../include/python/../numeric/../numeric/../numeric/integer.h \
  ../../include/python/../numeric/../numeric/../numeric/../declarations.h \
  ../../include/python/../numeric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/python/../numeric/../numeric/../numeric/../numeric/arithmetic.h \
  ../../include/python/../numeric/float64.h \
  ../../include/python/../numeric/../numeric/integer.h \
  ../../include/python/../numeric/mpfloat.h \
  ../../include/python/../numeric/dyadic.h \
  ../../include/python/../numeric/../base/dyadic.h \
  ../../include/python/../numeric/rational.h \
  ../../include/python/../declarations.h

export_map.o: export_map.cc ../../include/system/map.h \
  ../../include/system/../declarations.h ../../include/python/typedefs.h \
  ../../include/python/../real_typedef.h \
  ../../include/python/../numeric/numerical_types.h \
  ../../include/python/../numeric/mpfloat.h /home/collins/include/mpfr.h \
  ../../include/python/../numeric/../declarations.h \
  ../../include/python/../numeric/../numeric/numerical_traits.h \
  ../../include/python/../numeric/../numeric/approximation.h \
  ../../include/python/../numeric/../numeric/function.h \
  ../../include/python/../numeric/../numeric/rational.h \
  ../../include/python/../numeric/../numeric/../declarations.h \
  ../../include/python/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/python/../numeric/../numeric/../numeric/arithmetic.h \
  ../../include/python/../numeric/../numeric/../numeric/function.h \
  ../../include/python/../numeric/../numeric/../numeric/integer.h \
  ../../include/python/../numeric/../numeric/../numeric/../declarations.h \
  ../../include/python/../numeric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/python/../numeric/../numeric/../numeric/../numeric/arithmetic.h \
  ../../include/python/../numeric/float64.h \
  ../../include/python/../numeric/../numeric/integer.h \
  ../../include/python/../numeric/mpfloat.h \
  ../../include/python/../numeric/dyadic.h \
  ../../include/python/../numeric/../base/dyadic.h \
  ../../include/python/../numeric/rational.h \
  ../../include/python/../declarations.h

export_affine_map.o: export_affine_map.cc ../../include/geometry/point.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../base/array.h \
  ../../include/geometry/../base/../declarations.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../linear_algebra/../numeric/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/approximation.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../../include/geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/rounding.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/function.h \
  ../../include/geometry/../linear_algebra/matrix.h \
  ../../include/geometry/rectangle.h \
  ../../include/geometry/../numeric/arithmetic.h \
  ../../include/geometry/../numeric/function.h \
  ../../include/geometry/../numeric/interval.h \
  ../../include/geometry/../geometry/ppl_polyhedron.h \
  ../../include/geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../numeric/rational.h \
  ../../include/geometry/../geometry/../numeric/../declarations.h \
  ../../include/geometry/../geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../geometry/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../numeric/../numeric/function.h \
  ../../include/geometry/../geometry/../numeric/../numeric/integer.h \
  ../../include/geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/../geometry/rectangle_expression.h \
  ../../include/geometry/parallelotope.h \
  ../../include/geometry/../geometry/zonotope.h \
  ../../include/geometry/../geometry/../numeric/interval.h \
  ../../include/geometry/../geometry/../geometry/ppl_polyhedron.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/simplex.h \
  ../../include/geometry/../geometry/point_list.h \
  ../../include/geometry/../geometry/polytope.h \
  ../../include/geometry/../geometry/../geometry/point_list.h \
  ../../include/geometry/../geometry/../geometry/rectangle.h \
  ../../include/geometry/zonotope.h ../../include/geometry/polytope.h \
  ../../include/system/affine_map.h \
  ../../include/system/../declarations.h \
  ../../include/system/../linear_algebra/vector.h \
  ../../include/system/../linear_algebra/matrix.h \
  ../../include/system/../system/map.h \
  ../../include/system/../system/../declarations.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../numeric/numerical_types.h \
  ../../include/python/../numeric/mpfloat.h /home/collins/include/mpfr.h \
  ../../include/python/../numeric/../declarations.h \
  ../../include/python/../numeric/../numeric/numerical_traits.h \
  ../../include/python/../numeric/../numeric/function.h \
  ../../include/python/../numeric/../numeric/rational.h \
  ../../include/python/../numeric/float64.h \
  ../../include/python/../numeric/../numeric/integer.h \
  ../../include/python/../numeric/mpfloat.h \
  ../../include/python/../numeric/dyadic.h \
  ../../include/python/../numeric/../base/dyadic.h \
  ../../include/python/../numeric/rational.h \
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
  ../../include/geometry/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/approximation.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../../include/geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/rounding.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/function.h \
  ../../include/geometry/../linear_algebra/matrix.h \
  ../../include/geometry/rectangle.h \
  ../../include/geometry/../numeric/arithmetic.h \
  ../../include/geometry/../numeric/function.h \
  ../../include/geometry/../numeric/interval.h \
  ../../include/geometry/../geometry/ppl_polyhedron.h \
  ../../include/geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../numeric/rational.h \
  ../../include/geometry/../geometry/../numeric/../declarations.h \
  ../../include/geometry/../geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../geometry/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../numeric/../numeric/function.h \
  ../../include/geometry/../geometry/../numeric/../numeric/integer.h \
  ../../include/geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/../geometry/rectangle_expression.h \
  ../../include/geometry/parallelotope.h \
  ../../include/geometry/../geometry/zonotope.h \
  ../../include/geometry/../geometry/../numeric/interval.h \
  ../../include/geometry/../geometry/../geometry/ppl_polyhedron.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/simplex.h \
  ../../include/geometry/../geometry/point_list.h \
  ../../include/geometry/../geometry/polytope.h \
  ../../include/geometry/../geometry/../geometry/point_list.h \
  ../../include/geometry/../geometry/../geometry/rectangle.h \
  ../../include/geometry/zonotope.h ../../include/geometry/polytope.h \
  ../../include/system/affine_multimap.h ../../include/system/multimap.h \
  ../../include/system/../declarations.h ../../include/python/typedefs.h \
  ../../include/python/../real_typedef.h \
  ../../include/python/../numeric/numerical_types.h \
  ../../include/python/../numeric/mpfloat.h /home/collins/include/mpfr.h \
  ../../include/python/../numeric/../declarations.h \
  ../../include/python/../numeric/../numeric/numerical_traits.h \
  ../../include/python/../numeric/../numeric/function.h \
  ../../include/python/../numeric/../numeric/rational.h \
  ../../include/python/../numeric/float64.h \
  ../../include/python/../numeric/../numeric/integer.h \
  ../../include/python/../numeric/mpfloat.h \
  ../../include/python/../numeric/dyadic.h \
  ../../include/python/../numeric/../base/dyadic.h \
  ../../include/python/../numeric/rational.h \
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
  ../../include/geometry/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/approximation.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../../include/geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/rounding.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/function.h \
  ../../include/geometry/../linear_algebra/matrix.h \
  ../../include/geometry/rectangle.h \
  ../../include/geometry/../numeric/arithmetic.h \
  ../../include/geometry/../numeric/function.h \
  ../../include/geometry/../numeric/interval.h \
  ../../include/geometry/../geometry/ppl_polyhedron.h \
  ../../include/geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../numeric/rational.h \
  ../../include/geometry/../geometry/../numeric/../declarations.h \
  ../../include/geometry/../geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../geometry/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../numeric/../numeric/function.h \
  ../../include/geometry/../geometry/../numeric/../numeric/integer.h \
  ../../include/geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/../geometry/rectangle_expression.h \
  ../../include/system/polynomial_map.h \
  ../../include/system/../base/array.h \
  ../../include/system/../linear_algebra/matrix.h \
  ../../include/system/../system/map.h \
  ../../include/system/../system/../declarations.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../numeric/numerical_types.h \
  ../../include/python/../numeric/mpfloat.h /home/collins/include/mpfr.h \
  ../../include/python/../numeric/../declarations.h \
  ../../include/python/../numeric/../numeric/numerical_traits.h \
  ../../include/python/../numeric/../numeric/function.h \
  ../../include/python/../numeric/../numeric/rational.h \
  ../../include/python/../numeric/float64.h \
  ../../include/python/../numeric/../numeric/integer.h \
  ../../include/python/../numeric/mpfloat.h \
  ../../include/python/../numeric/dyadic.h \
  ../../include/python/../numeric/../base/dyadic.h \
  ../../include/python/../numeric/rational.h \
  ../../include/python/../declarations.h

export_vector_field.o: export_vector_field.cc \
  ../../include/system/vector_field.h \
  ../../include/system/../declarations.h ../../include/python/typedefs.h \
  ../../include/python/../real_typedef.h \
  ../../include/python/../numeric/numerical_types.h \
  ../../include/python/../numeric/mpfloat.h /home/collins/include/mpfr.h \
  ../../include/python/../numeric/../declarations.h \
  ../../include/python/../numeric/../numeric/numerical_traits.h \
  ../../include/python/../numeric/../numeric/approximation.h \
  ../../include/python/../numeric/../numeric/function.h \
  ../../include/python/../numeric/../numeric/rational.h \
  ../../include/python/../numeric/../numeric/../declarations.h \
  ../../include/python/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/python/../numeric/../numeric/../numeric/arithmetic.h \
  ../../include/python/../numeric/../numeric/../numeric/function.h \
  ../../include/python/../numeric/../numeric/../numeric/integer.h \
  ../../include/python/../numeric/../numeric/../numeric/../declarations.h \
  ../../include/python/../numeric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/python/../numeric/../numeric/../numeric/../numeric/arithmetic.h \
  ../../include/python/../numeric/float64.h \
  ../../include/python/../numeric/../numeric/integer.h \
  ../../include/python/../numeric/mpfloat.h \
  ../../include/python/../numeric/dyadic.h \
  ../../include/python/../numeric/../base/dyadic.h \
  ../../include/python/../numeric/rational.h \
  ../../include/python/../declarations.h

export_affine_vector_field.o: export_affine_vector_field.cc \
  ../../include/numeric/interval.h \
  ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/approximation.h \
  ../../include/numeric/../numeric/arithmetic.h \
  ../../include/numeric/../numeric/rounding.h \
  ../../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/../numeric/function.h \
  ../../include/linear_algebra/matrix.h \
  ../../include/linear_algebra/../declarations.h \
  ../../include/linear_algebra/../numeric/integer.h \
  ../../include/linear_algebra/../numeric/../declarations.h \
  ../../include/linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/linear_algebra/../numeric/../numeric/arithmetic.h \
  ../../include/linear_algebra/../numeric/interval.h \
  ../../include/geometry/point.h ../../include/geometry/../declarations.h \
  ../../include/geometry/../base/array.h \
  ../../include/geometry/../base/../declarations.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../../include/geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../linear_algebra/matrix.h \
  ../../include/geometry/rectangle.h \
  ../../include/geometry/../numeric/arithmetic.h \
  ../../include/geometry/../numeric/function.h \
  ../../include/geometry/../numeric/interval.h \
  ../../include/geometry/../geometry/ppl_polyhedron.h \
  ../../include/geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../numeric/rational.h \
  ../../include/geometry/../geometry/../numeric/../declarations.h \
  ../../include/geometry/../geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../geometry/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../numeric/../numeric/function.h \
  ../../include/geometry/../geometry/../numeric/../numeric/integer.h \
  ../../include/geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/../geometry/rectangle_expression.h \
  ../../include/geometry/parallelotope.h \
  ../../include/geometry/../geometry/zonotope.h \
  ../../include/geometry/../geometry/../numeric/interval.h \
  ../../include/geometry/../geometry/../geometry/ppl_polyhedron.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/simplex.h \
  ../../include/geometry/../geometry/point_list.h \
  ../../include/geometry/../geometry/polytope.h \
  ../../include/geometry/../geometry/../geometry/point_list.h \
  ../../include/geometry/../geometry/../geometry/rectangle.h \
  ../../include/geometry/zonotope.h ../../include/geometry/polytope.h \
  ../../include/system/affine_vector_field.h \
  ../../include/system/../linear_algebra/vector.h \
  ../../include/system/../linear_algebra/matrix.h \
  ../../include/system/../system/vector_field.h \
  ../../include/system/../system/../declarations.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../numeric/numerical_types.h \
  ../../include/python/../numeric/mpfloat.h /home/collins/include/mpfr.h \
  ../../include/python/../numeric/../declarations.h \
  ../../include/python/../numeric/../numeric/numerical_traits.h \
  ../../include/python/../numeric/../numeric/function.h \
  ../../include/python/../numeric/../numeric/rational.h \
  ../../include/python/../numeric/float64.h \
  ../../include/python/../numeric/../numeric/integer.h \
  ../../include/python/../numeric/mpfloat.h \
  ../../include/python/../numeric/dyadic.h \
  ../../include/python/../numeric/../base/dyadic.h \
  ../../include/python/../numeric/rational.h \
  ../../include/python/../declarations.h \
  ../../include/python/python_utilities.h

export_henon_map.o: export_henon_map.cc ../../include/system/henon_map.h \
  ../../include/system/../linear_algebra/matrix.h \
  ../../include/system/../linear_algebra/../declarations.h \
  ../../include/system/../linear_algebra/../numeric/integer.h \
  ../../include/system/../linear_algebra/../numeric/../declarations.h \
  ../../include/system/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/system/../linear_algebra/../numeric/../numeric/approximation.h \
  ../../include/system/../linear_algebra/../numeric/../numeric/arithmetic.h \
  ../../include/system/../linear_algebra/../numeric/interval.h \
  ../../include/system/../linear_algebra/../numeric/../numeric/rounding.h \
  ../../include/system/../linear_algebra/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/system/../linear_algebra/../numeric/../numeric/../numeric/function.h \
  ../../include/system/../geometry/point.h \
  ../../include/system/../geometry/../declarations.h \
  ../../include/system/../geometry/../base/array.h \
  ../../include/system/../geometry/../base/../declarations.h \
  ../../include/system/../geometry/../linear_algebra/vector.h \
  ../../include/system/../geometry/../linear_algebra/../declarations.h \
  ../../include/system/../geometry/../linear_algebra/../numeric/integer.h \
  ../../include/system/../geometry/../linear_algebra/../numeric/arithmetic.h \
  ../../include/system/../geometry/../linear_algebra/../numeric/interval.h \
  ../../include/system/../geometry/../linear_algebra/matrix.h \
  ../../include/system/../geometry/rectangle.h \
  ../../include/system/../geometry/../numeric/arithmetic.h \
  ../../include/system/../geometry/../numeric/function.h \
  ../../include/system/../geometry/../numeric/interval.h \
  ../../include/system/../geometry/../geometry/ppl_polyhedron.h \
  ../../include/system/../geometry/../geometry/../declarations.h \
  ../../include/system/../geometry/../geometry/../numeric/rational.h \
  ../../include/system/../geometry/../geometry/../numeric/../declarations.h \
  ../../include/system/../geometry/../geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/system/../geometry/../geometry/../numeric/../numeric/arithmetic.h \
  ../../include/system/../geometry/../geometry/../numeric/../numeric/function.h \
  ../../include/system/../geometry/../geometry/../numeric/../numeric/integer.h \
  ../../include/system/../geometry/../geometry/../linear_algebra/vector.h \
  ../../include/system/../geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/system/../geometry/../geometry/point.h \
  ../../include/system/../geometry/../geometry/rectangle_expression.h \
  ../../include/system/../geometry/parallelotope.h \
  ../../include/system/../geometry/../geometry/zonotope.h \
  ../../include/system/../geometry/../geometry/../numeric/interval.h \
  ../../include/system/../geometry/../geometry/../geometry/ppl_polyhedron.h \
  ../../include/system/../geometry/../geometry/../geometry/point.h \
  ../../include/system/../system/map.h \
  ../../include/system/../system/../declarations.h \
  ../../include/system/../evaluation/apply.h \
  ../../include/system/../evaluation/../declarations.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../numeric/numerical_types.h \
  ../../include/python/../numeric/mpfloat.h /home/collins/include/mpfr.h \
  ../../include/python/../numeric/../declarations.h \
  ../../include/python/../numeric/../numeric/numerical_traits.h \
  ../../include/python/../numeric/../numeric/function.h \
  ../../include/python/../numeric/../numeric/rational.h \
  ../../include/python/../numeric/float64.h \
  ../../include/python/../numeric/../numeric/integer.h \
  ../../include/python/../numeric/mpfloat.h \
  ../../include/python/../numeric/dyadic.h \
  ../../include/python/../numeric/../base/dyadic.h \
  ../../include/python/../numeric/rational.h \
  ../../include/python/../declarations.h

export_lorenz_system.o: export_lorenz_system.cc \
  ../../include/system/lorenz_system.h \
  ../../include/system/../linear_algebra/vector.h \
  ../../include/system/../linear_algebra/../declarations.h \
  ../../include/system/../linear_algebra/../numeric/integer.h \
  ../../include/system/../linear_algebra/../numeric/../declarations.h \
  ../../include/system/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/system/../linear_algebra/../numeric/../numeric/approximation.h \
  ../../include/system/../linear_algebra/../numeric/../numeric/arithmetic.h \
  ../../include/system/../linear_algebra/../numeric/arithmetic.h \
  ../../include/system/../linear_algebra/../numeric/interval.h \
  ../../include/system/../linear_algebra/../numeric/../numeric/rounding.h \
  ../../include/system/../linear_algebra/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/system/../linear_algebra/../numeric/../numeric/../numeric/function.h \
  ../../include/system/../linear_algebra/matrix.h \
  ../../include/system/../geometry/point.h \
  ../../include/system/../geometry/../declarations.h \
  ../../include/system/../geometry/../base/array.h \
  ../../include/system/../geometry/../base/../declarations.h \
  ../../include/system/../geometry/../linear_algebra/vector.h \
  ../../include/system/../geometry/../linear_algebra/matrix.h \
  ../../include/system/../geometry/rectangle.h \
  ../../include/system/../geometry/../numeric/arithmetic.h \
  ../../include/system/../geometry/../numeric/function.h \
  ../../include/system/../geometry/../numeric/interval.h \
  ../../include/system/../geometry/../geometry/ppl_polyhedron.h \
  ../../include/system/../geometry/../geometry/../declarations.h \
  ../../include/system/../geometry/../geometry/../numeric/rational.h \
  ../../include/system/../geometry/../geometry/../numeric/../declarations.h \
  ../../include/system/../geometry/../geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/system/../geometry/../geometry/../numeric/../numeric/arithmetic.h \
  ../../include/system/../geometry/../geometry/../numeric/../numeric/function.h \
  ../../include/system/../geometry/../geometry/../numeric/../numeric/integer.h \
  ../../include/system/../geometry/../geometry/../linear_algebra/vector.h \
  ../../include/system/../geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/system/../geometry/../geometry/point.h \
  ../../include/system/../geometry/../geometry/rectangle_expression.h \
  ../../include/system/../geometry/parallelotope.h \
  ../../include/system/../geometry/../geometry/zonotope.h \
  ../../include/system/../geometry/../geometry/../numeric/interval.h \
  ../../include/system/../geometry/../geometry/../geometry/ppl_polyhedron.h \
  ../../include/system/../geometry/../geometry/../geometry/point.h \
  ../../include/system/../system/vector_field.h \
  ../../include/system/../system/../declarations.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../numeric/numerical_types.h \
  ../../include/python/../numeric/mpfloat.h /home/collins/include/mpfr.h \
  ../../include/python/../numeric/../declarations.h \
  ../../include/python/../numeric/../numeric/numerical_traits.h \
  ../../include/python/../numeric/../numeric/function.h \
  ../../include/python/../numeric/../numeric/rational.h \
  ../../include/python/../numeric/float64.h \
  ../../include/python/../numeric/../numeric/integer.h \
  ../../include/python/../numeric/mpfloat.h \
  ../../include/python/../numeric/dyadic.h \
  ../../include/python/../numeric/../base/dyadic.h \
  ../../include/python/../numeric/rational.h \
  ../../include/python/../declarations.h

linear_algebra_module.o: linear_algebra_module.cc

export_vector.o: export_vector.cc ../../include/numeric/numerical_types.h \
  ../../include/numeric/mpfloat.h /home/collins/include/mpfr.h \
  ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/approximation.h \
  ../../include/numeric/../numeric/function.h \
  ../../include/numeric/../numeric/rational.h \
  ../../include/numeric/../numeric/../declarations.h \
  ../../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/../numeric/arithmetic.h \
  ../../include/numeric/../numeric/../numeric/function.h \
  ../../include/numeric/../numeric/../numeric/integer.h \
  ../../include/numeric/../numeric/../numeric/../declarations.h \
  ../../include/numeric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/../numeric/../numeric/arithmetic.h \
  ../../include/numeric/interval.h \
  ../../include/numeric/../numeric/arithmetic.h \
  ../../include/numeric/../numeric/rounding.h \
  ../../include/linear_algebra/vector.h \
  ../../include/linear_algebra/../declarations.h \
  ../../include/linear_algebra/../numeric/integer.h \
  ../../include/linear_algebra/../numeric/arithmetic.h \
  ../../include/linear_algebra/../numeric/interval.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../numeric/numerical_types.h \
  ../../include/python/../numeric/float64.h \
  ../../include/python/../numeric/../declarations.h \
  ../../include/python/../numeric/../numeric/numerical_traits.h \
  ../../include/python/../numeric/../numeric/function.h \
  ../../include/python/../numeric/../numeric/integer.h \
  ../../include/python/../numeric/mpfloat.h \
  ../../include/python/../numeric/dyadic.h \
  ../../include/python/../numeric/../base/dyadic.h \
  ../../include/python/../numeric/rational.h \
  ../../include/python/../declarations.h

export_matrix.o: export_matrix.cc ../../include/numeric/interval.h \
  ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/approximation.h \
  ../../include/numeric/../numeric/arithmetic.h \
  ../../include/numeric/../numeric/rounding.h \
  ../../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/../numeric/function.h \
  ../../include/numeric/integer.h ../../include/linear_algebra/vector.h \
  ../../include/linear_algebra/../declarations.h \
  ../../include/linear_algebra/../numeric/integer.h \
  ../../include/linear_algebra/../numeric/arithmetic.h \
  ../../include/linear_algebra/../numeric/interval.h \
  ../../include/linear_algebra/matrix.h \
  ../../include/linear_algebra/matrix_function.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../numeric/numerical_types.h \
  ../../include/python/../numeric/mpfloat.h /home/collins/include/mpfr.h \
  ../../include/python/../numeric/../declarations.h \
  ../../include/python/../numeric/../numeric/numerical_traits.h \
  ../../include/python/../numeric/../numeric/function.h \
  ../../include/python/../numeric/../numeric/rational.h \
  ../../include/python/../numeric/../numeric/../declarations.h \
  ../../include/python/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/python/../numeric/../numeric/../numeric/arithmetic.h \
  ../../include/python/../numeric/../numeric/../numeric/function.h \
  ../../include/python/../numeric/../numeric/../numeric/integer.h \
  ../../include/python/../numeric/float64.h \
  ../../include/python/../numeric/../numeric/integer.h \
  ../../include/python/../numeric/mpfloat.h \
  ../../include/python/../numeric/dyadic.h \
  ../../include/python/../numeric/../base/dyadic.h \
  ../../include/python/../numeric/rational.h \
  ../../include/python/../declarations.h \
  ../../include/python/python_utilities.h

export_linear_program.o: export_linear_program.cc \
  ../../include/linear_algebra/vector.h \
  ../../include/linear_algebra/../declarations.h \
  ../../include/linear_algebra/../numeric/integer.h \
  ../../include/linear_algebra/../numeric/../declarations.h \
  ../../include/linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/linear_algebra/../numeric/../numeric/approximation.h \
  ../../include/linear_algebra/../numeric/../numeric/arithmetic.h \
  ../../include/linear_algebra/../numeric/arithmetic.h \
  ../../include/linear_algebra/../numeric/interval.h \
  ../../include/linear_algebra/../numeric/../numeric/rounding.h \
  ../../include/linear_algebra/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/linear_algebra/../numeric/../numeric/../numeric/function.h \
  ../../include/linear_algebra/linear_program.h \
  ../../include/linear_algebra/../linear_algebra/matrix.h \
  ../../include/linear_algebra/../linear_algebra/../declarations.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/integer.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/interval.h \
  ../../include/python/typedefs.h ../../include/python/../real_typedef.h \
  ../../include/python/../numeric/numerical_types.h \
  ../../include/python/../numeric/mpfloat.h /home/collins/include/mpfr.h \
  ../../include/python/../numeric/../declarations.h \
  ../../include/python/../numeric/../numeric/numerical_traits.h \
  ../../include/python/../numeric/../numeric/function.h \
  ../../include/python/../numeric/../numeric/rational.h \
  ../../include/python/../numeric/../numeric/../declarations.h \
  ../../include/python/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/python/../numeric/../numeric/../numeric/arithmetic.h \
  ../../include/python/../numeric/../numeric/../numeric/function.h \
  ../../include/python/../numeric/../numeric/../numeric/integer.h \
  ../../include/python/../numeric/float64.h \
  ../../include/python/../numeric/../numeric/integer.h \
  ../../include/python/../numeric/mpfloat.h \
  ../../include/python/../numeric/dyadic.h \
  ../../include/python/../numeric/../base/dyadic.h \
  ../../include/python/../numeric/rational.h \
  ../../include/python/../declarations.h

