include ../../config.mk

OBJ = map.o vector_field.o \
	affine_map.o affine_multimap.o affine_vector_field.o \
	polynomial_map.o \
	lattice_map.o

all: ${OBJ}

.cc.o: 
	${CXX} ${INCLUDE} ${CXXFLAGS} -c $<

clean:                                                                          
	rm -f *.bak *.o *~ *% core libariadne.*

dep:  
	@sed -e '/^### Do Not edit this line$$/,$$d' Makefile \
		> MMM.$$$$ && mv MMM.$$$$ Makefile
	@echo "### Do Not edit this line" >> Makefile
	@echo "### Everything below is auto-generated" >> Makefile
	@for ff in ${OBJ}; \
	do f=`basename $$ff .o`.cc ; \
	g++ -MM ${CXXFLAGS} $$f; echo ""; \
	done >> Makefile

### Do Not edit this line
### Everything below is auto-generated
map.o: map.cc ../../include/system/map.h \
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
  ../../include/system/map.tpl ../../include/system/../geometry/point.h \
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
  ../../include/system/../geometry/../linear_algebra/../declarations.h \
  ../../include/system/../geometry/../linear_algebra/../numeric/interval.h \
  ../../include/system/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/system/../geometry/../geometry/point.h \
  ../../include/system/../geometry/list_set.h \
  ../../include/system/../geometry/../base/utility.h \
  ../../include/system/../geometry/../base/../declarations.h \
  ../../include/system/../geometry/parallelotope.h \
  ../../include/system/../geometry/../linear_algebra/matrix.h \
  ../../include/system/../geometry/polyhedron.h \
  ../../include/system/../geometry/../linear_algebra/constraint_system.h \
  ../../include/system/../geometry/../linear_algebra/../numeric/numerical_types.h \
  ../../include/system/../geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/system/../geometry/../linear_algebra/generator_system.h \
  ../../include/system/../geometry/../geometry/rectangle.h \
  ../../include/system/../system/map.h ../../include/real_typedef.h \
  ../../include/numeric/numerical_types.h

vector_field.o: vector_field.cc ../../include/system/vector_field.h \
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
  ../../include/system/vector_field.tpl \
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
  ../../include/system/../system/vector_field.h \
  ../../include/real_typedef.h ../../include/numeric/numerical_types.h

affine_map.o: affine_map.cc ../../include/system/affine_map.h \
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
  ../../include/system/../geometry/point.h \
  ../../include/system/../geometry/../declarations.h \
  ../../include/system/../geometry/../base/array.h \
  ../../include/system/../geometry/../linear_algebra/vector.h \
  ../../include/system/../geometry/simplex.h \
  ../../include/system/../geometry/../geometry/point.h \
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
  ../../include/system/../geometry/../linear_algebra/../declarations.h \
  ../../include/system/../geometry/../linear_algebra/../numeric/interval.h \
  ../../include/system/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/system/../geometry/parallelotope.h \
  ../../include/system/../geometry/../linear_algebra/matrix.h \
  ../../include/system/../geometry/zonotope.h \
  ../../include/system/../geometry/../linear_algebra/transformation_system.h \
  ../../include/system/../geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/system/../geometry/../linear_algebra/../linear_algebra/interval_vector.h \
  ../../include/system/../geometry/../linear_algebra/../linear_algebra/interval_matrix.h \
  ../../include/system/../geometry/polyhedron.h \
  ../../include/system/../geometry/../linear_algebra/constraint_system.h \
  ../../include/system/../geometry/../linear_algebra/../numeric/numerical_types.h \
  ../../include/system/../geometry/../linear_algebra/generator_system.h \
  ../../include/system/../geometry/../geometry/rectangle.h \
  ../../include/system/../geometry/list_set.h \
  ../../include/system/../geometry/../base/utility.h \
  ../../include/system/../geometry/../base/../declarations.h \
  ../../include/system/../geometry/grid_set.h \
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
  ../../include/system/affine_map.tpl \
  ../../include/system/../system/affine_map.h \
  ../../include/real_typedef.h ../../include/numeric/numerical_types.h

affine_multimap.o: affine_multimap.cc \
  ../../include/system/affine_multimap.h \
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
  ../../include/system/../geometry/point.h \
  ../../include/system/../geometry/../declarations.h \
  ../../include/system/../geometry/../base/array.h \
  ../../include/system/../geometry/../linear_algebra/vector.h \
  ../../include/system/../geometry/simplex.h \
  ../../include/system/../geometry/../geometry/point.h \
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
  ../../include/system/../geometry/../linear_algebra/../declarations.h \
  ../../include/system/../geometry/../linear_algebra/../numeric/interval.h \
  ../../include/system/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/system/../geometry/parallelotope.h \
  ../../include/system/../geometry/../linear_algebra/matrix.h \
  ../../include/system/../geometry/zonotope.h \
  ../../include/system/../geometry/../linear_algebra/transformation_system.h \
  ../../include/system/../geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/system/../geometry/../linear_algebra/../linear_algebra/interval_vector.h \
  ../../include/system/../geometry/../linear_algebra/../linear_algebra/interval_matrix.h \
  ../../include/system/../geometry/polyhedron.h \
  ../../include/system/../geometry/../linear_algebra/constraint_system.h \
  ../../include/system/../geometry/../linear_algebra/../numeric/numerical_types.h \
  ../../include/system/../geometry/../linear_algebra/generator_system.h \
  ../../include/system/../geometry/../geometry/rectangle.h \
  ../../include/system/../geometry/list_set.h \
  ../../include/system/../geometry/../base/utility.h \
  ../../include/system/../geometry/../base/../declarations.h \
  ../../include/system/../geometry/grid_set.h \
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
  ../../include/system/affine_multimap.tpl \
  ../../include/system/../system/affine_multimap.h \
  ../../include/real_typedef.h ../../include/numeric/numerical_types.h

affine_vector_field.o: affine_vector_field.cc \
  ../../include/system/affine_vector_field.h \
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
  ../../include/system/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/system/../linear_algebra/../linear_algebra/vector.h \
  ../../include/system/../linear_algebra/interval_matrix.h \
  ../../include/system/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/system/../system/vector_field.h \
  ../../include/system/../system/../declarations.h \
  ../../include/system/../system/../numeric/interval.h \
  ../../include/system/../system/../linear_algebra/vector.h \
  ../../include/system/../system/../linear_algebra/matrix.h \
  ../../include/system/../system/../linear_algebra/interval_vector.h \
  ../../include/system/../system/../linear_algebra/interval_matrix.h \
  ../../include/system/../system/affine_map.h \
  ../../include/system/../system/../geometry/point.h \
  ../../include/system/../system/../geometry/../declarations.h \
  ../../include/system/../system/../geometry/../base/array.h \
  ../../include/system/../system/../geometry/../linear_algebra/vector.h \
  ../../include/system/../system/../geometry/simplex.h \
  ../../include/system/../system/../geometry/../geometry/point.h \
  ../../include/system/../system/../geometry/rectangle.h \
  ../../include/system/../system/../geometry/../numeric/arithmetic.h \
  ../../include/system/../system/../geometry/../numeric/../numeric/numerical_types.h \
  ../../include/system/../system/../geometry/../numeric/../numeric/../declarations.h \
  ../../include/system/../system/../geometry/../numeric/../numeric/../numeric/integer.h \
  ../../include/system/../system/../geometry/../numeric/../numeric/../numeric/dyadic.h \
  ../../include/system/../system/../geometry/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/system/../system/../geometry/../numeric/../numeric/../utility/stlio.h \
  ../../include/system/../system/../geometry/../numeric/../numeric/approximation.h \
  ../../include/system/../system/../geometry/../numeric/../numeric/../numeric/numerical_types.h \
  ../../include/system/../system/../geometry/../numeric/interval.h \
  ../../include/system/../system/../geometry/../linear_algebra/interval_vector.h \
  ../../include/system/../system/../geometry/parallelotope.h \
  ../../include/system/../system/../geometry/../linear_algebra/matrix.h \
  ../../include/system/../system/../geometry/zonotope.h \
  ../../include/system/../system/../geometry/../linear_algebra/transformation_system.h \
  ../../include/system/../system/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/system/../system/../geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/system/../system/../geometry/../linear_algebra/../linear_algebra/interval_vector.h \
  ../../include/system/../system/../geometry/../linear_algebra/../linear_algebra/interval_matrix.h \
  ../../include/system/../system/../geometry/polyhedron.h \
  ../../include/system/../system/../geometry/../linear_algebra/constraint_system.h \
  ../../include/system/../system/../geometry/../linear_algebra/../declarations.h \
  ../../include/system/../system/../geometry/../linear_algebra/../numeric/numerical_types.h \
  ../../include/system/../system/../geometry/../linear_algebra/generator_system.h \
  ../../include/system/../system/../geometry/../geometry/rectangle.h \
  ../../include/system/../system/../geometry/list_set.h \
  ../../include/system/../system/../geometry/../base/utility.h \
  ../../include/system/../system/../geometry/../base/../declarations.h \
  ../../include/system/../system/../geometry/grid_set.h \
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
  ../../include/system/affine_vector_field.tpl \
  ../../include/system/../system/affine_vector_field.h \
  ../../include/real_typedef.h ../../include/numeric/numerical_types.h

polynomial_map.o: polynomial_map.cc ../../include/system/polynomial_map.h \
  ../../include/system/../system/map.h \
  ../../include/system/../system/../declarations.h \
  ../../include/system/../system/../linear_algebra/vector.h \
  ../../include/system/../system/../linear_algebra/../declarations.h \
  ../../include/system/../system/../linear_algebra/../numeric/integer.h \
  ../../include/system/../system/../linear_algebra/../numeric/../declarations.h \
  ../../include/system/../system/../linear_algebra/../numeric/../utility/stlio.h \
  ../../include/system/../system/../linear_algebra/../numeric/../utility/../base/array.h \
  ../../include/system/../system/../linear_algebra/../numeric/../utility/../base/../declarations.h \
  ../../include/system/../system/../linear_algebra/matrix.h \
  ../../include/system/../system/../linear_algebra/interval_matrix.h \
  ../../include/system/../system/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/system/../system/../linear_algebra/../numeric/interval.h \
  ../../include/system/../system/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/system/../linear_algebra/matrix.h \
  ../../include/system/polynomial_map.tpl \
  ../../include/system/polynomial_map.h \
  ../../include/system/../linear_algebra/vector.h \
  ../../include/system/../utility/stlio.h \
  ../../include/system/../base/array.h \
  ../../include/system/../numeric/interval.h \
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
  ../../include/system/../geometry/../linear_algebra/../declarations.h \
  ../../include/system/../geometry/../linear_algebra/../numeric/interval.h \
  ../../include/system/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/system/../geometry/../geometry/point.h \
  ../../include/system/../geometry/polyhedron.h \
  ../../include/system/../geometry/../linear_algebra/constraint_system.h \
  ../../include/system/../geometry/../linear_algebra/../numeric/numerical_types.h \
  ../../include/system/../geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/system/../geometry/../linear_algebra/matrix.h \
  ../../include/system/../geometry/../linear_algebra/generator_system.h \
  ../../include/system/../geometry/../geometry/rectangle.h \
  ../../include/system/../geometry/list_set.h \
  ../../include/system/../geometry/../base/utility.h \
  ../../include/system/../geometry/../base/../declarations.h \
  ../../include/real_typedef.h ../../include/numeric/numerical_types.h

lattice_map.o: lattice_map.cc ../../include/system/lattice_map.h \
  ../../include/system/../geometry/lattice_set.h \
  ../../include/system/../geometry/../base/array.h \
  ../../include/system/../geometry/../base/../declarations.h \
  ../../include/system/../geometry/../base/iterator.h \
  ../../include/system/../geometry/../base/../base/array.h \
  ../../include/system/../geometry/../numeric/interval.h \
  ../../include/system/../geometry/../numeric/../declarations.h \
  ../../include/system/../geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/system/../geometry/../base/array_operations.h \
  ../../include/utility/stlio.h ../../include/utility/../base/array.h \
  ../../include/geometry/lattice_set.h ../../include/debug.h

