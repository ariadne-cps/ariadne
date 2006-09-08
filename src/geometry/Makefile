include ../../config.mk

OBJ = point.o \
	ppl_polyhedron.o \
  rectangle.o polyhedron.o simplex.o zonotope.o parallelotope.o \
	sphere.o ellipsoid.o \
  list_set.o \
	lattice_set.o grid.o grid_set.o \
	subdivision_tree_set.o partition_tree_set.o

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

depclean:
	@sed -e '/^### Do Not edit this line$$/,$$d' Makefile \
		> MMM.$$$$ && mv MMM.$$$$ Makefile
	@echo "### Do Not edit this line" >> Makefile

### Do Not edit this line
### Everything below is auto-generated
point.o: point.cc ../../include/geometry/point.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../base/array.h \
  ../../include/geometry/../base/../declarations.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../linear_algebra/../numeric/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/integer.h \
  ../../include/geometry/../linear_algebra/matrix.h \
  ../../include/geometry/point.tpl ../../include/geometry/point.h \
  ../../include/geometry/../numeric/numerical_types.h \
  ../../include/geometry/../numeric/mpfloat.h \
  ../../include/geometry/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../numeric/../numeric/function.h \
  ../../include/geometry/../numeric/../numeric/integer.h \
  ../../include/geometry/../numeric/../numeric/rational.h \
  ../../include/geometry/../numeric/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../numeric/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../numeric/../numeric/../numeric/function.h \
  ../../include/geometry/../numeric/../numeric/../numeric/integer.h \
  ../../include/geometry/../base/stlio.h \
  ../../include/geometry/../base/../base/array.h \
  ../../include/numeric/float64.h ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/function.h \
  ../../include/numeric/../numeric/integer.h \
  ../../include/numeric/mpfloat.h ../../include/numeric/rational.h

ppl_polyhedron.o: ppl_polyhedron.cc \
  ../../include/numeric/numerical_types.h ../../include/numeric/mpfloat.h \
  ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/function.h \
  ../../include/numeric/../numeric/integer.h \
  ../../include/numeric/../numeric/../declarations.h \
  ../../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/rational.h \
  ../../include/numeric/../numeric/../numeric/arithmetic.h \
  ../../include/numeric/../numeric/../numeric/../declarations.h \
  ../../include/numeric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/../numeric/../numeric/integer.h \
  ../../include/numeric/../numeric/../numeric/function.h \
  ../../include/numeric/../numeric/../numeric/integer.h \
  ../../include/numeric/interval.h \
  ../../include/numeric/../numeric/arithmetic.h \
  ../../include/numeric/../numeric/rounding.h \
  ../../include/geometry/ppl_polyhedron.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../numeric/rational.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../../include/geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../linear_algebra/interval_vector.h \
  ../../include/geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/vector.h

rectangle.o: rectangle.cc ../../include/geometry/rectangle.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../numeric/arithmetic.h \
  ../../include/geometry/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../numeric/../numeric/integer.h \
  ../../include/geometry/../numeric/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../numeric/interval.h \
  ../../include/geometry/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../numeric/../numeric/rounding.h \
  ../../include/geometry/../numeric/../numeric/../numeric/function.h \
  ../../include/geometry/../numeric/../numeric/../numeric/integer.h \
  ../../include/geometry/../linear_algebra/interval_vector.h \
  ../../include/geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../numeric/arithmetic.h \
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
  ../../include/geometry/../geometry/../linear_algebra/interval_vector.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../base/../declarations.h \
  ../../include/geometry/rectangle.tpl ../../include/geometry/rectangle.h \
  ../../include/geometry/../base/binary_word.h \
  ../../include/geometry/../base/../declarations.h \
  ../../include/geometry/../base/array.h \
  ../../include/geometry/../geometry/point_list.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/list_set.h \
  ../../include/geometry/../geometry/../base/utility.h \
  ../../include/geometry/../geometry/grid_set.h \
  ../../include/geometry/../geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/../geometry/grid.h \
  ../../include/geometry/../geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../geometry/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../geometry/../geometry/lattice_set.h \
  ../../include/geometry/../geometry/../geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../geometry/../geometry/../base/iterator.h \
  ../../include/geometry/../geometry/../geometry/../geometry/../base/../base/array.h \
  ../../include/geometry/../geometry/../geometry/../geometry/../numeric/interval.h \
  ../../include/geometry/../geometry/../geometry/../geometry/../base/array_operations.h \
  ../../include/geometry/../geometry/../geometry/../geometry/../base/../declarations.h \
  ../../include/geometry/../geometry/../geometry/lattice_set.h \
  ../../include/numeric/float64.h ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/function.h \
  ../../include/numeric/../numeric/integer.h \
  ../../include/numeric/rational.h ../../include/real_typedef.h \
  ../../include/numeric/numerical_types.h ../../include/numeric/mpfloat.h \
  ../../include/numeric/../numeric/rational.h

polyhedron.o: polyhedron.cc ../../include/geometry/polyhedron.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../linear_algebra/../numeric/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/integer.h \
  ../../include/geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../linear_algebra/interval_vector.h \
  ../../include/geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/rounding.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/function.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/integer.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/vector.h \
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
  ../../include/geometry/../geometry/../linear_algebra/interval_vector.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../base/../declarations.h \
  ../../include/geometry/../geometry/point_list.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../numeric/interval.h \
  ../../include/geometry/../geometry/../geometry/ppl_polyhedron.h \
  ../../include/geometry/polyhedron.tpl \
  ../../include/geometry/polyhedron.h \
  ../../include/geometry/../base/stlio.h \
  ../../include/geometry/../base/../base/array.h \
  ../../include/geometry/../numeric/interval.h \
  ../../include/geometry/../numeric/arithmetic.h \
  ../../include/real_typedef.h ../../include/numeric/numerical_types.h \
  ../../include/numeric/mpfloat.h ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/function.h \
  ../../include/numeric/../numeric/rational.h

simplex.o: simplex.cc ../../include/geometry/simplex.h \
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
  ../../include/geometry/../geometry/../linear_algebra/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/../numeric/integer.h \
  ../../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/point_list.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/polyhedron.h \
  ../../include/geometry/../geometry/../linear_algebra/interval_vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/../numeric/rounding.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/../numeric/../numeric/function.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/../numeric/../numeric/integer.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/vector.h \
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
  ../../include/geometry/../geometry/../geometry/../linear_algebra/interval_vector.h \
  ../../include/geometry/../geometry/../geometry/point_list.h \
  ../../include/geometry/../geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/../geometry/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../geometry/../numeric/interval.h \
  ../../include/geometry/../geometry/../geometry/../geometry/ppl_polyhedron.h \
  ../../include/geometry/../geometry/../geometry/../geometry/point.h \
  ../../include/geometry/simplex.tpl ../../include/geometry/simplex.h \
  ../../include/geometry/../base/stlio.h \
  ../../include/geometry/../base/../base/array.h \
  ../../include/real_typedef.h ../../include/numeric/numerical_types.h \
  ../../include/numeric/mpfloat.h ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/function.h \
  ../../include/numeric/../numeric/rational.h

zonotope.o: zonotope.cc ../../include/geometry/zonotope.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../linear_algebra/../numeric/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/integer.h \
  ../../include/geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../numeric/interval.h \
  ../../include/geometry/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../numeric/../numeric/rounding.h \
  ../../include/geometry/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../numeric/../numeric/../numeric/function.h \
  ../../include/geometry/../numeric/../numeric/../numeric/integer.h \
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
  ../../include/geometry/../geometry/../linear_algebra/interval_vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../base/../declarations.h \
  ../../include/geometry/zonotope.tpl ../../include/geometry/zonotope.h \
  ../../include/geometry/../base/array.h \
  ../../include/geometry/../numeric/approximation.h \
  ../../include/geometry/../linear_algebra/interval_vector.h \
  ../../include/geometry/../linear_algebra/interval_matrix.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../geometry/lattice_set.h \
  ../../include/geometry/../geometry/../base/iterator.h \
  ../../include/geometry/../geometry/../base/../base/array.h \
  ../../include/geometry/../geometry/../numeric/interval.h \
  ../../include/geometry/../geometry/../base/array_operations.h \
  ../../include/geometry/../geometry/point_list.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../geometry/ppl_polyhedron.h \
  ../../include/geometry/../geometry/polyhedron.h \
  ../../include/geometry/../geometry/../geometry/point_list.h \
  ../../include/geometry/../geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/parallelotope.h \
  ../../include/geometry/../geometry/../geometry/zonotope.h \
  ../../include/geometry/../geometry/list_set.h \
  ../../include/geometry/../geometry/../base/utility.h \
  ../../include/geometry/../linear_algebra/matrix.tpl \
  ../../include/geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../linear_algebra/../base/array.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/lu_matrix.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../base/array.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/real_typedef.h ../../include/numeric/numerical_types.h \
  ../../include/numeric/mpfloat.h ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/function.h \
  ../../include/numeric/../numeric/rational.h

parallelotope.o: parallelotope.cc ../../include/geometry/parallelotope.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../numeric/interval.h \
  ../../include/geometry/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../numeric/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../numeric/../numeric/../numeric/integer.h \
  ../../include/geometry/../numeric/../numeric/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../numeric/../numeric/rounding.h \
  ../../include/geometry/../numeric/../numeric/../numeric/function.h \
  ../../include/geometry/../numeric/../numeric/../numeric/integer.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../linear_algebra/../numeric/arithmetic.h \
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
  ../../include/geometry/../geometry/../geometry/../linear_algebra/interval_vector.h \
  ../../include/geometry/../geometry/../geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../geometry/../geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../geometry/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/parallelotope.tpl \
  ../../include/geometry/parallelotope.h \
  ../../include/geometry/../base/stlio.h \
  ../../include/geometry/../base/../base/array.h \
  ../../include/geometry/../linear_algebra/interval_vector.h \
  ../../include/geometry/../linear_algebra/interval_matrix.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../linear_algebra/linear_program.h \
  ../../include/geometry/../geometry/lattice_set.h \
  ../../include/geometry/../geometry/../base/iterator.h \
  ../../include/geometry/../geometry/../base/../base/array.h \
  ../../include/geometry/../geometry/../base/array_operations.h \
  ../../include/geometry/../geometry/point_list.h \
  ../../include/geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../linear_algebra/interval_vector.h \
  ../../include/geometry/../geometry/list_set.h \
  ../../include/geometry/../geometry/../base/utility.h \
  ../../include/geometry/../geometry/polyhedron.h \
  ../../include/geometry/../geometry/../geometry/point_list.h \
  ../../include/geometry/../geometry/../geometry/rectangle.h \
  ../../include/real_typedef.h ../../include/numeric/numerical_types.h \
  ../../include/numeric/mpfloat.h ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/function.h \
  ../../include/numeric/../numeric/rational.h

sphere.o: sphere.cc ../../include/geometry/sphere.h \
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
  ../../include/geometry/../geometry/../linear_algebra/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/../numeric/integer.h \
  ../../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/geometry/sphere.tpl ../../include/geometry/sphere.h \
  ../../include/real_typedef.h ../../include/numeric/numerical_types.h \
  ../../include/numeric/mpfloat.h ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/function.h \
  ../../include/numeric/../numeric/integer.h \
  ../../include/numeric/../numeric/rational.h \
  ../../include/numeric/../numeric/../declarations.h \
  ../../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/../numeric/arithmetic.h \
  ../../include/numeric/../numeric/../numeric/function.h \
  ../../include/numeric/../numeric/../numeric/integer.h

ellipsoid.o: ellipsoid.cc ../../include/geometry/ellipsoid.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../base/array.h \
  ../../include/geometry/../base/../declarations.h \
  ../../include/geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../linear_algebra/../numeric/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/../declarations.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/../numeric/integer.h \
  ../../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/geometry/ellipsoid.tpl ../../include/geometry/ellipsoid.h \
  ../../include/geometry/../geometry/sphere.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/real_typedef.h ../../include/numeric/numerical_types.h \
  ../../include/numeric/mpfloat.h ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/function.h \
  ../../include/numeric/../numeric/integer.h \
  ../../include/numeric/../numeric/rational.h \
  ../../include/numeric/../numeric/../declarations.h \
  ../../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/../numeric/arithmetic.h \
  ../../include/numeric/../numeric/../numeric/function.h \
  ../../include/numeric/../numeric/../numeric/integer.h

list_set.o: list_set.cc ../../include/geometry/rectangle.tpl \
  ../../include/geometry/rectangle.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../numeric/arithmetic.h \
  ../../include/geometry/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../numeric/../numeric/integer.h \
  ../../include/geometry/../numeric/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../numeric/interval.h \
  ../../include/geometry/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../numeric/../numeric/rounding.h \
  ../../include/geometry/../numeric/../numeric/../numeric/function.h \
  ../../include/geometry/../numeric/../numeric/../numeric/integer.h \
  ../../include/geometry/../linear_algebra/interval_vector.h \
  ../../include/geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../numeric/arithmetic.h \
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
  ../../include/geometry/../geometry/../linear_algebra/interval_vector.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../base/../declarations.h \
  ../../include/geometry/../base/binary_word.h \
  ../../include/geometry/../base/../declarations.h \
  ../../include/geometry/../base/array.h \
  ../../include/geometry/../geometry/point_list.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/list_set.h \
  ../../include/geometry/../geometry/../base/utility.h \
  ../../include/geometry/../geometry/grid_set.h \
  ../../include/geometry/../geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/../geometry/grid.h \
  ../../include/geometry/../geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../geometry/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../geometry/../geometry/lattice_set.h \
  ../../include/geometry/../geometry/../geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../geometry/../geometry/../base/iterator.h \
  ../../include/geometry/../geometry/../geometry/../geometry/../base/../base/array.h \
  ../../include/geometry/../geometry/../geometry/../geometry/../numeric/interval.h \
  ../../include/geometry/../geometry/../geometry/../geometry/../base/array_operations.h \
  ../../include/geometry/../geometry/../geometry/../geometry/../base/../declarations.h \
  ../../include/geometry/../geometry/../geometry/lattice_set.h \
  ../../include/geometry/parallelotope.tpl \
  ../../include/geometry/parallelotope.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/zonotope.h \
  ../../include/geometry/../geometry/../numeric/interval.h \
  ../../include/geometry/../geometry/../geometry/ppl_polyhedron.h \
  ../../include/geometry/../base/stlio.h \
  ../../include/geometry/../base/../base/array.h \
  ../../include/geometry/../linear_algebra/interval_matrix.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../linear_algebra/linear_program.h \
  ../../include/geometry/../geometry/lattice_set.h \
  ../../include/geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/polyhedron.h \
  ../../include/geometry/../geometry/../geometry/point_list.h \
  ../../include/geometry/zonotope.tpl ../../include/geometry/zonotope.h \
  ../../include/geometry/../numeric/approximation.h \
  ../../include/geometry/../geometry/parallelotope.h \
  ../../include/geometry/../linear_algebra/matrix.tpl \
  ../../include/geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../linear_algebra/../base/array.h \
  ../../include/geometry/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/lu_matrix.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../base/array.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/list_set.h ../../include/geometry/list_set.tpl \
  ../../include/geometry/list_set.h ../../include/real_typedef.h \
  ../../include/numeric/numerical_types.h ../../include/numeric/mpfloat.h \
  ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/function.h \
  ../../include/numeric/../numeric/rational.h

lattice_set.o: lattice_set.cc ../../include/base/stlio.h \
  ../../include/base/../base/array.h \
  ../../include/base/../base/../declarations.h \
  ../../include/base/array_operations.h \
  ../../include/base/../declarations.h \
  ../../include/geometry/lattice_set.h \
  ../../include/geometry/../base/array.h \
  ../../include/geometry/../base/iterator.h \
  ../../include/geometry/../base/../base/array.h \
  ../../include/geometry/../numeric/interval.h \
  ../../include/geometry/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../numeric/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../numeric/../numeric/../numeric/integer.h \
  ../../include/geometry/../numeric/../numeric/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../numeric/../numeric/rounding.h \
  ../../include/geometry/../numeric/../numeric/../numeric/function.h \
  ../../include/geometry/../numeric/../numeric/../numeric/integer.h \
  ../../include/geometry/../base/array_operations.h

grid.o: grid.cc ../../include/geometry/grid.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../base/array.h \
  ../../include/geometry/../base/../declarations.h \
  ../../include/geometry/../numeric/arithmetic.h \
  ../../include/geometry/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../numeric/../numeric/integer.h \
  ../../include/geometry/../numeric/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../geometry/lattice_set.h \
  ../../include/geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../base/iterator.h \
  ../../include/geometry/../geometry/../base/../base/array.h \
  ../../include/geometry/../geometry/../numeric/interval.h \
  ../../include/geometry/../geometry/../numeric/../declarations.h \
  ../../include/geometry/../geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../geometry/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../numeric/../numeric/rounding.h \
  ../../include/geometry/../geometry/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../geometry/../numeric/../numeric/../numeric/function.h \
  ../../include/geometry/../geometry/../numeric/../numeric/../numeric/integer.h \
  ../../include/geometry/../geometry/../base/array_operations.h \
  ../../include/geometry/../geometry/../base/../declarations.h \
  ../../include/geometry/grid.tpl ../../include/geometry/grid.h \
  ../../include/geometry/../base/stlio.h \
  ../../include/geometry/../base/../base/array.h \
  ../../include/geometry/../geometry/grid_set.h \
  ../../include/geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../geometry/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../geometry/../numeric/interval.h \
  ../../include/geometry/../geometry/../geometry/../linear_algebra/interval_vector.h \
  ../../include/geometry/../geometry/../geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../geometry/../geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../geometry/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../geometry/../linear_algebra/../linear_algebra/../declarations.h \
  ../../include/geometry/../geometry/../geometry/../linear_algebra/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../geometry/../geometry/../linear_algebra/../linear_algebra/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../geometry/../geometry/ppl_polyhedron.h \
  ../../include/geometry/../geometry/../geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../geometry/../geometry/../numeric/rational.h \
  ../../include/geometry/../geometry/../geometry/../geometry/../numeric/../declarations.h \
  ../../include/geometry/../geometry/../geometry/../geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../geometry/../geometry/../geometry/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../geometry/../geometry/../numeric/../numeric/function.h \
  ../../include/geometry/../geometry/../geometry/../geometry/../numeric/../numeric/integer.h \
  ../../include/geometry/../geometry/../geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../geometry/../geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../geometry/../geometry/../geometry/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../geometry/../geometry/../geometry/../linear_algebra/interval_vector.h \
  ../../include/geometry/../geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../geometry/grid.h \
  ../../include/geometry/../geometry/../geometry/lattice_set.h \
  ../../include/geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/list_set.h \
  ../../include/geometry/../geometry/../base/utility.h \
  ../../include/real_typedef.h ../../include/numeric/numerical_types.h \
  ../../include/numeric/mpfloat.h ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/function.h \
  ../../include/numeric/../numeric/rational.h

grid_set.o: grid_set.cc ../../include/geometry/grid_set.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../base/array.h \
  ../../include/geometry/../base/../declarations.h \
  ../../include/geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../numeric/../declarations.h \
  ../../include/geometry/../geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../geometry/../numeric/../numeric/integer.h \
  ../../include/geometry/../geometry/../numeric/../numeric/../declarations.h \
  ../../include/geometry/../geometry/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../geometry/../numeric/interval.h \
  ../../include/geometry/../geometry/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../numeric/../numeric/rounding.h \
  ../../include/geometry/../geometry/../numeric/../numeric/../numeric/function.h \
  ../../include/geometry/../geometry/../numeric/../numeric/../numeric/integer.h \
  ../../include/geometry/../geometry/../linear_algebra/interval_vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/../declarations.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/../numeric/arithmetic.h \
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
  ../../include/geometry/../geometry/../geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../geometry/../geometry/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../geometry/../geometry/../linear_algebra/interval_vector.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/grid.h \
  ../../include/geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../geometry/lattice_set.h \
  ../../include/geometry/../geometry/../geometry/../base/iterator.h \
  ../../include/geometry/../geometry/../geometry/../base/../base/array.h \
  ../../include/geometry/../geometry/../geometry/../numeric/interval.h \
  ../../include/geometry/../geometry/../geometry/../base/array_operations.h \
  ../../include/geometry/../geometry/../geometry/../base/../declarations.h \
  ../../include/geometry/../geometry/lattice_set.h \
  ../../include/geometry/grid_set.tpl ../../include/geometry/grid_set.h \
  ../../include/geometry/../base/stlio.h \
  ../../include/geometry/../base/../base/array.h \
  ../../include/geometry/../base/array_operations.h \
  ../../include/geometry/../geometry/zonotope.h \
  ../../include/geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/parallelotope.h \
  ../../include/geometry/../geometry/../geometry/zonotope.h \
  ../../include/geometry/../geometry/polyhedron.h \
  ../../include/geometry/../geometry/../geometry/point_list.h \
  ../../include/geometry/../geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/list_set.h \
  ../../include/geometry/../geometry/../base/utility.h \
  ../../include/geometry/../geometry/../base/../declarations.h \
  ../../include/geometry/../geometry/partition_tree_set.h \
  ../../include/geometry/../geometry/../base/binary_word.h \
  ../../include/geometry/../geometry/../base/binary_tree.h \
  ../../include/geometry/../geometry/../base/../base/array.h \
  ../../include/geometry/../geometry/../base/../base/binary_word.h \
  ../../include/geometry/../geometry/../base/iterator.h \
  ../../include/geometry/../geometry/../geometry/subdivision_tree_set.h \
  ../../include/geometry/../geometry/../geometry/../base/sequence.h \
  ../../include/geometry/../geometry/../geometry/../base/binary_word.h \
  ../../include/geometry/../geometry/../geometry/../base/binary_tree.h \
  ../../include/geometry/../geometry/../geometry/../geometry/rectangle.h \
  ../../include/real_typedef.h ../../include/numeric/numerical_types.h \
  ../../include/numeric/mpfloat.h ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/function.h \
  ../../include/numeric/../numeric/rational.h

subdivision_tree_set.o: subdivision_tree_set.cc \
  ../../include/numeric/arithmetic.h \
  ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/integer.h \
  ../../include/numeric/../numeric/../declarations.h \
  ../../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/rectangle.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../numeric/arithmetic.h \
  ../../include/geometry/../numeric/interval.h \
  ../../include/geometry/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../numeric/../numeric/rounding.h \
  ../../include/geometry/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../numeric/../numeric/../numeric/function.h \
  ../../include/geometry/../numeric/../numeric/../numeric/integer.h \
  ../../include/geometry/../linear_algebra/interval_vector.h \
  ../../include/geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../numeric/arithmetic.h \
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
  ../../include/geometry/../geometry/../linear_algebra/interval_vector.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../base/../declarations.h \
  ../../include/geometry/lattice_set.h \
  ../../include/geometry/../base/array.h \
  ../../include/geometry/../base/iterator.h \
  ../../include/geometry/../base/../base/array.h \
  ../../include/geometry/../base/array_operations.h \
  ../../include/geometry/../base/../declarations.h \
  ../../include/geometry/subdivision_tree_set.h \
  ../../include/geometry/../base/sequence.h \
  ../../include/geometry/../base/binary_word.h \
  ../../include/geometry/../base/binary_tree.h \
  ../../include/geometry/../base/../base/binary_word.h \
  ../../include/geometry/../geometry/rectangle.h

partition_tree_set.o: partition_tree_set.cc \
  ../../include/geometry/parallelotope.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../numeric/interval.h \
  ../../include/geometry/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../numeric/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../numeric/../numeric/../numeric/integer.h \
  ../../include/geometry/../numeric/../numeric/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../numeric/../numeric/rounding.h \
  ../../include/geometry/../numeric/../numeric/../numeric/function.h \
  ../../include/geometry/../numeric/../numeric/../numeric/integer.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../linear_algebra/../numeric/arithmetic.h \
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
  ../../include/geometry/../geometry/../geometry/../linear_algebra/interval_vector.h \
  ../../include/geometry/../geometry/../geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../geometry/../geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../geometry/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/partition_tree_set.h \
  ../../include/geometry/../base/binary_word.h \
  ../../include/geometry/../base/../declarations.h \
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
  ../../include/geometry/../geometry/../geometry/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../geometry/../numeric/interval.h \
  ../../include/geometry/../geometry/../geometry/../geometry/ppl_polyhedron.h \
  ../../include/geometry/../geometry/../geometry/../geometry/point.h \
  ../../include/geometry/partition_tree_set.tpl \
  ../../include/geometry/partition_tree_set.h \
  ../../include/geometry/../base/stlio.h \
  ../../include/geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/list_set.h \
  ../../include/geometry/../geometry/../base/utility.h \
  ../../include/geometry/../geometry/grid_set.h \
  ../../include/geometry/../geometry/../geometry/grid.h \
  ../../include/geometry/../geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../geometry/../geometry/lattice_set.h \
  ../../include/geometry/../geometry/../geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../geometry/../geometry/../base/iterator.h \
  ../../include/geometry/../geometry/../geometry/../geometry/../numeric/interval.h \
  ../../include/geometry/../geometry/../geometry/../geometry/../base/array_operations.h \
  ../../include/geometry/../geometry/../geometry/../geometry/../base/../declarations.h \
  ../../include/geometry/../geometry/../geometry/lattice_set.h \
  ../../include/real_typedef.h ../../include/numeric/numerical_types.h \
  ../../include/numeric/mpfloat.h ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/function.h \
  ../../include/numeric/../numeric/rational.h

