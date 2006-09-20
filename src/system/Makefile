include ../../config.mk

OBJ = map.o multimap.o vector_field.o \
	affine_map.o affine_multimap.o affine_vector_field.o \
	polynomial_map.o \

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
map.o: map.cc ../../include/system/map.h \
  ../../include/system/../declarations.h ../../include/system/map.tpl \
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
  ../../include/system/../geometry/polytope.h \
  ../../include/system/../geometry/../geometry/point_list.h \
  ../../include/system/../geometry/../geometry/../geometry/point.h \
  ../../include/system/../geometry/../geometry/rectangle.h \
  ../../include/system/../geometry/simplex.h \
  ../../include/system/../geometry/../geometry/polytope.h \
  ../../include/system/../geometry/zonotope.h \
  ../../include/system/../geometry/parallelotope.h \
  ../../include/system/../geometry/../geometry/zonotope.h \
  ../../include/system/../system/map.h ../../include/real_typedef.h \
  ../../include/numeric/numerical_types.h ../../include/numeric/mpfloat.h \
  /home/collins/include/mpfr.h ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/function.h \
  ../../include/numeric/../numeric/rational.h

multimap.o: multimap.cc ../../include/system/multimap.h \
  ../../include/system/../declarations.h \
  ../../include/system/multimap.tpl \
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
  ../../include/system/../geometry/zonotope.h \
  ../../include/system/../geometry/simplex.h \
  ../../include/system/../geometry/../geometry/point_list.h \
  ../../include/system/../geometry/../geometry/polytope.h \
  ../../include/system/../geometry/../geometry/../geometry/point_list.h \
  ../../include/system/../geometry/../geometry/../geometry/rectangle.h \
  ../../include/system/../geometry/polytope.h \
  ../../include/system/../geometry/list_set.h \
  ../../include/system/../geometry/../base/utility.h \
  ../../include/system/../system/multimap.h ../../include/real_typedef.h \
  ../../include/numeric/numerical_types.h ../../include/numeric/mpfloat.h \
  /home/collins/include/mpfr.h ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/function.h \
  ../../include/numeric/../numeric/rational.h

vector_field.o: vector_field.cc ../../include/system/vector_field.h \
  ../../include/system/../declarations.h \
  ../../include/system/vector_field.tpl \
  ../../include/system/vector_field.h \
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
  ../../include/real_typedef.h ../../include/numeric/numerical_types.h \
  ../../include/numeric/mpfloat.h /home/collins/include/mpfr.h \
  ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/function.h \
  ../../include/numeric/../numeric/rational.h

affine_map.o: affine_map.cc ../../include/system/affine_map.h \
  ../../include/system/../declarations.h \
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
  ../../include/system/../system/map.h \
  ../../include/system/../system/../declarations.h \
  ../../include/system/affine_map.tpl ../../include/system/affine_map.h \
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
  ../../include/system/../geometry/zonotope.h \
  ../../include/real_typedef.h ../../include/numeric/numerical_types.h \
  ../../include/numeric/mpfloat.h /home/collins/include/mpfr.h \
  ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/function.h \
  ../../include/numeric/../numeric/rational.h

affine_multimap.o: affine_multimap.cc \
  ../../include/system/affine_multimap.h ../../include/system/multimap.h \
  ../../include/system/../declarations.h \
  ../../include/system/affine_multimap.tpl \
  ../../include/system/affine_multimap.h \
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
  ../../include/system/../geometry/list_set.h \
  ../../include/system/../geometry/../base/utility.h \
  ../../include/system/../geometry/parallelotope.h \
  ../../include/system/../geometry/../geometry/zonotope.h \
  ../../include/system/../geometry/../geometry/../numeric/interval.h \
  ../../include/system/../geometry/../geometry/../geometry/ppl_polyhedron.h \
  ../../include/system/../geometry/../geometry/../geometry/point.h \
  ../../include/system/../geometry/zonotope.h \
  ../../include/system/../geometry/simplex.h \
  ../../include/system/../geometry/../geometry/point_list.h \
  ../../include/system/../geometry/../geometry/polytope.h \
  ../../include/system/../geometry/../geometry/../geometry/point_list.h \
  ../../include/system/../geometry/../geometry/../geometry/rectangle.h \
  ../../include/system/../geometry/polyhedron.h \
  ../../include/system/../geometry/../geometry/rectangle.h \
  ../../include/system/../system/affine_map.h \
  ../../include/system/../system/../declarations.h \
  ../../include/system/../system/../linear_algebra/vector.h \
  ../../include/system/../system/../linear_algebra/matrix.h \
  ../../include/system/../system/../system/map.h \
  ../../include/system/../system/../system/../declarations.h \
  ../../include/real_typedef.h ../../include/numeric/numerical_types.h \
  ../../include/numeric/mpfloat.h /home/collins/include/mpfr.h \
  ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/function.h \
  ../../include/numeric/../numeric/rational.h

affine_vector_field.o: affine_vector_field.cc \
  ../../include/system/affine_vector_field.h \
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
  ../../include/system/../system/vector_field.h \
  ../../include/system/../system/../declarations.h \
  ../../include/system/affine_vector_field.tpl \
  ../../include/system/affine_vector_field.h \
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
  ../../include/system/../geometry/zonotope.h \
  ../../include/system/../geometry/parallelotope.h \
  ../../include/system/../geometry/../geometry/zonotope.h \
  ../../include/system/../geometry/simplex.h \
  ../../include/system/../geometry/../geometry/point_list.h \
  ../../include/system/../geometry/../geometry/../geometry/point.h \
  ../../include/system/../geometry/../geometry/polytope.h \
  ../../include/system/../geometry/../geometry/../geometry/ppl_polyhedron.h \
  ../../include/system/../geometry/../geometry/../geometry/point_list.h \
  ../../include/system/../geometry/../geometry/../geometry/rectangle.h \
  ../../include/system/../geometry/polyhedron.h \
  ../../include/system/../geometry/../geometry/rectangle.h \
  ../../include/real_typedef.h ../../include/numeric/numerical_types.h \
  ../../include/numeric/mpfloat.h /home/collins/include/mpfr.h \
  ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/function.h \
  ../../include/numeric/../numeric/rational.h

polynomial_map.o: polynomial_map.cc ../../include/system/polynomial_map.h \
  ../../include/system/../base/array.h \
  ../../include/system/../base/../declarations.h \
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
  ../../include/system/../system/map.h \
  ../../include/system/../system/../declarations.h \
  ../../include/system/polynomial_map.tpl \
  ../../include/system/polynomial_map.h \
  ../../include/system/../base/stlio.h \
  ../../include/system/../base/../base/array.h \
  ../../include/system/../numeric/interval.h \
  ../../include/system/../linear_algebra/vector.h \
  ../../include/system/../linear_algebra/../numeric/arithmetic.h \
  ../../include/system/../geometry/point.h \
  ../../include/system/../geometry/../declarations.h \
  ../../include/system/../geometry/../base/array.h \
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
  ../../include/real_typedef.h ../../include/numeric/numerical_types.h \
  ../../include/numeric/mpfloat.h /home/collins/include/mpfr.h \
  ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/function.h \
  ../../include/numeric/../numeric/rational.h

