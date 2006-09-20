include ../../config.mk

OBJ = newton.o apply.o integrator.o lohner_integrator.o

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
newton.o: newton.cc ../../include/evaluation/newton.h \
  ../../include/evaluation/../declarations.h \
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
  ../../include/linear_algebra/matrix.h \
  ../../include/geometry/rectangle.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../numeric/arithmetic.h \
  ../../include/geometry/../numeric/function.h \
  ../../include/geometry/../numeric/interval.h \
  ../../include/geometry/../linear_algebra/vector.h \
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
  ../../include/geometry/../geometry/rectangle_expression.h \
  ../../include/system/map.h ../../include/system/../declarations.h \
  ../../include/system/vector_field.h ../../include/evaluation/newton.tpl \
  ../../include/evaluation/newton.h \
  ../../include/evaluation/../geometry/point.h \
  ../../include/evaluation/../geometry/rectangle.h \
  ../../include/evaluation/../system/vector_field.h \
  ../../include/real_typedef.h ../../include/numeric/numerical_types.h \
  ../../include/numeric/mpfloat.h /home/collins/include/mpfr.h \
  ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/function.h \
  ../../include/numeric/../numeric/rational.h

apply.o: apply.cc ../../include/evaluation/apply.h \
  ../../include/evaluation/../declarations.h \
  ../../include/evaluation/apply.tpl ../../include/evaluation/apply.h \
  ../../include/evaluation/../numeric/interval.h \
  ../../include/evaluation/../numeric/../declarations.h \
  ../../include/evaluation/../numeric/../numeric/numerical_traits.h \
  ../../include/evaluation/../numeric/../numeric/approximation.h \
  ../../include/evaluation/../numeric/../numeric/arithmetic.h \
  ../../include/evaluation/../numeric/../numeric/rounding.h \
  ../../include/evaluation/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/evaluation/../numeric/../numeric/../numeric/function.h \
  ../../include/evaluation/../linear_algebra/vector.h \
  ../../include/evaluation/../linear_algebra/../declarations.h \
  ../../include/evaluation/../linear_algebra/../numeric/integer.h \
  ../../include/evaluation/../linear_algebra/../numeric/../declarations.h \
  ../../include/evaluation/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/evaluation/../linear_algebra/../numeric/../numeric/arithmetic.h \
  ../../include/evaluation/../linear_algebra/../numeric/arithmetic.h \
  ../../include/evaluation/../linear_algebra/../numeric/interval.h \
  ../../include/evaluation/../linear_algebra/matrix.h \
  ../../include/evaluation/../combinatoric/lattice_set.h \
  ../../include/evaluation/../combinatoric/../base/array.h \
  ../../include/evaluation/../combinatoric/../base/../declarations.h \
  ../../include/evaluation/../combinatoric/../base/iterator.h \
  ../../include/evaluation/../combinatoric/../base/../base/array.h \
  ../../include/evaluation/../combinatoric/../base/array_operations.h \
  ../../include/evaluation/../combinatoric/../numeric/interval.h \
  ../../include/evaluation/../geometry/rectangle.h \
  ../../include/evaluation/../geometry/../declarations.h \
  ../../include/evaluation/../geometry/../numeric/arithmetic.h \
  ../../include/evaluation/../geometry/../numeric/function.h \
  ../../include/evaluation/../geometry/../numeric/interval.h \
  ../../include/evaluation/../geometry/../linear_algebra/vector.h \
  ../../include/evaluation/../geometry/../geometry/ppl_polyhedron.h \
  ../../include/evaluation/../geometry/../geometry/../declarations.h \
  ../../include/evaluation/../geometry/../geometry/../numeric/rational.h \
  ../../include/evaluation/../geometry/../geometry/../numeric/../declarations.h \
  ../../include/evaluation/../geometry/../geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/evaluation/../geometry/../geometry/../numeric/../numeric/arithmetic.h \
  ../../include/evaluation/../geometry/../geometry/../numeric/../numeric/function.h \
  ../../include/evaluation/../geometry/../geometry/../numeric/../numeric/integer.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/vector.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/evaluation/../geometry/../geometry/point.h \
  ../../include/evaluation/../geometry/../geometry/../base/array.h \
  ../../include/evaluation/../geometry/../geometry/rectangle_expression.h \
  ../../include/evaluation/../geometry/parallelotope.h \
  ../../include/evaluation/../geometry/../linear_algebra/matrix.h \
  ../../include/evaluation/../geometry/../geometry/zonotope.h \
  ../../include/evaluation/../geometry/../geometry/../numeric/interval.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/ppl_polyhedron.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/point.h \
  ../../include/evaluation/../geometry/list_set.h \
  ../../include/evaluation/../geometry/../base/utility.h \
  ../../include/evaluation/../geometry/../base/../declarations.h \
  ../../include/evaluation/../geometry/grid_set.h \
  ../../include/evaluation/../geometry/../base/array.h \
  ../../include/evaluation/../geometry/../combinatoric/lattice_set.h \
  ../../include/evaluation/../geometry/../geometry/rectangle.h \
  ../../include/evaluation/../geometry/../geometry/grid.h \
  ../../include/evaluation/../geometry/../geometry/../numeric/arithmetic.h \
  ../../include/evaluation/../geometry/../geometry/../combinatoric/lattice_set.h \
  ../../include/evaluation/../system/map.h \
  ../../include/evaluation/../system/../declarations.h \
  ../../include/real_typedef.h ../../include/numeric/numerical_types.h \
  ../../include/numeric/mpfloat.h /home/collins/include/mpfr.h \
  ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/function.h \
  ../../include/numeric/../numeric/rational.h

integrator.o: integrator.cc ../../include/evaluation/integrator.h \
  ../../include/evaluation/../declarations.h \
  ../../include/evaluation/integrator.tpl \
  ../../include/evaluation/integrator.h \
  ../../include/evaluation/../base/array.h \
  ../../include/evaluation/../base/../declarations.h \
  ../../include/evaluation/../numeric/arithmetic.h \
  ../../include/evaluation/../numeric/interval.h \
  ../../include/evaluation/../numeric/../declarations.h \
  ../../include/evaluation/../numeric/../numeric/numerical_traits.h \
  ../../include/evaluation/../numeric/../numeric/approximation.h \
  ../../include/evaluation/../numeric/../numeric/arithmetic.h \
  ../../include/evaluation/../numeric/../numeric/rounding.h \
  ../../include/evaluation/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/evaluation/../numeric/../numeric/../numeric/function.h \
  ../../include/evaluation/../linear_algebra/vector.h \
  ../../include/evaluation/../linear_algebra/../declarations.h \
  ../../include/evaluation/../linear_algebra/../numeric/integer.h \
  ../../include/evaluation/../linear_algebra/../numeric/../declarations.h \
  ../../include/evaluation/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/evaluation/../linear_algebra/../numeric/../numeric/arithmetic.h \
  ../../include/evaluation/../linear_algebra/../numeric/arithmetic.h \
  ../../include/evaluation/../linear_algebra/../numeric/interval.h \
  ../../include/evaluation/../linear_algebra/matrix.h \
  ../../include/evaluation/../geometry/rectangle.h \
  ../../include/evaluation/../geometry/../declarations.h \
  ../../include/evaluation/../geometry/../numeric/arithmetic.h \
  ../../include/evaluation/../geometry/../numeric/function.h \
  ../../include/evaluation/../geometry/../numeric/interval.h \
  ../../include/evaluation/../geometry/../linear_algebra/vector.h \
  ../../include/evaluation/../geometry/../geometry/ppl_polyhedron.h \
  ../../include/evaluation/../geometry/../geometry/../declarations.h \
  ../../include/evaluation/../geometry/../geometry/../numeric/rational.h \
  ../../include/evaluation/../geometry/../geometry/../numeric/../declarations.h \
  ../../include/evaluation/../geometry/../geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/evaluation/../geometry/../geometry/../numeric/../numeric/arithmetic.h \
  ../../include/evaluation/../geometry/../geometry/../numeric/../numeric/function.h \
  ../../include/evaluation/../geometry/../geometry/../numeric/../numeric/integer.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/vector.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/evaluation/../geometry/../geometry/point.h \
  ../../include/evaluation/../geometry/../geometry/../base/array.h \
  ../../include/evaluation/../geometry/../geometry/rectangle_expression.h \
  ../../include/evaluation/../geometry/parallelotope.h \
  ../../include/evaluation/../geometry/../linear_algebra/matrix.h \
  ../../include/evaluation/../geometry/../geometry/zonotope.h \
  ../../include/evaluation/../geometry/../geometry/../numeric/interval.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/ppl_polyhedron.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/point.h \
  ../../include/evaluation/../geometry/zonotope.h \
  ../../include/evaluation/../geometry/list_set.h \
  ../../include/evaluation/../geometry/../base/utility.h \
  ../../include/evaluation/../geometry/../base/../declarations.h \
  ../../include/evaluation/../geometry/grid.h \
  ../../include/evaluation/../geometry/../base/array.h \
  ../../include/evaluation/../geometry/../combinatoric/lattice_set.h \
  ../../include/evaluation/../geometry/../combinatoric/../base/array.h \
  ../../include/evaluation/../geometry/../combinatoric/../base/iterator.h \
  ../../include/evaluation/../geometry/../combinatoric/../base/../base/array.h \
  ../../include/evaluation/../geometry/../combinatoric/../base/array_operations.h \
  ../../include/evaluation/../geometry/../combinatoric/../base/../declarations.h \
  ../../include/evaluation/../geometry/../combinatoric/../numeric/interval.h \
  ../../include/evaluation/../geometry/grid_set.h \
  ../../include/evaluation/../geometry/../geometry/rectangle.h \
  ../../include/evaluation/../geometry/../geometry/grid.h \
  ../../include/evaluation/../system/vector_field.h \
  ../../include/evaluation/../system/../declarations.h \
  ../../include/evaluation/../system/affine_vector_field.h \
  ../../include/evaluation/../system/../linear_algebra/vector.h \
  ../../include/evaluation/../system/../linear_algebra/matrix.h \
  ../../include/evaluation/../system/../system/vector_field.h \
  ../../include/real_typedef.h ../../include/numeric/numerical_types.h \
  ../../include/numeric/mpfloat.h /home/collins/include/mpfr.h \
  ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/function.h \
  ../../include/numeric/../numeric/rational.h

lohner_integrator.o: lohner_integrator.cc \
  ../../include/evaluation/integrator.h \
  ../../include/evaluation/../declarations.h \
  ../../include/evaluation/lohner_integrator.h \
  ../../include/evaluation/../evaluation/integrator.h \
  ../../include/evaluation/lohner_integrator.tpl \
  ../../include/evaluation/lohner_integrator.h \
  ../../include/evaluation/../base/array.h \
  ../../include/evaluation/../base/../declarations.h \
  ../../include/evaluation/../numeric/arithmetic.h \
  ../../include/evaluation/../numeric/interval.h \
  ../../include/evaluation/../numeric/../declarations.h \
  ../../include/evaluation/../numeric/../numeric/numerical_traits.h \
  ../../include/evaluation/../numeric/../numeric/approximation.h \
  ../../include/evaluation/../numeric/../numeric/arithmetic.h \
  ../../include/evaluation/../numeric/../numeric/rounding.h \
  ../../include/evaluation/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/evaluation/../numeric/../numeric/../numeric/function.h \
  ../../include/evaluation/../linear_algebra/vector.h \
  ../../include/evaluation/../linear_algebra/../declarations.h \
  ../../include/evaluation/../linear_algebra/../numeric/integer.h \
  ../../include/evaluation/../linear_algebra/../numeric/../declarations.h \
  ../../include/evaluation/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/evaluation/../linear_algebra/../numeric/../numeric/arithmetic.h \
  ../../include/evaluation/../linear_algebra/../numeric/arithmetic.h \
  ../../include/evaluation/../linear_algebra/../numeric/interval.h \
  ../../include/evaluation/../linear_algebra/matrix.h \
  ../../include/evaluation/../linear_algebra/matrix_function.h \
  ../../include/evaluation/../geometry/rectangle.h \
  ../../include/evaluation/../geometry/../declarations.h \
  ../../include/evaluation/../geometry/../numeric/arithmetic.h \
  ../../include/evaluation/../geometry/../numeric/function.h \
  ../../include/evaluation/../geometry/../numeric/interval.h \
  ../../include/evaluation/../geometry/../linear_algebra/vector.h \
  ../../include/evaluation/../geometry/../geometry/ppl_polyhedron.h \
  ../../include/evaluation/../geometry/../geometry/../declarations.h \
  ../../include/evaluation/../geometry/../geometry/../numeric/rational.h \
  ../../include/evaluation/../geometry/../geometry/../numeric/../declarations.h \
  ../../include/evaluation/../geometry/../geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/evaluation/../geometry/../geometry/../numeric/../numeric/arithmetic.h \
  ../../include/evaluation/../geometry/../geometry/../numeric/../numeric/function.h \
  ../../include/evaluation/../geometry/../geometry/../numeric/../numeric/integer.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/vector.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/evaluation/../geometry/../geometry/point.h \
  ../../include/evaluation/../geometry/../geometry/../base/array.h \
  ../../include/evaluation/../geometry/../geometry/rectangle_expression.h \
  ../../include/evaluation/../geometry/parallelotope.h \
  ../../include/evaluation/../geometry/../linear_algebra/matrix.h \
  ../../include/evaluation/../geometry/../geometry/zonotope.h \
  ../../include/evaluation/../geometry/../geometry/../numeric/interval.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/ppl_polyhedron.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/point.h \
  ../../include/evaluation/../geometry/zonotope.h \
  ../../include/evaluation/../geometry/list_set.h \
  ../../include/evaluation/../geometry/../base/utility.h \
  ../../include/evaluation/../geometry/../base/../declarations.h \
  ../../include/evaluation/../geometry/grid.h \
  ../../include/evaluation/../geometry/../base/array.h \
  ../../include/evaluation/../geometry/../combinatoric/lattice_set.h \
  ../../include/evaluation/../geometry/../combinatoric/../base/array.h \
  ../../include/evaluation/../geometry/../combinatoric/../base/iterator.h \
  ../../include/evaluation/../geometry/../combinatoric/../base/../base/array.h \
  ../../include/evaluation/../geometry/../combinatoric/../base/array_operations.h \
  ../../include/evaluation/../geometry/../combinatoric/../base/../declarations.h \
  ../../include/evaluation/../geometry/../combinatoric/../numeric/interval.h \
  ../../include/evaluation/../geometry/grid_set.h \
  ../../include/evaluation/../geometry/../geometry/rectangle.h \
  ../../include/evaluation/../geometry/../geometry/grid.h \
  ../../include/evaluation/../system/vector_field.h \
  ../../include/evaluation/../system/../declarations.h \
  ../../include/evaluation/../system/affine_vector_field.h \
  ../../include/evaluation/../system/../linear_algebra/vector.h \
  ../../include/evaluation/../system/../linear_algebra/matrix.h \
  ../../include/evaluation/../system/../system/vector_field.h \
  ../../include/real_typedef.h ../../include/numeric/numerical_types.h \
  ../../include/numeric/mpfloat.h /home/collins/include/mpfr.h \
  ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/function.h \
  ../../include/numeric/../numeric/rational.h

