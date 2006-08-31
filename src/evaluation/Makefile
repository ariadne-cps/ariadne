include ../../config.mk

OBJ = newton.o apply.o integrator.o

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
  ../../include/linear_algebra/interval_vector.h \
  ../../include/linear_algebra/../declarations.h \
  ../../include/linear_algebra/../numeric/interval.h \
  ../../include/linear_algebra/../numeric/../declarations.h \
  ../../include/linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/linear_algebra/../numeric/../numeric/arithmetic.h \
  ../../include/linear_algebra/../numeric/../numeric/../numeric/numerical_types.h \
  ../../include/linear_algebra/../numeric/../numeric/../numeric/../declarations.h \
  ../../include/linear_algebra/../numeric/../numeric/../numeric/../numeric/integer.h \
  ../../include/linear_algebra/../numeric/../numeric/../numeric/../numeric/../declarations.h \
  ../../include/linear_algebra/../numeric/../numeric/../numeric/../numeric/../utility/stlio.h \
  ../../include/linear_algebra/../numeric/../numeric/../numeric/../numeric/../utility/../base/array.h \
  ../../include/linear_algebra/../numeric/../numeric/../numeric/../numeric/../utility/../base/../declarations.h \
  ../../include/linear_algebra/../numeric/../numeric/../numeric/../numeric/dyadic.h \
  ../../include/linear_algebra/../numeric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/linear_algebra/../numeric/../numeric/../numeric/../utility/stlio.h \
  ../../include/linear_algebra/../numeric/../numeric/../numeric/approximation.h \
  ../../include/linear_algebra/../numeric/../numeric/../numeric/../numeric/numerical_types.h \
  ../../include/linear_algebra/../linear_algebra/vector.h \
  ../../include/linear_algebra/../linear_algebra/../declarations.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/integer.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/arithmetic.h \
  ../../include/linear_algebra/interval_matrix.h \
  ../../include/linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/rectangle.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../numeric/arithmetic.h \
  ../../include/geometry/../numeric/interval.h \
  ../../include/geometry/../linear_algebra/interval_vector.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../linear_algebra/vector.h \
  ../../include/system/map.h ../../include/system/../declarations.h \
  ../../include/system/vector_field.h ../../include/evaluation/newton.tpl \
  ../../include/evaluation/newton.h \
  ../../include/evaluation/../linear_algebra/matrix.h \
  ../../include/evaluation/../linear_algebra/interval_vector.h \
  ../../include/evaluation/../linear_algebra/interval_matrix.h \
  ../../include/evaluation/../geometry/point.h \
  ../../include/evaluation/../geometry/rectangle.h \
  ../../include/evaluation/../system/vector_field.h \
  ../../include/real_typedef.h ../../include/numeric/numerical_types.h

apply.o: apply.cc ../../include/evaluation/apply.h \
  ../../include/evaluation/../declarations.h \
  ../../include/evaluation/apply.tpl ../../include/evaluation/apply.h \
  ../../include/evaluation/../numeric/interval.h \
  ../../include/evaluation/../numeric/../declarations.h \
  ../../include/evaluation/../numeric/../numeric/numerical_traits.h \
  ../../include/evaluation/../numeric/../numeric/arithmetic.h \
  ../../include/evaluation/../numeric/../numeric/../numeric/numerical_types.h \
  ../../include/evaluation/../numeric/../numeric/../numeric/../declarations.h \
  ../../include/evaluation/../numeric/../numeric/../numeric/../numeric/integer.h \
  ../../include/evaluation/../numeric/../numeric/../numeric/../numeric/../declarations.h \
  ../../include/evaluation/../numeric/../numeric/../numeric/../numeric/../utility/stlio.h \
  ../../include/evaluation/../numeric/../numeric/../numeric/../numeric/../utility/../base/array.h \
  ../../include/evaluation/../numeric/../numeric/../numeric/../numeric/../utility/../base/../declarations.h \
  ../../include/evaluation/../numeric/../numeric/../numeric/../numeric/dyadic.h \
  ../../include/evaluation/../numeric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/evaluation/../numeric/../numeric/../numeric/../utility/stlio.h \
  ../../include/evaluation/../numeric/../numeric/../numeric/approximation.h \
  ../../include/evaluation/../numeric/../numeric/../numeric/../numeric/numerical_types.h \
  ../../include/evaluation/../linear_algebra/vector.h \
  ../../include/evaluation/../linear_algebra/../declarations.h \
  ../../include/evaluation/../linear_algebra/../numeric/integer.h \
  ../../include/evaluation/../linear_algebra/../numeric/arithmetic.h \
  ../../include/evaluation/../linear_algebra/matrix.h \
  ../../include/evaluation/../linear_algebra/interval_vector.h \
  ../../include/evaluation/../linear_algebra/../numeric/interval.h \
  ../../include/evaluation/../linear_algebra/../linear_algebra/vector.h \
  ../../include/evaluation/../linear_algebra/interval_matrix.h \
  ../../include/evaluation/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/evaluation/../geometry/rectangle.h \
  ../../include/evaluation/../geometry/../declarations.h \
  ../../include/evaluation/../geometry/../numeric/arithmetic.h \
  ../../include/evaluation/../geometry/../numeric/interval.h \
  ../../include/evaluation/../geometry/../linear_algebra/interval_vector.h \
  ../../include/evaluation/../geometry/../geometry/point.h \
  ../../include/evaluation/../geometry/../geometry/../declarations.h \
  ../../include/evaluation/../geometry/../geometry/../base/array.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/vector.h \
  ../../include/evaluation/../geometry/parallelotope.h \
  ../../include/evaluation/../geometry/../linear_algebra/vector.h \
  ../../include/evaluation/../geometry/../linear_algebra/matrix.h \
  ../../include/evaluation/../geometry/../geometry/polyhedron.h \
  /ufs/collins/include/ppl.hh \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/constraint_system.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/../declarations.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/../numeric/numerical_types.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/interval_vector.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/generator_system.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/point.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/rectangle.h \
  ../../include/evaluation/../geometry/list_set.h \
  ../../include/evaluation/../geometry/../base/utility.h \
  ../../include/evaluation/../geometry/../base/../declarations.h \
  ../../include/evaluation/../geometry/grid_set.h \
  ../../include/evaluation/../geometry/../base/array.h \
  ../../include/evaluation/../geometry/../geometry/rectangle.h \
  ../../include/evaluation/../geometry/../geometry/grid.h \
  ../../include/evaluation/../geometry/../geometry/../numeric/arithmetic.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/lattice_set.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/../base/array.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/../base/iterator.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/../base/../base/array.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/../numeric/interval.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/../base/array_operations.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/../base/../declarations.h \
  ../../include/evaluation/../geometry/../geometry/lattice_set.h \
  ../../include/evaluation/../geometry/lattice_set.h \
  ../../include/evaluation/../system/map.h \
  ../../include/evaluation/../system/../declarations.h \
  ../../include/real_typedef.h ../../include/numeric/numerical_types.h

integrator.o: integrator.cc ../../include/evaluation/integrator.h \
  ../../include/evaluation/../declarations.h \
  ../../include/evaluation/lohner_integrator.h \
  ../../include/evaluation/../evaluation/integrator.h \
  ../../include/evaluation/integrator.tpl \
  ../../include/evaluation/integrator.h \
  ../../include/evaluation/../utility/stlio.h \
  ../../include/evaluation/../utility/../base/array.h \
  ../../include/evaluation/../utility/../base/../declarations.h \
  ../../include/evaluation/../base/array.h \
  ../../include/evaluation/../numeric/arithmetic.h \
  ../../include/evaluation/../numeric/../numeric/numerical_types.h \
  ../../include/evaluation/../numeric/../numeric/../declarations.h \
  ../../include/evaluation/../numeric/../numeric/../numeric/integer.h \
  ../../include/evaluation/../numeric/../numeric/../numeric/../declarations.h \
  ../../include/evaluation/../numeric/../numeric/../numeric/../utility/stlio.h \
  ../../include/evaluation/../numeric/../numeric/../numeric/dyadic.h \
  ../../include/evaluation/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/evaluation/../numeric/../numeric/../utility/stlio.h \
  ../../include/evaluation/../numeric/../numeric/approximation.h \
  ../../include/evaluation/../numeric/../numeric/../numeric/numerical_types.h \
  ../../include/evaluation/../numeric/function.h \
  /ufs/collins/include/mpfr.h \
  ../../include/evaluation/../numeric/interval.h \
  ../../include/evaluation/../numeric/../declarations.h \
  ../../include/evaluation/../numeric/../numeric/numerical_traits.h \
  ../../include/evaluation/../numeric/../numeric/arithmetic.h \
  ../../include/evaluation/../linear_algebra/vector.h \
  ../../include/evaluation/../linear_algebra/../declarations.h \
  ../../include/evaluation/../linear_algebra/../numeric/integer.h \
  ../../include/evaluation/../linear_algebra/../numeric/arithmetic.h \
  ../../include/evaluation/../linear_algebra/interval_vector.h \
  ../../include/evaluation/../linear_algebra/../numeric/interval.h \
  ../../include/evaluation/../linear_algebra/../linear_algebra/vector.h \
  ../../include/evaluation/../linear_algebra/transformation_system.h \
  ../../include/evaluation/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/evaluation/../linear_algebra/../linear_algebra/../declarations.h \
  ../../include/evaluation/../linear_algebra/../linear_algebra/../numeric/integer.h \
  ../../include/evaluation/../linear_algebra/../linear_algebra/interval_vector.h \
  ../../include/evaluation/../linear_algebra/../linear_algebra/interval_matrix.h \
  ../../include/evaluation/../linear_algebra/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/evaluation/../linear_algebra/../linear_algebra/../numeric/interval.h \
  ../../include/evaluation/../linear_algebra/matrix.h \
  ../../include/evaluation/../linear_algebra/interval_matrix.h \
  ../../include/evaluation/../geometry/rectangle.h \
  ../../include/evaluation/../geometry/../declarations.h \
  ../../include/evaluation/../geometry/../numeric/arithmetic.h \
  ../../include/evaluation/../geometry/../numeric/interval.h \
  ../../include/evaluation/../geometry/../linear_algebra/interval_vector.h \
  ../../include/evaluation/../geometry/../geometry/point.h \
  ../../include/evaluation/../geometry/../geometry/../declarations.h \
  ../../include/evaluation/../geometry/../geometry/../base/array.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/vector.h \
  ../../include/evaluation/../geometry/parallelotope.h \
  ../../include/evaluation/../geometry/../linear_algebra/vector.h \
  ../../include/evaluation/../geometry/../linear_algebra/matrix.h \
  ../../include/evaluation/../geometry/../geometry/polyhedron.h \
  /ufs/collins/include/ppl.hh \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/constraint_system.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/../declarations.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/../numeric/numerical_types.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/interval_vector.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/generator_system.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/point.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/rectangle.h \
  ../../include/evaluation/../geometry/zonotope.h \
  ../../include/evaluation/../geometry/../linear_algebra/transformation_system.h \
  ../../include/evaluation/../geometry/list_set.h \
  ../../include/evaluation/../geometry/../base/utility.h \
  ../../include/evaluation/../geometry/../base/../declarations.h \
  ../../include/evaluation/../geometry/grid.h \
  ../../include/evaluation/../geometry/../base/array.h \
  ../../include/evaluation/../geometry/../geometry/lattice_set.h \
  ../../include/evaluation/../geometry/../geometry/../base/iterator.h \
  ../../include/evaluation/../geometry/../geometry/../base/../base/array.h \
  ../../include/evaluation/../geometry/../geometry/../numeric/interval.h \
  ../../include/evaluation/../geometry/../geometry/../base/array_operations.h \
  ../../include/evaluation/../geometry/../geometry/../base/../declarations.h \
  ../../include/evaluation/../geometry/grid_set.h \
  ../../include/evaluation/../geometry/../geometry/rectangle.h \
  ../../include/evaluation/../geometry/../geometry/grid.h \
  ../../include/evaluation/../system/vector_field.h \
  ../../include/evaluation/../system/../declarations.h \
  ../../include/evaluation/../system/affine_vector_field.h \
  ../../include/evaluation/../system/../linear_algebra/vector.h \
  ../../include/evaluation/../system/../linear_algebra/matrix.h \
  ../../include/evaluation/../system/../linear_algebra/interval_vector.h \
  ../../include/evaluation/../system/../linear_algebra/interval_matrix.h \
  ../../include/evaluation/../system/../system/vector_field.h \
  ../../include/evaluation/../system/../system/affine_map.h \
  ../../include/evaluation/../system/../system/../numeric/numerical_types.h \
  ../../include/evaluation/../system/../system/../linear_algebra/vector.h \
  ../../include/evaluation/../system/../system/../linear_algebra/matrix.h \
  ../../include/evaluation/../system/../system/../linear_algebra/interval_matrix.h \
  ../../include/evaluation/../system/../system/../geometry/point.h \
  ../../include/evaluation/../system/../system/../geometry/simplex.h \
  ../../include/evaluation/../system/../system/../geometry/../declarations.h \
  ../../include/evaluation/../system/../system/../geometry/../base/array.h \
  ../../include/evaluation/../system/../system/../geometry/../geometry/point.h \
  ../../include/evaluation/../system/../system/../geometry/rectangle.h \
  ../../include/evaluation/../system/../system/../geometry/parallelotope.h \
  ../../include/evaluation/../system/../system/../geometry/zonotope.h \
  ../../include/evaluation/../system/../system/../geometry/polyhedron.h \
  ../../include/evaluation/../system/../system/../geometry/list_set.h \
  ../../include/evaluation/../system/../system/../geometry/grid_set.h \
  ../../include/evaluation/../system/../system/../system/map.h \
  ../../include/evaluation/../system/../system/../system/../declarations.h \
  ../../include/evaluation/lohner_integrator.tpl \
  ../../include/evaluation/lohner_integrator.h \
  ../../include/real_typedef.h ../../include/numeric/numerical_types.h

