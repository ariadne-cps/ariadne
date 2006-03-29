include ../../config.mk

OBJ = map.o vector_field.o \
	apply.o integrate.o

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
map.o: map.cc ../../include/evaluation/map.h \
  ../../include/evaluation/../geometry/geometry_declarations.h \
  ../../include/evaluation/../linear_algebra/vector.h \
  ../../include/evaluation/../linear_algebra/../base/basic_type.h \
  ../../include/evaluation/../linear_algebra/../base/../base/array.h \
  ../../include/evaluation/../linear_algebra/../base/../base/sequence.h \
  ../../include/evaluation/../linear_algebra/../base/numerical_type.h \
  ../../include/evaluation/../linear_algebra/../base/../base/dyadic.h \
  ../../include/evaluation/../linear_algebra/../base/interval.h \
  ../../include/evaluation/../linear_algebra/matrix.h \
  ../../include/evaluation/../linear_algebra/../base/utility.h \
  ../../include/evaluation/map.tpl \
  ../../include/evaluation/../geometry/point.h \
  ../../include/evaluation/../geometry/../base/utility.h \
  ../../include/evaluation/../geometry/../linear_algebra/vector.h \
  ../../include/evaluation/../geometry/../geometry/geometry_declarations.h \
  ../../include/evaluation/../geometry/rectangle.h \
  ../../include/evaluation/../geometry/../base/interval.h \
  ../../include/evaluation/../geometry/../base/binary_word.h \
  ../../include/evaluation/../geometry/../base/../base/basic_type.h \
  ../../include/evaluation/../geometry/../geometry/point.h \
  ../../include/evaluation/../geometry/list_set.h \
  ../../include/evaluation/../geometry/parallelotope.h \
  ../../include/evaluation/../geometry/../linear_algebra/matrix.h \
  ../../include/evaluation/../geometry/../geometry/rectangle.h \
  ../../include/evaluation/../geometry/../geometry/list_set.h \
  ../../include/evaluation/../geometry/../geometry/polyhedron.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/constraint.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/vector.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/generator.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/point.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/rectangle.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/geometry_declarations.h \
  ../../include/evaluation/../geometry/../geometry/zonotope.h \
  ../../include/evaluation/../geometry/../geometry/../base/utility.h \
  ../../include/evaluation/../geometry/../geometry/../base/interval.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/list_set.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/parallelotope.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/polyhedron.h \
  ../../include/evaluation/../geometry/polyhedron.h \
  ../../include/evaluation/../evaluation/map.h

vector_field.o: vector_field.cc ../../include/evaluation/vector_field.h \
  ../../include/evaluation/../base/interval.h \
  ../../include/evaluation/../linear_algebra/vector.h \
  ../../include/evaluation/../linear_algebra/../base/basic_type.h \
  ../../include/evaluation/../linear_algebra/../base/../base/array.h \
  ../../include/evaluation/../linear_algebra/../base/../base/sequence.h \
  ../../include/evaluation/../linear_algebra/../base/numerical_type.h \
  ../../include/evaluation/../linear_algebra/../base/../base/dyadic.h \
  ../../include/evaluation/../linear_algebra/../base/interval.h \
  ../../include/evaluation/../linear_algebra/matrix.h \
  ../../include/evaluation/../linear_algebra/../base/utility.h \
  ../../include/evaluation/../linear_algebra/interval_vector.h \
  ../../include/evaluation/../linear_algebra/interval_matrix.h \
  ../../include/evaluation/../geometry/geometry_declarations.h \
  ../../include/evaluation/../evaluation/evaluation_declarations.h \
  ../../include/evaluation/vector_field.tpl \
  ../../include/evaluation/../geometry/point.h \
  ../../include/evaluation/../geometry/../base/utility.h \
  ../../include/evaluation/../geometry/../linear_algebra/vector.h \
  ../../include/evaluation/../geometry/../geometry/geometry_declarations.h \
  ../../include/evaluation/../geometry/rectangle.h \
  ../../include/evaluation/../geometry/../base/interval.h \
  ../../include/evaluation/../geometry/../base/binary_word.h \
  ../../include/evaluation/../geometry/../base/../base/basic_type.h \
  ../../include/evaluation/../geometry/../geometry/point.h \
  ../../include/evaluation/../geometry/parallelotope.h \
  ../../include/evaluation/../geometry/../linear_algebra/matrix.h \
  ../../include/evaluation/../geometry/../geometry/rectangle.h \
  ../../include/evaluation/../geometry/../geometry/list_set.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/geometry_declarations.h \
  ../../include/evaluation/../geometry/../geometry/polyhedron.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/constraint.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/vector.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/generator.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/point.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/rectangle.h \
  ../../include/evaluation/../geometry/../geometry/zonotope.h \
  ../../include/evaluation/../geometry/../geometry/../base/utility.h \
  ../../include/evaluation/../geometry/../geometry/../base/interval.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/list_set.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/parallelotope.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/polyhedron.h \
  ../../include/evaluation/../evaluation/vector_field.h

apply.o: apply.cc ../../include/evaluation/apply.h \
  ../../include/evaluation/../geometry/geometry_declarations.h \
  ../../include/evaluation/../evaluation/evaluation_declarations.h \
  ../../include/evaluation/../evaluation/map.h \
  ../../include/evaluation/../evaluation/../geometry/geometry_declarations.h \
  ../../include/evaluation/../evaluation/../linear_algebra/vector.h \
  ../../include/evaluation/../evaluation/../linear_algebra/../base/basic_type.h \
  ../../include/evaluation/../evaluation/../linear_algebra/../base/../base/array.h \
  ../../include/evaluation/../evaluation/../linear_algebra/../base/../base/sequence.h \
  ../../include/evaluation/../evaluation/../linear_algebra/../base/numerical_type.h \
  ../../include/evaluation/../evaluation/../linear_algebra/../base/../base/dyadic.h \
  ../../include/evaluation/../evaluation/../linear_algebra/../base/interval.h \
  ../../include/evaluation/../evaluation/../linear_algebra/matrix.h \
  ../../include/evaluation/../evaluation/../linear_algebra/../base/utility.h \
  ../../include/evaluation/apply.tpl \
  ../../include/evaluation/../base/interval.h \
  ../../include/evaluation/../linear_algebra/vector.h \
  ../../include/evaluation/../linear_algebra/matrix.h \
  ../../include/evaluation/../linear_algebra/interval_vector.h \
  ../../include/evaluation/../linear_algebra/../base/basic_type.h \
  ../../include/evaluation/../linear_algebra/../base/numerical_type.h \
  ../../include/evaluation/../linear_algebra/../base/interval.h \
  ../../include/evaluation/../linear_algebra/interval_matrix.h \
  ../../include/evaluation/../geometry/rectangle.h \
  ../../include/evaluation/../geometry/../base/utility.h \
  ../../include/evaluation/../geometry/../base/interval.h \
  ../../include/evaluation/../geometry/../base/binary_word.h \
  ../../include/evaluation/../geometry/../base/../base/basic_type.h \
  ../../include/evaluation/../geometry/../geometry/point.h \
  ../../include/evaluation/../geometry/../geometry/../base/utility.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/vector.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/geometry_declarations.h \
  ../../include/evaluation/../geometry/../geometry/geometry_declarations.h \
  ../../include/evaluation/../geometry/parallelotope.h \
  ../../include/evaluation/../geometry/../linear_algebra/vector.h \
  ../../include/evaluation/../geometry/../linear_algebra/matrix.h \
  ../../include/evaluation/../geometry/../geometry/rectangle.h \
  ../../include/evaluation/../geometry/../geometry/list_set.h \
  ../../include/evaluation/../geometry/../geometry/polyhedron.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/constraint.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/generator.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/point.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/rectangle.h \
  ../../include/evaluation/../geometry/../geometry/zonotope.h \
  ../../include/evaluation/../geometry/../geometry/../base/interval.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/list_set.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/parallelotope.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/polyhedron.h \
  ../../include/evaluation/../geometry/list_set.h \
  ../../include/evaluation/../geometry/grid_set.h \
  ../../include/evaluation/../geometry/../base/array.h \
  ../../include/evaluation/../geometry/../geometry/grid.h \
  ../../include/evaluation/../geometry/../geometry/../base/array.h \
  ../../include/evaluation/../geometry/../geometry/../base/binary_word.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/lattice_set.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/../base/array.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/../base/iterator.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/../base/../base/array.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/../base/basic_type.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/../base/interval.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/../base/array_operations.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/../base/../base/basic_type.h \
  ../../include/evaluation/../geometry/../geometry/lattice_set.h \
  ../../include/evaluation/../evaluation/apply.h

integrate.o: integrate.cc ../../include/evaluation/integrate.h \
  ../../include/evaluation/../geometry/geometry_declarations.h \
  ../../include/evaluation/../linear_algebra/linear_algebra_declarations.h \
  ../../include/evaluation/../evaluation/evaluation_declarations.h \
  ../../include/evaluation/../evaluation/vector_field.h \
  ../../include/evaluation/../evaluation/../base/interval.h \
  ../../include/evaluation/../evaluation/../linear_algebra/vector.h \
  ../../include/evaluation/../evaluation/../linear_algebra/../base/basic_type.h \
  ../../include/evaluation/../evaluation/../linear_algebra/../base/../base/array.h \
  ../../include/evaluation/../evaluation/../linear_algebra/../base/../base/sequence.h \
  ../../include/evaluation/../evaluation/../linear_algebra/../base/numerical_type.h \
  ../../include/evaluation/../evaluation/../linear_algebra/../base/../base/dyadic.h \
  ../../include/evaluation/../evaluation/../linear_algebra/../base/interval.h \
  ../../include/evaluation/../evaluation/../linear_algebra/matrix.h \
  ../../include/evaluation/../evaluation/../linear_algebra/../base/utility.h \
  ../../include/evaluation/../evaluation/../linear_algebra/interval_vector.h \
  ../../include/evaluation/../evaluation/../linear_algebra/interval_matrix.h \
  ../../include/evaluation/../evaluation/../geometry/geometry_declarations.h \
  ../../include/evaluation/../evaluation/../evaluation/evaluation_declarations.h \
  ../../include/evaluation/integrate.tpl \
  ../../include/evaluation/../base/utility.h \
  ../../include/evaluation/../base/array.h \
  ../../include/evaluation/../base/arithmetic.h \
  ../../include/evaluation/../base/../base/numerical_type.h \
  ../../include/evaluation/../base/interval.h \
  ../../include/evaluation/../linear_algebra/vector.h \
  ../../include/evaluation/../linear_algebra/matrix.h \
  ../../include/evaluation/../linear_algebra/interval_vector.h \
  ../../include/evaluation/../linear_algebra/interval_matrix.h \
  ../../include/evaluation/../geometry/rectangle.h \
  ../../include/evaluation/../geometry/../base/utility.h \
  ../../include/evaluation/../geometry/../base/interval.h \
  ../../include/evaluation/../geometry/../base/binary_word.h \
  ../../include/evaluation/../geometry/../base/../base/basic_type.h \
  ../../include/evaluation/../geometry/../geometry/point.h \
  ../../include/evaluation/../geometry/../geometry/../base/utility.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/vector.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/geometry_declarations.h \
  ../../include/evaluation/../geometry/../geometry/geometry_declarations.h \
  ../../include/evaluation/../geometry/parallelotope.h \
  ../../include/evaluation/../geometry/../linear_algebra/vector.h \
  ../../include/evaluation/../geometry/../linear_algebra/matrix.h \
  ../../include/evaluation/../geometry/../geometry/rectangle.h \
  ../../include/evaluation/../geometry/../geometry/list_set.h \
  ../../include/evaluation/../geometry/../geometry/polyhedron.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/constraint.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/evaluation/../geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/generator.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/point.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/rectangle.h \
  ../../include/evaluation/../geometry/../geometry/zonotope.h \
  ../../include/evaluation/../geometry/../geometry/../base/interval.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/list_set.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/parallelotope.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/polyhedron.h \
  ../../include/evaluation/../geometry/list_set.h \
  ../../include/evaluation/../geometry/grid_set.h \
  ../../include/evaluation/../geometry/../base/array.h \
  ../../include/evaluation/../geometry/../geometry/grid.h \
  ../../include/evaluation/../geometry/../geometry/../base/array.h \
  ../../include/evaluation/../geometry/../geometry/../base/binary_word.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/lattice_set.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/../base/array.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/../base/iterator.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/../base/../base/array.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/../base/basic_type.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/../base/interval.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/../base/array_operations.h \
  ../../include/evaluation/../geometry/../geometry/../geometry/../base/../base/basic_type.h \
  ../../include/evaluation/../geometry/../geometry/lattice_set.h

