include ../../config.mk

OBJ = parallelotope.o zonotope.o simplex.o polyhedron.o \
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

### Do Not edit this line
### Everything below is auto-generated
parallelotope.o: parallelotope.cc ../../include/geometry/parallelotope.h \
  ../../include/geometry/../base/utility.h \
  ../../include/geometry/../base/../base/array.h \
  ../../include/geometry/../base/interval.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/../base/basic_type.h \
  ../../include/geometry/../linear_algebra/../base/../base/array.h \
  ../../include/geometry/../linear_algebra/../base/../base/sequence.h \
  ../../include/geometry/../linear_algebra/../base/numerical_type.h \
  ../../include/geometry/../linear_algebra/../base/../base/dyadic.h \
  ../../include/geometry/../linear_algebra/../base/interval.h \
  ../../include/geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../linear_algebra/../base/utility.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../base/utility.h \
  ../../include/geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../geometry/geometry_declarations.h \
  ../../include/geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/../base/interval.h \
  ../../include/geometry/../geometry/../base/binary_word.h \
  ../../include/geometry/../geometry/../base/../base/basic_type.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/list_set.h \
  ../../include/geometry/../geometry/polyhedron.h \
  ../../include/geometry/../geometry/../linear_algebra/constraint.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../geometry/generator.h \
  ../../include/geometry/../geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/zonotope.h \
  ../../include/geometry/../geometry/../geometry/list_set.h \
  ../../include/geometry/../geometry/../geometry/parallelotope.h \
  ../../include/geometry/../geometry/../geometry/polyhedron.h \
  ../../include/geometry/../geometry/geometry_declarations.h \
  ../../include/geometry/parallelotope.tpl \
  ../../include/geometry/../linear_algebra/linear_program.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/lattice_set.h \
  ../../include/geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../base/iterator.h \
  ../../include/geometry/../geometry/../base/../base/array.h \
  ../../include/geometry/../geometry/../base/basic_type.h \
  ../../include/geometry/../geometry/../base/array_operations.h \
  ../../include/geometry/../geometry/parallelotope.h

zonotope.o: zonotope.cc ../../include/geometry/zonotope.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/../base/basic_type.h \
  ../../include/geometry/../linear_algebra/../base/../base/array.h \
  ../../include/geometry/../linear_algebra/../base/../base/sequence.h \
  ../../include/geometry/../linear_algebra/../base/numerical_type.h \
  ../../include/geometry/../linear_algebra/../base/../base/dyadic.h \
  ../../include/geometry/../linear_algebra/../base/interval.h \
  ../../include/geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../linear_algebra/../base/utility.h \
  ../../include/geometry/../base/utility.h \
  ../../include/geometry/../base/interval.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../base/utility.h \
  ../../include/geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../geometry/geometry_declarations.h \
  ../../include/geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/../base/interval.h \
  ../../include/geometry/../geometry/../base/binary_word.h \
  ../../include/geometry/../geometry/../base/../base/basic_type.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/list_set.h \
  ../../include/geometry/../geometry/parallelotope.h \
  ../../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/../geometry/list_set.h \
  ../../include/geometry/../geometry/../geometry/polyhedron.h \
  ../../include/geometry/../geometry/../geometry/../linear_algebra/constraint.h \
  ../../include/geometry/../geometry/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../geometry/../geometry/generator.h \
  ../../include/geometry/../geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/../geometry/../geometry/geometry_declarations.h \
  ../../include/geometry/../geometry/../geometry/zonotope.h \
  ../../include/geometry/../geometry/polyhedron.h \
  ../../include/geometry/zonotope.tpl \
  ../../include/geometry/../geometry/zonotope.h

simplex.o: simplex.cc ../../include/geometry/simplex.h \
  ../../include/geometry/../base/array.h \
  ../../include/geometry/../base/utility.h \
  ../../include/geometry/../base/../base/array.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/../base/basic_type.h \
  ../../include/geometry/../linear_algebra/../base/../base/array.h \
  ../../include/geometry/../linear_algebra/../base/../base/sequence.h \
  ../../include/geometry/../linear_algebra/../base/numerical_type.h \
  ../../include/geometry/../linear_algebra/../base/../base/dyadic.h \
  ../../include/geometry/../linear_algebra/../base/interval.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../base/utility.h \
  ../../include/geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../geometry/geometry_declarations.h \
  ../../include/geometry/../geometry/polyhedron.h \
  ../../include/geometry/../geometry/../linear_algebra/constraint.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/../base/basic_type.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/../base/numerical_type.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/../base/interval.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/../base/utility.h \
  ../../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../geometry/generator.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/../geometry/../base/utility.h \
  ../../include/geometry/../geometry/../geometry/../base/interval.h \
  ../../include/geometry/../geometry/../geometry/../base/binary_word.h \
  ../../include/geometry/../geometry/../geometry/../base/../base/basic_type.h \
  ../../include/geometry/../geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../geometry/../geometry/geometry_declarations.h \
  ../../include/geometry/../geometry/list_set.h \
  ../../include/geometry/../geometry/geometry_declarations.h \
  ../../include/geometry/simplex.tpl \
  ../../include/geometry/../geometry/simplex.h

polyhedron.o: polyhedron.cc ../../include/geometry/polyhedron.h \
  ../../include/geometry/../linear_algebra/constraint.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../base/basic_type.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../base/../base/array.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../base/../base/sequence.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../base/numerical_type.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../base/../base/dyadic.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../base/interval.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../base/utility.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/generator.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../base/utility.h \
  ../../include/geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../geometry/geometry_declarations.h \
  ../../include/geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/../base/interval.h \
  ../../include/geometry/../geometry/../base/binary_word.h \
  ../../include/geometry/../geometry/../base/../base/basic_type.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/geometry_declarations.h \
  ../../include/geometry/polyhedron.tpl \
  ../../include/geometry/../base/utility.h \
  ../../include/geometry/../base/interval.h \
  ../../include/geometry/../geometry/list_set.h \
  ../../include/geometry/../geometry/parallelotope.h \
  ../../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/../geometry/list_set.h \
  ../../include/geometry/../geometry/../geometry/polyhedron.h \
  ../../include/geometry/../geometry/../geometry/zonotope.h \
  ../../include/geometry/../geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../geometry/../base/utility.h \
  ../../include/geometry/../geometry/../geometry/../base/interval.h \
  ../../include/geometry/../geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/../geometry/../geometry/list_set.h \
  ../../include/geometry/../geometry/../geometry/../geometry/parallelotope.h \
  ../../include/geometry/../geometry/../geometry/../geometry/polyhedron.h \
  ../../include/geometry/../geometry/polyhedron.h

lattice_set.o: lattice_set.cc ../../include/base/array_operations.h \
  ../../include/base/../base/array.h \
  ../../include/base/../base/basic_type.h \
  ../../include/base/../base/../base/array.h \
  ../../include/base/../base/../base/sequence.h \
  ../../include/geometry/lattice_set.h \
  ../../include/geometry/../base/array.h \
  ../../include/geometry/../base/iterator.h \
  ../../include/geometry/../base/../base/array.h \
  ../../include/geometry/../base/basic_type.h \
  ../../include/geometry/../base/interval.h \
  ../../include/geometry/../base/array_operations.h

grid.o: grid.cc ../../include/geometry/grid.h \
  ../../include/geometry/../base/array.h \
  ../../include/geometry/../base/interval.h \
  ../../include/geometry/../base/binary_word.h \
  ../../include/geometry/../base/../base/basic_type.h \
  ../../include/geometry/../base/../base/../base/array.h \
  ../../include/geometry/../base/../base/../base/sequence.h \
  ../../include/geometry/../geometry/geometry_declarations.h \
  ../../include/geometry/../geometry/lattice_set.h \
  ../../include/geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../base/iterator.h \
  ../../include/geometry/../geometry/../base/../base/array.h \
  ../../include/geometry/../geometry/../base/basic_type.h \
  ../../include/geometry/../geometry/../base/interval.h \
  ../../include/geometry/../geometry/../base/array_operations.h \
  ../../include/geometry/../geometry/../base/../base/basic_type.h \
  ../../include/geometry/grid.tpl \
  ../../include/geometry/../base/arithmetic.h \
  ../../include/geometry/../base/../base/numerical_type.h \
  ../../include/geometry/../base/../base/../base/dyadic.h \
  ../../include/geometry/../geometry/grid.h \
  ../../include/geometry/../geometry/grid_set.h \
  ../../include/geometry/../geometry/../base/binary_word.h \
  ../../include/geometry/../geometry/../geometry/geometry_declarations.h \
  ../../include/geometry/../geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/../geometry/../base/utility.h \
  ../../include/geometry/../geometry/../geometry/../base/../base/array.h \
  ../../include/geometry/../geometry/../geometry/../base/interval.h \
  ../../include/geometry/../geometry/../geometry/../base/binary_word.h \
  ../../include/geometry/../geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../geometry/../geometry/../base/utility.h \
  ../../include/geometry/../geometry/../geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../geometry/../geometry/../linear_algebra/../base/basic_type.h \
  ../../include/geometry/../geometry/../geometry/../geometry/../linear_algebra/../base/numerical_type.h \
  ../../include/geometry/../geometry/../geometry/../geometry/../linear_algebra/../base/interval.h \
  ../../include/geometry/../geometry/../geometry/../geometry/../geometry/geometry_declarations.h \
  ../../include/geometry/../geometry/../geometry/../geometry/geometry_declarations.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../geometry/grid.h \
  ../../include/geometry/../geometry/../geometry/lattice_set.h \
  ../../include/geometry/../geometry/../geometry/list_set.h \
  ../../include/geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/list_set.h

grid_set.o: grid_set.cc ../../include/geometry/grid_set.h \
  ../../include/geometry/../base/array.h \
  ../../include/geometry/../base/interval.h \
  ../../include/geometry/../base/binary_word.h \
  ../../include/geometry/../base/../base/basic_type.h \
  ../../include/geometry/../base/../base/../base/array.h \
  ../../include/geometry/../base/../base/../base/sequence.h \
  ../../include/geometry/../geometry/geometry_declarations.h \
  ../../include/geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/../base/utility.h \
  ../../include/geometry/../geometry/../base/../base/array.h \
  ../../include/geometry/../geometry/../base/interval.h \
  ../../include/geometry/../geometry/../base/binary_word.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../geometry/../base/utility.h \
  ../../include/geometry/../geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../geometry/../linear_algebra/../base/basic_type.h \
  ../../include/geometry/../geometry/../geometry/../linear_algebra/../base/numerical_type.h \
  ../../include/geometry/../geometry/../geometry/../linear_algebra/../base/../base/dyadic.h \
  ../../include/geometry/../geometry/../geometry/../linear_algebra/../base/interval.h \
  ../../include/geometry/../geometry/../geometry/../geometry/geometry_declarations.h \
  ../../include/geometry/../geometry/../geometry/geometry_declarations.h \
  ../../include/geometry/../geometry/point.h \
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
  ../../include/geometry/../geometry/list_set.h \
  ../../include/geometry/grid_set.tpl \
  ../../include/geometry/../base/arithmetic.h \
  ../../include/geometry/../base/../base/numerical_type.h \
  ../../include/geometry/../geometry/parallelotope.h \
  ../../include/geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../linear_algebra/../base/basic_type.h \
  ../../include/geometry/../geometry/../linear_algebra/../base/numerical_type.h \
  ../../include/geometry/../geometry/../linear_algebra/../base/interval.h \
  ../../include/geometry/../geometry/../linear_algebra/../base/utility.h \
  ../../include/geometry/../geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/../geometry/list_set.h \
  ../../include/geometry/../geometry/../geometry/polyhedron.h \
  ../../include/geometry/../geometry/../geometry/../linear_algebra/constraint.h \
  ../../include/geometry/../geometry/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../geometry/../geometry/generator.h \
  ../../include/geometry/../geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/../geometry/zonotope.h \
  ../../include/geometry/../geometry/../geometry/../geometry/list_set.h \
  ../../include/geometry/../geometry/../geometry/../geometry/parallelotope.h \
  ../../include/geometry/../geometry/../geometry/../geometry/polyhedron.h \
  ../../include/geometry/../geometry/partition_tree_set.h \
  ../../include/geometry/../geometry/../base/basic_type.h \
  ../../include/geometry/../geometry/../base/binary_tree.h \
  ../../include/geometry/../geometry/../base/../base/basic_type.h \
  ../../include/geometry/../geometry/../base/../base/binary_word.h \
  ../../include/geometry/../geometry/../base/iterator.h \
  ../../include/geometry/../geometry/../geometry/subdivision_tree_set.h \
  ../../include/geometry/../geometry/../geometry/../base/sequence.h \
  ../../include/geometry/../geometry/../geometry/../base/binary_word.h \
  ../../include/geometry/../geometry/../geometry/../base/binary_tree.h \
  ../../include/geometry/../geometry/grid_set.h

subdivision_tree_set.o: subdivision_tree_set.cc \
  ../../include/base/arithmetic.h \
  ../../include/base/../base/numerical_type.h \
  ../../include/base/../base/../base/dyadic.h \
  ../../include/geometry/rectangle.h \
  ../../include/geometry/../base/utility.h \
  ../../include/geometry/../base/../base/array.h \
  ../../include/geometry/../base/interval.h \
  ../../include/geometry/../base/binary_word.h \
  ../../include/geometry/../base/../base/basic_type.h \
  ../../include/geometry/../base/../base/../base/array.h \
  ../../include/geometry/../base/../base/../base/sequence.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../base/utility.h \
  ../../include/geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../base/basic_type.h \
  ../../include/geometry/../geometry/../linear_algebra/../base/numerical_type.h \
  ../../include/geometry/../geometry/../linear_algebra/../base/interval.h \
  ../../include/geometry/../geometry/../geometry/geometry_declarations.h \
  ../../include/geometry/../geometry/geometry_declarations.h \
  ../../include/geometry/lattice_set.h \
  ../../include/geometry/../base/array.h \
  ../../include/geometry/../base/iterator.h \
  ../../include/geometry/../base/basic_type.h \
  ../../include/geometry/../base/array_operations.h \
  ../../include/geometry/subdivision_tree_set.h \
  ../../include/geometry/../base/sequence.h \
  ../../include/geometry/../base/binary_tree.h \
  ../../include/geometry/../base/../base/binary_word.h \
  ../../include/geometry/../geometry/rectangle.h

partition_tree_set.o: partition_tree_set.cc \
  ../../include/geometry/partition_tree_set.h \
  ../../include/geometry/../base/basic_type.h \
  ../../include/geometry/../base/../base/array.h \
  ../../include/geometry/../base/../base/sequence.h \
  ../../include/geometry/../base/binary_word.h \
  ../../include/geometry/../base/../base/basic_type.h \
  ../../include/geometry/../base/binary_tree.h \
  ../../include/geometry/../base/../base/binary_word.h \
  ../../include/geometry/../base/iterator.h \
  ../../include/geometry/../base/utility.h \
  ../../include/geometry/../geometry/geometry_declarations.h \
  ../../include/geometry/../geometry/subdivision_tree_set.h \
  ../../include/geometry/../geometry/../base/basic_type.h \
  ../../include/geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../base/sequence.h \
  ../../include/geometry/../geometry/../base/interval.h \
  ../../include/geometry/../geometry/../base/binary_word.h \
  ../../include/geometry/../geometry/../base/binary_tree.h \
  ../../include/geometry/../geometry/../base/iterator.h \
  ../../include/geometry/../geometry/../base/utility.h \
  ../../include/geometry/../geometry/../geometry/geometry_declarations.h \
  ../../include/geometry/../geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/../geometry/../base/utility.h \
  ../../include/geometry/../geometry/../geometry/../base/interval.h \
  ../../include/geometry/../geometry/../geometry/../base/binary_word.h \
  ../../include/geometry/../geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../geometry/../geometry/../base/utility.h \
  ../../include/geometry/../geometry/../geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../geometry/../geometry/../linear_algebra/../base/basic_type.h \
  ../../include/geometry/../geometry/../geometry/../geometry/../linear_algebra/../base/numerical_type.h \
  ../../include/geometry/../geometry/../geometry/../geometry/../linear_algebra/../base/../base/dyadic.h \
  ../../include/geometry/../geometry/../geometry/../geometry/../linear_algebra/../base/interval.h \
  ../../include/geometry/../geometry/../geometry/../geometry/../geometry/geometry_declarations.h \
  ../../include/geometry/../geometry/../geometry/../geometry/geometry_declarations.h \
  ../../include/geometry/partition_tree_set.tpl \
  ../../include/geometry/../geometry/partition_tree_set.h \
  ../../include/geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/list_set.h \
  ../../include/geometry/../geometry/grid_set.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../geometry/grid.h \
  ../../include/geometry/../geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../geometry/../geometry/lattice_set.h \
  ../../include/geometry/../geometry/../geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../geometry/../geometry/../base/iterator.h \
  ../../include/geometry/../geometry/../geometry/../geometry/../base/basic_type.h \
  ../../include/geometry/../geometry/../geometry/../geometry/../base/interval.h \
  ../../include/geometry/../geometry/../geometry/../geometry/../base/array_operations.h \
  ../../include/geometry/../geometry/../geometry/../geometry/../base/../base/array.h \
  ../../include/geometry/../geometry/../geometry/../geometry/../base/../base/basic_type.h \
  ../../include/geometry/../geometry/../geometry/lattice_set.h \
  ../../include/geometry/../geometry/../geometry/list_set.h

