include ../../config.mk

OBJ = parallelopiped.o \
	grid_operations.o unit_grid_set.o grid.o grid_set.o \
	partition_tree_operations.o unit_partition_tree_set.o partition_tree_set.o

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
parallelopiped.o: parallelopiped.cc \
  ../../include/geometry/parallelopiped.h ../../include/base/utility.h \
  ../../include/base/array.h ../../include/base/interval.h \
  ../../include/linear_algebra/vector.h ../../include/base/basic_type.h \
  ../../include/base/sequence.h ../../include/base/numerical_type.h \
  ../../include/base/dyadic.h ../../include/base/interval.h \
  ../../include/linear_algebra/matrix.h \
  ../../include/linear_algebra/constraint.h \
  ../../include/linear_algebra/vector.h \
  ../../include/linear_algebra/matrix.h ../../include/geometry/point.h \
  ../../include/geometry/geometry_declarations.h \
  ../../include/geometry/rectangle.h ../../include/base/binary_word.h \
  ../../include/base/utility.h ../../include/base/basic_type.h \
  ../../include/geometry/list_set.h ../../include/geometry/polyhedron.h \
  ../../include/geometry/generator.h \
  ../../include/geometry/parallelopiped.tpl \
  ../../include/geometry/unit_grid_set.h ../../include/base/array.h \
  ../../include/base/iterator.h ../../include/base/../base/array.h \
  ../../include/base/../base/basic_type.h ../../include/base/basic_type.h \
  ../../include/geometry/grid_operations.h \
  ../../include/geometry/parallelopiped.h

grid_operations.o: grid_operations.cc \
  ../../include/geometry/grid_operations.h ../../include/base/array.h \
  ../../include/base/basic_type.h ../../include/base/array.h \
  ../../include/base/sequence.h ../../include/geometry/rectangle.h \
  ../../include/base/utility.h ../../include/base/interval.h \
  ../../include/base/binary_word.h ../../include/base/utility.h \
  ../../include/base/basic_type.h ../../include/geometry/point.h \
  ../../include/linear_algebra/vector.h ../../include/base/basic_type.h \
  ../../include/base/numerical_type.h ../../include/base/dyadic.h \
  ../../include/base/interval.h \
  ../../include/geometry/geometry_declarations.h \
  ../../include/geometry/unit_grid_set.h ../../include/base/iterator.h \
  ../../include/base/../base/array.h \
  ../../include/base/../base/basic_type.h \
  ../../include/geometry/grid_operations.h

unit_grid_set.o: unit_grid_set.cc ../../include/geometry/unit_grid_set.h \
  ../../include/base/array.h ../../include/base/iterator.h \
  ../../include/base/../base/array.h \
  ../../include/base/../base/basic_type.h ../../include/base/array.h \
  ../../include/base/sequence.h ../../include/base/basic_type.h \
  ../../include/base/interval.h ../../include/geometry/grid_operations.h \
  ../../include/geometry/rectangle.h ../../include/base/utility.h \
  ../../include/base/array.h ../../include/base/binary_word.h \
  ../../include/base/utility.h ../../include/base/basic_type.h \
  ../../include/geometry/point.h ../../include/linear_algebra/vector.h \
  ../../include/base/basic_type.h ../../include/base/numerical_type.h \
  ../../include/base/dyadic.h ../../include/base/interval.h \
  ../../include/geometry/geometry_declarations.h \
  ../../include/geometry/grid_operations.h

grid.o: grid.cc ../../include/geometry/grid.h ../../include/base/array.h \
  ../../include/base/interval.h ../../include/base/binary_word.h \
  ../../include/base/utility.h ../../include/base/array.h \
  ../../include/base/basic_type.h ../../include/base/sequence.h \
  ../../include/geometry/geometry_declarations.h \
  ../../include/geometry/rectangle.h ../../include/base/utility.h \
  ../../include/geometry/point.h ../../include/linear_algebra/vector.h \
  ../../include/base/basic_type.h ../../include/base/numerical_type.h \
  ../../include/base/dyadic.h ../../include/base/interval.h \
  ../../include/geometry/unit_grid_set.h ../../include/base/iterator.h \
  ../../include/base/../base/array.h \
  ../../include/base/../base/basic_type.h ../../include/base/basic_type.h \
  ../../include/geometry/grid_operations.h \
  ../../include/geometry/grid.tpl \
  ../../include/geometry/../base/arithmetic.h \
  ../../include/base/numerical_type.h \
  ../../include/geometry/../geometry/grid.h \
  ../../include/geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/list_set.h \
  ../../include/geometry/geometry_declarations.h

grid_set.o: grid_set.cc ../../include/geometry/grid_set.h \
  ../../include/base/array.h ../../include/base/interval.h \
  ../../include/base/binary_word.h ../../include/base/utility.h \
  ../../include/base/array.h ../../include/base/basic_type.h \
  ../../include/base/sequence.h \
  ../../include/geometry/geometry_declarations.h \
  ../../include/geometry/rectangle.h ../../include/base/utility.h \
  ../../include/geometry/point.h ../../include/linear_algebra/vector.h \
  ../../include/base/basic_type.h ../../include/base/numerical_type.h \
  ../../include/base/dyadic.h ../../include/base/interval.h \
  ../../include/geometry/grid.h ../../include/geometry/unit_grid_set.h \
  ../../include/base/iterator.h ../../include/base/../base/array.h \
  ../../include/base/../base/basic_type.h ../../include/base/basic_type.h \
  ../../include/geometry/grid_operations.h \
  ../../include/geometry/unit_grid_set.h \
  ../../include/geometry/list_set.h ../../include/geometry/grid_set.tpl \
  ../../include/geometry/../base/arithmetic.h \
  ../../include/base/numerical_type.h \
  ../../include/geometry/../geometry/parallelopiped.h \
  ../../include/base/utility.h ../../include/base/interval.h \
  ../../include/linear_algebra/vector.h \
  ../../include/linear_algebra/matrix.h \
  ../../include/linear_algebra/constraint.h \
  ../../include/linear_algebra/vector.h \
  ../../include/linear_algebra/matrix.h ../../include/geometry/point.h \
  ../../include/geometry/rectangle.h ../../include/geometry/list_set.h \
  ../../include/geometry/polyhedron.h \
  ../../include/linear_algebra/constraint.h \
  ../../include/linear_algebra/matrix.h \
  ../../include/geometry/generator.h \
  ../../include/geometry/geometry_declarations.h \
  ../../include/geometry/../geometry/partition_tree_set.h \
  ../../include/base/basic_type.h ../../include/base/binary_word.h \
  ../../include/base/binary_tree.h ../../include/base/binary_word.h \
  ../../include/base/iterator.h \
  ../../include/geometry/unit_partition_tree_set.h \
  ../../include/base/sequence.h ../../include/base/binary_tree.h \
  ../../include/geometry/../geometry/grid_set.h

partition_tree_operations.o: partition_tree_operations.cc \
  ../../include/geometry/partition_tree_operations.h \
  ../../include/base/basic_type.h ../../include/base/array.h \
  ../../include/base/sequence.h \
  ../../include/geometry/geometry_declarations.h

unit_partition_tree_set.o: unit_partition_tree_set.cc \
  ../../include/base/arithmetic.h ../../include/base/numerical_type.h \
  ../../include/base/dyadic.h ../../include/geometry/rectangle.h \
  ../../include/base/utility.h ../../include/base/array.h \
  ../../include/base/interval.h ../../include/base/binary_word.h \
  ../../include/base/utility.h ../../include/base/basic_type.h \
  ../../include/base/sequence.h ../../include/geometry/point.h \
  ../../include/linear_algebra/vector.h ../../include/base/basic_type.h \
  ../../include/base/numerical_type.h ../../include/base/interval.h \
  ../../include/geometry/geometry_declarations.h \
  ../../include/geometry/unit_grid_set.h ../../include/base/array.h \
  ../../include/base/iterator.h ../../include/base/../base/array.h \
  ../../include/base/../base/basic_type.h ../../include/base/basic_type.h \
  ../../include/geometry/grid_operations.h \
  ../../include/geometry/rectangle.h \
  ../../include/geometry/unit_partition_tree_set.h \
  ../../include/base/sequence.h ../../include/base/binary_tree.h \
  ../../include/base/binary_word.h \
  ../../include/geometry/partition_tree_operations.h

partition_tree_set.o: partition_tree_set.cc \
  ../../include/geometry/partition_tree_set.h \
  ../../include/base/basic_type.h ../../include/base/array.h \
  ../../include/base/sequence.h ../../include/base/binary_word.h \
  ../../include/base/utility.h ../../include/base/basic_type.h \
  ../../include/base/binary_tree.h ../../include/base/binary_word.h \
  ../../include/base/iterator.h ../../include/base/../base/array.h \
  ../../include/base/../base/basic_type.h ../../include/base/utility.h \
  ../../include/geometry/geometry_declarations.h \
  ../../include/geometry/unit_partition_tree_set.h \
  ../../include/base/array.h ../../include/base/sequence.h \
  ../../include/base/interval.h ../../include/geometry/rectangle.h \
  ../../include/geometry/point.h ../../include/linear_algebra/vector.h \
  ../../include/base/basic_type.h ../../include/base/numerical_type.h \
  ../../include/base/dyadic.h ../../include/base/interval.h \
  ../../include/geometry/partition_tree_set.tpl \
  ../../include/geometry/partition_tree_set.h \
  ../../include/geometry/partition_tree_operations.h \
  ../../include/geometry/list_set.h ../../include/geometry/grid_set.h \
  ../../include/geometry/grid.h ../../include/geometry/unit_grid_set.h \
  ../../include/geometry/grid_operations.h \
  ../../include/geometry/unit_grid_set.h

