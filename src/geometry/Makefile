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

