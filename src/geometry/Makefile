include ../../config.mk

OBJ = point.o \
	rectangle.o parallelotope.o zonotope.o simplex.o polyhedron.o \
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

### Do Not edit this line
### Everything below is auto-generated
point.o: point.cc ../../include/geometry/point.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../linear_algebra/../numeric/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/../base/utility.h \
  ../../include/geometry/../linear_algebra/../numeric/../base/../declarations.h \
  ../../include/geometry/point.tpl ../../include/geometry/point.h \
  ../../include/real_typedef.h ../../include/numeric/numerical_types.h \
  ../../include/numeric/../declarations.h \
  ../../include/numeric/../base/utility.h \
  ../../include/numeric/../numeric/dyadic.h

rectangle.o: rectangle.cc ../../include/geometry/rectangle.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../numeric/arithmetic.h \
  ../../include/geometry/../numeric/../numeric/numerical_types.h \
  ../../include/geometry/../numeric/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/../base/utility.h \
  ../../include/geometry/../numeric/../numeric/../base/../declarations.h \
  ../../include/geometry/../numeric/../numeric/../numeric/dyadic.h \
  ../../include/geometry/../numeric/../numeric/approximation.h \
  ../../include/geometry/../numeric/../numeric/../numeric/numerical_types.h \
  ../../include/geometry/../numeric/interval.h \
  ../../include/geometry/../numeric/../declarations.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/../declarations.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/../base/utility.h \
  ../../include/geometry/rectangle.tpl ../../include/geometry/rectangle.h \
  ../../include/geometry/../base/binary_word.h \
  ../../include/geometry/../base/../declarations.h \
  ../../include/geometry/../linear_algebra/interval_vector.h \
  ../../include/geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/lattice_set.h \
  ../../include/geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../base/../declarations.h \
  ../../include/geometry/../geometry/../base/iterator.h \
  ../../include/geometry/../geometry/../base/../base/array.h \
  ../../include/geometry/../geometry/../numeric/interval.h \
  ../../include/geometry/../geometry/../base/array_operations.h \
  ../../include/geometry/../geometry/list_set.h \
  ../../include/geometry/../geometry/../base/utility.h \
  ../../include/geometry/../geometry/grid_set.h \
  ../../include/geometry/../geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/../geometry/grid.h \
  ../../include/geometry/../geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../geometry/../geometry/lattice_set.h \
  ../../include/geometry/../geometry/../geometry/lattice_set.h \
  ../../include/real_typedef.h ../../include/numeric/numerical_types.h

parallelotope.o: parallelotope.cc ../../include/geometry/parallelotope.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../numeric/interval.h \
  ../../include/geometry/../numeric/../declarations.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../linear_algebra/../numeric/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/../base/utility.h \
  ../../include/geometry/../linear_algebra/../numeric/../base/../declarations.h \
  ../../include/geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/parallelotope.tpl \
  ../../include/geometry/parallelotope.h \
  ../../include/geometry/../utility/stlio.h \
  ../../include/geometry/../utility/../base/array.h \
  ../../include/geometry/../utility/../base/../declarations.h \
  ../../include/geometry/../base/utility.h \
  ../../include/geometry/../linear_algebra/linear_program.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/lattice_set.h \
  ../../include/geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../base/iterator.h \
  ../../include/geometry/../geometry/../base/../base/array.h \
  ../../include/geometry/../geometry/../numeric/interval.h \
  ../../include/geometry/../geometry/../base/array_operations.h \
  ../../include/geometry/../geometry/../base/../declarations.h \
  ../../include/geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../numeric/../numeric/numerical_types.h \
  ../../include/geometry/../geometry/../numeric/../numeric/../declarations.h \
  ../../include/geometry/../geometry/../numeric/../numeric/../base/utility.h \
  ../../include/geometry/../geometry/../numeric/../numeric/../numeric/dyadic.h \
  ../../include/geometry/../geometry/../numeric/../numeric/approximation.h \
  ../../include/geometry/../geometry/../numeric/../numeric/../numeric/numerical_types.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/list_set.h \
  ../../include/geometry/../geometry/../base/utility.h \
  ../../include/geometry/../geometry/zonotope.h \
  ../../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../linear_algebra/zonotopic_vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/interval_vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/../declarations.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/interval_matrix.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/polyhedron.h \
  /ufs/collins/include/ppl.hh \
  ../../include/geometry/../geometry/../linear_algebra/constraint.h \
  ../../include/geometry/../geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/numerical_types.h \
  ../../include/geometry/../geometry/../geometry/generator.h \
  ../../include/geometry/../geometry/../geometry/rectangle.h \
  ../../include/real_typedef.h ../../include/numeric/numerical_types.h

zonotope.o: zonotope.cc ../../include/real_typedef.h \
  ../../include/numeric/numerical_types.h \
  ../../include/numeric/../declarations.h \
  ../../include/numeric/../base/utility.h \
  ../../include/numeric/../base/../declarations.h \
  ../../include/numeric/../numeric/dyadic.h \
  ../../include/geometry/zonotope.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../linear_algebra/../numeric/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/../base/utility.h \
  ../../include/geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../linear_algebra/zonotopic_vector.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/interval_vector.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../numeric/../declarations.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/interval_matrix.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../numeric/interval.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/zonotope.tpl ../../include/geometry/zonotope.h \
  ../../include/geometry/../base/array.h \
  ../../include/geometry/../base/../declarations.h \
  ../../include/geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../numeric/../numeric/numerical_types.h \
  ../../include/geometry/../geometry/../numeric/../numeric/approximation.h \
  ../../include/geometry/../geometry/../numeric/../numeric/../declarations.h \
  ../../include/geometry/../geometry/../numeric/../numeric/../numeric/numerical_types.h \
  ../../include/geometry/../geometry/../numeric/interval.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/parallelotope.h \
  ../../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/polyhedron.h \
  /ufs/collins/include/ppl.hh \
  ../../include/geometry/../geometry/../linear_algebra/constraint.h \
  ../../include/geometry/../geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/numerical_types.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../geometry/generator.h \
  ../../include/geometry/../geometry/../geometry/rectangle.h

simplex.o: simplex.cc ../../include/geometry/simplex.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../base/array.h \
  ../../include/geometry/../base/../declarations.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/../declarations.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/../base/utility.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/../base/../declarations.h \
  ../../include/geometry/simplex.tpl ../../include/geometry/simplex.h \
  ../../include/geometry/../utility/stlio.h \
  ../../include/geometry/../utility/../base/array.h \
  ../../include/geometry/../geometry/parallelotope.h \
  ../../include/geometry/../geometry/../numeric/interval.h \
  ../../include/geometry/../geometry/../numeric/../declarations.h \
  ../../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/polyhedron.h \
  /ufs/collins/include/ppl.hh \
  ../../include/geometry/../geometry/../linear_algebra/constraint.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/numerical_types.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/../numeric/dyadic.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/../geometry/generator.h \
  ../../include/geometry/../geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../geometry/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../geometry/../numeric/../numeric/numerical_types.h \
  ../../include/geometry/../geometry/../geometry/../numeric/../numeric/approximation.h \
  ../../include/geometry/../geometry/../geometry/../numeric/../numeric/../declarations.h \
  ../../include/geometry/../geometry/../geometry/../numeric/../numeric/../numeric/numerical_types.h \
  ../../include/geometry/../geometry/../geometry/../numeric/interval.h \
  ../../include/geometry/../geometry/../geometry/../geometry/point.h \
  ../../include/real_typedef.h ../../include/numeric/numerical_types.h

