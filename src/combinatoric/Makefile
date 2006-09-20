include ../../config.mk

OBJ = binary_tree.o lattice_set.o lattice_map.o subdivision_tree_set.o

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
binary_tree.o: binary_tree.cc ../../include/combinatoric/binary_tree.h \
  ../../include/combinatoric/../declarations.h \
  ../../include/combinatoric/../base/array.h \
  ../../include/combinatoric/../base/../declarations.h \
  ../../include/combinatoric/../combinatoric/binary_word.h \
  ../../include/combinatoric/../combinatoric/../declarations.h \
  ../../include/base/array.h ../../include/base/stlio.h \
  ../../include/base/../base/array.h

lattice_set.o: lattice_set.cc ../../include/base/stlio.h \
  ../../include/base/../base/array.h \
  ../../include/base/../base/../declarations.h \
  ../../include/base/array_operations.h \
  ../../include/base/../declarations.h \
  ../../include/combinatoric/lattice_set.h \
  ../../include/combinatoric/../base/array.h \
  ../../include/combinatoric/../base/iterator.h \
  ../../include/combinatoric/../base/../base/array.h \
  ../../include/combinatoric/../base/array_operations.h \
  ../../include/combinatoric/../numeric/interval.h \
  ../../include/combinatoric/../numeric/../declarations.h \
  ../../include/combinatoric/../numeric/../numeric/numerical_traits.h \
  ../../include/combinatoric/../numeric/../numeric/approximation.h \
  ../../include/combinatoric/../numeric/../numeric/arithmetic.h \
  ../../include/combinatoric/../numeric/../numeric/rounding.h \
  ../../include/combinatoric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/combinatoric/../numeric/../numeric/../numeric/function.h

lattice_map.o: lattice_map.cc ../../include/debug.h \
  ../../include/base/stlio.h ../../include/base/../base/array.h \
  ../../include/base/../base/../declarations.h \
  ../../include/combinatoric/lattice_set.h \
  ../../include/combinatoric/../base/array.h \
  ../../include/combinatoric/../base/iterator.h \
  ../../include/combinatoric/../base/../base/array.h \
  ../../include/combinatoric/../base/array_operations.h \
  ../../include/combinatoric/../base/../declarations.h \
  ../../include/combinatoric/../numeric/interval.h \
  ../../include/combinatoric/../numeric/../declarations.h \
  ../../include/combinatoric/../numeric/../numeric/numerical_traits.h \
  ../../include/combinatoric/../numeric/../numeric/approximation.h \
  ../../include/combinatoric/../numeric/../numeric/arithmetic.h \
  ../../include/combinatoric/../numeric/../numeric/rounding.h \
  ../../include/combinatoric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/combinatoric/../numeric/../numeric/../numeric/function.h \
  ../../include/combinatoric/lattice_map.h \
  ../../include/combinatoric/../combinatoric/lattice_set.h

subdivision_tree_set.o: subdivision_tree_set.cc \
  ../../include/numeric/arithmetic.h \
  ../../include/combinatoric/lattice_set.h \
  ../../include/combinatoric/../base/array.h \
  ../../include/combinatoric/../base/../declarations.h \
  ../../include/combinatoric/../base/iterator.h \
  ../../include/combinatoric/../base/../base/array.h \
  ../../include/combinatoric/../base/array_operations.h \
  ../../include/combinatoric/../numeric/interval.h \
  ../../include/combinatoric/../numeric/../declarations.h \
  ../../include/combinatoric/../numeric/../numeric/numerical_traits.h \
  ../../include/combinatoric/../numeric/../numeric/approximation.h \
  ../../include/combinatoric/../numeric/../numeric/arithmetic.h \
  ../../include/combinatoric/../numeric/../numeric/rounding.h \
  ../../include/combinatoric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/combinatoric/../numeric/../numeric/../numeric/function.h \
  ../../include/combinatoric/subdivision_tree_set.h \
  ../../include/combinatoric/../declarations.h \
  ../../include/combinatoric/../base/sequence.h \
  ../../include/combinatoric/../combinatoric/binary_word.h \
  ../../include/combinatoric/../combinatoric/../declarations.h \
  ../../include/combinatoric/../combinatoric/binary_tree.h \
  ../../include/combinatoric/../combinatoric/../base/array.h \
  ../../include/combinatoric/../combinatoric/../combinatoric/binary_word.h \
  ../../include/combinatoric/../geometry/rectangle.h \
  ../../include/combinatoric/../geometry/../declarations.h \
  ../../include/combinatoric/../geometry/../numeric/arithmetic.h \
  ../../include/combinatoric/../geometry/../numeric/function.h \
  ../../include/combinatoric/../geometry/../numeric/interval.h \
  ../../include/combinatoric/../geometry/../linear_algebra/vector.h \
  ../../include/combinatoric/../geometry/../linear_algebra/../declarations.h \
  ../../include/combinatoric/../geometry/../linear_algebra/../numeric/integer.h \
  ../../include/combinatoric/../geometry/../linear_algebra/../numeric/../declarations.h \
  ../../include/combinatoric/../geometry/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/combinatoric/../geometry/../linear_algebra/../numeric/../numeric/arithmetic.h \
  ../../include/combinatoric/../geometry/../linear_algebra/../numeric/arithmetic.h \
  ../../include/combinatoric/../geometry/../linear_algebra/../numeric/interval.h \
  ../../include/combinatoric/../geometry/../geometry/ppl_polyhedron.h \
  ../../include/combinatoric/../geometry/../geometry/../declarations.h \
  ../../include/combinatoric/../geometry/../geometry/../numeric/rational.h \
  ../../include/combinatoric/../geometry/../geometry/../numeric/../declarations.h \
  ../../include/combinatoric/../geometry/../geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/combinatoric/../geometry/../geometry/../numeric/../numeric/arithmetic.h \
  ../../include/combinatoric/../geometry/../geometry/../numeric/../numeric/function.h \
  ../../include/combinatoric/../geometry/../geometry/../numeric/../numeric/integer.h \
  ../../include/combinatoric/../geometry/../geometry/../linear_algebra/vector.h \
  ../../include/combinatoric/../geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/combinatoric/../geometry/../geometry/../linear_algebra/../declarations.h \
  ../../include/combinatoric/../geometry/../geometry/../linear_algebra/../numeric/integer.h \
  ../../include/combinatoric/../geometry/../geometry/../linear_algebra/../numeric/interval.h \
  ../../include/combinatoric/../geometry/../geometry/point.h \
  ../../include/combinatoric/../geometry/../geometry/../base/array.h \
  ../../include/combinatoric/../geometry/../geometry/rectangle_expression.h \
  ../../include/geometry/rectangle.h

