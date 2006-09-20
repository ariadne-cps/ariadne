include ../../config.mk

OBJ = \
	ppl_polyhedron.o \
	point.o point_list.o \
	rectangle.o polyhedron.o polytope.o simplex.o zonotope.o parallelotope.o \
	sphere.o ellipsoid.o \
	list_set.o \
	grid.o grid_set.o \
	partition_tree_set.o

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
ppl_polyhedron.o: ppl_polyhedron.cc \
  ../../include/numeric/numerical_types.h ../../include/numeric/mpfloat.h \
  /home/collins/include/mpfr.h ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/approximation.h \
  ../../include/numeric/../numeric/function.h \
  ../../include/numeric/../numeric/rational.h \
  ../../include/numeric/../numeric/../declarations.h \
  ../../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/../numeric/arithmetic.h \
  ../../include/numeric/../numeric/../numeric/function.h \
  ../../include/numeric/../numeric/../numeric/integer.h \
  ../../include/numeric/../numeric/../numeric/../declarations.h \
  ../../include/numeric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/../numeric/../numeric/arithmetic.h \
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
  ../../include/geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../linear_algebra/matrix.h \
  /home/collins/include/ppl.hh

point.o: point.cc ../../include/geometry/point.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../base/array.h \
  ../../include/geometry/../base/../declarations.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../linear_algebra/../numeric/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/approximation.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../../include/geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/rounding.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/function.h \
  ../../include/geometry/../linear_algebra/matrix.h \
  ../../include/geometry/point.tpl ../../include/geometry/point.h \
  ../../include/geometry/../numeric/numerical_types.h \
  ../../include/geometry/../numeric/mpfloat.h \
  /home/collins/include/mpfr.h \
  ../../include/geometry/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../numeric/../numeric/function.h \
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

point_list.o: point_list.cc ../../include/geometry/point_list.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../linear_algebra/../numeric/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/approximation.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../../include/geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/rounding.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/function.h \
  ../../include/geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../base/../declarations.h \
  ../../include/geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/geometry/point_list.tpl \
  ../../include/geometry/point_list.h ../../include/numeric/float64.h \
  /home/collins/include/mpfr.h ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/function.h \
  ../../include/numeric/../numeric/integer.h \
  ../../include/numeric/mpfloat.h \
  ../../include/numeric/../numeric/rational.h \
  ../../include/numeric/../numeric/../declarations.h \
  ../../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/../numeric/arithmetic.h \
  ../../include/numeric/../numeric/../numeric/function.h \
  ../../include/numeric/../numeric/../numeric/integer.h \
  ../../include/numeric/rational.h

rectangle.o: rectangle.cc ../../include/geometry/rectangle.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../numeric/arithmetic.h \
  ../../include/geometry/../numeric/function.h \
  ../../include/geometry/../numeric/interval.h \
  ../../include/geometry/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../numeric/../numeric/approximation.h \
  ../../include/geometry/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../numeric/../numeric/rounding.h \
  ../../include/geometry/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../numeric/../numeric/../numeric/function.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../linear_algebra/../numeric/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../../include/geometry/../linear_algebra/../numeric/interval.h \
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
  ../../include/geometry/../geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../base/../declarations.h \
  ../../include/geometry/../geometry/rectangle_expression.h \
  ../../include/geometry/rectangle.tpl /home/collins/include/ppl.hh \
  ../../include/geometry/rectangle.h \
  ../../include/geometry/../combinatoric/binary_word.h \
  ../../include/geometry/../combinatoric/../declarations.h \
  ../../include/geometry/../base/array.h \
  ../../include/geometry/../geometry/point_list.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/list_set.h \
  ../../include/geometry/../geometry/../base/utility.h \
  ../../include/geometry/../geometry/grid_set.h \
  ../../include/geometry/../geometry/../combinatoric/lattice_set.h \
  ../../include/geometry/../geometry/../combinatoric/../base/array.h \
  ../../include/geometry/../geometry/../combinatoric/../base/iterator.h \
  ../../include/geometry/../geometry/../combinatoric/../base/../base/array.h \
  ../../include/geometry/../geometry/../combinatoric/../base/array_operations.h \
  ../../include/geometry/../geometry/../combinatoric/../base/../declarations.h \
  ../../include/geometry/../geometry/../combinatoric/../numeric/interval.h \
  ../../include/geometry/../geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/../geometry/grid.h \
  ../../include/geometry/../geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../geometry/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../geometry/../combinatoric/lattice_set.h \
  ../../include/numeric/float64.h /home/collins/include/mpfr.h \
  ../../include/numeric/../declarations.h \
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
  ../../include/geometry/../linear_algebra/../numeric/../numeric/approximation.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../../include/geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/rounding.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/function.h \
  ../../include/geometry/../linear_algebra/matrix.h \
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
  ../../include/geometry/../geometry/point_list.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../numeric/function.h \
  ../../include/geometry/../geometry/../numeric/interval.h \
  ../../include/geometry/../geometry/../geometry/ppl_polyhedron.h \
  ../../include/geometry/../geometry/../geometry/rectangle_expression.h \
  ../../include/geometry/../geometry/../geometry/../declarations.h \
  ../../include/geometry/polyhedron.tpl /home/collins/include/ppl.hh \
  ../../include/geometry/polyhedron.h \
  ../../include/geometry/../base/stlio.h \
  ../../include/geometry/../base/../base/array.h \
  ../../include/geometry/../numeric/interval.h \
  ../../include/geometry/../numeric/arithmetic.h \
  ../../include/real_typedef.h ../../include/numeric/numerical_types.h \
  ../../include/numeric/mpfloat.h /home/collins/include/mpfr.h \
  ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/function.h \
  ../../include/numeric/../numeric/rational.h

polytope.o: polytope.cc ../../include/geometry/polytope.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../linear_algebra/../numeric/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/approximation.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../../include/geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/rounding.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/function.h \
  ../../include/geometry/../linear_algebra/matrix.h \
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
  ../../include/geometry/../geometry/point_list.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../numeric/function.h \
  ../../include/geometry/../geometry/../numeric/interval.h \
  ../../include/geometry/../geometry/../geometry/ppl_polyhedron.h \
  ../../include/geometry/../geometry/../geometry/rectangle_expression.h \
  ../../include/geometry/../geometry/../geometry/../declarations.h \
  ../../include/geometry/polytope.tpl /home/collins/include/ppl.hh \
  ../../include/geometry/polytope.h \
  ../../include/geometry/../base/stlio.h \
  ../../include/geometry/../base/../base/array.h \
  ../../include/geometry/../numeric/interval.h \
  ../../include/geometry/../numeric/arithmetic.h \
  ../../include/real_typedef.h ../../include/numeric/numerical_types.h \
  ../../include/numeric/mpfloat.h /home/collins/include/mpfr.h \
  ../../include/numeric/../declarations.h \
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
  ../../include/geometry/../geometry/../linear_algebra/../numeric/../numeric/approximation.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/../numeric/rounding.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/../numeric/../numeric/function.h \
  ../../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/point_list.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/polytope.h \
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
  ../../include/geometry/../geometry/../geometry/point_list.h \
  ../../include/geometry/../geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/../geometry/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../geometry/../numeric/function.h \
  ../../include/geometry/../geometry/../geometry/../numeric/interval.h \
  ../../include/geometry/../geometry/../geometry/../geometry/ppl_polyhedron.h \
  ../../include/geometry/../geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../geometry/../geometry/rectangle_expression.h \
  ../../include/geometry/../geometry/../geometry/../geometry/../declarations.h \
  ../../include/geometry/simplex.tpl ../../include/geometry/simplex.h \
  /home/collins/include/ppl.hh ../../include/geometry/../base/stlio.h \
  ../../include/geometry/../base/../base/array.h \
  ../../include/real_typedef.h ../../include/numeric/numerical_types.h \
  ../../include/numeric/mpfloat.h /home/collins/include/mpfr.h \
  ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/function.h \
  ../../include/numeric/../numeric/rational.h

zonotope.o: zonotope.cc ../../include/geometry/zonotope.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../linear_algebra/../numeric/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/approximation.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../../include/geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/rounding.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/function.h \
  ../../include/geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../numeric/interval.h \
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
  ../../include/geometry/zonotope.tpl /home/collins/include/ppl.hh \
  ../../include/geometry/zonotope.h \
  ../../include/geometry/../base/array.h \
  ../../include/geometry/../numeric/approximation.h \
  ../../include/geometry/../combinatoric/lattice_set.h \
  ../../include/geometry/../combinatoric/../base/array.h \
  ../../include/geometry/../combinatoric/../base/iterator.h \
  ../../include/geometry/../combinatoric/../base/../base/array.h \
  ../../include/geometry/../combinatoric/../base/array_operations.h \
  ../../include/geometry/../combinatoric/../base/../declarations.h \
  ../../include/geometry/../combinatoric/../numeric/interval.h \
  ../../include/geometry/../geometry/point_list.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../numeric/function.h \
  ../../include/geometry/../geometry/../numeric/interval.h \
  ../../include/geometry/../geometry/../geometry/ppl_polyhedron.h \
  ../../include/geometry/../geometry/../geometry/rectangle_expression.h \
  ../../include/geometry/../geometry/../geometry/../declarations.h \
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
  ../../include/geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/lu_matrix.h \
  /home/collins/include/tblas/geset.hpp \
  /home/collins/include/tblas/tblas.hpp \
  /home/collins/include/tblas/fblas.hpp \
  /home/collins/include/tblas/amax.hpp \
  /home/collins/include/tblas/asum.hpp \
  /home/collins/include/tblas/axpy.hpp \
  /home/collins/include/tblas/copy.hpp \
  /home/collins/include/tblas/dot.hpp \
  /home/collins/include/tblas/gbmv.hpp \
  /home/collins/include/tblas/gecpy.hpp \
  /home/collins/include/tblas/gemm.hpp \
  /home/collins/include/tblas/gemv.hpp \
  /home/collins/include/tblas/ger.hpp \
  /home/collins/include/tblas/geset.hpp \
  /home/collins/include/tblas/iamax.hpp \
  /home/collins/include/tblas/nrm2.hpp \
  /home/collins/include/tblas/output.hpp \
  /home/collins/include/tblas/rot.hpp \
  /home/collins/include/tblas/rotg.hpp \
  /home/collins/include/tblas/rotm.hpp \
  /home/collins/include/tblas/rotmg.hpp \
  /home/collins/include/tblas/sbmv.hpp \
  /home/collins/include/tblas/scal.hpp \
  /home/collins/include/tblas/set.hpp \
  /home/collins/include/tblas/setm.hpp \
  /home/collins/include/tblas/setv.hpp \
  /home/collins/include/tblas/spmv.hpp \
  /home/collins/include/tblas/spr.hpp \
  /home/collins/include/tblas/spr2.hpp \
  /home/collins/include/tblas/swap.hpp \
  /home/collins/include/tblas/symm.hpp \
  /home/collins/include/tblas/symv.hpp \
  /home/collins/include/tblas/syr.hpp \
  /home/collins/include/tblas/syr2.hpp \
  /home/collins/include/tblas/syr2k.hpp \
  /home/collins/include/tblas/syrk.hpp \
  /home/collins/include/tblas/tbmv.hpp \
  /home/collins/include/tblas/tbsv.hpp \
  /home/collins/include/tblas/tpmv.hpp \
  /home/collins/include/tblas/tpsv.hpp \
  /home/collins/include/tblas/trcpy.hpp \
  /home/collins/include/tblas/trmm.hpp \
  /home/collins/include/tblas/trmv.hpp \
  /home/collins/include/tblas/trset.hpp \
  /home/collins/include/tblas/trsm.hpp \
  /home/collins/include/tblas/trsv.hpp \
  /home/collins/include/tblas/zdotc.hpp \
  /home/collins/include/tblas/zdotu.hpp \
  /home/collins/include/tblas/zgerc.hpp \
  /home/collins/include/tblas/zgeru.hpp \
  /home/collins/include/tblas/zhbmv.hpp \
  /home/collins/include/tblas/zhemm.hpp \
  /home/collins/include/tblas/zhemv.hpp \
  /home/collins/include/tblas/zher.hpp \
  /home/collins/include/tblas/zher2.hpp \
  /home/collins/include/tblas/zher2k.hpp \
  /home/collins/include/tblas/zherk.hpp \
  /home/collins/include/tblas/zhpmv.hpp \
  /home/collins/include/tblas/zhpr.hpp \
  /home/collins/include/tblas/zhpr2.hpp \
  /home/collins/include/tblas/trcpy.hpp \
  /home/collins/include/tblas/gemm.hpp \
  /home/collins/include/tblas/trmm.hpp \
  /home/collins/include/tlapack/laswp.hpp \
  /home/collins/include/tlapack/tlapack.hpp \
  /home/collins/include/tblas/tblas.hpp \
  /home/collins/include/tlapack/geqrf.hpp \
  /home/collins/include/tlapack/larf.hpp \
  /home/collins/include/tblas/ger.hpp \
  /home/collins/include/tlapack/larfg.hpp \
  /home/collins/include/tblas/nrm2.hpp \
  /home/collins/include/tblas/scal.hpp \
  /home/collins/include/tlapack/gesv.hpp \
  /home/collins/include/tlapack/getrf.hpp \
  /home/collins/include/tblas/iamax.hpp \
  /home/collins/include/tblas/swap.hpp \
  /home/collins/include/tlapack/getrs.hpp \
  /home/collins/include/tlapack/gesvd.hpp \
  /home/collins/include/tlapack/laset.hpp \
  /home/collins/include/tlapack/laswp.hpp \
  /home/collins/include/tlapack/orgqr.hpp \
  /home/collins/include/tlapack/getrf.hpp \
  /home/collins/include/tlapack/getrs.hpp \
  ../../include/geometry/../linear_algebra/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../base/array.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/real_typedef.h ../../include/numeric/numerical_types.h \
  ../../include/numeric/mpfloat.h /home/collins/include/mpfr.h \
  ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/function.h \
  ../../include/numeric/../numeric/rational.h

parallelotope.o: parallelotope.cc ../../include/geometry/parallelotope.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../numeric/interval.h \
  ../../include/geometry/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../numeric/../numeric/approximation.h \
  ../../include/geometry/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../numeric/../numeric/rounding.h \
  ../../include/geometry/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../numeric/../numeric/../numeric/function.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../linear_algebra/../numeric/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../../include/geometry/../linear_algebra/../numeric/interval.h \
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
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/parallelotope.tpl \
  ../../include/geometry/parallelotope.h \
  ../../include/geometry/../base/stlio.h \
  ../../include/geometry/../base/../base/array.h \
  ../../include/geometry/../linear_algebra/linear_program.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../combinatoric/lattice_set.h \
  ../../include/geometry/../combinatoric/../base/array.h \
  ../../include/geometry/../combinatoric/../base/iterator.h \
  ../../include/geometry/../combinatoric/../base/../base/array.h \
  ../../include/geometry/../combinatoric/../base/array_operations.h \
  ../../include/geometry/../combinatoric/../base/../declarations.h \
  ../../include/geometry/../combinatoric/../numeric/interval.h \
  ../../include/geometry/../geometry/point_list.h \
  ../../include/geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../numeric/function.h \
  ../../include/geometry/../geometry/../geometry/rectangle_expression.h \
  ../../include/geometry/../geometry/list_set.h \
  ../../include/geometry/../geometry/../base/utility.h \
  ../../include/geometry/../geometry/polyhedron.h \
  ../../include/geometry/../geometry/../geometry/point_list.h \
  ../../include/geometry/../geometry/../geometry/rectangle.h \
  ../../include/real_typedef.h ../../include/numeric/numerical_types.h \
  ../../include/numeric/mpfloat.h /home/collins/include/mpfr.h \
  ../../include/numeric/../declarations.h \
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
  ../../include/geometry/../geometry/../linear_algebra/../numeric/../numeric/approximation.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/../numeric/rounding.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/../numeric/../numeric/function.h \
  ../../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/geometry/sphere.tpl ../../include/geometry/sphere.h \
  ../../include/real_typedef.h ../../include/numeric/numerical_types.h \
  ../../include/numeric/mpfloat.h /home/collins/include/mpfr.h \
  ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/function.h \
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
  ../../include/geometry/../linear_algebra/../numeric/../numeric/approximation.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/rounding.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/function.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/geometry/ellipsoid.tpl ../../include/geometry/ellipsoid.h \
  ../../include/geometry/../geometry/sphere.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/real_typedef.h ../../include/numeric/numerical_types.h \
  ../../include/numeric/mpfloat.h /home/collins/include/mpfr.h \
  ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/function.h \
  ../../include/numeric/../numeric/rational.h \
  ../../include/numeric/../numeric/../declarations.h \
  ../../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/../numeric/arithmetic.h \
  ../../include/numeric/../numeric/../numeric/function.h \
  ../../include/numeric/../numeric/../numeric/integer.h

list_set.o: list_set.cc ../../include/geometry/rectangle.tpl \
  /home/collins/include/ppl.hh ../../include/geometry/rectangle.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../numeric/arithmetic.h \
  ../../include/geometry/../numeric/function.h \
  ../../include/geometry/../numeric/interval.h \
  ../../include/geometry/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../numeric/../numeric/approximation.h \
  ../../include/geometry/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../numeric/../numeric/rounding.h \
  ../../include/geometry/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../numeric/../numeric/../numeric/function.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../linear_algebra/../numeric/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../../include/geometry/../linear_algebra/../numeric/interval.h \
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
  ../../include/geometry/../geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../base/../declarations.h \
  ../../include/geometry/../geometry/rectangle_expression.h \
  ../../include/geometry/../combinatoric/binary_word.h \
  ../../include/geometry/../combinatoric/../declarations.h \
  ../../include/geometry/../base/array.h \
  ../../include/geometry/../geometry/point_list.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/list_set.h \
  ../../include/geometry/../geometry/../base/utility.h \
  ../../include/geometry/../geometry/grid_set.h \
  ../../include/geometry/../geometry/../combinatoric/lattice_set.h \
  ../../include/geometry/../geometry/../combinatoric/../base/array.h \
  ../../include/geometry/../geometry/../combinatoric/../base/iterator.h \
  ../../include/geometry/../geometry/../combinatoric/../base/../base/array.h \
  ../../include/geometry/../geometry/../combinatoric/../base/array_operations.h \
  ../../include/geometry/../geometry/../combinatoric/../base/../declarations.h \
  ../../include/geometry/../geometry/../combinatoric/../numeric/interval.h \
  ../../include/geometry/../geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/../geometry/grid.h \
  ../../include/geometry/../geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../geometry/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../geometry/../combinatoric/lattice_set.h \
  ../../include/geometry/parallelotope.tpl \
  ../../include/geometry/parallelotope.h \
  ../../include/geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/zonotope.h \
  ../../include/geometry/../geometry/../numeric/interval.h \
  ../../include/geometry/../geometry/../geometry/ppl_polyhedron.h \
  ../../include/geometry/../base/stlio.h \
  ../../include/geometry/../base/../base/array.h \
  ../../include/geometry/../linear_algebra/linear_program.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/../combinatoric/lattice_set.h \
  ../../include/geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/polyhedron.h \
  ../../include/geometry/../geometry/../geometry/point_list.h \
  ../../include/geometry/zonotope.tpl ../../include/geometry/zonotope.h \
  ../../include/geometry/../numeric/approximation.h \
  ../../include/geometry/../geometry/parallelotope.h \
  ../../include/geometry/../linear_algebra/matrix.tpl \
  ../../include/geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../linear_algebra/../base/array.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/lu_matrix.h \
  /home/collins/include/tblas/geset.hpp \
  /home/collins/include/tblas/tblas.hpp \
  /home/collins/include/tblas/fblas.hpp \
  /home/collins/include/tblas/amax.hpp \
  /home/collins/include/tblas/asum.hpp \
  /home/collins/include/tblas/axpy.hpp \
  /home/collins/include/tblas/copy.hpp \
  /home/collins/include/tblas/dot.hpp \
  /home/collins/include/tblas/gbmv.hpp \
  /home/collins/include/tblas/gecpy.hpp \
  /home/collins/include/tblas/gemm.hpp \
  /home/collins/include/tblas/gemv.hpp \
  /home/collins/include/tblas/ger.hpp \
  /home/collins/include/tblas/geset.hpp \
  /home/collins/include/tblas/iamax.hpp \
  /home/collins/include/tblas/nrm2.hpp \
  /home/collins/include/tblas/output.hpp \
  /home/collins/include/tblas/rot.hpp \
  /home/collins/include/tblas/rotg.hpp \
  /home/collins/include/tblas/rotm.hpp \
  /home/collins/include/tblas/rotmg.hpp \
  /home/collins/include/tblas/sbmv.hpp \
  /home/collins/include/tblas/scal.hpp \
  /home/collins/include/tblas/set.hpp \
  /home/collins/include/tblas/setm.hpp \
  /home/collins/include/tblas/setv.hpp \
  /home/collins/include/tblas/spmv.hpp \
  /home/collins/include/tblas/spr.hpp \
  /home/collins/include/tblas/spr2.hpp \
  /home/collins/include/tblas/swap.hpp \
  /home/collins/include/tblas/symm.hpp \
  /home/collins/include/tblas/symv.hpp \
  /home/collins/include/tblas/syr.hpp \
  /home/collins/include/tblas/syr2.hpp \
  /home/collins/include/tblas/syr2k.hpp \
  /home/collins/include/tblas/syrk.hpp \
  /home/collins/include/tblas/tbmv.hpp \
  /home/collins/include/tblas/tbsv.hpp \
  /home/collins/include/tblas/tpmv.hpp \
  /home/collins/include/tblas/tpsv.hpp \
  /home/collins/include/tblas/trcpy.hpp \
  /home/collins/include/tblas/trmm.hpp \
  /home/collins/include/tblas/trmv.hpp \
  /home/collins/include/tblas/trset.hpp \
  /home/collins/include/tblas/trsm.hpp \
  /home/collins/include/tblas/trsv.hpp \
  /home/collins/include/tblas/zdotc.hpp \
  /home/collins/include/tblas/zdotu.hpp \
  /home/collins/include/tblas/zgerc.hpp \
  /home/collins/include/tblas/zgeru.hpp \
  /home/collins/include/tblas/zhbmv.hpp \
  /home/collins/include/tblas/zhemm.hpp \
  /home/collins/include/tblas/zhemv.hpp \
  /home/collins/include/tblas/zher.hpp \
  /home/collins/include/tblas/zher2.hpp \
  /home/collins/include/tblas/zher2k.hpp \
  /home/collins/include/tblas/zherk.hpp \
  /home/collins/include/tblas/zhpmv.hpp \
  /home/collins/include/tblas/zhpr.hpp \
  /home/collins/include/tblas/zhpr2.hpp \
  /home/collins/include/tblas/trcpy.hpp \
  /home/collins/include/tblas/gemm.hpp \
  /home/collins/include/tblas/trmm.hpp \
  /home/collins/include/tlapack/laswp.hpp \
  /home/collins/include/tlapack/tlapack.hpp \
  /home/collins/include/tblas/tblas.hpp \
  /home/collins/include/tlapack/geqrf.hpp \
  /home/collins/include/tlapack/larf.hpp \
  /home/collins/include/tblas/ger.hpp \
  /home/collins/include/tlapack/larfg.hpp \
  /home/collins/include/tblas/nrm2.hpp \
  /home/collins/include/tblas/scal.hpp \
  /home/collins/include/tlapack/gesv.hpp \
  /home/collins/include/tlapack/getrf.hpp \
  /home/collins/include/tblas/iamax.hpp \
  /home/collins/include/tblas/swap.hpp \
  /home/collins/include/tlapack/getrs.hpp \
  /home/collins/include/tlapack/gesvd.hpp \
  /home/collins/include/tlapack/laset.hpp \
  /home/collins/include/tlapack/laswp.hpp \
  /home/collins/include/tlapack/orgqr.hpp \
  /home/collins/include/tlapack/getrf.hpp \
  /home/collins/include/tlapack/getrs.hpp \
  ../../include/geometry/../linear_algebra/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../base/array.h \
  ../../include/geometry/../linear_algebra/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/geometry/list_set.h ../../include/geometry/list_set.tpl \
  ../../include/geometry/list_set.h ../../include/real_typedef.h \
  ../../include/numeric/numerical_types.h ../../include/numeric/mpfloat.h \
  /home/collins/include/mpfr.h ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/function.h \
  ../../include/numeric/../numeric/rational.h

grid.o: grid.cc ../../include/geometry/grid.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../base/array.h \
  ../../include/geometry/../base/../declarations.h \
  ../../include/geometry/../numeric/arithmetic.h \
  ../../include/geometry/../combinatoric/lattice_set.h \
  ../../include/geometry/../combinatoric/../base/array.h \
  ../../include/geometry/../combinatoric/../base/iterator.h \
  ../../include/geometry/../combinatoric/../base/../base/array.h \
  ../../include/geometry/../combinatoric/../base/array_operations.h \
  ../../include/geometry/../combinatoric/../base/../declarations.h \
  ../../include/geometry/../combinatoric/../numeric/interval.h \
  ../../include/geometry/../combinatoric/../numeric/../declarations.h \
  ../../include/geometry/../combinatoric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../combinatoric/../numeric/../numeric/approximation.h \
  ../../include/geometry/../combinatoric/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../combinatoric/../numeric/../numeric/rounding.h \
  ../../include/geometry/../combinatoric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../combinatoric/../numeric/../numeric/../numeric/function.h \
  ../../include/geometry/grid.tpl ../../include/geometry/grid.h \
  ../../include/geometry/../base/stlio.h \
  ../../include/geometry/../base/../base/array.h \
  ../../include/geometry/../geometry/grid_set.h \
  ../../include/geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../combinatoric/lattice_set.h \
  ../../include/geometry/../geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../geometry/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../geometry/../numeric/function.h \
  ../../include/geometry/../geometry/../geometry/../numeric/interval.h \
  ../../include/geometry/../geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../geometry/../geometry/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../geometry/../geometry/../linear_algebra/../numeric/../declarations.h \
  ../../include/geometry/../geometry/../geometry/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../geometry/../geometry/../linear_algebra/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../geometry/../linear_algebra/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../geometry/../linear_algebra/../numeric/interval.h \
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
  ../../include/geometry/../geometry/../geometry/../geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../geometry/../geometry/rectangle_expression.h \
  ../../include/geometry/../geometry/../geometry/grid.h \
  ../../include/geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/list_set.h \
  ../../include/geometry/../geometry/../base/utility.h \
  ../../include/geometry/../geometry/../base/../declarations.h \
  ../../include/real_typedef.h ../../include/numeric/numerical_types.h \
  ../../include/numeric/mpfloat.h /home/collins/include/mpfr.h \
  ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/function.h \
  ../../include/numeric/../numeric/rational.h

grid_set.o: grid_set.cc ../../include/geometry/grid_set.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../base/array.h \
  ../../include/geometry/../base/../declarations.h \
  ../../include/geometry/../combinatoric/lattice_set.h \
  ../../include/geometry/../combinatoric/../base/array.h \
  ../../include/geometry/../combinatoric/../base/iterator.h \
  ../../include/geometry/../combinatoric/../base/../base/array.h \
  ../../include/geometry/../combinatoric/../base/array_operations.h \
  ../../include/geometry/../combinatoric/../base/../declarations.h \
  ../../include/geometry/../combinatoric/../numeric/interval.h \
  ../../include/geometry/../combinatoric/../numeric/../declarations.h \
  ../../include/geometry/../combinatoric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../combinatoric/../numeric/../numeric/approximation.h \
  ../../include/geometry/../combinatoric/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../combinatoric/../numeric/../numeric/rounding.h \
  ../../include/geometry/../combinatoric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../combinatoric/../numeric/../numeric/../numeric/function.h \
  ../../include/geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/../declarations.h \
  ../../include/geometry/../geometry/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../numeric/function.h \
  ../../include/geometry/../geometry/../numeric/interval.h \
  ../../include/geometry/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/../geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/../declarations.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../linear_algebra/../numeric/interval.h \
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
  ../../include/geometry/../geometry/../geometry/../linear_algebra/../numeric/interval.h \
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../geometry/rectangle_expression.h \
  ../../include/geometry/../geometry/grid.h \
  ../../include/geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../combinatoric/lattice_set.h \
  ../../include/geometry/grid_set.tpl ../../include/geometry/grid_set.h \
  ../../include/geometry/../base/stlio.h \
  ../../include/geometry/../base/../base/array.h \
  ../../include/geometry/../base/array_operations.h \
  ../../include/geometry/../geometry/zonotope.h \
  ../../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../../include/geometry/../geometry/parallelotope.h \
  ../../include/geometry/../geometry/../geometry/zonotope.h \
  ../../include/geometry/../geometry/polytope.h \
  ../../include/geometry/../geometry/../geometry/point_list.h \
  ../../include/geometry/../geometry/../geometry/../geometry/point.h \
  ../../include/geometry/../geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/list_set.h \
  ../../include/geometry/../geometry/../base/utility.h \
  ../../include/geometry/../geometry/../base/../declarations.h \
  ../../include/geometry/../geometry/partition_tree_set.h \
  ../../include/geometry/../geometry/../base/iterator.h \
  ../../include/geometry/../geometry/../combinatoric/binary_word.h \
  ../../include/geometry/../geometry/../combinatoric/../declarations.h \
  ../../include/geometry/../geometry/../combinatoric/binary_tree.h \
  ../../include/geometry/../geometry/../combinatoric/../base/array.h \
  ../../include/geometry/../geometry/../combinatoric/../combinatoric/binary_word.h \
  ../../include/geometry/../geometry/../combinatoric/subdivision_tree_set.h \
  ../../include/geometry/../geometry/../combinatoric/../base/sequence.h \
  ../../include/geometry/../geometry/../combinatoric/../base/iterator.h \
  ../../include/geometry/../geometry/../combinatoric/../combinatoric/binary_tree.h \
  ../../include/geometry/../geometry/../combinatoric/../geometry/rectangle.h \
  ../../include/real_typedef.h ../../include/numeric/numerical_types.h \
  ../../include/numeric/mpfloat.h /home/collins/include/mpfr.h \
  ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/function.h \
  ../../include/numeric/../numeric/rational.h

partition_tree_set.o: partition_tree_set.cc \
  ../../include/geometry/parallelotope.h \
  ../../include/geometry/../declarations.h \
  ../../include/geometry/../numeric/interval.h \
  ../../include/geometry/../numeric/../declarations.h \
  ../../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../numeric/../numeric/approximation.h \
  ../../include/geometry/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../numeric/../numeric/rounding.h \
  ../../include/geometry/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../numeric/../numeric/../numeric/function.h \
  ../../include/geometry/../linear_algebra/vector.h \
  ../../include/geometry/../linear_algebra/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/integer.h \
  ../../include/geometry/../linear_algebra/../numeric/../declarations.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/geometry/../linear_algebra/../numeric/../numeric/arithmetic.h \
  ../../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../../include/geometry/../linear_algebra/../numeric/interval.h \
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
  ../../include/geometry/../geometry/../geometry/point.h \
  ../../include/geometry/partition_tree_set.h \
  ../../include/geometry/../base/iterator.h \
  ../../include/geometry/../base/../base/array.h \
  ../../include/geometry/../combinatoric/binary_word.h \
  ../../include/geometry/../combinatoric/../declarations.h \
  ../../include/geometry/../combinatoric/binary_tree.h \
  ../../include/geometry/../combinatoric/../base/array.h \
  ../../include/geometry/../combinatoric/../combinatoric/binary_word.h \
  ../../include/geometry/../combinatoric/subdivision_tree_set.h \
  ../../include/geometry/../combinatoric/../base/sequence.h \
  ../../include/geometry/../combinatoric/../base/iterator.h \
  ../../include/geometry/../combinatoric/../combinatoric/binary_tree.h \
  ../../include/geometry/../combinatoric/../geometry/rectangle.h \
  ../../include/geometry/../combinatoric/../geometry/../declarations.h \
  ../../include/geometry/../combinatoric/../geometry/../numeric/arithmetic.h \
  ../../include/geometry/../combinatoric/../geometry/../numeric/function.h \
  ../../include/geometry/../combinatoric/../geometry/../numeric/interval.h \
  ../../include/geometry/../combinatoric/../geometry/../linear_algebra/vector.h \
  ../../include/geometry/../combinatoric/../geometry/../geometry/ppl_polyhedron.h \
  ../../include/geometry/../combinatoric/../geometry/../geometry/point.h \
  ../../include/geometry/../combinatoric/../geometry/../geometry/rectangle_expression.h \
  ../../include/geometry/../combinatoric/../geometry/../geometry/../declarations.h \
  ../../include/geometry/partition_tree_set.tpl \
  ../../include/geometry/partition_tree_set.h \
  ../../include/geometry/../base/stlio.h \
  ../../include/geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/list_set.h \
  ../../include/geometry/../geometry/../base/utility.h \
  ../../include/geometry/../geometry/grid_set.h \
  ../../include/geometry/../geometry/../combinatoric/lattice_set.h \
  ../../include/geometry/../geometry/../combinatoric/../base/array.h \
  ../../include/geometry/../geometry/../combinatoric/../base/iterator.h \
  ../../include/geometry/../geometry/../combinatoric/../base/array_operations.h \
  ../../include/geometry/../geometry/../combinatoric/../base/../declarations.h \
  ../../include/geometry/../geometry/../combinatoric/../numeric/interval.h \
  ../../include/geometry/../geometry/../geometry/rectangle.h \
  ../../include/geometry/../geometry/../geometry/grid.h \
  ../../include/geometry/../geometry/../geometry/../base/array.h \
  ../../include/geometry/../geometry/../geometry/../numeric/arithmetic.h \
  ../../include/geometry/../geometry/../geometry/../combinatoric/lattice_set.h \
  ../../include/real_typedef.h ../../include/numeric/numerical_types.h \
  ../../include/numeric/mpfloat.h /home/collins/include/mpfr.h \
  ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/function.h \
  ../../include/numeric/../numeric/rational.h

