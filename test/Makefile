include ../config.mk

TESTS = \
  test_array test_binary_word \
  test_lattice \
  test_arithmetic test_interval \
  test_vector test_matrix test_linear_algebra \
  test_point \
  test_rectangle test_parallelotope test_zonotope \
  test_simplex test_polyhedron \
  test_list_set test_grid_set test_partition_tree_set \
  test_epsfstream \
  test_interval_newton \
  test_chainreach \
  test_integration_step 

CXXFLAGS += -Wall -I../include/ 
LDXXFLAGS += -L../src/ -lariadne -lgmp -lgmpxx -lppl #-lblas

all: ${TESTS}

.cc.o:
	${CXX} ${INCLUDE} ${CXXFLAGS} -c $<

clean:
	-rm -rf core *.log *.eps *.o ${TESTS} 

dep:  
	@sed -e '/^### Do Not edit this line$$/,$$d' Makefile \
		> MMM.$$$$ && mv MMM.$$$$ Makefile
	@echo "### Do Not edit this line" >> Makefile
	@echo "### Everything below is auto-generated" >> Makefile
	@echo "" >> Makefile 
	@echo "all: $$""{TESTS}" >> Makefile
	@for ff in ${TESTS}; do echo "	-./$$ff"; done >> Makefile
	@echo "" >> Makefile 
	@for ff in ${TESTS}; do g++ -MM ${CXXFLAGS} $$ff.cc; \
	echo "	-$$""{CXX} $$""{CXXFLAGS} -c $$ff.cc "; \
	echo " "; \
	echo "$$ff: $$ff.o"; \
	echo "	-$$""{CXX} $$""{LDXXFLAGS} $$ff.o -o $$ff "; \
	echo ""; done >> Makefile


### Do Not edit this line
### Everything below is auto-generated

all: ${TESTS}
	-./test_array
	-./test_binary_word
	-./test_lattice
	-./test_arithmetic
	-./test_interval
	-./test_vector
	-./test_matrix
	-./test_linear_algebra
	-./test_point
	-./test_rectangle
	-./test_parallelotope
	-./test_zonotope
	-./test_simplex
	-./test_polyhedron
	-./test_list_set
	-./test_grid_set
	-./test_partition_tree_set
	-./test_epsfstream
	-./test_interval_newton
	-./test_chainreach
	-./test_integration_step

test_array.o: test_array.cc ../include/ariadne.h ../include/base/array.h \
  ../include/base/../declarations.h ../include/numeric/numerical_types.h \
  ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../utility/stlio.h \
  ../include/numeric/../numeric/../utility/../base/array.h \
  ../include/numeric/../numeric/dyadic.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../utility/stlio.h test.h
	-${CXX} ${CXXFLAGS} -c test_array.cc 
 
test_array: test_array.o
	-${CXX} ${LDXXFLAGS} test_array.o -o test_array 

test_binary_word.o: test_binary_word.cc ../include/ariadne.h \
  ../include/base/exception.h ../include/base/binary_word.h \
  ../include/base/../declarations.h ../include/base/binary_tree.h \
  ../include/base/../base/array.h \
  ../include/base/../base/../declarations.h \
  ../include/base/../base/binary_word.h \
  ../include/numeric/numerical_types.h \
  ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../utility/stlio.h \
  ../include/numeric/../numeric/../utility/../base/array.h \
  ../include/numeric/../numeric/dyadic.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../utility/stlio.h ../include/base/utility.h \
  ../include/geometry/point.h ../include/geometry/../declarations.h \
  ../include/geometry/../base/array.h \
  ../include/geometry/../linear_algebra/vector.h \
  ../include/geometry/../linear_algebra/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/integer.h \
  ../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/numerical_types.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/approximation.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/numerical_types.h \
  ../include/geometry/rectangle.h \
  ../include/geometry/../numeric/arithmetic.h \
  ../include/geometry/../numeric/interval.h \
  ../include/geometry/../numeric/../declarations.h \
  ../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../numeric/../numeric/arithmetic.h \
  ../include/geometry/../linear_algebra/interval_vector.h \
  ../include/geometry/../linear_algebra/../numeric/interval.h \
  ../include/geometry/../linear_algebra/../linear_algebra/vector.h \
  ../include/geometry/../geometry/point.h test.h
	-${CXX} ${CXXFLAGS} -c test_binary_word.cc 
 
test_binary_word: test_binary_word.o
	-${CXX} ${LDXXFLAGS} test_binary_word.o -o test_binary_word 

test_lattice.o: test_lattice.cc ../include/declarations.h \
  ../include/real_typedef.h ../include/numeric/numerical_types.h \
  ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../utility/stlio.h \
  ../include/numeric/../numeric/../utility/../base/array.h \
  ../include/numeric/../numeric/../utility/../base/../declarations.h \
  ../include/numeric/../numeric/dyadic.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../utility/stlio.h ../include/geometry/lattice_set.h \
  ../include/geometry/../base/array.h \
  ../include/geometry/../base/iterator.h \
  ../include/geometry/../base/../base/array.h \
  ../include/geometry/../numeric/interval.h \
  ../include/geometry/../numeric/../declarations.h \
  ../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../numeric/../numeric/arithmetic.h \
  ../include/geometry/../numeric/../numeric/../numeric/numerical_types.h \
  ../include/geometry/../numeric/../numeric/../numeric/approximation.h \
  ../include/geometry/../numeric/../numeric/../numeric/../declarations.h \
  ../include/geometry/../numeric/../numeric/../numeric/../numeric/numerical_types.h \
  ../include/geometry/../base/array_operations.h \
  ../include/geometry/../base/../declarations.h \
  ../include/system/lattice_map.h \
  ../include/system/../geometry/lattice_set.h
	-${CXX} ${CXXFLAGS} -c test_lattice.cc 
 
test_lattice: test_lattice.o
	-${CXX} ${LDXXFLAGS} test_lattice.o -o test_lattice 

test_arithmetic.o: test_arithmetic.cc ../include/numeric/arithmetic.h \
  ../include/numeric/../numeric/numerical_types.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/../utility/stlio.h \
  ../include/numeric/../numeric/../numeric/../utility/../base/array.h \
  ../include/numeric/../numeric/../numeric/../utility/../base/../declarations.h \
  ../include/numeric/../numeric/../numeric/dyadic.h \
  ../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../utility/stlio.h \
  ../include/numeric/../numeric/approximation.h \
  ../include/numeric/../numeric/../numeric/numerical_types.h test.h
	-${CXX} ${CXXFLAGS} -c test_arithmetic.cc 
 
test_arithmetic: test_arithmetic.o
	-${CXX} ${LDXXFLAGS} test_arithmetic.o -o test_arithmetic 

test_interval.o: test_interval.cc ../include/ariadne.h \
  ../include/real_typedef.h ../include/numeric/numerical_types.h \
  ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../utility/stlio.h \
  ../include/numeric/../numeric/../utility/../base/array.h \
  ../include/numeric/../numeric/../utility/../base/../declarations.h \
  ../include/numeric/../numeric/dyadic.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../utility/stlio.h \
  ../include/numeric/numerical_types.h ../include/numeric/interval.h \
  ../include/numeric/../numeric/arithmetic.h \
  ../include/numeric/../numeric/../numeric/numerical_types.h \
  ../include/numeric/../numeric/../numeric/approximation.h \
  ../include/numeric/../numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/../numeric/numerical_types.h \
  test.h
	-${CXX} ${CXXFLAGS} -c test_interval.cc 
 
test_interval: test_interval.o
	-${CXX} ${LDXXFLAGS} test_interval.o -o test_interval 

test_vector.o: test_vector.cc ../include/real_typedef.h \
  ../include/numeric/numerical_types.h \
  ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../utility/stlio.h \
  ../include/numeric/../numeric/../utility/../base/array.h \
  ../include/numeric/../numeric/../utility/../base/../declarations.h \
  ../include/numeric/../numeric/dyadic.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../utility/stlio.h \
  ../include/numeric/numerical_types.h ../include/linear_algebra/vector.h \
  ../include/linear_algebra/../declarations.h \
  ../include/linear_algebra/../numeric/integer.h \
  ../include/linear_algebra/../numeric/arithmetic.h \
  ../include/linear_algebra/../numeric/../numeric/numerical_types.h \
  ../include/linear_algebra/../numeric/../numeric/approximation.h \
  ../include/linear_algebra/../numeric/../numeric/../declarations.h \
  ../include/linear_algebra/../numeric/../numeric/../numeric/numerical_types.h \
  ../include/linear_algebra/interval_vector.h \
  ../include/linear_algebra/../numeric/interval.h \
  ../include/linear_algebra/../numeric/../declarations.h \
  ../include/linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../include/linear_algebra/../numeric/../numeric/arithmetic.h \
  ../include/linear_algebra/../linear_algebra/vector.h \
  ../include/linear_algebra/vector.tpl ../include/linear_algebra/vector.h \
  ../include/linear_algebra/../utility/stlio.h
	-${CXX} ${CXXFLAGS} -c test_vector.cc 
 
test_vector: test_vector.o
	-${CXX} ${LDXXFLAGS} test_vector.o -o test_vector 

test_matrix.o: test_matrix.cc ../include/declarations.h \
  ../include/real_typedef.h ../include/numeric/numerical_types.h \
  ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../utility/stlio.h \
  ../include/numeric/../numeric/../utility/../base/array.h \
  ../include/numeric/../numeric/../utility/../base/../declarations.h \
  ../include/numeric/../numeric/dyadic.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../utility/stlio.h \
  ../include/numeric/numerical_types.h ../include/linear_algebra/vector.h \
  ../include/linear_algebra/../declarations.h \
  ../include/linear_algebra/../numeric/integer.h \
  ../include/linear_algebra/../numeric/arithmetic.h \
  ../include/linear_algebra/../numeric/../numeric/numerical_types.h \
  ../include/linear_algebra/../numeric/../numeric/approximation.h \
  ../include/linear_algebra/../numeric/../numeric/../declarations.h \
  ../include/linear_algebra/../numeric/../numeric/../numeric/numerical_types.h \
  ../include/linear_algebra/vector.tpl ../include/linear_algebra/vector.h \
  ../include/linear_algebra/../utility/stlio.h \
  ../include/linear_algebra/matrix.h ../include/linear_algebra/matrix.tpl \
  ../include/linear_algebra/matrix.h \
  ../include/linear_algebra/../base/array.h \
  ../include/linear_algebra/../numeric/function.h \
  /ufs/collins/include/mpfr.h \
  ../include/linear_algebra/../linear_algebra/vector.h \
  ../include/linear_algebra/../linear_algebra/matrix.h \
  ../include/linear_algebra/../linear_algebra/lu_matrix.h \
  /ufs/collins/include/tblas/geset.hpp \
  /ufs/collins/include/tblas/tblas.hpp \
  /ufs/collins/include/tblas/fblas.hpp \
  /ufs/collins/include/tblas/amax.hpp /ufs/collins/include/tblas/asum.hpp \
  /ufs/collins/include/tblas/axpy.hpp /ufs/collins/include/tblas/copy.hpp \
  /ufs/collins/include/tblas/dot.hpp /ufs/collins/include/tblas/gbmv.hpp \
  /ufs/collins/include/tblas/gecpy.hpp \
  /ufs/collins/include/tblas/gemm.hpp /ufs/collins/include/tblas/gemv.hpp \
  /ufs/collins/include/tblas/ger.hpp /ufs/collins/include/tblas/geset.hpp \
  /ufs/collins/include/tblas/iamax.hpp \
  /ufs/collins/include/tblas/nrm2.hpp \
  /ufs/collins/include/tblas/output.hpp \
  /ufs/collins/include/tblas/rot.hpp /ufs/collins/include/tblas/rotg.hpp \
  /ufs/collins/include/tblas/rotm.hpp \
  /ufs/collins/include/tblas/rotmg.hpp \
  /ufs/collins/include/tblas/sbmv.hpp /ufs/collins/include/tblas/scal.hpp \
  /ufs/collins/include/tblas/set.hpp /ufs/collins/include/tblas/setm.hpp \
  /ufs/collins/include/tblas/setv.hpp /ufs/collins/include/tblas/spmv.hpp \
  /ufs/collins/include/tblas/spr.hpp /ufs/collins/include/tblas/spr2.hpp \
  /ufs/collins/include/tblas/swap.hpp /ufs/collins/include/tblas/symm.hpp \
  /ufs/collins/include/tblas/symv.hpp /ufs/collins/include/tblas/syr.hpp \
  /ufs/collins/include/tblas/syr2.hpp \
  /ufs/collins/include/tblas/syr2k.hpp \
  /ufs/collins/include/tblas/syrk.hpp /ufs/collins/include/tblas/tbmv.hpp \
  /ufs/collins/include/tblas/tbsv.hpp /ufs/collins/include/tblas/tpmv.hpp \
  /ufs/collins/include/tblas/tpsv.hpp \
  /ufs/collins/include/tblas/trcpy.hpp \
  /ufs/collins/include/tblas/trmm.hpp /ufs/collins/include/tblas/trmv.hpp \
  /ufs/collins/include/tblas/trset.hpp \
  /ufs/collins/include/tblas/trsm.hpp /ufs/collins/include/tblas/trsv.hpp \
  /ufs/collins/include/tblas/zdotc.hpp \
  /ufs/collins/include/tblas/zdotu.hpp \
  /ufs/collins/include/tblas/zgerc.hpp \
  /ufs/collins/include/tblas/zgeru.hpp \
  /ufs/collins/include/tblas/zhbmv.hpp \
  /ufs/collins/include/tblas/zhemm.hpp \
  /ufs/collins/include/tblas/zhemv.hpp \
  /ufs/collins/include/tblas/zher.hpp \
  /ufs/collins/include/tblas/zher2.hpp \
  /ufs/collins/include/tblas/zher2k.hpp \
  /ufs/collins/include/tblas/zherk.hpp \
  /ufs/collins/include/tblas/zhpmv.hpp \
  /ufs/collins/include/tblas/zhpr.hpp \
  /ufs/collins/include/tblas/zhpr2.hpp \
  /ufs/collins/include/tblas/trcpy.hpp \
  /ufs/collins/include/tblas/gemm.hpp /ufs/collins/include/tblas/trmm.hpp \
  /ufs/collins/include/tlapack/laswp.hpp \
  /ufs/collins/include/tlapack/tlapack.hpp \
  /ufs/collins/include/tblas/tblas.hpp \
  /ufs/collins/include/tlapack/geqrf.hpp \
  /ufs/collins/include/tlapack/larf.hpp \
  /ufs/collins/include/tblas/ger.hpp \
  /ufs/collins/include/tlapack/larfg.hpp \
  /ufs/collins/include/tblas/nrm2.hpp /ufs/collins/include/tblas/scal.hpp \
  /ufs/collins/include/tlapack/gesv.hpp \
  /ufs/collins/include/tlapack/getrf.hpp \
  /ufs/collins/include/tblas/iamax.hpp \
  /ufs/collins/include/tblas/swap.hpp \
  /ufs/collins/include/tlapack/getrs.hpp \
  /ufs/collins/include/tlapack/gesvd.hpp \
  /ufs/collins/include/tlapack/laset.hpp \
  /ufs/collins/include/tlapack/laswp.hpp \
  /ufs/collins/include/tlapack/orgqr.hpp \
  /ufs/collins/include/tlapack/getrf.hpp \
  /ufs/collins/include/tlapack/getrs.hpp \
  ../include/linear_algebra/../linear_algebra/../declarations.h \
  ../include/linear_algebra/../linear_algebra/../base/array.h \
  ../include/linear_algebra/../linear_algebra/../linear_algebra/matrix.h \
  ../include/linear_algebra/interval_matrix.h \
  ../include/linear_algebra/../numeric/interval.h \
  ../include/linear_algebra/../numeric/../declarations.h \
  ../include/linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../include/linear_algebra/../numeric/../numeric/arithmetic.h \
  ../include/linear_algebra/interval_matrix.tpl ../include/declarations.h \
  ../include/linear_algebra/interval_matrix.h \
  ../include/linear_algebra/../linear_algebra/interval_vector.h \
  ../include/linear_algebra/../linear_algebra/../numeric/interval.h \
  ../include/linear_algebra/../linear_algebra/../linear_algebra/vector.h
	-${CXX} ${CXXFLAGS} -c test_matrix.cc 
 
test_matrix: test_matrix.o
	-${CXX} ${LDXXFLAGS} test_matrix.o -o test_matrix 

test_linear_algebra.o: test_linear_algebra.cc ../include/real_typedef.h \
  ../include/numeric/numerical_types.h \
  ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../utility/stlio.h \
  ../include/numeric/../numeric/../utility/../base/array.h \
  ../include/numeric/../numeric/../utility/../base/../declarations.h \
  ../include/numeric/../numeric/dyadic.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../utility/stlio.h \
  ../include/numeric/numerical_types.h ../include/linear_algebra/vector.h \
  ../include/linear_algebra/../declarations.h \
  ../include/linear_algebra/../numeric/integer.h \
  ../include/linear_algebra/../numeric/arithmetic.h \
  ../include/linear_algebra/../numeric/../numeric/numerical_types.h \
  ../include/linear_algebra/../numeric/../numeric/approximation.h \
  ../include/linear_algebra/../numeric/../numeric/../declarations.h \
  ../include/linear_algebra/../numeric/../numeric/../numeric/numerical_types.h \
  ../include/linear_algebra/matrix.h \
  ../include/linear_algebra/lu_matrix.h \
  /ufs/collins/include/tblas/geset.hpp \
  /ufs/collins/include/tblas/tblas.hpp \
  /ufs/collins/include/tblas/fblas.hpp \
  /ufs/collins/include/tblas/amax.hpp /ufs/collins/include/tblas/asum.hpp \
  /ufs/collins/include/tblas/axpy.hpp /ufs/collins/include/tblas/copy.hpp \
  /ufs/collins/include/tblas/dot.hpp /ufs/collins/include/tblas/gbmv.hpp \
  /ufs/collins/include/tblas/gecpy.hpp \
  /ufs/collins/include/tblas/gemm.hpp /ufs/collins/include/tblas/gemv.hpp \
  /ufs/collins/include/tblas/ger.hpp /ufs/collins/include/tblas/geset.hpp \
  /ufs/collins/include/tblas/iamax.hpp \
  /ufs/collins/include/tblas/nrm2.hpp \
  /ufs/collins/include/tblas/output.hpp \
  /ufs/collins/include/tblas/rot.hpp /ufs/collins/include/tblas/rotg.hpp \
  /ufs/collins/include/tblas/rotm.hpp \
  /ufs/collins/include/tblas/rotmg.hpp \
  /ufs/collins/include/tblas/sbmv.hpp /ufs/collins/include/tblas/scal.hpp \
  /ufs/collins/include/tblas/set.hpp /ufs/collins/include/tblas/setm.hpp \
  /ufs/collins/include/tblas/setv.hpp /ufs/collins/include/tblas/spmv.hpp \
  /ufs/collins/include/tblas/spr.hpp /ufs/collins/include/tblas/spr2.hpp \
  /ufs/collins/include/tblas/swap.hpp /ufs/collins/include/tblas/symm.hpp \
  /ufs/collins/include/tblas/symv.hpp /ufs/collins/include/tblas/syr.hpp \
  /ufs/collins/include/tblas/syr2.hpp \
  /ufs/collins/include/tblas/syr2k.hpp \
  /ufs/collins/include/tblas/syrk.hpp /ufs/collins/include/tblas/tbmv.hpp \
  /ufs/collins/include/tblas/tbsv.hpp /ufs/collins/include/tblas/tpmv.hpp \
  /ufs/collins/include/tblas/tpsv.hpp \
  /ufs/collins/include/tblas/trcpy.hpp \
  /ufs/collins/include/tblas/trmm.hpp /ufs/collins/include/tblas/trmv.hpp \
  /ufs/collins/include/tblas/trset.hpp \
  /ufs/collins/include/tblas/trsm.hpp /ufs/collins/include/tblas/trsv.hpp \
  /ufs/collins/include/tblas/zdotc.hpp \
  /ufs/collins/include/tblas/zdotu.hpp \
  /ufs/collins/include/tblas/zgerc.hpp \
  /ufs/collins/include/tblas/zgeru.hpp \
  /ufs/collins/include/tblas/zhbmv.hpp \
  /ufs/collins/include/tblas/zhemm.hpp \
  /ufs/collins/include/tblas/zhemv.hpp \
  /ufs/collins/include/tblas/zher.hpp \
  /ufs/collins/include/tblas/zher2.hpp \
  /ufs/collins/include/tblas/zher2k.hpp \
  /ufs/collins/include/tblas/zherk.hpp \
  /ufs/collins/include/tblas/zhpmv.hpp \
  /ufs/collins/include/tblas/zhpr.hpp \
  /ufs/collins/include/tblas/zhpr2.hpp \
  /ufs/collins/include/tblas/trcpy.hpp \
  /ufs/collins/include/tblas/gemm.hpp /ufs/collins/include/tblas/trmm.hpp \
  /ufs/collins/include/tlapack/laswp.hpp \
  /ufs/collins/include/tlapack/tlapack.hpp \
  /ufs/collins/include/tblas/tblas.hpp \
  /ufs/collins/include/tlapack/geqrf.hpp \
  /ufs/collins/include/tlapack/larf.hpp \
  /ufs/collins/include/tblas/ger.hpp \
  /ufs/collins/include/tlapack/larfg.hpp \
  /ufs/collins/include/tblas/nrm2.hpp /ufs/collins/include/tblas/scal.hpp \
  /ufs/collins/include/tlapack/gesv.hpp \
  /ufs/collins/include/tlapack/getrf.hpp \
  /ufs/collins/include/tblas/iamax.hpp \
  /ufs/collins/include/tblas/swap.hpp \
  /ufs/collins/include/tlapack/getrs.hpp \
  /ufs/collins/include/tlapack/gesvd.hpp \
  /ufs/collins/include/tlapack/laset.hpp \
  /ufs/collins/include/tlapack/laswp.hpp \
  /ufs/collins/include/tlapack/orgqr.hpp \
  /ufs/collins/include/tlapack/getrf.hpp \
  /ufs/collins/include/tlapack/getrs.hpp \
  ../include/linear_algebra/../base/array.h \
  ../include/linear_algebra/../linear_algebra/matrix.h \
  ../include/linear_algebra/qr_matrix.h \
  /ufs/collins/include/tblas/gecpy.hpp \
  /ufs/collins/include/tblas/trset.hpp \
  /ufs/collins/include/tlapack/larf.hpp \
  /ufs/collins/include/tlapack/larfg.hpp \
  /ufs/collins/include/tlapack/geqrf.hpp \
  /ufs/collins/include/tlapack/orgqr.hpp \
  ../include/linear_algebra/svd_matrix.h \
  /ufs/collins/include/tblas/copy.hpp \
  /ufs/collins/include/tlapack/gesvd.hpp \
  ../include/linear_algebra/../linear_algebra/vector.h test.h
	-${CXX} ${CXXFLAGS} -c test_linear_algebra.cc 
 
test_linear_algebra: test_linear_algebra.o
	-${CXX} ${LDXXFLAGS} test_linear_algebra.o -o test_linear_algebra 

test_point.o: test_point.cc ../include/real_typedef.h \
  ../include/numeric/numerical_types.h \
  ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../utility/stlio.h \
  ../include/numeric/../numeric/../utility/../base/array.h \
  ../include/numeric/../numeric/../utility/../base/../declarations.h \
  ../include/numeric/../numeric/dyadic.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../utility/stlio.h ../include/geometry/point.h \
  ../include/geometry/../declarations.h \
  ../include/geometry/../base/array.h \
  ../include/geometry/../linear_algebra/vector.h \
  ../include/geometry/../linear_algebra/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/integer.h \
  ../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/numerical_types.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/approximation.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/numerical_types.h \
  test.h
	-${CXX} ${CXXFLAGS} -c test_point.cc 
 
test_point: test_point.o
	-${CXX} ${LDXXFLAGS} test_point.o -o test_point 

test_rectangle.o: test_rectangle.cc ../include/ariadne.h \
  ../include/real_typedef.h ../include/numeric/numerical_types.h \
  ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../utility/stlio.h \
  ../include/numeric/../numeric/../utility/../base/array.h \
  ../include/numeric/../numeric/../utility/../base/../declarations.h \
  ../include/numeric/../numeric/dyadic.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../utility/stlio.h ../include/base/exception.h \
  ../include/base/utility.h ../include/base/../declarations.h \
  ../include/geometry/rectangle.h ../include/geometry/../declarations.h \
  ../include/geometry/../numeric/arithmetic.h \
  ../include/geometry/../numeric/../numeric/numerical_types.h \
  ../include/geometry/../numeric/../numeric/approximation.h \
  ../include/geometry/../numeric/../numeric/../declarations.h \
  ../include/geometry/../numeric/../numeric/../numeric/numerical_types.h \
  ../include/geometry/../numeric/interval.h \
  ../include/geometry/../numeric/../declarations.h \
  ../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../numeric/../numeric/arithmetic.h \
  ../include/geometry/../linear_algebra/interval_vector.h \
  ../include/geometry/../linear_algebra/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/interval.h \
  ../include/geometry/../linear_algebra/../linear_algebra/vector.h \
  ../include/geometry/../linear_algebra/../linear_algebra/../declarations.h \
  ../include/geometry/../linear_algebra/../linear_algebra/../numeric/integer.h \
  ../include/geometry/../linear_algebra/../linear_algebra/../numeric/arithmetic.h \
  ../include/geometry/../geometry/point.h \
  ../include/geometry/../geometry/../declarations.h \
  ../include/geometry/../geometry/../base/array.h \
  ../include/geometry/../geometry/../linear_algebra/vector.h \
  ../include/geometry/list_set.h ../include/geometry/../base/utility.h \
  test.h
	-${CXX} ${CXXFLAGS} -c test_rectangle.cc 
 
test_rectangle: test_rectangle.o
	-${CXX} ${LDXXFLAGS} test_rectangle.o -o test_rectangle 

test_parallelotope.o: test_parallelotope.cc ../include/ariadne.h \
  ../include/real_typedef.h ../include/numeric/numerical_types.h \
  ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../utility/stlio.h \
  ../include/numeric/../numeric/../utility/../base/array.h \
  ../include/numeric/../numeric/../utility/../base/../declarations.h \
  ../include/numeric/../numeric/dyadic.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../utility/stlio.h ../include/base/exception.h \
  ../include/base/utility.h ../include/base/../declarations.h \
  ../include/numeric/numerical_types.h ../include/geometry/point.h \
  ../include/geometry/../declarations.h \
  ../include/geometry/../base/array.h \
  ../include/geometry/../linear_algebra/vector.h \
  ../include/geometry/../linear_algebra/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/integer.h \
  ../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/numerical_types.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/approximation.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/numerical_types.h \
  ../include/geometry/rectangle.h \
  ../include/geometry/../numeric/arithmetic.h \
  ../include/geometry/../numeric/interval.h \
  ../include/geometry/../numeric/../declarations.h \
  ../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../numeric/../numeric/arithmetic.h \
  ../include/geometry/../linear_algebra/interval_vector.h \
  ../include/geometry/../linear_algebra/../numeric/interval.h \
  ../include/geometry/../linear_algebra/../linear_algebra/vector.h \
  ../include/geometry/../geometry/point.h \
  ../include/geometry/parallelotope.h \
  ../include/geometry/../linear_algebra/matrix.h \
  ../include/geometry/../geometry/polyhedron.h \
  /ufs/collins/include/ppl.hh \
  ../include/geometry/../geometry/../declarations.h \
  ../include/geometry/../geometry/../linear_algebra/constraint_system.h \
  ../include/geometry/../geometry/../linear_algebra/../declarations.h \
  ../include/geometry/../geometry/../linear_algebra/../numeric/numerical_types.h \
  ../include/geometry/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../include/geometry/../geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../include/geometry/../geometry/../linear_algebra/vector.h \
  ../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../include/geometry/../geometry/../linear_algebra/interval_vector.h \
  ../include/geometry/../geometry/../linear_algebra/generator_system.h \
  ../include/geometry/../geometry/../geometry/point.h \
  ../include/geometry/../geometry/../geometry/rectangle.h
	-${CXX} ${CXXFLAGS} -c test_parallelotope.cc 
 
test_parallelotope: test_parallelotope.o
	-${CXX} ${LDXXFLAGS} test_parallelotope.o -o test_parallelotope 

test_zonotope.o: test_zonotope.cc ../include/ariadne.h \
  ../include/real_typedef.h ../include/numeric/numerical_types.h \
  ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../utility/stlio.h \
  ../include/numeric/../numeric/../utility/../base/array.h \
  ../include/numeric/../numeric/../utility/../base/../declarations.h \
  ../include/numeric/../numeric/dyadic.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../utility/stlio.h ../include/base/exception.h \
  ../include/base/utility.h ../include/base/../declarations.h \
  ../include/numeric/numerical_types.h ../include/geometry/point.h \
  ../include/geometry/../declarations.h \
  ../include/geometry/../base/array.h \
  ../include/geometry/../linear_algebra/vector.h \
  ../include/geometry/../linear_algebra/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/integer.h \
  ../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/numerical_types.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/approximation.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/numerical_types.h \
  ../include/geometry/rectangle.h \
  ../include/geometry/../numeric/arithmetic.h \
  ../include/geometry/../numeric/interval.h \
  ../include/geometry/../numeric/../declarations.h \
  ../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../numeric/../numeric/arithmetic.h \
  ../include/geometry/../linear_algebra/interval_vector.h \
  ../include/geometry/../linear_algebra/../numeric/interval.h \
  ../include/geometry/../linear_algebra/../linear_algebra/vector.h \
  ../include/geometry/../geometry/point.h \
  ../include/geometry/parallelotope.h \
  ../include/geometry/../linear_algebra/matrix.h \
  ../include/geometry/../geometry/polyhedron.h \
  /ufs/collins/include/ppl.hh \
  ../include/geometry/../geometry/../declarations.h \
  ../include/geometry/../geometry/../linear_algebra/constraint_system.h \
  ../include/geometry/../geometry/../linear_algebra/../declarations.h \
  ../include/geometry/../geometry/../linear_algebra/../numeric/numerical_types.h \
  ../include/geometry/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../include/geometry/../geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../include/geometry/../geometry/../linear_algebra/vector.h \
  ../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../include/geometry/../geometry/../linear_algebra/interval_vector.h \
  ../include/geometry/../geometry/../linear_algebra/generator_system.h \
  ../include/geometry/../geometry/../geometry/point.h \
  ../include/geometry/../geometry/../geometry/rectangle.h \
  ../include/geometry/zonotope.h \
  ../include/geometry/../linear_algebra/transformation_system.h \
  ../include/geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../include/geometry/../linear_algebra/../linear_algebra/interval_vector.h \
  ../include/geometry/../linear_algebra/../linear_algebra/interval_matrix.h \
  ../include/geometry/../linear_algebra/../linear_algebra/../linear_algebra/matrix.h \
  ../include/geometry/../linear_algebra/../linear_algebra/../declarations.h \
  ../include/geometry/../linear_algebra/../linear_algebra/../numeric/interval.h \
  test.h
	-${CXX} ${CXXFLAGS} -c test_zonotope.cc 
 
test_zonotope: test_zonotope.o
	-${CXX} ${LDXXFLAGS} test_zonotope.o -o test_zonotope 

test_simplex.o: test_simplex.cc ../include/ariadne.h \
  ../include/real_typedef.h ../include/numeric/numerical_types.h \
  ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../utility/stlio.h \
  ../include/numeric/../numeric/../utility/../base/array.h \
  ../include/numeric/../numeric/../utility/../base/../declarations.h \
  ../include/numeric/../numeric/dyadic.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../utility/stlio.h ../include/base/exception.h \
  ../include/base/utility.h ../include/base/../declarations.h \
  ../include/numeric/numerical_types.h ../include/geometry/point.h \
  ../include/geometry/../declarations.h \
  ../include/geometry/../base/array.h \
  ../include/geometry/../linear_algebra/vector.h \
  ../include/geometry/../linear_algebra/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/integer.h \
  ../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/numerical_types.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/approximation.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/numerical_types.h \
  ../include/geometry/simplex.h ../include/geometry/../geometry/point.h \
  ../include/geometry/polyhedron.h /ufs/collins/include/ppl.hh \
  ../include/geometry/../linear_algebra/constraint_system.h \
  ../include/geometry/../linear_algebra/../numeric/numerical_types.h \
  ../include/geometry/../linear_algebra/../linear_algebra/vector.h \
  ../include/geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../include/geometry/../linear_algebra/../linear_algebra/../declarations.h \
  ../include/geometry/../linear_algebra/../linear_algebra/../numeric/integer.h \
  ../include/geometry/../linear_algebra/matrix.h \
  ../include/geometry/../linear_algebra/interval_vector.h \
  ../include/geometry/../linear_algebra/../numeric/interval.h \
  ../include/geometry/../linear_algebra/../numeric/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/arithmetic.h \
  ../include/geometry/../linear_algebra/generator_system.h \
  ../include/geometry/../geometry/rectangle.h \
  ../include/geometry/../geometry/../declarations.h \
  ../include/geometry/../geometry/../numeric/arithmetic.h \
  ../include/geometry/../geometry/../numeric/interval.h \
  ../include/geometry/../geometry/../linear_algebra/interval_vector.h \
  ../include/geometry/../geometry/../geometry/point.h test.h
	-${CXX} ${CXXFLAGS} -c test_simplex.cc 
 
test_simplex: test_simplex.o
	-${CXX} ${LDXXFLAGS} test_simplex.o -o test_simplex 

test_polyhedron.o: test_polyhedron.cc ../include/ariadne.h \
  ../include/real_typedef.h ../include/numeric/numerical_types.h \
  ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../utility/stlio.h \
  ../include/numeric/../numeric/../utility/../base/array.h \
  ../include/numeric/../numeric/../utility/../base/../declarations.h \
  ../include/numeric/../numeric/dyadic.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../utility/stlio.h ../include/base/exception.h \
  ../include/base/utility.h ../include/base/../declarations.h \
  ../include/numeric/numerical_types.h ../include/geometry/point.h \
  ../include/geometry/../declarations.h \
  ../include/geometry/../base/array.h \
  ../include/geometry/../linear_algebra/vector.h \
  ../include/geometry/../linear_algebra/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/integer.h \
  ../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/numerical_types.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/approximation.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/numerical_types.h \
  ../include/geometry/polyhedron.h /ufs/collins/include/ppl.hh \
  ../include/geometry/../linear_algebra/constraint_system.h \
  ../include/geometry/../linear_algebra/../numeric/numerical_types.h \
  ../include/geometry/../linear_algebra/../linear_algebra/vector.h \
  ../include/geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../include/geometry/../linear_algebra/../linear_algebra/../declarations.h \
  ../include/geometry/../linear_algebra/../linear_algebra/../numeric/integer.h \
  ../include/geometry/../linear_algebra/matrix.h \
  ../include/geometry/../linear_algebra/interval_vector.h \
  ../include/geometry/../linear_algebra/../numeric/interval.h \
  ../include/geometry/../linear_algebra/../numeric/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/arithmetic.h \
  ../include/geometry/../linear_algebra/generator_system.h \
  ../include/geometry/../geometry/point.h \
  ../include/geometry/../geometry/rectangle.h \
  ../include/geometry/../geometry/../declarations.h \
  ../include/geometry/../geometry/../numeric/arithmetic.h \
  ../include/geometry/../geometry/../numeric/interval.h \
  ../include/geometry/../geometry/../linear_algebra/interval_vector.h \
  ../include/geometry/../geometry/../geometry/point.h test.h
	-${CXX} ${CXXFLAGS} -c test_polyhedron.cc 
 
test_polyhedron: test_polyhedron.o
	-${CXX} ${LDXXFLAGS} test_polyhedron.o -o test_polyhedron 

test_list_set.o: test_list_set.cc ../include/ariadne.h \
  ../include/real_typedef.h ../include/numeric/numerical_types.h \
  ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../utility/stlio.h \
  ../include/numeric/../numeric/../utility/../base/array.h \
  ../include/numeric/../numeric/../utility/../base/../declarations.h \
  ../include/numeric/../numeric/dyadic.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../utility/stlio.h ../include/base/exception.h \
  ../include/base/utility.h ../include/base/../declarations.h \
  ../include/numeric/numerical_types.h ../include/geometry/point.h \
  ../include/geometry/../declarations.h \
  ../include/geometry/../base/array.h \
  ../include/geometry/../linear_algebra/vector.h \
  ../include/geometry/../linear_algebra/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/integer.h \
  ../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/numerical_types.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/approximation.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/numerical_types.h \
  ../include/geometry/rectangle.h \
  ../include/geometry/../numeric/arithmetic.h \
  ../include/geometry/../numeric/interval.h \
  ../include/geometry/../numeric/../declarations.h \
  ../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../numeric/../numeric/arithmetic.h \
  ../include/geometry/../linear_algebra/interval_vector.h \
  ../include/geometry/../linear_algebra/../numeric/interval.h \
  ../include/geometry/../linear_algebra/../linear_algebra/vector.h \
  ../include/geometry/../geometry/point.h ../include/geometry/list_set.h \
  ../include/geometry/../base/utility.h test.h
	-${CXX} ${CXXFLAGS} -c test_list_set.cc 
 
test_list_set: test_list_set.o
	-${CXX} ${LDXXFLAGS} test_list_set.o -o test_list_set 

test_grid_set.o: test_grid_set.cc ../include/ariadne.h \
  ../include/real_typedef.h ../include/numeric/numerical_types.h \
  ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../utility/stlio.h \
  ../include/numeric/../numeric/../utility/../base/array.h \
  ../include/numeric/../numeric/../utility/../base/../declarations.h \
  ../include/numeric/../numeric/dyadic.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../utility/stlio.h ../include/base/exception.h \
  ../include/base/utility.h ../include/base/../declarations.h \
  ../include/base/binary_word.h ../include/numeric/numerical_types.h \
  ../include/geometry/point.h ../include/geometry/../declarations.h \
  ../include/geometry/../base/array.h \
  ../include/geometry/../linear_algebra/vector.h \
  ../include/geometry/../linear_algebra/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/integer.h \
  ../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/numerical_types.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/approximation.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/numerical_types.h \
  ../include/geometry/rectangle.h \
  ../include/geometry/../numeric/arithmetic.h \
  ../include/geometry/../numeric/interval.h \
  ../include/geometry/../numeric/../declarations.h \
  ../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../numeric/../numeric/arithmetic.h \
  ../include/geometry/../linear_algebra/interval_vector.h \
  ../include/geometry/../linear_algebra/../numeric/interval.h \
  ../include/geometry/../linear_algebra/../linear_algebra/vector.h \
  ../include/geometry/../geometry/point.h ../include/geometry/grid_set.h \
  ../include/geometry/../geometry/rectangle.h \
  ../include/geometry/../geometry/grid.h \
  ../include/geometry/../geometry/../declarations.h \
  ../include/geometry/../geometry/../base/array.h \
  ../include/geometry/../geometry/../numeric/arithmetic.h \
  ../include/geometry/../geometry/../geometry/lattice_set.h \
  ../include/geometry/../geometry/../geometry/../base/array.h \
  ../include/geometry/../geometry/../geometry/../base/iterator.h \
  ../include/geometry/../geometry/../geometry/../base/../base/array.h \
  ../include/geometry/../geometry/../geometry/../numeric/interval.h \
  ../include/geometry/../geometry/../geometry/../base/array_operations.h \
  ../include/geometry/../geometry/../geometry/../base/../declarations.h \
  ../include/geometry/../geometry/lattice_set.h \
  ../include/geometry/list_set.h ../include/geometry/../base/utility.h \
  test.h
	-${CXX} ${CXXFLAGS} -c test_grid_set.cc 
 
test_grid_set: test_grid_set.o
	-${CXX} ${LDXXFLAGS} test_grid_set.o -o test_grid_set 

test_partition_tree_set.o: test_partition_tree_set.cc \
  ../include/ariadne.h ../include/real_typedef.h \
  ../include/numeric/numerical_types.h \
  ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../utility/stlio.h \
  ../include/numeric/../numeric/../utility/../base/array.h \
  ../include/numeric/../numeric/../utility/../base/../declarations.h \
  ../include/numeric/../numeric/dyadic.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../utility/stlio.h ../include/base/exception.h \
  ../include/base/utility.h ../include/base/../declarations.h \
  ../include/base/binary_word.h ../include/numeric/numerical_types.h \
  ../include/linear_algebra/matrix.h \
  ../include/linear_algebra/../declarations.h \
  ../include/linear_algebra/../numeric/integer.h \
  ../include/linear_algebra/vector.h \
  ../include/linear_algebra/../numeric/arithmetic.h \
  ../include/linear_algebra/../numeric/../numeric/numerical_types.h \
  ../include/linear_algebra/../numeric/../numeric/approximation.h \
  ../include/linear_algebra/../numeric/../numeric/../declarations.h \
  ../include/linear_algebra/../numeric/../numeric/../numeric/numerical_types.h \
  ../include/geometry/point.h ../include/geometry/../declarations.h \
  ../include/geometry/../base/array.h \
  ../include/geometry/../linear_algebra/vector.h \
  ../include/geometry/rectangle.h \
  ../include/geometry/../numeric/arithmetic.h \
  ../include/geometry/../numeric/interval.h \
  ../include/geometry/../numeric/../declarations.h \
  ../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../numeric/../numeric/arithmetic.h \
  ../include/geometry/../linear_algebra/interval_vector.h \
  ../include/geometry/../linear_algebra/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/interval.h \
  ../include/geometry/../linear_algebra/../linear_algebra/vector.h \
  ../include/geometry/../geometry/point.h \
  ../include/geometry/parallelotope.h \
  ../include/geometry/../linear_algebra/matrix.h \
  ../include/geometry/../geometry/polyhedron.h \
  /ufs/collins/include/ppl.hh \
  ../include/geometry/../geometry/../declarations.h \
  ../include/geometry/../geometry/../linear_algebra/constraint_system.h \
  ../include/geometry/../geometry/../linear_algebra/../declarations.h \
  ../include/geometry/../geometry/../linear_algebra/../numeric/numerical_types.h \
  ../include/geometry/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../include/geometry/../geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../include/geometry/../geometry/../linear_algebra/vector.h \
  ../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../include/geometry/../geometry/../linear_algebra/interval_vector.h \
  ../include/geometry/../geometry/../linear_algebra/generator_system.h \
  ../include/geometry/../geometry/../geometry/point.h \
  ../include/geometry/../geometry/../geometry/rectangle.h \
  ../include/geometry/list_set.h ../include/geometry/../base/utility.h \
  ../include/geometry/partition_tree_set.h \
  ../include/geometry/../base/binary_word.h \
  ../include/geometry/../base/binary_tree.h \
  ../include/geometry/../base/../declarations.h \
  ../include/geometry/../base/../base/array.h \
  ../include/geometry/../base/../base/binary_word.h \
  ../include/geometry/../base/iterator.h \
  ../include/geometry/../geometry/subdivision_tree_set.h \
  ../include/geometry/../geometry/../base/array.h \
  ../include/geometry/../geometry/../base/sequence.h \
  ../include/geometry/../geometry/../base/binary_word.h \
  ../include/geometry/../geometry/../base/binary_tree.h \
  ../include/geometry/../geometry/../base/iterator.h \
  ../include/utility/epsfstream.h ../include/utility/../geometry/point.h \
  ../include/utility/../geometry/rectangle.h \
  ../include/utility/../geometry/list_set.h \
  ../include/utility/../geometry/grid_set.h \
  ../include/utility/../geometry/../declarations.h \
  ../include/utility/../geometry/../base/array.h \
  ../include/utility/../geometry/../geometry/rectangle.h \
  ../include/utility/../geometry/../geometry/grid.h \
  ../include/utility/../geometry/../geometry/../declarations.h \
  ../include/utility/../geometry/../geometry/../base/array.h \
  ../include/utility/../geometry/../geometry/../numeric/arithmetic.h \
  ../include/utility/../geometry/../geometry/../geometry/lattice_set.h \
  ../include/utility/../geometry/../geometry/../geometry/../base/array.h \
  ../include/utility/../geometry/../geometry/../geometry/../base/iterator.h \
  ../include/utility/../geometry/../geometry/../geometry/../numeric/interval.h \
  ../include/utility/../geometry/../geometry/../geometry/../base/array_operations.h \
  ../include/utility/../geometry/../geometry/../geometry/../base/../declarations.h \
  ../include/utility/../geometry/../geometry/lattice_set.h \
  ../include/utility/../geometry/parallelotope.h \
  ../include/utility/../geometry/zonotope.h \
  ../include/utility/../geometry/../linear_algebra/vector.h \
  ../include/utility/../geometry/../linear_algebra/matrix.h \
  ../include/utility/../geometry/../linear_algebra/transformation_system.h \
  ../include/utility/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../include/utility/../geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../include/utility/../geometry/../linear_algebra/../linear_algebra/interval_vector.h \
  ../include/utility/../geometry/../linear_algebra/../linear_algebra/interval_matrix.h \
  ../include/utility/../geometry/../linear_algebra/../linear_algebra/../linear_algebra/matrix.h \
  ../include/utility/../geometry/../linear_algebra/../linear_algebra/../declarations.h \
  ../include/utility/../geometry/../linear_algebra/../linear_algebra/../numeric/interval.h \
  ../include/utility/../geometry/../numeric/interval.h \
  ../include/utility/../geometry/../geometry/point.h \
  ../include/utility/../geometry/polyhedron.h \
  ../include/utility/../geometry/partition_tree_set.h \
  ../include/utility/../linear_algebra/matrix.h \
  ../include/utility/../system/affine_map.h \
  ../include/utility/../system/../numeric/numerical_types.h \
  ../include/utility/../system/../linear_algebra/vector.h \
  ../include/utility/../system/../linear_algebra/matrix.h \
  ../include/utility/../system/../linear_algebra/interval_matrix.h \
  ../include/utility/../system/../geometry/point.h \
  ../include/utility/../system/../geometry/simplex.h \
  ../include/utility/../system/../geometry/../declarations.h \
  ../include/utility/../system/../geometry/../base/array.h \
  ../include/utility/../system/../geometry/../geometry/point.h \
  ../include/utility/../system/../geometry/rectangle.h \
  ../include/utility/../system/../geometry/parallelotope.h \
  ../include/utility/../system/../geometry/zonotope.h \
  ../include/utility/../system/../geometry/polyhedron.h \
  ../include/utility/../system/../geometry/list_set.h \
  ../include/utility/../system/../geometry/grid_set.h \
  ../include/utility/../system/../system/map.h \
  ../include/utility/../system/../system/../declarations.h test.h
	-${CXX} ${CXXFLAGS} -c test_partition_tree_set.cc 
 
test_partition_tree_set: test_partition_tree_set.o
	-${CXX} ${LDXXFLAGS} test_partition_tree_set.o -o test_partition_tree_set 

test_epsfstream.o: test_epsfstream.cc ../include/ariadne.h \
  ../include/real_typedef.h ../include/numeric/numerical_types.h \
  ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../utility/stlio.h \
  ../include/numeric/../numeric/../utility/../base/array.h \
  ../include/numeric/../numeric/../utility/../base/../declarations.h \
  ../include/numeric/../numeric/dyadic.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../utility/stlio.h ../include/base/exception.h \
  ../include/base/utility.h ../include/base/../declarations.h \
  ../include/numeric/numerical_types.h ../include/geometry/point.h \
  ../include/geometry/../declarations.h \
  ../include/geometry/../base/array.h \
  ../include/geometry/../linear_algebra/vector.h \
  ../include/geometry/../linear_algebra/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/integer.h \
  ../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/numerical_types.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/approximation.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/numerical_types.h \
  ../include/geometry/rectangle.h \
  ../include/geometry/../numeric/arithmetic.h \
  ../include/geometry/../numeric/interval.h \
  ../include/geometry/../numeric/../declarations.h \
  ../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../numeric/../numeric/arithmetic.h \
  ../include/geometry/../linear_algebra/interval_vector.h \
  ../include/geometry/../linear_algebra/../numeric/interval.h \
  ../include/geometry/../linear_algebra/../linear_algebra/vector.h \
  ../include/geometry/../geometry/point.h ../include/geometry/list_set.h \
  ../include/geometry/../base/utility.h ../include/utility/epsfstream.h \
  ../include/utility/../geometry/point.h \
  ../include/utility/../geometry/rectangle.h \
  ../include/utility/../geometry/list_set.h \
  ../include/utility/../geometry/grid_set.h \
  ../include/utility/../geometry/../declarations.h \
  ../include/utility/../geometry/../base/array.h \
  ../include/utility/../geometry/../geometry/rectangle.h \
  ../include/utility/../geometry/../geometry/grid.h \
  ../include/utility/../geometry/../geometry/../declarations.h \
  ../include/utility/../geometry/../geometry/../base/array.h \
  ../include/utility/../geometry/../geometry/../numeric/arithmetic.h \
  ../include/utility/../geometry/../geometry/../geometry/lattice_set.h \
  ../include/utility/../geometry/../geometry/../geometry/../base/array.h \
  ../include/utility/../geometry/../geometry/../geometry/../base/iterator.h \
  ../include/utility/../geometry/../geometry/../geometry/../base/../base/array.h \
  ../include/utility/../geometry/../geometry/../geometry/../numeric/interval.h \
  ../include/utility/../geometry/../geometry/../geometry/../base/array_operations.h \
  ../include/utility/../geometry/../geometry/../geometry/../base/../declarations.h \
  ../include/utility/../geometry/../geometry/lattice_set.h \
  ../include/utility/../geometry/parallelotope.h \
  ../include/utility/../geometry/../numeric/interval.h \
  ../include/utility/../geometry/../linear_algebra/vector.h \
  ../include/utility/../geometry/../linear_algebra/matrix.h \
  ../include/utility/../geometry/../linear_algebra/../declarations.h \
  ../include/utility/../geometry/../linear_algebra/../numeric/integer.h \
  ../include/utility/../geometry/../geometry/point.h \
  ../include/utility/../geometry/../geometry/polyhedron.h \
  /ufs/collins/include/ppl.hh \
  ../include/utility/../geometry/../geometry/../linear_algebra/constraint_system.h \
  ../include/utility/../geometry/../geometry/../linear_algebra/../declarations.h \
  ../include/utility/../geometry/../geometry/../linear_algebra/../numeric/numerical_types.h \
  ../include/utility/../geometry/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../include/utility/../geometry/../geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../include/utility/../geometry/../geometry/../linear_algebra/vector.h \
  ../include/utility/../geometry/../geometry/../linear_algebra/matrix.h \
  ../include/utility/../geometry/../geometry/../linear_algebra/interval_vector.h \
  ../include/utility/../geometry/../geometry/../linear_algebra/generator_system.h \
  ../include/utility/../geometry/../geometry/../geometry/point.h \
  ../include/utility/../geometry/../geometry/../geometry/rectangle.h \
  ../include/utility/../geometry/zonotope.h \
  ../include/utility/../geometry/../linear_algebra/transformation_system.h \
  ../include/utility/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../include/utility/../geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../include/utility/../geometry/../linear_algebra/../linear_algebra/interval_vector.h \
  ../include/utility/../geometry/../linear_algebra/../linear_algebra/interval_matrix.h \
  ../include/utility/../geometry/../linear_algebra/../linear_algebra/../linear_algebra/matrix.h \
  ../include/utility/../geometry/../linear_algebra/../linear_algebra/../declarations.h \
  ../include/utility/../geometry/../linear_algebra/../linear_algebra/../numeric/interval.h \
  ../include/utility/../geometry/polyhedron.h \
  ../include/utility/../geometry/partition_tree_set.h \
  ../include/utility/../geometry/../base/binary_word.h \
  ../include/utility/../geometry/../base/../declarations.h \
  ../include/utility/../geometry/../base/binary_tree.h \
  ../include/utility/../geometry/../base/../base/array.h \
  ../include/utility/../geometry/../base/../base/binary_word.h \
  ../include/utility/../geometry/../base/iterator.h \
  ../include/utility/../geometry/../geometry/subdivision_tree_set.h \
  ../include/utility/../geometry/../geometry/../base/sequence.h \
  ../include/utility/../geometry/../geometry/../base/binary_word.h \
  ../include/utility/../geometry/../geometry/../base/binary_tree.h \
  ../include/utility/../geometry/../geometry/../base/iterator.h \
  ../include/utility/../linear_algebra/matrix.h \
  ../include/utility/../system/affine_map.h \
  ../include/utility/../system/../numeric/numerical_types.h \
  ../include/utility/../system/../linear_algebra/vector.h \
  ../include/utility/../system/../linear_algebra/matrix.h \
  ../include/utility/../system/../linear_algebra/interval_matrix.h \
  ../include/utility/../system/../geometry/point.h \
  ../include/utility/../system/../geometry/simplex.h \
  ../include/utility/../system/../geometry/../declarations.h \
  ../include/utility/../system/../geometry/../base/array.h \
  ../include/utility/../system/../geometry/../geometry/point.h \
  ../include/utility/../system/../geometry/rectangle.h \
  ../include/utility/../system/../geometry/parallelotope.h \
  ../include/utility/../system/../geometry/zonotope.h \
  ../include/utility/../system/../geometry/polyhedron.h \
  ../include/utility/../system/../geometry/list_set.h \
  ../include/utility/../system/../geometry/grid_set.h \
  ../include/utility/../system/../system/map.h \
  ../include/utility/../system/../system/../declarations.h test.h
	-${CXX} ${CXXFLAGS} -c test_epsfstream.cc 
 
test_epsfstream: test_epsfstream.o
	-${CXX} ${LDXXFLAGS} test_epsfstream.o -o test_epsfstream 

test_interval_newton.o: test_interval_newton.cc \
  ../include/geometry/point.h ../include/geometry/../declarations.h \
  ../include/geometry/../base/array.h \
  ../include/geometry/../base/../declarations.h \
  ../include/geometry/../linear_algebra/vector.h \
  ../include/geometry/../linear_algebra/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/integer.h \
  ../include/geometry/../linear_algebra/../numeric/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/../utility/stlio.h \
  ../include/geometry/../linear_algebra/../numeric/../utility/../base/array.h \
  ../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/numerical_types.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/integer.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/dyadic.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/../utility/stlio.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/approximation.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/numerical_types.h \
  ../include/geometry/rectangle.h \
  ../include/geometry/../numeric/arithmetic.h \
  ../include/geometry/../numeric/interval.h \
  ../include/geometry/../numeric/../declarations.h \
  ../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../numeric/../numeric/arithmetic.h \
  ../include/geometry/../linear_algebra/interval_vector.h \
  ../include/geometry/../linear_algebra/../numeric/interval.h \
  ../include/geometry/../linear_algebra/../linear_algebra/vector.h \
  ../include/geometry/../geometry/point.h ../include/system/henon_map.h \
  ../include/system/../linear_algebra/matrix.h \
  ../include/system/../linear_algebra/../declarations.h \
  ../include/system/../linear_algebra/../numeric/integer.h \
  ../include/system/../linear_algebra/interval_matrix.h \
  ../include/system/../linear_algebra/../linear_algebra/matrix.h \
  ../include/system/../linear_algebra/../numeric/interval.h \
  ../include/system/../geometry/point.h \
  ../include/system/../geometry/rectangle.h \
  ../include/system/../geometry/parallelotope.h \
  ../include/system/../geometry/../declarations.h \
  ../include/system/../geometry/../numeric/interval.h \
  ../include/system/../geometry/../linear_algebra/vector.h \
  ../include/system/../geometry/../linear_algebra/matrix.h \
  ../include/system/../geometry/../geometry/point.h \
  ../include/system/../geometry/../geometry/polyhedron.h \
  /ufs/collins/include/ppl.hh \
  ../include/system/../geometry/../geometry/../declarations.h \
  ../include/system/../geometry/../geometry/../linear_algebra/constraint_system.h \
  ../include/system/../geometry/../geometry/../linear_algebra/../declarations.h \
  ../include/system/../geometry/../geometry/../linear_algebra/../numeric/numerical_types.h \
  ../include/system/../geometry/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../include/system/../geometry/../geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../include/system/../geometry/../geometry/../linear_algebra/vector.h \
  ../include/system/../geometry/../geometry/../linear_algebra/matrix.h \
  ../include/system/../geometry/../geometry/../linear_algebra/interval_vector.h \
  ../include/system/../geometry/../geometry/../linear_algebra/generator_system.h \
  ../include/system/../geometry/../geometry/../geometry/point.h \
  ../include/system/../geometry/../geometry/../geometry/rectangle.h \
  ../include/system/../system/map.h \
  ../include/system/../system/../declarations.h \
  ../include/system/../evaluation/apply.h \
  ../include/system/../evaluation/../declarations.h \
  ../include/evaluation/newton.h ../include/evaluation/../declarations.h \
  ../include/linear_algebra/interval_vector.h \
  ../include/linear_algebra/interval_matrix.h \
  ../include/geometry/rectangle.h ../include/system/map.h \
  ../include/system/vector_field.h ../include/system/../declarations.h \
  test.h ../include/real_typedef.h ../include/numeric/numerical_types.h
	-${CXX} ${CXXFLAGS} -c test_interval_newton.cc 
 
test_interval_newton: test_interval_newton.o
	-${CXX} ${LDXXFLAGS} test_interval_newton.o -o test_interval_newton 

test_chainreach.o: test_chainreach.cc ../include/real_typedef.h \
  ../include/numeric/numerical_types.h \
  ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../utility/stlio.h \
  ../include/numeric/../numeric/../utility/../base/array.h \
  ../include/numeric/../numeric/../utility/../base/../declarations.h \
  ../include/numeric/../numeric/dyadic.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../utility/stlio.h ../include/utility/epsfstream.h \
  ../include/utility/../geometry/point.h \
  ../include/utility/../geometry/../declarations.h \
  ../include/utility/../geometry/../base/array.h \
  ../include/utility/../geometry/../linear_algebra/vector.h \
  ../include/utility/../geometry/../linear_algebra/../declarations.h \
  ../include/utility/../geometry/../linear_algebra/../numeric/integer.h \
  ../include/utility/../geometry/../linear_algebra/../numeric/arithmetic.h \
  ../include/utility/../geometry/../linear_algebra/../numeric/../numeric/numerical_types.h \
  ../include/utility/../geometry/../linear_algebra/../numeric/../numeric/approximation.h \
  ../include/utility/../geometry/../linear_algebra/../numeric/../numeric/../declarations.h \
  ../include/utility/../geometry/../linear_algebra/../numeric/../numeric/../numeric/numerical_types.h \
  ../include/utility/../geometry/rectangle.h \
  ../include/utility/../geometry/../numeric/arithmetic.h \
  ../include/utility/../geometry/../numeric/interval.h \
  ../include/utility/../geometry/../numeric/../declarations.h \
  ../include/utility/../geometry/../numeric/../numeric/numerical_traits.h \
  ../include/utility/../geometry/../numeric/../numeric/arithmetic.h \
  ../include/utility/../geometry/../linear_algebra/interval_vector.h \
  ../include/utility/../geometry/../linear_algebra/../numeric/interval.h \
  ../include/utility/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../include/utility/../geometry/../geometry/point.h \
  ../include/utility/../geometry/list_set.h \
  ../include/utility/../geometry/../base/utility.h \
  ../include/utility/../geometry/../base/../declarations.h \
  ../include/utility/../geometry/grid_set.h \
  ../include/utility/../geometry/../geometry/rectangle.h \
  ../include/utility/../geometry/../geometry/grid.h \
  ../include/utility/../geometry/../geometry/../declarations.h \
  ../include/utility/../geometry/../geometry/../base/array.h \
  ../include/utility/../geometry/../geometry/../numeric/arithmetic.h \
  ../include/utility/../geometry/../geometry/../geometry/lattice_set.h \
  ../include/utility/../geometry/../geometry/../geometry/../base/array.h \
  ../include/utility/../geometry/../geometry/../geometry/../base/iterator.h \
  ../include/utility/../geometry/../geometry/../geometry/../base/../base/array.h \
  ../include/utility/../geometry/../geometry/../geometry/../numeric/interval.h \
  ../include/utility/../geometry/../geometry/../geometry/../base/array_operations.h \
  ../include/utility/../geometry/../geometry/../geometry/../base/../declarations.h \
  ../include/utility/../geometry/../geometry/lattice_set.h \
  ../include/utility/../geometry/parallelotope.h \
  ../include/utility/../geometry/../linear_algebra/matrix.h \
  ../include/utility/../geometry/../geometry/polyhedron.h \
  /ufs/collins/include/ppl.hh \
  ../include/utility/../geometry/../geometry/../linear_algebra/constraint_system.h \
  ../include/utility/../geometry/../geometry/../linear_algebra/../declarations.h \
  ../include/utility/../geometry/../geometry/../linear_algebra/../numeric/numerical_types.h \
  ../include/utility/../geometry/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../include/utility/../geometry/../geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../include/utility/../geometry/../geometry/../linear_algebra/vector.h \
  ../include/utility/../geometry/../geometry/../linear_algebra/matrix.h \
  ../include/utility/../geometry/../geometry/../linear_algebra/interval_vector.h \
  ../include/utility/../geometry/../geometry/../linear_algebra/generator_system.h \
  ../include/utility/../geometry/../geometry/../geometry/point.h \
  ../include/utility/../geometry/../geometry/../geometry/rectangle.h \
  ../include/utility/../geometry/zonotope.h \
  ../include/utility/../geometry/../linear_algebra/transformation_system.h \
  ../include/utility/../geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../include/utility/../geometry/../linear_algebra/../linear_algebra/interval_vector.h \
  ../include/utility/../geometry/../linear_algebra/../linear_algebra/interval_matrix.h \
  ../include/utility/../geometry/../linear_algebra/../linear_algebra/../linear_algebra/matrix.h \
  ../include/utility/../geometry/../linear_algebra/../linear_algebra/../declarations.h \
  ../include/utility/../geometry/../linear_algebra/../linear_algebra/../numeric/interval.h \
  ../include/utility/../geometry/polyhedron.h \
  ../include/utility/../geometry/partition_tree_set.h \
  ../include/utility/../geometry/../base/binary_word.h \
  ../include/utility/../geometry/../base/binary_tree.h \
  ../include/utility/../geometry/../base/../base/array.h \
  ../include/utility/../geometry/../base/../base/binary_word.h \
  ../include/utility/../geometry/../base/iterator.h \
  ../include/utility/../geometry/../geometry/subdivision_tree_set.h \
  ../include/utility/../geometry/../geometry/../base/sequence.h \
  ../include/utility/../geometry/../geometry/../base/binary_word.h \
  ../include/utility/../geometry/../geometry/../base/binary_tree.h \
  ../include/utility/../geometry/../geometry/../base/iterator.h \
  ../include/utility/../linear_algebra/matrix.h \
  ../include/utility/../system/affine_map.h \
  ../include/utility/../system/../numeric/numerical_types.h \
  ../include/utility/../system/../linear_algebra/vector.h \
  ../include/utility/../system/../linear_algebra/matrix.h \
  ../include/utility/../system/../linear_algebra/interval_matrix.h \
  ../include/utility/../system/../geometry/point.h \
  ../include/utility/../system/../geometry/simplex.h \
  ../include/utility/../system/../geometry/../declarations.h \
  ../include/utility/../system/../geometry/../base/array.h \
  ../include/utility/../system/../geometry/../geometry/point.h \
  ../include/utility/../system/../geometry/rectangle.h \
  ../include/utility/../system/../geometry/parallelotope.h \
  ../include/utility/../system/../geometry/zonotope.h \
  ../include/utility/../system/../geometry/polyhedron.h \
  ../include/utility/../system/../geometry/list_set.h \
  ../include/utility/../system/../geometry/grid_set.h \
  ../include/utility/../system/../system/map.h \
  ../include/utility/../system/../system/../declarations.h \
  ../include/geometry/point.h ../include/geometry/rectangle.h \
  ../include/geometry/grid.h ../include/geometry/grid_set.h \
  ../include/geometry/partition_tree_set.h ../include/system/henon_map.h \
  ../include/system/../linear_algebra/matrix.h \
  ../include/system/../linear_algebra/interval_matrix.h \
  ../include/system/../geometry/point.h \
  ../include/system/../geometry/rectangle.h \
  ../include/system/../geometry/parallelotope.h \
  ../include/system/../system/map.h \
  ../include/system/../evaluation/apply.h \
  ../include/system/../evaluation/../declarations.h \
  ../include/evaluation/apply.h test.h
	-${CXX} ${CXXFLAGS} -c test_chainreach.cc 
 
test_chainreach: test_chainreach.o
	-${CXX} ${LDXXFLAGS} test_chainreach.o -o test_chainreach 

test_integration_step.o: test_integration_step.cc ../include/ariadne.h \
  ../include/real_typedef.h ../include/numeric/numerical_types.h \
  ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../utility/stlio.h \
  ../include/numeric/../numeric/../utility/../base/array.h \
  ../include/numeric/../numeric/../utility/../base/../declarations.h \
  ../include/numeric/../numeric/dyadic.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../utility/stlio.h \
  ../include/numeric/numerical_types.h ../include/linear_algebra/vector.h \
  ../include/linear_algebra/../declarations.h \
  ../include/linear_algebra/../numeric/integer.h \
  ../include/linear_algebra/../numeric/arithmetic.h \
  ../include/linear_algebra/../numeric/../numeric/numerical_types.h \
  ../include/linear_algebra/../numeric/../numeric/approximation.h \
  ../include/linear_algebra/../numeric/../numeric/../declarations.h \
  ../include/linear_algebra/../numeric/../numeric/../numeric/numerical_types.h \
  ../include/linear_algebra/matrix.h ../include/geometry/rectangle.h \
  ../include/geometry/../declarations.h \
  ../include/geometry/../numeric/arithmetic.h \
  ../include/geometry/../numeric/interval.h \
  ../include/geometry/../numeric/../declarations.h \
  ../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../numeric/../numeric/arithmetic.h \
  ../include/geometry/../linear_algebra/interval_vector.h \
  ../include/geometry/../linear_algebra/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/interval.h \
  ../include/geometry/../linear_algebra/../linear_algebra/vector.h \
  ../include/geometry/../geometry/point.h \
  ../include/geometry/../geometry/../declarations.h \
  ../include/geometry/../geometry/../base/array.h \
  ../include/geometry/../geometry/../linear_algebra/vector.h \
  ../include/geometry/parallelotope.h \
  ../include/geometry/../linear_algebra/vector.h \
  ../include/geometry/../linear_algebra/matrix.h \
  ../include/geometry/../geometry/polyhedron.h \
  /ufs/collins/include/ppl.hh \
  ../include/geometry/../geometry/../linear_algebra/constraint_system.h \
  ../include/geometry/../geometry/../linear_algebra/../declarations.h \
  ../include/geometry/../geometry/../linear_algebra/../numeric/numerical_types.h \
  ../include/geometry/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../include/geometry/../geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../include/geometry/../geometry/../linear_algebra/interval_vector.h \
  ../include/geometry/../geometry/../linear_algebra/generator_system.h \
  ../include/geometry/../geometry/../geometry/point.h \
  ../include/geometry/../geometry/../geometry/rectangle.h \
  ../include/system/affine_vector_field.h \
  ../include/system/../linear_algebra/vector.h \
  ../include/system/../linear_algebra/matrix.h \
  ../include/system/../linear_algebra/interval_vector.h \
  ../include/system/../linear_algebra/interval_matrix.h \
  ../include/system/../linear_algebra/../linear_algebra/matrix.h \
  ../include/system/../linear_algebra/../declarations.h \
  ../include/system/../linear_algebra/../numeric/interval.h \
  ../include/system/../system/vector_field.h \
  ../include/system/../system/../declarations.h \
  ../include/system/../system/affine_map.h \
  ../include/system/../system/../numeric/numerical_types.h \
  ../include/system/../system/../linear_algebra/vector.h \
  ../include/system/../system/../linear_algebra/matrix.h \
  ../include/system/../system/../linear_algebra/interval_matrix.h \
  ../include/system/../system/../geometry/point.h \
  ../include/system/../system/../geometry/simplex.h \
  ../include/system/../system/../geometry/../declarations.h \
  ../include/system/../system/../geometry/../base/array.h \
  ../include/system/../system/../geometry/../geometry/point.h \
  ../include/system/../system/../geometry/rectangle.h \
  ../include/system/../system/../geometry/parallelotope.h \
  ../include/system/../system/../geometry/zonotope.h \
  ../include/system/../system/../geometry/../linear_algebra/vector.h \
  ../include/system/../system/../geometry/../linear_algebra/matrix.h \
  ../include/system/../system/../geometry/../linear_algebra/transformation_system.h \
  ../include/system/../system/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../include/system/../system/../geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../include/system/../system/../geometry/../linear_algebra/../linear_algebra/interval_vector.h \
  ../include/system/../system/../geometry/../linear_algebra/../linear_algebra/interval_matrix.h \
  ../include/system/../system/../geometry/../numeric/interval.h \
  ../include/system/../system/../geometry/polyhedron.h \
  ../include/system/../system/../geometry/list_set.h \
  ../include/system/../system/../geometry/../base/utility.h \
  ../include/system/../system/../geometry/../base/../declarations.h \
  ../include/system/../system/../geometry/grid_set.h \
  ../include/system/../system/../geometry/../geometry/rectangle.h \
  ../include/system/../system/../geometry/../geometry/grid.h \
  ../include/system/../system/../geometry/../geometry/../declarations.h \
  ../include/system/../system/../geometry/../geometry/../base/array.h \
  ../include/system/../system/../geometry/../geometry/../numeric/arithmetic.h \
  ../include/system/../system/../geometry/../geometry/../geometry/lattice_set.h \
  ../include/system/../system/../geometry/../geometry/../geometry/../base/array.h \
  ../include/system/../system/../geometry/../geometry/../geometry/../base/iterator.h \
  ../include/system/../system/../geometry/../geometry/../geometry/../base/../base/array.h \
  ../include/system/../system/../geometry/../geometry/../geometry/../numeric/interval.h \
  ../include/system/../system/../geometry/../geometry/../geometry/../base/array_operations.h \
  ../include/system/../system/../geometry/../geometry/../geometry/../base/../declarations.h \
  ../include/system/../system/../geometry/../geometry/lattice_set.h \
  ../include/system/../system/../system/map.h \
  ../include/system/../system/../system/../declarations.h \
  ../include/evaluation/lohner_integrator.h \
  ../include/evaluation/../declarations.h \
  ../include/evaluation/../evaluation/integrator.h \
  ../include/evaluation/../evaluation/../declarations.h test.h
	-${CXX} ${CXXFLAGS} -c test_integration_step.cc 
 
test_integration_step: test_integration_step.o
	-${CXX} ${LDXXFLAGS} test_integration_step.o -o test_integration_step 

