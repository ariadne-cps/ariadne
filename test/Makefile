include ../config.mk

TESTS = \
  test_array test_binary_word \
  test_lattice \
  test_arithmetic test_function test_interval \
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
#LDXXFLAGS += .-L../src/ -lariadne -lmpfr -lgmp -lgmpxx -lppl #-lblas
LDXXFLAGS += ../src/libariadne.so -lmpfr -lgmp -lgmpxx -lppl #-lblas

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
	echo "	-$$""{CXX} $$ff.o $$""{LDXXFLAGS} -o $$ff "; \
	echo ""; done >> Makefile

depclean:
	@sed -e '/^### Do Not edit this line$$/,$$d' Makefile \
		> MMM.$$$$ && mv MMM.$$$$ Makefile
	@echo "### Do Not edit this line" >> Makefile

### Do Not edit this line
### Everything below is auto-generated

all: ${TESTS}
	-./test_array
	-./test_binary_word
	-./test_lattice
	-./test_arithmetic
	-./test_function
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
	-${CXX} test_array.o ${LDXXFLAGS} -o test_array 

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
	-${CXX} test_binary_word.o ${LDXXFLAGS} -o test_binary_word 

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
	-${CXX} test_lattice.o ${LDXXFLAGS} -o test_lattice 

test_arithmetic.o: test_arithmetic.cc /ufs/collins/include/mpfr.h \
  ../include/numeric/arithmetic.h \
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
	-${CXX} test_arithmetic.o ${LDXXFLAGS} -o test_arithmetic 

test_function.o: test_function.cc ../include/real_typedef.h \
  ../include/numeric/numerical_types.h \
  ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../utility/stlio.h \
  ../include/numeric/../numeric/../utility/../base/array.h \
  ../include/numeric/../numeric/../utility/../base/../declarations.h \
  ../include/numeric/../numeric/dyadic.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../utility/stlio.h ../include/numeric/function.h \
  /ufs/collins/include/mpfr.h \
  ../include/numeric/../numeric/numerical_types.h test.h
	-${CXX} ${CXXFLAGS} -c test_function.cc 
 
test_function: test_function.o
	-${CXX} test_function.o ${LDXXFLAGS} -o test_function 

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
	-${CXX} test_interval.o ${LDXXFLAGS} -o test_interval 

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
	-${CXX} test_vector.o ${LDXXFLAGS} -o test_vector 

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
	-${CXX} test_matrix.o ${LDXXFLAGS} -o test_matrix 

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
  ../include/linear_algebra/interval_matrix.h \
  ../include/linear_algebra/../linear_algebra/matrix.h \
  ../include/linear_algebra/../numeric/interval.h \
  ../include/linear_algebra/../numeric/../declarations.h \
  ../include/linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../include/linear_algebra/../numeric/../numeric/arithmetic.h \
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
	-${CXX} test_linear_algebra.o ${LDXXFLAGS} -o test_linear_algebra 

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
	-${CXX} test_point.o ${LDXXFLAGS} -o test_point 

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
	-${CXX} test_rectangle.o ${LDXXFLAGS} -o test_rectangle 

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
	-${CXX} test_parallelotope.o ${LDXXFLAGS} -o test_parallelotope 

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
	-${CXX} test_zonotope.o ${LDXXFLAGS} -o test_zonotope 

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
	-${CXX} test_simplex.o ${LDXXFLAGS} -o test_simplex 

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
	-${CXX} test_polyhedron.o ${LDXXFLAGS} -o test_polyhedron 

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
	-${CXX} test_list_set.o ${LDXXFLAGS} -o test_list_set 

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
	-${CXX} test_grid_set.o ${LDXXFLAGS} -o test_grid_set 

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
  ../include/output/epsfstream.h ../include/output/../geometry/point.h \
  ../include/output/../geometry/rectangle.h \
  ../include/output/../geometry/list_set.h \
  ../include/output/../geometry/grid_set.h \
  ../include/output/../geometry/../declarations.h \
  ../include/output/../geometry/../base/array.h \
  ../include/output/../geometry/../geometry/rectangle.h \
  ../include/output/../geometry/../geometry/grid.h \
  ../include/output/../geometry/../geometry/../declarations.h \
  ../include/output/../geometry/../geometry/../base/array.h \
  ../include/output/../geometry/../geometry/../numeric/arithmetic.h \
  ../include/output/../geometry/../geometry/../geometry/lattice_set.h \
  ../include/output/../geometry/../geometry/../geometry/../base/array.h \
  ../include/output/../geometry/../geometry/../geometry/../base/iterator.h \
  ../include/output/../geometry/../geometry/../geometry/../numeric/interval.h \
  ../include/output/../geometry/../geometry/../geometry/../base/array_operations.h \
  ../include/output/../geometry/../geometry/../geometry/../base/../declarations.h \
  ../include/output/../geometry/../geometry/lattice_set.h \
  ../include/output/../geometry/parallelotope.h \
  ../include/output/../geometry/zonotope.h \
  ../include/output/../geometry/../linear_algebra/vector.h \
  ../include/output/../geometry/../linear_algebra/matrix.h \
  ../include/output/../geometry/../linear_algebra/transformation_system.h \
  ../include/output/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../include/output/../geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../include/output/../geometry/../linear_algebra/../linear_algebra/interval_vector.h \
  ../include/output/../geometry/../linear_algebra/../linear_algebra/interval_matrix.h \
  ../include/output/../geometry/../linear_algebra/../linear_algebra/../linear_algebra/matrix.h \
  ../include/output/../geometry/../linear_algebra/../linear_algebra/../declarations.h \
  ../include/output/../geometry/../linear_algebra/../linear_algebra/../numeric/interval.h \
  ../include/output/../geometry/../numeric/interval.h \
  ../include/output/../geometry/../geometry/point.h \
  ../include/output/../geometry/polyhedron.h \
  ../include/output/../geometry/partition_tree_set.h \
  ../include/output/../linear_algebra/matrix.h \
  ../include/output/../system/affine_map.h \
  ../include/output/../system/../numeric/numerical_types.h \
  ../include/output/../system/../linear_algebra/vector.h \
  ../include/output/../system/../linear_algebra/matrix.h \
  ../include/output/../system/../linear_algebra/interval_matrix.h \
  ../include/output/../system/../geometry/point.h \
  ../include/output/../system/../geometry/simplex.h \
  ../include/output/../system/../geometry/../declarations.h \
  ../include/output/../system/../geometry/../base/array.h \
  ../include/output/../system/../geometry/../geometry/point.h \
  ../include/output/../system/../geometry/rectangle.h \
  ../include/output/../system/../geometry/parallelotope.h \
  ../include/output/../system/../geometry/zonotope.h \
  ../include/output/../system/../geometry/polyhedron.h \
  ../include/output/../system/../geometry/list_set.h \
  ../include/output/../system/../geometry/grid_set.h \
  ../include/output/../system/../system/map.h \
  ../include/output/../system/../system/../declarations.h test.h
	-${CXX} ${CXXFLAGS} -c test_partition_tree_set.cc 
 
test_partition_tree_set: test_partition_tree_set.o
	-${CXX} test_partition_tree_set.o ${LDXXFLAGS} -o test_partition_tree_set 

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
  ../include/geometry/../base/utility.h ../include/output/epsfstream.h \
  ../include/output/../geometry/point.h \
  ../include/output/../geometry/rectangle.h \
  ../include/output/../geometry/list_set.h \
  ../include/output/../geometry/grid_set.h \
  ../include/output/../geometry/../declarations.h \
  ../include/output/../geometry/../base/array.h \
  ../include/output/../geometry/../geometry/rectangle.h \
  ../include/output/../geometry/../geometry/grid.h \
  ../include/output/../geometry/../geometry/../declarations.h \
  ../include/output/../geometry/../geometry/../base/array.h \
  ../include/output/../geometry/../geometry/../numeric/arithmetic.h \
  ../include/output/../geometry/../geometry/../geometry/lattice_set.h \
  ../include/output/../geometry/../geometry/../geometry/../base/array.h \
  ../include/output/../geometry/../geometry/../geometry/../base/iterator.h \
  ../include/output/../geometry/../geometry/../geometry/../base/../base/array.h \
  ../include/output/../geometry/../geometry/../geometry/../numeric/interval.h \
  ../include/output/../geometry/../geometry/../geometry/../base/array_operations.h \
  ../include/output/../geometry/../geometry/../geometry/../base/../declarations.h \
  ../include/output/../geometry/../geometry/lattice_set.h \
  ../include/output/../geometry/parallelotope.h \
  ../include/output/../geometry/../numeric/interval.h \
  ../include/output/../geometry/../linear_algebra/vector.h \
  ../include/output/../geometry/../linear_algebra/matrix.h \
  ../include/output/../geometry/../linear_algebra/../declarations.h \
  ../include/output/../geometry/../linear_algebra/../numeric/integer.h \
  ../include/output/../geometry/../geometry/point.h \
  ../include/output/../geometry/../geometry/polyhedron.h \
  /ufs/collins/include/ppl.hh \
  ../include/output/../geometry/../geometry/../linear_algebra/constraint_system.h \
  ../include/output/../geometry/../geometry/../linear_algebra/../declarations.h \
  ../include/output/../geometry/../geometry/../linear_algebra/../numeric/numerical_types.h \
  ../include/output/../geometry/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../include/output/../geometry/../geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../include/output/../geometry/../geometry/../linear_algebra/vector.h \
  ../include/output/../geometry/../geometry/../linear_algebra/matrix.h \
  ../include/output/../geometry/../geometry/../linear_algebra/interval_vector.h \
  ../include/output/../geometry/../geometry/../linear_algebra/generator_system.h \
  ../include/output/../geometry/../geometry/../geometry/point.h \
  ../include/output/../geometry/../geometry/../geometry/rectangle.h \
  ../include/output/../geometry/zonotope.h \
  ../include/output/../geometry/../linear_algebra/transformation_system.h \
  ../include/output/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../include/output/../geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../include/output/../geometry/../linear_algebra/../linear_algebra/interval_vector.h \
  ../include/output/../geometry/../linear_algebra/../linear_algebra/interval_matrix.h \
  ../include/output/../geometry/../linear_algebra/../linear_algebra/../linear_algebra/matrix.h \
  ../include/output/../geometry/../linear_algebra/../linear_algebra/../declarations.h \
  ../include/output/../geometry/../linear_algebra/../linear_algebra/../numeric/interval.h \
  ../include/output/../geometry/polyhedron.h \
  ../include/output/../geometry/partition_tree_set.h \
  ../include/output/../geometry/../base/binary_word.h \
  ../include/output/../geometry/../base/../declarations.h \
  ../include/output/../geometry/../base/binary_tree.h \
  ../include/output/../geometry/../base/../base/array.h \
  ../include/output/../geometry/../base/../base/binary_word.h \
  ../include/output/../geometry/../base/iterator.h \
  ../include/output/../geometry/../geometry/subdivision_tree_set.h \
  ../include/output/../geometry/../geometry/../base/sequence.h \
  ../include/output/../geometry/../geometry/../base/binary_word.h \
  ../include/output/../geometry/../geometry/../base/binary_tree.h \
  ../include/output/../geometry/../geometry/../base/iterator.h \
  ../include/output/../linear_algebra/matrix.h \
  ../include/output/../system/affine_map.h \
  ../include/output/../system/../numeric/numerical_types.h \
  ../include/output/../system/../linear_algebra/vector.h \
  ../include/output/../system/../linear_algebra/matrix.h \
  ../include/output/../system/../linear_algebra/interval_matrix.h \
  ../include/output/../system/../geometry/point.h \
  ../include/output/../system/../geometry/simplex.h \
  ../include/output/../system/../geometry/../declarations.h \
  ../include/output/../system/../geometry/../base/array.h \
  ../include/output/../system/../geometry/../geometry/point.h \
  ../include/output/../system/../geometry/rectangle.h \
  ../include/output/../system/../geometry/parallelotope.h \
  ../include/output/../system/../geometry/zonotope.h \
  ../include/output/../system/../geometry/polyhedron.h \
  ../include/output/../system/../geometry/list_set.h \
  ../include/output/../system/../geometry/grid_set.h \
  ../include/output/../system/../system/map.h \
  ../include/output/../system/../system/../declarations.h test.h
	-${CXX} ${CXXFLAGS} -c test_epsfstream.cc 
 
test_epsfstream: test_epsfstream.o
	-${CXX} test_epsfstream.o ${LDXXFLAGS} -o test_epsfstream 

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
	-${CXX} test_interval_newton.o ${LDXXFLAGS} -o test_interval_newton 

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
  ../include/geometry/rectangle.h \
  ../include/geometry/../numeric/arithmetic.h \
  ../include/geometry/../numeric/interval.h \
  ../include/geometry/../numeric/../declarations.h \
  ../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../numeric/../numeric/arithmetic.h \
  ../include/geometry/../linear_algebra/interval_vector.h \
  ../include/geometry/../linear_algebra/../numeric/interval.h \
  ../include/geometry/../linear_algebra/../linear_algebra/vector.h \
  ../include/geometry/../geometry/point.h ../include/geometry/grid.h \
  ../include/geometry/../geometry/lattice_set.h \
  ../include/geometry/../geometry/../base/array.h \
  ../include/geometry/../geometry/../base/iterator.h \
  ../include/geometry/../geometry/../base/../base/array.h \
  ../include/geometry/../geometry/../numeric/interval.h \
  ../include/geometry/../geometry/../base/array_operations.h \
  ../include/geometry/../geometry/../base/../declarations.h \
  ../include/geometry/grid_set.h \
  ../include/geometry/../geometry/rectangle.h \
  ../include/geometry/../geometry/grid.h \
  ../include/geometry/partition_tree_set.h \
  ../include/geometry/../base/binary_word.h \
  ../include/geometry/../base/../declarations.h \
  ../include/geometry/../base/binary_tree.h \
  ../include/geometry/../base/../base/array.h \
  ../include/geometry/../base/../base/binary_word.h \
  ../include/geometry/../base/iterator.h \
  ../include/geometry/../geometry/subdivision_tree_set.h \
  ../include/geometry/../geometry/../declarations.h \
  ../include/geometry/../geometry/../base/sequence.h \
  ../include/geometry/../geometry/../base/binary_word.h \
  ../include/geometry/../geometry/../base/binary_tree.h \
  ../include/geometry/../geometry/../geometry/rectangle.h \
  ../include/system/henon_map.h \
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
  ../include/evaluation/apply.h ../include/output/epsfstream.h \
  ../include/output/../geometry/point.h \
  ../include/output/../geometry/rectangle.h \
  ../include/output/../geometry/list_set.h \
  ../include/output/../geometry/../declarations.h \
  ../include/output/../geometry/../base/utility.h \
  ../include/output/../geometry/../base/../declarations.h \
  ../include/output/../geometry/grid_set.h \
  ../include/output/../geometry/parallelotope.h \
  ../include/output/../geometry/zonotope.h \
  ../include/output/../geometry/../linear_algebra/vector.h \
  ../include/output/../geometry/../linear_algebra/matrix.h \
  ../include/output/../geometry/../linear_algebra/transformation_system.h \
  ../include/output/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../include/output/../geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../include/output/../geometry/../linear_algebra/../linear_algebra/interval_vector.h \
  ../include/output/../geometry/../linear_algebra/../linear_algebra/interval_matrix.h \
  ../include/output/../geometry/../numeric/interval.h \
  ../include/output/../geometry/../geometry/point.h \
  ../include/output/../geometry/polyhedron.h \
  ../include/output/../geometry/partition_tree_set.h \
  ../include/output/../linear_algebra/matrix.h \
  ../include/output/../system/affine_map.h \
  ../include/output/../system/../numeric/numerical_types.h \
  ../include/output/../system/../linear_algebra/vector.h \
  ../include/output/../system/../linear_algebra/matrix.h \
  ../include/output/../system/../linear_algebra/interval_matrix.h \
  ../include/output/../system/../geometry/point.h \
  ../include/output/../system/../geometry/simplex.h \
  ../include/output/../system/../geometry/../declarations.h \
  ../include/output/../system/../geometry/../base/array.h \
  ../include/output/../system/../geometry/../geometry/point.h \
  ../include/output/../system/../geometry/rectangle.h \
  ../include/output/../system/../geometry/parallelotope.h \
  ../include/output/../system/../geometry/zonotope.h \
  ../include/output/../system/../geometry/polyhedron.h \
  ../include/output/../system/../geometry/list_set.h \
  ../include/output/../system/../geometry/grid_set.h \
  ../include/output/../system/../system/map.h test.h
	-${CXX} ${CXXFLAGS} -c test_chainreach.cc 
 
test_chainreach: test_chainreach.o
	-${CXX} test_chainreach.o ${LDXXFLAGS} -o test_chainreach 

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
	-${CXX} test_integration_step.o ${LDXXFLAGS} -o test_integration_step 

