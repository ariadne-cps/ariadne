include ../config.mk

TESTS = \
  test_array test_binary_word \
  test_lattice \
  test_mpfr test_arithmetic test_function test_interval \
  test_vector test_matrix test_linear_algebra \
  test_point \
  test_rectangle test_parallelotope test_zonotope \
  test_simplex test_polytope test_polyhedron \
  test_list_set test_grid_set test_partition_tree_set \
  test_epsfstream \
  test_interval_newton \
  test_chainreach \
  test_integration_step 

OBJS=$(addsuffix .o,$(TESTS))

CXXFLAGS += -Wall -I../include/ 
LDLIBS = ../src/libariadne.so -lmpfr -lgmp -lgmpxx -lppl #-lblas

all: ${TESTS}

check: $(TESTS)
	@pass=0; fail=0; \
	for test in $(TESTS); \
	do ./$$test > $$test.log; \
	  result=$$?; \
	  if test $$result -eq 0; then \
	    pass=`expr $$pass + 1`; echo "PASS: $$test"; \
	  else \
	    fail=`expr $$fail + 1`; echo "FAIL ($$result): $$test"; \
	  fi; \
	done; \
	if test $$fail -eq 0; then \
	  echo "Passed all of $$pass tests."; \
	else \
	  echo "Failed $$fail of" `expr $$pass + $$fail` "tests."; \
	fi


$(TESTS): %: %.o ../src/libariadne.so
	$(CXX) $@.o $(LDXXFLAGS) $(LDLIBS) -o $@

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
	@for ff in ${TESTS}; do g++ -MM ${CXXFLAGS} $$ff.cc; \
	echo "	-$$""{CXX} $$""{CXXFLAGS} -c $$ff.cc "; \
	echo ""; done >> Makefile

depclean:
	@sed -e '/^### Do Not edit this line$$/,$$d' Makefile \
		> MMM.$$$$ && mv MMM.$$$$ Makefile
	@echo "### Do Not edit this line" >> Makefile

### Do Not edit this line
### Everything below is auto-generated

test_array.o: test_array.cc ../include/ariadne.h ../include/base/stlio.h \
  ../include/base/../base/array.h \
  ../include/base/../base/../declarations.h ../include/base/array.h \
  ../include/numeric/numerical_types.h ../include/numeric/mpfloat.h \
  /home/collins/include/mpfr.h ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/approximation.h \
  ../include/numeric/../numeric/function.h \
  ../include/numeric/../numeric/rational.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../numeric/arithmetic.h \
  ../include/numeric/../numeric/../numeric/function.h \
  ../include/numeric/../numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../numeric/../numeric/arithmetic.h test.h
	-${CXX} ${CXXFLAGS} -c test_array.cc 

test_binary_word.o: test_binary_word.cc ../include/ariadne.h \
  ../include/base/stlio.h ../include/base/../base/array.h \
  ../include/base/../base/../declarations.h ../include/base/exception.h \
  ../include/combinatoric/binary_word.h \
  ../include/combinatoric/../declarations.h \
  ../include/combinatoric/binary_tree.h \
  ../include/combinatoric/../base/array.h \
  ../include/combinatoric/../combinatoric/binary_word.h \
  ../include/numeric/numerical_types.h ../include/numeric/mpfloat.h \
  /home/collins/include/mpfr.h ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/approximation.h \
  ../include/numeric/../numeric/function.h \
  ../include/numeric/../numeric/rational.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../numeric/arithmetic.h \
  ../include/numeric/../numeric/../numeric/function.h \
  ../include/numeric/../numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../numeric/../numeric/arithmetic.h \
  ../include/base/utility.h ../include/base/../declarations.h \
  ../include/geometry/point.h ../include/geometry/../declarations.h \
  ../include/geometry/../base/array.h \
  ../include/geometry/../linear_algebra/vector.h \
  ../include/geometry/../linear_algebra/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/integer.h \
  ../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../include/geometry/../linear_algebra/../numeric/interval.h \
  ../include/geometry/../linear_algebra/../numeric/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/arithmetic.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/rounding.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/function.h \
  ../include/geometry/../linear_algebra/matrix.h \
  ../include/geometry/rectangle.h \
  ../include/geometry/../numeric/arithmetic.h \
  ../include/geometry/../numeric/function.h \
  ../include/geometry/../numeric/interval.h \
  ../include/geometry/../geometry/ppl_polyhedron.h \
  ../include/geometry/../geometry/../declarations.h \
  ../include/geometry/../geometry/../numeric/rational.h \
  ../include/geometry/../geometry/../linear_algebra/vector.h \
  ../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../include/geometry/../geometry/point.h \
  ../include/geometry/../geometry/rectangle_expression.h test.h
	-${CXX} ${CXXFLAGS} -c test_binary_word.cc 

test_lattice.o: test_lattice.cc ../include/declarations.h \
  ../include/real_typedef.h ../include/numeric/numerical_types.h \
  ../include/numeric/mpfloat.h /home/collins/include/mpfr.h \
  ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/approximation.h \
  ../include/numeric/../numeric/function.h \
  ../include/numeric/../numeric/rational.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../numeric/arithmetic.h \
  ../include/numeric/../numeric/../numeric/function.h \
  ../include/numeric/../numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../numeric/../numeric/arithmetic.h \
  ../include/combinatoric/lattice_set.h \
  ../include/combinatoric/../base/array.h \
  ../include/combinatoric/../base/../declarations.h \
  ../include/combinatoric/../base/iterator.h \
  ../include/combinatoric/../base/../base/array.h \
  ../include/combinatoric/../base/array_operations.h \
  ../include/combinatoric/../numeric/interval.h \
  ../include/combinatoric/../numeric/../declarations.h \
  ../include/combinatoric/../numeric/../numeric/numerical_traits.h \
  ../include/combinatoric/../numeric/../numeric/arithmetic.h \
  ../include/combinatoric/../numeric/../numeric/rounding.h \
  ../include/combinatoric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/combinatoric/../numeric/../numeric/../numeric/function.h \
  ../include/combinatoric/lattice_map.h \
  ../include/combinatoric/../combinatoric/lattice_set.h
	-${CXX} ${CXXFLAGS} -c test_lattice.cc 

test_mpfr.o: test_mpfr.cc /home/collins/include/mpfr.h
	-${CXX} ${CXXFLAGS} -c test_mpfr.cc 

test_arithmetic.o: test_arithmetic.cc /home/collins/include/mpfr.h \
  ../include/numeric/arithmetic.h ../include/numeric/interval.h \
  ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/approximation.h \
  ../include/numeric/../numeric/arithmetic.h \
  ../include/numeric/../numeric/rounding.h \
  ../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../numeric/function.h \
  ../include/numeric/float64.h ../include/numeric/../numeric/function.h \
  ../include/numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/arithmetic.h \
  ../include/numeric/mpfloat.h ../include/numeric/../numeric/rational.h \
  ../include/numeric/../numeric/../numeric/integer.h \
  ../include/numeric/rational.h test.h
	-${CXX} ${CXXFLAGS} -c test_arithmetic.cc 

test_function.o: test_function.cc ../include/numeric/function.h \
  ../include/numeric/float64.h /home/collins/include/mpfr.h \
  ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/approximation.h \
  ../include/numeric/../numeric/function.h \
  ../include/numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../numeric/arithmetic.h \
  ../include/numeric/mpfloat.h ../include/numeric/../numeric/rational.h \
  ../include/numeric/../numeric/../numeric/function.h \
  ../include/numeric/../numeric/../numeric/integer.h test.h
	-${CXX} ${CXXFLAGS} -c test_function.cc 

test_interval.o: test_interval.cc ../include/numeric/interval.h \
  ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/approximation.h \
  ../include/numeric/../numeric/arithmetic.h \
  ../include/numeric/../numeric/rounding.h \
  ../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../numeric/function.h \
  ../include/numeric/float64.h /home/collins/include/mpfr.h \
  ../include/numeric/../numeric/function.h \
  ../include/numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/arithmetic.h \
  ../include/numeric/mpfloat.h ../include/numeric/../numeric/rational.h \
  ../include/numeric/../numeric/../numeric/integer.h \
  ../include/numeric/rational.h test.h
	-${CXX} ${CXXFLAGS} -c test_interval.cc 

test_vector.o: test_vector.cc ../include/real_typedef.h \
  ../include/numeric/numerical_types.h ../include/numeric/mpfloat.h \
  /home/collins/include/mpfr.h ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/approximation.h \
  ../include/numeric/../numeric/function.h \
  ../include/numeric/../numeric/rational.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../numeric/arithmetic.h \
  ../include/numeric/../numeric/../numeric/function.h \
  ../include/numeric/../numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../numeric/../numeric/arithmetic.h \
  ../include/numeric/numerical_types.h ../include/numeric/interval.h \
  ../include/numeric/../numeric/arithmetic.h \
  ../include/numeric/../numeric/rounding.h \
  ../include/linear_algebra/vector.h \
  ../include/linear_algebra/../declarations.h \
  ../include/linear_algebra/../numeric/integer.h \
  ../include/linear_algebra/../numeric/arithmetic.h \
  ../include/linear_algebra/../numeric/interval.h \
  ../include/linear_algebra/vector.tpl ../include/linear_algebra/vector.h \
  ../include/linear_algebra/../base/stlio.h \
  ../include/linear_algebra/../base/../base/array.h \
  ../include/linear_algebra/../base/../base/../declarations.h
	-${CXX} ${CXXFLAGS} -c test_vector.cc 

test_matrix.o: test_matrix.cc ../include/declarations.h \
  ../include/real_typedef.h ../include/numeric/numerical_types.h \
  ../include/numeric/mpfloat.h /home/collins/include/mpfr.h \
  ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/approximation.h \
  ../include/numeric/../numeric/function.h \
  ../include/numeric/../numeric/rational.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../numeric/arithmetic.h \
  ../include/numeric/../numeric/../numeric/function.h \
  ../include/numeric/../numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../numeric/../numeric/arithmetic.h \
  ../include/numeric/numerical_types.h ../include/linear_algebra/vector.h \
  ../include/linear_algebra/../declarations.h \
  ../include/linear_algebra/../numeric/integer.h \
  ../include/linear_algebra/../numeric/arithmetic.h \
  ../include/linear_algebra/../numeric/interval.h \
  ../include/linear_algebra/../numeric/../declarations.h \
  ../include/linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../include/linear_algebra/../numeric/../numeric/arithmetic.h \
  ../include/linear_algebra/../numeric/../numeric/rounding.h \
  ../include/linear_algebra/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/linear_algebra/../numeric/../numeric/../numeric/function.h \
  ../include/linear_algebra/vector.tpl ../include/linear_algebra/vector.h \
  ../include/linear_algebra/../base/stlio.h \
  ../include/linear_algebra/../base/../base/array.h \
  ../include/linear_algebra/../base/../base/../declarations.h \
  ../include/linear_algebra/matrix.h ../include/linear_algebra/matrix.tpl \
  ../include/linear_algebra/matrix.h \
  ../include/linear_algebra/../base/array.h \
  ../include/linear_algebra/../linear_algebra/vector.h \
  ../include/linear_algebra/../linear_algebra/matrix.h \
  ../include/linear_algebra/../linear_algebra/lu_matrix.h \
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
  ../include/linear_algebra/../linear_algebra/../declarations.h \
  ../include/linear_algebra/../linear_algebra/../base/array.h \
  ../include/linear_algebra/../linear_algebra/../linear_algebra/matrix.h
	-${CXX} ${CXXFLAGS} -c test_matrix.cc 

test_linear_algebra.o: test_linear_algebra.cc ../include/real_typedef.h \
  ../include/numeric/numerical_types.h ../include/numeric/mpfloat.h \
  /home/collins/include/mpfr.h ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/approximation.h \
  ../include/numeric/../numeric/function.h \
  ../include/numeric/../numeric/rational.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../numeric/arithmetic.h \
  ../include/numeric/../numeric/../numeric/function.h \
  ../include/numeric/../numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../numeric/../numeric/arithmetic.h \
  ../include/numeric/numerical_types.h ../include/linear_algebra/vector.h \
  ../include/linear_algebra/../declarations.h \
  ../include/linear_algebra/../numeric/integer.h \
  ../include/linear_algebra/../numeric/arithmetic.h \
  ../include/linear_algebra/../numeric/interval.h \
  ../include/linear_algebra/../numeric/../declarations.h \
  ../include/linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../include/linear_algebra/../numeric/../numeric/arithmetic.h \
  ../include/linear_algebra/../numeric/../numeric/rounding.h \
  ../include/linear_algebra/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/linear_algebra/../numeric/../numeric/../numeric/function.h \
  ../include/linear_algebra/matrix.h \
  ../include/linear_algebra/lu_matrix.h \
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
  ../include/linear_algebra/../base/array.h \
  ../include/linear_algebra/../base/../declarations.h \
  ../include/linear_algebra/../linear_algebra/matrix.h \
  ../include/linear_algebra/qr_matrix.h \
  /home/collins/include/tblas/gecpy.hpp \
  /home/collins/include/tblas/trset.hpp \
  /home/collins/include/tlapack/larf.hpp \
  /home/collins/include/tlapack/larfg.hpp \
  /home/collins/include/tlapack/geqrf.hpp \
  /home/collins/include/tlapack/orgqr.hpp \
  ../include/linear_algebra/svd_matrix.h \
  /home/collins/include/tblas/copy.hpp \
  /home/collins/include/tlapack/gesvd.hpp \
  ../include/linear_algebra/../linear_algebra/vector.h test.h
	-${CXX} ${CXXFLAGS} -c test_linear_algebra.cc 

test_point.o: test_point.cc ../include/real_typedef.h \
  ../include/numeric/numerical_types.h ../include/numeric/mpfloat.h \
  /home/collins/include/mpfr.h ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/approximation.h \
  ../include/numeric/../numeric/function.h \
  ../include/numeric/../numeric/rational.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../numeric/arithmetic.h \
  ../include/numeric/../numeric/../numeric/function.h \
  ../include/numeric/../numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../numeric/../numeric/arithmetic.h \
  ../include/geometry/point.h ../include/geometry/../declarations.h \
  ../include/geometry/../base/array.h \
  ../include/geometry/../base/../declarations.h \
  ../include/geometry/../linear_algebra/vector.h \
  ../include/geometry/../linear_algebra/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/integer.h \
  ../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../include/geometry/../linear_algebra/../numeric/interval.h \
  ../include/geometry/../linear_algebra/../numeric/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/arithmetic.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/rounding.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/function.h \
  ../include/geometry/../linear_algebra/matrix.h test.h
	-${CXX} ${CXXFLAGS} -c test_point.cc 

test_rectangle.o: test_rectangle.cc ../include/ariadne.h \
  ../include/real_typedef.h ../include/numeric/numerical_types.h \
  ../include/numeric/mpfloat.h /home/collins/include/mpfr.h \
  ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/approximation.h \
  ../include/numeric/../numeric/function.h \
  ../include/numeric/../numeric/rational.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../numeric/arithmetic.h \
  ../include/numeric/../numeric/../numeric/function.h \
  ../include/numeric/../numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../numeric/../numeric/arithmetic.h \
  ../include/base/exception.h ../include/base/utility.h \
  ../include/base/../declarations.h ../include/geometry/rectangle.h \
  ../include/geometry/../declarations.h \
  ../include/geometry/../numeric/arithmetic.h \
  ../include/geometry/../numeric/function.h \
  ../include/geometry/../numeric/interval.h \
  ../include/geometry/../numeric/../declarations.h \
  ../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../numeric/../numeric/arithmetic.h \
  ../include/geometry/../numeric/../numeric/rounding.h \
  ../include/geometry/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../numeric/../numeric/../numeric/function.h \
  ../include/geometry/../linear_algebra/vector.h \
  ../include/geometry/../linear_algebra/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/integer.h \
  ../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../include/geometry/../linear_algebra/../numeric/interval.h \
  ../include/geometry/../geometry/ppl_polyhedron.h \
  ../include/geometry/../geometry/../declarations.h \
  ../include/geometry/../geometry/../numeric/rational.h \
  ../include/geometry/../geometry/../linear_algebra/vector.h \
  ../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../include/geometry/../geometry/../linear_algebra/../declarations.h \
  ../include/geometry/../geometry/../linear_algebra/../numeric/integer.h \
  ../include/geometry/../geometry/../linear_algebra/../numeric/interval.h \
  ../include/geometry/../geometry/point.h \
  ../include/geometry/../geometry/../base/array.h \
  ../include/geometry/../geometry/../base/../declarations.h \
  ../include/geometry/../geometry/rectangle_expression.h \
  ../include/geometry/list_set.h ../include/geometry/../base/utility.h \
  test.h
	-${CXX} ${CXXFLAGS} -c test_rectangle.cc 

test_parallelotope.o: test_parallelotope.cc ../include/ariadne.h \
  ../include/real_typedef.h ../include/numeric/numerical_types.h \
  ../include/numeric/mpfloat.h /home/collins/include/mpfr.h \
  ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/approximation.h \
  ../include/numeric/../numeric/function.h \
  ../include/numeric/../numeric/rational.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../numeric/arithmetic.h \
  ../include/numeric/../numeric/../numeric/function.h \
  ../include/numeric/../numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../numeric/../numeric/arithmetic.h \
  ../include/base/exception.h ../include/base/utility.h \
  ../include/base/../declarations.h ../include/numeric/numerical_types.h \
  ../include/geometry/point.h ../include/geometry/../declarations.h \
  ../include/geometry/../base/array.h \
  ../include/geometry/../base/../declarations.h \
  ../include/geometry/../linear_algebra/vector.h \
  ../include/geometry/../linear_algebra/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/integer.h \
  ../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../include/geometry/../linear_algebra/../numeric/interval.h \
  ../include/geometry/../linear_algebra/../numeric/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/arithmetic.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/rounding.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/function.h \
  ../include/geometry/../linear_algebra/matrix.h \
  ../include/geometry/rectangle.h \
  ../include/geometry/../numeric/arithmetic.h \
  ../include/geometry/../numeric/function.h \
  ../include/geometry/../numeric/interval.h \
  ../include/geometry/../geometry/ppl_polyhedron.h \
  ../include/geometry/../geometry/../declarations.h \
  ../include/geometry/../geometry/../numeric/rational.h \
  ../include/geometry/../geometry/../linear_algebra/vector.h \
  ../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../include/geometry/../geometry/point.h \
  ../include/geometry/../geometry/rectangle_expression.h \
  ../include/geometry/parallelotope.h \
  ../include/geometry/../geometry/zonotope.h \
  ../include/geometry/../geometry/../numeric/interval.h \
  ../include/geometry/../geometry/../geometry/ppl_polyhedron.h \
  ../include/geometry/../geometry/../geometry/point.h
	-${CXX} ${CXXFLAGS} -c test_parallelotope.cc 

test_zonotope.o: test_zonotope.cc ../include/ariadne.h \
  ../include/real_typedef.h ../include/numeric/numerical_types.h \
  ../include/numeric/mpfloat.h /home/collins/include/mpfr.h \
  ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/approximation.h \
  ../include/numeric/../numeric/function.h \
  ../include/numeric/../numeric/rational.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../numeric/arithmetic.h \
  ../include/numeric/../numeric/../numeric/function.h \
  ../include/numeric/../numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../numeric/../numeric/arithmetic.h \
  ../include/base/exception.h ../include/base/utility.h \
  ../include/base/../declarations.h ../include/numeric/numerical_types.h \
  ../include/geometry/point.h ../include/geometry/../declarations.h \
  ../include/geometry/../base/array.h \
  ../include/geometry/../base/../declarations.h \
  ../include/geometry/../linear_algebra/vector.h \
  ../include/geometry/../linear_algebra/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/integer.h \
  ../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../include/geometry/../linear_algebra/../numeric/interval.h \
  ../include/geometry/../linear_algebra/../numeric/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/arithmetic.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/rounding.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/function.h \
  ../include/geometry/../linear_algebra/matrix.h \
  ../include/geometry/rectangle.h \
  ../include/geometry/../numeric/arithmetic.h \
  ../include/geometry/../numeric/function.h \
  ../include/geometry/../numeric/interval.h \
  ../include/geometry/../geometry/ppl_polyhedron.h \
  ../include/geometry/../geometry/../declarations.h \
  ../include/geometry/../geometry/../numeric/rational.h \
  ../include/geometry/../geometry/../linear_algebra/vector.h \
  ../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../include/geometry/../geometry/point.h \
  ../include/geometry/../geometry/rectangle_expression.h \
  ../include/geometry/parallelotope.h \
  ../include/geometry/../geometry/zonotope.h \
  ../include/geometry/../geometry/../numeric/interval.h \
  ../include/geometry/../geometry/../geometry/ppl_polyhedron.h \
  ../include/geometry/../geometry/../geometry/point.h \
  ../include/geometry/zonotope.h ../include/output/epsfstream.h \
  ../include/output/../numeric/approximation.h \
  ../include/output/../numeric/float64.h \
  ../include/output/../numeric/../declarations.h \
  ../include/output/../numeric/../numeric/numerical_traits.h \
  ../include/output/../numeric/../numeric/function.h \
  ../include/output/../numeric/../numeric/integer.h \
  ../include/output/../geometry/point.h \
  ../include/output/../geometry/rectangle.h \
  ../include/output/../geometry/list_set.h \
  ../include/output/../geometry/../declarations.h \
  ../include/output/../geometry/../base/utility.h \
  ../include/output/../geometry/grid_set.h \
  ../include/output/../geometry/../base/array.h \
  ../include/output/../geometry/../combinatoric/lattice_set.h \
  ../include/output/../geometry/../combinatoric/../base/array.h \
  ../include/output/../geometry/../combinatoric/../base/iterator.h \
  ../include/output/../geometry/../combinatoric/../base/../base/array.h \
  ../include/output/../geometry/../combinatoric/../base/array_operations.h \
  ../include/output/../geometry/../combinatoric/../base/../declarations.h \
  ../include/output/../geometry/../combinatoric/../numeric/interval.h \
  ../include/output/../geometry/../geometry/rectangle.h \
  ../include/output/../geometry/../geometry/grid.h \
  ../include/output/../geometry/../geometry/../declarations.h \
  ../include/output/../geometry/../geometry/../base/array.h \
  ../include/output/../geometry/../geometry/../numeric/arithmetic.h \
  ../include/output/../geometry/../geometry/../combinatoric/lattice_set.h \
  ../include/output/../geometry/parallelotope.h \
  ../include/output/../geometry/zonotope.h \
  ../include/output/../geometry/polytope.h \
  ../include/output/../geometry/../linear_algebra/vector.h \
  ../include/output/../geometry/../linear_algebra/matrix.h \
  ../include/output/../geometry/../geometry/ppl_polyhedron.h \
  ../include/output/../geometry/../geometry/point.h \
  ../include/output/../geometry/../geometry/point_list.h \
  ../include/output/../geometry/../geometry/../linear_algebra/vector.h \
  ../include/output/../geometry/../geometry/../linear_algebra/matrix.h \
  ../include/output/../geometry/../geometry/../geometry/point.h \
  ../include/output/../geometry/partition_tree_set.h \
  ../include/output/../geometry/../base/iterator.h \
  ../include/output/../geometry/../combinatoric/binary_word.h \
  ../include/output/../geometry/../combinatoric/../declarations.h \
  ../include/output/../geometry/../combinatoric/binary_tree.h \
  ../include/output/../geometry/../combinatoric/../combinatoric/binary_word.h \
  ../include/output/../geometry/../combinatoric/subdivision_tree_set.h \
  ../include/output/../geometry/../combinatoric/../base/sequence.h \
  ../include/output/../geometry/../combinatoric/../combinatoric/binary_tree.h \
  ../include/output/../geometry/../combinatoric/../geometry/rectangle.h \
  ../include/output/../linear_algebra/matrix.h \
  ../include/output/../system/affine_map.h \
  ../include/output/../system/../declarations.h \
  ../include/output/../system/../linear_algebra/vector.h \
  ../include/output/../system/../linear_algebra/matrix.h \
  ../include/output/../system/../system/map.h \
  ../include/output/../system/../system/../declarations.h test.h
	-${CXX} ${CXXFLAGS} -c test_zonotope.cc 

test_simplex.o: test_simplex.cc ../include/ariadne.h \
  ../include/real_typedef.h ../include/numeric/numerical_types.h \
  ../include/numeric/mpfloat.h /home/collins/include/mpfr.h \
  ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/approximation.h \
  ../include/numeric/../numeric/function.h \
  ../include/numeric/../numeric/rational.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../numeric/arithmetic.h \
  ../include/numeric/../numeric/../numeric/function.h \
  ../include/numeric/../numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../numeric/../numeric/arithmetic.h \
  ../include/base/exception.h ../include/base/utility.h \
  ../include/base/../declarations.h ../include/numeric/numerical_types.h \
  ../include/geometry/point.h ../include/geometry/../declarations.h \
  ../include/geometry/../base/array.h \
  ../include/geometry/../base/../declarations.h \
  ../include/geometry/../linear_algebra/vector.h \
  ../include/geometry/../linear_algebra/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/integer.h \
  ../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../include/geometry/../linear_algebra/../numeric/interval.h \
  ../include/geometry/../linear_algebra/../numeric/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/arithmetic.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/rounding.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/function.h \
  ../include/geometry/../linear_algebra/matrix.h \
  ../include/geometry/simplex.h ../include/geometry/../geometry/point.h \
  ../include/geometry/../geometry/point_list.h \
  ../include/geometry/../geometry/../declarations.h \
  ../include/geometry/../geometry/../linear_algebra/vector.h \
  ../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../include/geometry/../geometry/../geometry/point.h \
  ../include/geometry/../geometry/polytope.h \
  ../include/geometry/../geometry/../geometry/ppl_polyhedron.h \
  ../include/geometry/../geometry/../geometry/../declarations.h \
  ../include/geometry/../geometry/../geometry/../numeric/rational.h \
  ../include/geometry/../geometry/../geometry/../linear_algebra/vector.h \
  ../include/geometry/../geometry/../geometry/../linear_algebra/matrix.h \
  ../include/geometry/../geometry/../geometry/point_list.h \
  ../include/geometry/../geometry/../geometry/rectangle.h \
  ../include/geometry/../geometry/../geometry/../numeric/arithmetic.h \
  ../include/geometry/../geometry/../geometry/../numeric/function.h \
  ../include/geometry/../geometry/../geometry/../numeric/interval.h \
  ../include/geometry/../geometry/../geometry/../geometry/ppl_polyhedron.h \
  ../include/geometry/../geometry/../geometry/../geometry/point.h \
  ../include/geometry/../geometry/../geometry/../geometry/rectangle_expression.h \
  ../include/geometry/../geometry/../geometry/../geometry/../declarations.h \
  ../include/geometry/polyhedron.h \
  ../include/geometry/../geometry/ppl_polyhedron.h \
  ../include/geometry/../geometry/rectangle.h test.h
	-${CXX} ${CXXFLAGS} -c test_simplex.cc 

test_polytope.o: test_polytope.cc /home/collins/include/ppl.hh \
  ../include/ariadne.h ../include/base/exception.h \
  ../include/base/utility.h ../include/base/../declarations.h \
  ../include/numeric/float64.h /home/collins/include/mpfr.h \
  ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/approximation.h \
  ../include/numeric/../numeric/function.h \
  ../include/numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../numeric/arithmetic.h \
  ../include/numeric/mpfloat.h ../include/numeric/../numeric/rational.h \
  ../include/numeric/../numeric/../numeric/function.h \
  ../include/numeric/../numeric/../numeric/integer.h \
  ../include/numeric/rational.h ../include/geometry/point.h \
  ../include/geometry/../declarations.h \
  ../include/geometry/../base/array.h \
  ../include/geometry/../base/../declarations.h \
  ../include/geometry/../linear_algebra/vector.h \
  ../include/geometry/../linear_algebra/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/integer.h \
  ../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../include/geometry/../linear_algebra/../numeric/interval.h \
  ../include/geometry/../linear_algebra/../numeric/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/arithmetic.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/rounding.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/function.h \
  ../include/geometry/../linear_algebra/matrix.h \
  ../include/geometry/point_list.h \
  ../include/geometry/../geometry/point.h ../include/geometry/polytope.h \
  ../include/geometry/../geometry/ppl_polyhedron.h \
  ../include/geometry/../geometry/../declarations.h \
  ../include/geometry/../geometry/../numeric/rational.h \
  ../include/geometry/../geometry/../linear_algebra/vector.h \
  ../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../include/geometry/../geometry/point_list.h \
  ../include/geometry/../geometry/rectangle.h \
  ../include/geometry/../geometry/../numeric/arithmetic.h \
  ../include/geometry/../geometry/../numeric/function.h \
  ../include/geometry/../geometry/../numeric/interval.h \
  ../include/geometry/../geometry/../geometry/ppl_polyhedron.h \
  ../include/geometry/../geometry/../geometry/point.h \
  ../include/geometry/../geometry/../geometry/rectangle_expression.h \
  ../include/geometry/../geometry/../geometry/../declarations.h \
  ../include/geometry/polytope.tpl ../include/geometry/polytope.h \
  ../include/geometry/../base/stlio.h \
  ../include/geometry/../base/../base/array.h \
  ../include/geometry/../numeric/interval.h \
  ../include/geometry/../numeric/arithmetic.h test.h
	-${CXX} ${CXXFLAGS} -c test_polytope.cc 

test_polyhedron.o: test_polyhedron.cc /home/collins/include/ppl.hh \
  ../include/ariadne.h ../include/base/exception.h \
  ../include/base/utility.h ../include/base/../declarations.h \
  ../include/numeric/float64.h /home/collins/include/mpfr.h \
  ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/approximation.h \
  ../include/numeric/../numeric/function.h \
  ../include/numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../numeric/arithmetic.h \
  ../include/numeric/mpfloat.h ../include/numeric/../numeric/rational.h \
  ../include/numeric/../numeric/../numeric/function.h \
  ../include/numeric/../numeric/../numeric/integer.h \
  ../include/numeric/rational.h ../include/geometry/point.h \
  ../include/geometry/../declarations.h \
  ../include/geometry/../base/array.h \
  ../include/geometry/../base/../declarations.h \
  ../include/geometry/../linear_algebra/vector.h \
  ../include/geometry/../linear_algebra/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/integer.h \
  ../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../include/geometry/../linear_algebra/../numeric/interval.h \
  ../include/geometry/../linear_algebra/../numeric/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/arithmetic.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/rounding.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/function.h \
  ../include/geometry/../linear_algebra/matrix.h \
  ../include/geometry/point_list.h \
  ../include/geometry/../geometry/point.h \
  ../include/geometry/polyhedron.h \
  ../include/geometry/../geometry/ppl_polyhedron.h \
  ../include/geometry/../geometry/../declarations.h \
  ../include/geometry/../geometry/../numeric/rational.h \
  ../include/geometry/../geometry/../linear_algebra/vector.h \
  ../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../include/geometry/../geometry/point_list.h \
  ../include/geometry/../geometry/rectangle.h \
  ../include/geometry/../geometry/../numeric/arithmetic.h \
  ../include/geometry/../geometry/../numeric/function.h \
  ../include/geometry/../geometry/../numeric/interval.h \
  ../include/geometry/../geometry/../geometry/ppl_polyhedron.h \
  ../include/geometry/../geometry/../geometry/point.h \
  ../include/geometry/../geometry/../geometry/rectangle_expression.h \
  ../include/geometry/../geometry/../geometry/../declarations.h \
  ../include/geometry/polyhedron.tpl ../include/geometry/polyhedron.h \
  ../include/geometry/../base/stlio.h \
  ../include/geometry/../base/../base/array.h \
  ../include/geometry/../numeric/interval.h \
  ../include/geometry/../numeric/arithmetic.h test.h
	-${CXX} ${CXXFLAGS} -c test_polyhedron.cc 

test_list_set.o: test_list_set.cc ../include/ariadne.h \
  ../include/real_typedef.h ../include/numeric/numerical_types.h \
  ../include/numeric/mpfloat.h /home/collins/include/mpfr.h \
  ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/approximation.h \
  ../include/numeric/../numeric/function.h \
  ../include/numeric/../numeric/rational.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../numeric/arithmetic.h \
  ../include/numeric/../numeric/../numeric/function.h \
  ../include/numeric/../numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../numeric/../numeric/arithmetic.h \
  ../include/base/exception.h ../include/base/utility.h \
  ../include/base/../declarations.h ../include/numeric/numerical_types.h \
  ../include/geometry/point.h ../include/geometry/../declarations.h \
  ../include/geometry/../base/array.h \
  ../include/geometry/../base/../declarations.h \
  ../include/geometry/../linear_algebra/vector.h \
  ../include/geometry/../linear_algebra/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/integer.h \
  ../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../include/geometry/../linear_algebra/../numeric/interval.h \
  ../include/geometry/../linear_algebra/../numeric/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/arithmetic.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/rounding.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/function.h \
  ../include/geometry/../linear_algebra/matrix.h \
  ../include/geometry/rectangle.h \
  ../include/geometry/../numeric/arithmetic.h \
  ../include/geometry/../numeric/function.h \
  ../include/geometry/../numeric/interval.h \
  ../include/geometry/../geometry/ppl_polyhedron.h \
  ../include/geometry/../geometry/../declarations.h \
  ../include/geometry/../geometry/../numeric/rational.h \
  ../include/geometry/../geometry/../linear_algebra/vector.h \
  ../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../include/geometry/../geometry/point.h \
  ../include/geometry/../geometry/rectangle_expression.h \
  ../include/geometry/list_set.h ../include/geometry/../base/utility.h \
  test.h
	-${CXX} ${CXXFLAGS} -c test_list_set.cc 

test_grid_set.o: test_grid_set.cc ../include/ariadne.h \
  ../include/real_typedef.h ../include/numeric/numerical_types.h \
  ../include/numeric/mpfloat.h /home/collins/include/mpfr.h \
  ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/approximation.h \
  ../include/numeric/../numeric/function.h \
  ../include/numeric/../numeric/rational.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../numeric/arithmetic.h \
  ../include/numeric/../numeric/../numeric/function.h \
  ../include/numeric/../numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../numeric/../numeric/arithmetic.h \
  ../include/base/exception.h ../include/base/utility.h \
  ../include/base/../declarations.h ../include/combinatoric/binary_word.h \
  ../include/combinatoric/../declarations.h \
  ../include/numeric/numerical_types.h ../include/geometry/point.h \
  ../include/geometry/../declarations.h \
  ../include/geometry/../base/array.h \
  ../include/geometry/../base/../declarations.h \
  ../include/geometry/../linear_algebra/vector.h \
  ../include/geometry/../linear_algebra/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/integer.h \
  ../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../include/geometry/../linear_algebra/../numeric/interval.h \
  ../include/geometry/../linear_algebra/../numeric/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/arithmetic.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/rounding.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/function.h \
  ../include/geometry/../linear_algebra/matrix.h \
  ../include/geometry/rectangle.h \
  ../include/geometry/../numeric/arithmetic.h \
  ../include/geometry/../numeric/function.h \
  ../include/geometry/../numeric/interval.h \
  ../include/geometry/../geometry/ppl_polyhedron.h \
  ../include/geometry/../geometry/../declarations.h \
  ../include/geometry/../geometry/../numeric/rational.h \
  ../include/geometry/../geometry/../linear_algebra/vector.h \
  ../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../include/geometry/../geometry/point.h \
  ../include/geometry/../geometry/rectangle_expression.h \
  ../include/geometry/grid_set.h \
  ../include/geometry/../combinatoric/lattice_set.h \
  ../include/geometry/../combinatoric/../base/array.h \
  ../include/geometry/../combinatoric/../base/iterator.h \
  ../include/geometry/../combinatoric/../base/../base/array.h \
  ../include/geometry/../combinatoric/../base/array_operations.h \
  ../include/geometry/../combinatoric/../base/../declarations.h \
  ../include/geometry/../combinatoric/../numeric/interval.h \
  ../include/geometry/../geometry/rectangle.h \
  ../include/geometry/../geometry/grid.h \
  ../include/geometry/../geometry/../base/array.h \
  ../include/geometry/../geometry/../numeric/arithmetic.h \
  ../include/geometry/../geometry/../combinatoric/lattice_set.h \
  ../include/geometry/list_set.h ../include/geometry/../base/utility.h \
  test.h
	-${CXX} ${CXXFLAGS} -c test_grid_set.cc 

test_partition_tree_set.o: test_partition_tree_set.cc \
  ../include/ariadne.h ../include/real_typedef.h \
  ../include/numeric/numerical_types.h ../include/numeric/mpfloat.h \
  /home/collins/include/mpfr.h ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/approximation.h \
  ../include/numeric/../numeric/function.h \
  ../include/numeric/../numeric/rational.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../numeric/arithmetic.h \
  ../include/numeric/../numeric/../numeric/function.h \
  ../include/numeric/../numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../numeric/../numeric/arithmetic.h \
  ../include/base/exception.h ../include/base/stlio.h \
  ../include/base/../base/array.h \
  ../include/base/../base/../declarations.h ../include/base/utility.h \
  ../include/base/../declarations.h ../include/combinatoric/binary_word.h \
  ../include/combinatoric/../declarations.h \
  ../include/numeric/numerical_types.h ../include/linear_algebra/matrix.h \
  ../include/linear_algebra/../declarations.h \
  ../include/linear_algebra/../numeric/integer.h \
  ../include/linear_algebra/../numeric/interval.h \
  ../include/linear_algebra/../numeric/../declarations.h \
  ../include/linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../include/linear_algebra/../numeric/../numeric/arithmetic.h \
  ../include/linear_algebra/../numeric/../numeric/rounding.h \
  ../include/linear_algebra/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/linear_algebra/../numeric/../numeric/../numeric/function.h \
  ../include/linear_algebra/vector.h \
  ../include/linear_algebra/../numeric/arithmetic.h \
  ../include/geometry/point.h ../include/geometry/../declarations.h \
  ../include/geometry/../base/array.h \
  ../include/geometry/../linear_algebra/vector.h \
  ../include/geometry/../linear_algebra/matrix.h \
  ../include/geometry/rectangle.h \
  ../include/geometry/../numeric/arithmetic.h \
  ../include/geometry/../numeric/function.h \
  ../include/geometry/../numeric/interval.h \
  ../include/geometry/../geometry/ppl_polyhedron.h \
  ../include/geometry/../geometry/../declarations.h \
  ../include/geometry/../geometry/../numeric/rational.h \
  ../include/geometry/../geometry/../linear_algebra/vector.h \
  ../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../include/geometry/../geometry/point.h \
  ../include/geometry/../geometry/rectangle_expression.h \
  ../include/geometry/parallelotope.h \
  ../include/geometry/../geometry/zonotope.h \
  ../include/geometry/../geometry/../numeric/interval.h \
  ../include/geometry/../geometry/../geometry/ppl_polyhedron.h \
  ../include/geometry/../geometry/../geometry/point.h \
  ../include/geometry/list_set.h ../include/geometry/../base/utility.h \
  ../include/geometry/partition_tree_set.h \
  ../include/geometry/../base/iterator.h \
  ../include/geometry/../base/../base/array.h \
  ../include/geometry/../combinatoric/binary_word.h \
  ../include/geometry/../combinatoric/binary_tree.h \
  ../include/geometry/../combinatoric/../declarations.h \
  ../include/geometry/../combinatoric/../base/array.h \
  ../include/geometry/../combinatoric/../combinatoric/binary_word.h \
  ../include/geometry/../combinatoric/subdivision_tree_set.h \
  ../include/geometry/../combinatoric/../base/sequence.h \
  ../include/geometry/../combinatoric/../base/iterator.h \
  ../include/geometry/../combinatoric/../combinatoric/binary_tree.h \
  ../include/geometry/../combinatoric/../geometry/rectangle.h \
  ../include/output/epsfstream.h \
  ../include/output/../numeric/approximation.h \
  ../include/output/../numeric/float64.h \
  ../include/output/../numeric/../declarations.h \
  ../include/output/../numeric/../numeric/numerical_traits.h \
  ../include/output/../numeric/../numeric/function.h \
  ../include/output/../numeric/../numeric/integer.h \
  ../include/output/../geometry/point.h \
  ../include/output/../geometry/rectangle.h \
  ../include/output/../geometry/list_set.h \
  ../include/output/../geometry/grid_set.h \
  ../include/output/../geometry/../declarations.h \
  ../include/output/../geometry/../base/array.h \
  ../include/output/../geometry/../combinatoric/lattice_set.h \
  ../include/output/../geometry/../combinatoric/../base/array.h \
  ../include/output/../geometry/../combinatoric/../base/iterator.h \
  ../include/output/../geometry/../combinatoric/../base/array_operations.h \
  ../include/output/../geometry/../combinatoric/../base/../declarations.h \
  ../include/output/../geometry/../combinatoric/../numeric/interval.h \
  ../include/output/../geometry/../geometry/rectangle.h \
  ../include/output/../geometry/../geometry/grid.h \
  ../include/output/../geometry/../geometry/../declarations.h \
  ../include/output/../geometry/../geometry/../base/array.h \
  ../include/output/../geometry/../geometry/../numeric/arithmetic.h \
  ../include/output/../geometry/../geometry/../combinatoric/lattice_set.h \
  ../include/output/../geometry/parallelotope.h \
  ../include/output/../geometry/zonotope.h \
  ../include/output/../geometry/polytope.h \
  ../include/output/../geometry/../linear_algebra/vector.h \
  ../include/output/../geometry/../linear_algebra/matrix.h \
  ../include/output/../geometry/../geometry/ppl_polyhedron.h \
  ../include/output/../geometry/../geometry/point.h \
  ../include/output/../geometry/../geometry/point_list.h \
  ../include/output/../geometry/../geometry/../linear_algebra/vector.h \
  ../include/output/../geometry/../geometry/../linear_algebra/matrix.h \
  ../include/output/../geometry/../geometry/../geometry/point.h \
  ../include/output/../geometry/partition_tree_set.h \
  ../include/output/../linear_algebra/matrix.h \
  ../include/output/../system/affine_map.h \
  ../include/output/../system/../declarations.h \
  ../include/output/../system/../linear_algebra/vector.h \
  ../include/output/../system/../linear_algebra/matrix.h \
  ../include/output/../system/../system/map.h \
  ../include/output/../system/../system/../declarations.h test.h
	-${CXX} ${CXXFLAGS} -c test_partition_tree_set.cc 

test_epsfstream.o: test_epsfstream.cc ../include/ariadne.h \
  ../include/real_typedef.h ../include/numeric/numerical_types.h \
  ../include/numeric/mpfloat.h /home/collins/include/mpfr.h \
  ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/approximation.h \
  ../include/numeric/../numeric/function.h \
  ../include/numeric/../numeric/rational.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../numeric/arithmetic.h \
  ../include/numeric/../numeric/../numeric/function.h \
  ../include/numeric/../numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../numeric/../numeric/arithmetic.h \
  ../include/base/exception.h ../include/base/utility.h \
  ../include/base/../declarations.h ../include/numeric/numerical_types.h \
  ../include/geometry/point.h ../include/geometry/../declarations.h \
  ../include/geometry/../base/array.h \
  ../include/geometry/../base/../declarations.h \
  ../include/geometry/../linear_algebra/vector.h \
  ../include/geometry/../linear_algebra/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/integer.h \
  ../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../include/geometry/../linear_algebra/../numeric/interval.h \
  ../include/geometry/../linear_algebra/../numeric/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/arithmetic.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/rounding.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/function.h \
  ../include/geometry/../linear_algebra/matrix.h \
  ../include/geometry/rectangle.h \
  ../include/geometry/../numeric/arithmetic.h \
  ../include/geometry/../numeric/function.h \
  ../include/geometry/../numeric/interval.h \
  ../include/geometry/../geometry/ppl_polyhedron.h \
  ../include/geometry/../geometry/../declarations.h \
  ../include/geometry/../geometry/../numeric/rational.h \
  ../include/geometry/../geometry/../linear_algebra/vector.h \
  ../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../include/geometry/../geometry/point.h \
  ../include/geometry/../geometry/rectangle_expression.h \
  ../include/geometry/list_set.h ../include/geometry/../base/utility.h \
  ../include/output/epsfstream.h \
  ../include/output/../numeric/approximation.h \
  ../include/output/../numeric/float64.h \
  ../include/output/../numeric/../declarations.h \
  ../include/output/../numeric/../numeric/numerical_traits.h \
  ../include/output/../numeric/../numeric/function.h \
  ../include/output/../numeric/../numeric/integer.h \
  ../include/output/../geometry/point.h \
  ../include/output/../geometry/rectangle.h \
  ../include/output/../geometry/list_set.h \
  ../include/output/../geometry/grid_set.h \
  ../include/output/../geometry/../declarations.h \
  ../include/output/../geometry/../base/array.h \
  ../include/output/../geometry/../combinatoric/lattice_set.h \
  ../include/output/../geometry/../combinatoric/../base/array.h \
  ../include/output/../geometry/../combinatoric/../base/iterator.h \
  ../include/output/../geometry/../combinatoric/../base/../base/array.h \
  ../include/output/../geometry/../combinatoric/../base/array_operations.h \
  ../include/output/../geometry/../combinatoric/../base/../declarations.h \
  ../include/output/../geometry/../combinatoric/../numeric/interval.h \
  ../include/output/../geometry/../geometry/rectangle.h \
  ../include/output/../geometry/../geometry/grid.h \
  ../include/output/../geometry/../geometry/../declarations.h \
  ../include/output/../geometry/../geometry/../base/array.h \
  ../include/output/../geometry/../geometry/../numeric/arithmetic.h \
  ../include/output/../geometry/../geometry/../combinatoric/lattice_set.h \
  ../include/output/../geometry/parallelotope.h \
  ../include/output/../geometry/../numeric/interval.h \
  ../include/output/../geometry/../linear_algebra/vector.h \
  ../include/output/../geometry/../linear_algebra/matrix.h \
  ../include/output/../geometry/../geometry/point.h \
  ../include/output/../geometry/../geometry/zonotope.h \
  ../include/output/../geometry/../geometry/../linear_algebra/vector.h \
  ../include/output/../geometry/../geometry/../linear_algebra/matrix.h \
  ../include/output/../geometry/../geometry/../numeric/interval.h \
  ../include/output/../geometry/../geometry/../geometry/ppl_polyhedron.h \
  ../include/output/../geometry/../geometry/../geometry/point.h \
  ../include/output/../geometry/zonotope.h \
  ../include/output/../geometry/polytope.h \
  ../include/output/../geometry/../geometry/ppl_polyhedron.h \
  ../include/output/../geometry/../geometry/point_list.h \
  ../include/output/../geometry/partition_tree_set.h \
  ../include/output/../geometry/../base/iterator.h \
  ../include/output/../geometry/../combinatoric/binary_word.h \
  ../include/output/../geometry/../combinatoric/../declarations.h \
  ../include/output/../geometry/../combinatoric/binary_tree.h \
  ../include/output/../geometry/../combinatoric/../combinatoric/binary_word.h \
  ../include/output/../geometry/../combinatoric/subdivision_tree_set.h \
  ../include/output/../geometry/../combinatoric/../base/sequence.h \
  ../include/output/../geometry/../combinatoric/../combinatoric/binary_tree.h \
  ../include/output/../geometry/../combinatoric/../geometry/rectangle.h \
  ../include/output/../linear_algebra/matrix.h \
  ../include/output/../system/affine_map.h \
  ../include/output/../system/../declarations.h \
  ../include/output/../system/../linear_algebra/vector.h \
  ../include/output/../system/../linear_algebra/matrix.h \
  ../include/output/../system/../system/map.h \
  ../include/output/../system/../system/../declarations.h test.h
	-${CXX} ${CXXFLAGS} -c test_epsfstream.cc 

test_interval_newton.o: test_interval_newton.cc \
  ../include/geometry/point.h ../include/geometry/../declarations.h \
  ../include/geometry/../base/array.h \
  ../include/geometry/../base/../declarations.h \
  ../include/geometry/../linear_algebra/vector.h \
  ../include/geometry/../linear_algebra/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/integer.h \
  ../include/geometry/../linear_algebra/../numeric/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/approximation.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/arithmetic.h \
  ../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../include/geometry/../linear_algebra/../numeric/interval.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/rounding.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/function.h \
  ../include/geometry/../linear_algebra/matrix.h \
  ../include/geometry/rectangle.h \
  ../include/geometry/../numeric/arithmetic.h \
  ../include/geometry/../numeric/function.h \
  ../include/geometry/../numeric/interval.h \
  ../include/geometry/../geometry/ppl_polyhedron.h \
  ../include/geometry/../geometry/../declarations.h \
  ../include/geometry/../geometry/../numeric/rational.h \
  ../include/geometry/../geometry/../numeric/../declarations.h \
  ../include/geometry/../geometry/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../geometry/../numeric/../numeric/arithmetic.h \
  ../include/geometry/../geometry/../numeric/../numeric/function.h \
  ../include/geometry/../geometry/../numeric/../numeric/integer.h \
  ../include/geometry/../geometry/../linear_algebra/vector.h \
  ../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../include/geometry/../geometry/point.h \
  ../include/geometry/../geometry/rectangle_expression.h \
  ../include/system/henon_map.h \
  ../include/system/../linear_algebra/matrix.h \
  ../include/system/../geometry/point.h \
  ../include/system/../geometry/rectangle.h \
  ../include/system/../geometry/parallelotope.h \
  ../include/system/../geometry/../declarations.h \
  ../include/system/../geometry/../numeric/interval.h \
  ../include/system/../geometry/../linear_algebra/vector.h \
  ../include/system/../geometry/../linear_algebra/matrix.h \
  ../include/system/../geometry/../geometry/point.h \
  ../include/system/../geometry/../geometry/zonotope.h \
  ../include/system/../geometry/../geometry/../linear_algebra/vector.h \
  ../include/system/../geometry/../geometry/../linear_algebra/matrix.h \
  ../include/system/../geometry/../geometry/../numeric/interval.h \
  ../include/system/../geometry/../geometry/../geometry/ppl_polyhedron.h \
  ../include/system/../geometry/../geometry/../geometry/point.h \
  ../include/system/../system/map.h \
  ../include/system/../system/../declarations.h \
  ../include/system/../evaluation/apply.h \
  ../include/system/../evaluation/../declarations.h \
  ../include/evaluation/newton.h ../include/evaluation/../declarations.h \
  ../include/linear_algebra/vector.h ../include/linear_algebra/matrix.h \
  ../include/geometry/rectangle.h ../include/system/map.h \
  ../include/system/vector_field.h ../include/system/../declarations.h \
  test.h ../include/real_typedef.h ../include/numeric/numerical_types.h \
  ../include/numeric/mpfloat.h /home/collins/include/mpfr.h \
  ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/function.h \
  ../include/numeric/../numeric/rational.h
	-${CXX} ${CXXFLAGS} -c test_interval_newton.cc 

test_chainreach.o: test_chainreach.cc ../include/real_typedef.h \
  ../include/numeric/numerical_types.h ../include/numeric/mpfloat.h \
  /home/collins/include/mpfr.h ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/approximation.h \
  ../include/numeric/../numeric/function.h \
  ../include/numeric/../numeric/rational.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../numeric/arithmetic.h \
  ../include/numeric/../numeric/../numeric/function.h \
  ../include/numeric/../numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../numeric/../numeric/arithmetic.h \
  ../include/geometry/point.h ../include/geometry/../declarations.h \
  ../include/geometry/../base/array.h \
  ../include/geometry/../base/../declarations.h \
  ../include/geometry/../linear_algebra/vector.h \
  ../include/geometry/../linear_algebra/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/integer.h \
  ../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../include/geometry/../linear_algebra/../numeric/interval.h \
  ../include/geometry/../linear_algebra/../numeric/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/arithmetic.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/rounding.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/function.h \
  ../include/geometry/../linear_algebra/matrix.h \
  ../include/geometry/rectangle.h \
  ../include/geometry/../numeric/arithmetic.h \
  ../include/geometry/../numeric/function.h \
  ../include/geometry/../numeric/interval.h \
  ../include/geometry/../geometry/ppl_polyhedron.h \
  ../include/geometry/../geometry/../declarations.h \
  ../include/geometry/../geometry/../numeric/rational.h \
  ../include/geometry/../geometry/../linear_algebra/vector.h \
  ../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../include/geometry/../geometry/point.h \
  ../include/geometry/../geometry/rectangle_expression.h \
  ../include/geometry/grid.h \
  ../include/geometry/../combinatoric/lattice_set.h \
  ../include/geometry/../combinatoric/../base/array.h \
  ../include/geometry/../combinatoric/../base/iterator.h \
  ../include/geometry/../combinatoric/../base/../base/array.h \
  ../include/geometry/../combinatoric/../base/array_operations.h \
  ../include/geometry/../combinatoric/../base/../declarations.h \
  ../include/geometry/../combinatoric/../numeric/interval.h \
  ../include/geometry/grid_set.h \
  ../include/geometry/../geometry/rectangle.h \
  ../include/geometry/../geometry/grid.h \
  ../include/geometry/partition_tree_set.h \
  ../include/geometry/../base/iterator.h \
  ../include/geometry/../combinatoric/binary_word.h \
  ../include/geometry/../combinatoric/../declarations.h \
  ../include/geometry/../combinatoric/binary_tree.h \
  ../include/geometry/../combinatoric/../combinatoric/binary_word.h \
  ../include/geometry/../combinatoric/subdivision_tree_set.h \
  ../include/geometry/../combinatoric/../base/sequence.h \
  ../include/geometry/../combinatoric/../combinatoric/binary_tree.h \
  ../include/geometry/../combinatoric/../geometry/rectangle.h \
  ../include/system/henon_map.h \
  ../include/system/../linear_algebra/matrix.h \
  ../include/system/../geometry/point.h \
  ../include/system/../geometry/rectangle.h \
  ../include/system/../geometry/parallelotope.h \
  ../include/system/../geometry/../declarations.h \
  ../include/system/../geometry/../numeric/interval.h \
  ../include/system/../geometry/../linear_algebra/vector.h \
  ../include/system/../geometry/../linear_algebra/matrix.h \
  ../include/system/../geometry/../geometry/point.h \
  ../include/system/../geometry/../geometry/zonotope.h \
  ../include/system/../geometry/../geometry/../linear_algebra/vector.h \
  ../include/system/../geometry/../geometry/../linear_algebra/matrix.h \
  ../include/system/../geometry/../geometry/../numeric/interval.h \
  ../include/system/../geometry/../geometry/../geometry/ppl_polyhedron.h \
  ../include/system/../geometry/../geometry/../geometry/point.h \
  ../include/system/../system/map.h \
  ../include/system/../system/../declarations.h \
  ../include/system/../evaluation/apply.h \
  ../include/system/../evaluation/../declarations.h \
  ../include/evaluation/apply.h ../include/output/epsfstream.h \
  ../include/output/../numeric/approximation.h \
  ../include/output/../numeric/float64.h \
  ../include/output/../numeric/../declarations.h \
  ../include/output/../numeric/../numeric/numerical_traits.h \
  ../include/output/../numeric/../numeric/function.h \
  ../include/output/../numeric/../numeric/integer.h \
  ../include/output/../geometry/point.h \
  ../include/output/../geometry/rectangle.h \
  ../include/output/../geometry/list_set.h \
  ../include/output/../geometry/../declarations.h \
  ../include/output/../geometry/../base/utility.h \
  ../include/output/../geometry/../base/../declarations.h \
  ../include/output/../geometry/grid_set.h \
  ../include/output/../geometry/parallelotope.h \
  ../include/output/../geometry/zonotope.h \
  ../include/output/../geometry/polytope.h \
  ../include/output/../geometry/../linear_algebra/vector.h \
  ../include/output/../geometry/../linear_algebra/matrix.h \
  ../include/output/../geometry/../geometry/ppl_polyhedron.h \
  ../include/output/../geometry/../geometry/point.h \
  ../include/output/../geometry/../geometry/point_list.h \
  ../include/output/../geometry/../geometry/../declarations.h \
  ../include/output/../geometry/../geometry/../linear_algebra/vector.h \
  ../include/output/../geometry/../geometry/../linear_algebra/matrix.h \
  ../include/output/../geometry/../geometry/../geometry/point.h \
  ../include/output/../geometry/../geometry/rectangle.h \
  ../include/output/../geometry/partition_tree_set.h \
  ../include/output/../linear_algebra/matrix.h \
  ../include/output/../system/affine_map.h \
  ../include/output/../system/../declarations.h \
  ../include/output/../system/../linear_algebra/vector.h \
  ../include/output/../system/../linear_algebra/matrix.h \
  ../include/output/../system/../system/map.h test.h
	-${CXX} ${CXXFLAGS} -c test_chainreach.cc 

test_integration_step.o: test_integration_step.cc ../include/ariadne.h \
  ../include/real_typedef.h ../include/numeric/numerical_types.h \
  ../include/numeric/mpfloat.h /home/collins/include/mpfr.h \
  ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/approximation.h \
  ../include/numeric/../numeric/function.h \
  ../include/numeric/../numeric/rational.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../numeric/arithmetic.h \
  ../include/numeric/../numeric/../numeric/function.h \
  ../include/numeric/../numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../numeric/../numeric/arithmetic.h \
  ../include/numeric/numerical_types.h ../include/linear_algebra/vector.h \
  ../include/linear_algebra/../declarations.h \
  ../include/linear_algebra/../numeric/integer.h \
  ../include/linear_algebra/../numeric/arithmetic.h \
  ../include/linear_algebra/../numeric/interval.h \
  ../include/linear_algebra/../numeric/../declarations.h \
  ../include/linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../include/linear_algebra/../numeric/../numeric/arithmetic.h \
  ../include/linear_algebra/../numeric/../numeric/rounding.h \
  ../include/linear_algebra/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/linear_algebra/../numeric/../numeric/../numeric/function.h \
  ../include/linear_algebra/matrix.h ../include/geometry/rectangle.h \
  ../include/geometry/../declarations.h \
  ../include/geometry/../numeric/arithmetic.h \
  ../include/geometry/../numeric/function.h \
  ../include/geometry/../numeric/interval.h \
  ../include/geometry/../linear_algebra/vector.h \
  ../include/geometry/../geometry/ppl_polyhedron.h \
  ../include/geometry/../geometry/../declarations.h \
  ../include/geometry/../geometry/../numeric/rational.h \
  ../include/geometry/../geometry/../linear_algebra/vector.h \
  ../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../include/geometry/../geometry/point.h \
  ../include/geometry/../geometry/../base/array.h \
  ../include/geometry/../geometry/../base/../declarations.h \
  ../include/geometry/../geometry/rectangle_expression.h \
  ../include/geometry/parallelotope.h \
  ../include/geometry/../linear_algebra/matrix.h \
  ../include/geometry/../geometry/zonotope.h \
  ../include/geometry/../geometry/../numeric/interval.h \
  ../include/geometry/../geometry/../geometry/ppl_polyhedron.h \
  ../include/geometry/../geometry/../geometry/point.h \
  ../include/system/affine_vector_field.h \
  ../include/system/../linear_algebra/vector.h \
  ../include/system/../linear_algebra/matrix.h \
  ../include/system/../system/vector_field.h \
  ../include/system/../system/../declarations.h \
  ../include/evaluation/lohner_integrator.h \
  ../include/evaluation/../declarations.h \
  ../include/evaluation/../evaluation/integrator.h \
  ../include/evaluation/../evaluation/../declarations.h test.h
	-${CXX} ${CXXFLAGS} -c test_integration_step.cc 

