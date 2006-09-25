include ../../config.mk

OBJ =   vector.o matrix.o tensor.o \
	linear_program.o \
	matrix_function.o
  
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
vector.o: vector.cc ../../include/numeric/interval.h \
  ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/approximation.h \
  ../../include/numeric/../numeric/arithmetic.h \
  ../../include/numeric/../numeric/rounding.h \
  ../../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/../numeric/function.h \
  ../../include/linear_algebra/vector.h \
  ../../include/linear_algebra/../declarations.h \
  ../../include/linear_algebra/../numeric/integer.h \
  ../../include/linear_algebra/../numeric/../declarations.h \
  ../../include/linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/linear_algebra/../numeric/../numeric/arithmetic.h \
  ../../include/linear_algebra/../numeric/arithmetic.h \
  ../../include/linear_algebra/../numeric/interval.h \
  ../../include/linear_algebra/vector.tpl \
  ../../include/linear_algebra/vector.h \
  ../../include/linear_algebra/../base/stlio.h \
  ../../include/linear_algebra/../base/../base/array.h \
  ../../include/linear_algebra/../base/../base/../declarations.h \
  ../../include/numeric/float64.h /home/collins/include/mpfr.h \
  ../../include/numeric/../numeric/function.h \
  ../../include/numeric/../numeric/integer.h \
  ../../include/numeric/mpfloat.h \
  ../../include/numeric/../numeric/rational.h \
  ../../include/numeric/../numeric/../declarations.h \
  ../../include/numeric/../numeric/../numeric/arithmetic.h \
  ../../include/numeric/../numeric/../numeric/integer.h \
  ../../include/numeric/rational.h

matrix.o: matrix.cc ../../include/numeric/interval.h \
  ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/approximation.h \
  ../../include/numeric/../numeric/arithmetic.h \
  ../../include/numeric/../numeric/rounding.h \
  ../../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/../numeric/function.h \
  ../../include/linear_algebra/matrix.h \
  ../../include/linear_algebra/../declarations.h \
  ../../include/linear_algebra/../numeric/integer.h \
  ../../include/linear_algebra/../numeric/../declarations.h \
  ../../include/linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/linear_algebra/../numeric/../numeric/arithmetic.h \
  ../../include/linear_algebra/../numeric/interval.h \
  ../../include/linear_algebra/matrix.tpl \
  ../../include/linear_algebra/matrix.h \
  ../../include/linear_algebra/../base/array.h \
  ../../include/linear_algebra/../base/../declarations.h \
  ../../include/linear_algebra/../numeric/arithmetic.h \
  ../../include/linear_algebra/../linear_algebra/vector.h \
  ../../include/linear_algebra/../linear_algebra/../declarations.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/integer.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/arithmetic.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/interval.h \
  ../../include/linear_algebra/../linear_algebra/matrix.h \
  ../../include/linear_algebra/../linear_algebra/lu_matrix.h \
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
  ../../include/linear_algebra/../linear_algebra/../base/array.h \
  ../../include/linear_algebra/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/numeric/float64.h /home/collins/include/mpfr.h \
  ../../include/numeric/../numeric/function.h \
  ../../include/numeric/../numeric/integer.h \
  ../../include/numeric/mpfloat.h \
  ../../include/numeric/../numeric/rational.h \
  ../../include/numeric/../numeric/../declarations.h \
  ../../include/numeric/../numeric/../numeric/arithmetic.h \
  ../../include/numeric/../numeric/../numeric/integer.h \
  ../../include/numeric/rational.h ../../include/real_typedef.h \
  ../../include/numeric/numerical_types.h ../../include/numeric/mpfloat.h

tensor.o: tensor.cc ../../include/linear_algebra/tensor.h \
  ../../include/linear_algebra/../declarations.h \
  ../../include/linear_algebra/../base/array.h \
  ../../include/linear_algebra/../base/../declarations.h \
  ../../include/linear_algebra/tensor.tpl \
  ../../include/linear_algebra/tensor.h \
  ../../include/linear_algebra/../linear_algebra/vector.h \
  ../../include/linear_algebra/../linear_algebra/../declarations.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/integer.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/../declarations.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/../numeric/approximation.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/../numeric/arithmetic.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/arithmetic.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/interval.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/../numeric/rounding.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/../numeric/../numeric/function.h \
  ../../include/linear_algebra/../linear_algebra/matrix.h \
  ../../include/numeric/float64.h /home/collins/include/mpfr.h \
  ../../include/numeric/../declarations.h \
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

linear_program.o: linear_program.cc \
  ../../include/linear_algebra/linear_program.h \
  ../../include/linear_algebra/../linear_algebra/matrix.h \
  ../../include/linear_algebra/../linear_algebra/../declarations.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/integer.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/../declarations.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/../numeric/approximation.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/../numeric/arithmetic.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/interval.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/../numeric/rounding.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/../numeric/../numeric/function.h \
  ../../include/linear_algebra/linear_program.tpl \
  ../../include/linear_algebra/linear_program.h \
  ../../include/linear_algebra/../linear_algebra/vector.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/arithmetic.h \
  ../../include/linear_algebra/../base/stlio.h \
  ../../include/linear_algebra/../base/../base/array.h \
  ../../include/linear_algebra/../base/../base/../declarations.h \
  ../../include/numeric/rational.h \
  ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/arithmetic.h \
  ../../include/numeric/../numeric/function.h \
  ../../include/numeric/../numeric/integer.h

matrix_function.o: matrix_function.cc \
  ../../include/linear_algebra/matrix_function.h \
  ../../include/linear_algebra/../declarations.h \
  ../../include/linear_algebra/matrix_function.tpl \
  ../../include/linear_algebra/matrix_function.h \
  ../../include/linear_algebra/../numeric/interval.h \
  ../../include/linear_algebra/../numeric/../declarations.h \
  ../../include/linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/linear_algebra/../numeric/../numeric/approximation.h \
  ../../include/linear_algebra/../numeric/../numeric/arithmetic.h \
  ../../include/linear_algebra/../numeric/../numeric/rounding.h \
  ../../include/linear_algebra/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/linear_algebra/../numeric/../numeric/../numeric/function.h \
  ../../include/linear_algebra/../linear_algebra/matrix.h \
  ../../include/linear_algebra/../linear_algebra/../declarations.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/integer.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/../declarations.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/../numeric/arithmetic.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/interval.h \
  ../../include/numeric/float64.h /home/collins/include/mpfr.h \
  ../../include/numeric/../declarations.h \
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

