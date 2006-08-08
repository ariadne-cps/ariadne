include ../../config.mk

OBJ =   vector.o matrix.o tensor.o \
	interval_vector.o interval_matrix.o \
	constraint_system.o generator_system.o transformation_system.o \
	linear_program.o
  
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
vector.o: vector.cc ../../include/real_typedef.h \
  ../../include/numeric/numerical_types.h \
  ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/integer.h \
  ../../include/numeric/../numeric/../declarations.h \
  ../../include/numeric/../numeric/../utility/stlio.h \
  ../../include/numeric/../numeric/../utility/../base/array.h \
  ../../include/numeric/../numeric/../utility/../base/../declarations.h \
  ../../include/numeric/../numeric/dyadic.h \
  ../../include/numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../utility/stlio.h \
  ../../include/linear_algebra/vector.h \
  ../../include/linear_algebra/../declarations.h \
  ../../include/linear_algebra/../numeric/integer.h \
  ../../include/linear_algebra/vector.tpl \
  ../../include/linear_algebra/vector.h \
  ../../include/linear_algebra/../utility/stlio.h

matrix.o: matrix.cc ../../include/real_typedef.h \
  ../../include/numeric/numerical_types.h \
  ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/integer.h \
  ../../include/numeric/../numeric/../declarations.h \
  ../../include/numeric/../numeric/../utility/stlio.h \
  ../../include/numeric/../numeric/../utility/../base/array.h \
  ../../include/numeric/../numeric/../utility/../base/../declarations.h \
  ../../include/numeric/../numeric/dyadic.h \
  ../../include/numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../utility/stlio.h \
  ../../include/linear_algebra/matrix.h \
  ../../include/linear_algebra/../declarations.h \
  ../../include/linear_algebra/../numeric/integer.h \
  ../../include/linear_algebra/matrix.tpl \
  ../../include/linear_algebra/matrix.h \
  ../../include/linear_algebra/../base/array.h \
  ../../include/linear_algebra/../numeric/arithmetic.h \
  ../../include/linear_algebra/../numeric/../numeric/numerical_types.h \
  ../../include/linear_algebra/../numeric/../numeric/approximation.h \
  ../../include/linear_algebra/../numeric/../numeric/../declarations.h \
  ../../include/linear_algebra/../numeric/../numeric/../numeric/numerical_types.h \
  ../../include/linear_algebra/../numeric/function.h \
  /ufs/collins/include/mpfr.h \
  ../../include/linear_algebra/../linear_algebra/vector.h \
  ../../include/linear_algebra/../linear_algebra/../declarations.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/integer.h \
  ../../include/linear_algebra/../linear_algebra/matrix.h \
  ../../include/linear_algebra/../linear_algebra/lu_matrix.h \
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
  ../../include/linear_algebra/../linear_algebra/../base/array.h \
  ../../include/linear_algebra/../linear_algebra/../linear_algebra/matrix.h

tensor.o: tensor.cc ../../include/real_typedef.h \
  ../../include/numeric/numerical_types.h \
  ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/integer.h \
  ../../include/numeric/../numeric/../declarations.h \
  ../../include/numeric/../numeric/../utility/stlio.h \
  ../../include/numeric/../numeric/../utility/../base/array.h \
  ../../include/numeric/../numeric/../utility/../base/../declarations.h \
  ../../include/numeric/../numeric/dyadic.h \
  ../../include/numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../utility/stlio.h \
  ../../include/linear_algebra/tensor.h \
  ../../include/linear_algebra/../declarations.h \
  ../../include/linear_algebra/../base/array.h \
  ../../include/linear_algebra/tensor.tpl \
  ../../include/linear_algebra/tensor.h \
  ../../include/linear_algebra/../linear_algebra/vector.h \
  ../../include/linear_algebra/../linear_algebra/../declarations.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/integer.h \
  ../../include/linear_algebra/../linear_algebra/matrix.h

interval_vector.o: interval_vector.cc ../../include/numeric/arithmetic.h \
  ../../include/numeric/../numeric/numerical_types.h \
  ../../include/numeric/../numeric/../declarations.h \
  ../../include/numeric/../numeric/../numeric/integer.h \
  ../../include/numeric/../numeric/../numeric/../declarations.h \
  ../../include/numeric/../numeric/../numeric/../utility/stlio.h \
  ../../include/numeric/../numeric/../numeric/../utility/../base/array.h \
  ../../include/numeric/../numeric/../numeric/../utility/../base/../declarations.h \
  ../../include/numeric/../numeric/../numeric/dyadic.h \
  ../../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/../utility/stlio.h \
  ../../include/numeric/../numeric/approximation.h \
  ../../include/numeric/../numeric/../numeric/numerical_types.h \
  ../../include/real_typedef.h ../../include/numeric/numerical_types.h \
  ../../include/linear_algebra/interval_vector.h \
  ../../include/linear_algebra/../declarations.h \
  ../../include/linear_algebra/../numeric/interval.h \
  ../../include/linear_algebra/../numeric/../declarations.h \
  ../../include/linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/linear_algebra/../linear_algebra/vector.h \
  ../../include/linear_algebra/../linear_algebra/../declarations.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/integer.h \
  ../../include/linear_algebra/interval_vector.tpl \
  ../../include/linear_algebra/interval_vector.h

interval_matrix.o: interval_matrix.cc ../../include/real_typedef.h \
  ../../include/numeric/numerical_types.h \
  ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/integer.h \
  ../../include/numeric/../numeric/../declarations.h \
  ../../include/numeric/../numeric/../utility/stlio.h \
  ../../include/numeric/../numeric/../utility/../base/array.h \
  ../../include/numeric/../numeric/../utility/../base/../declarations.h \
  ../../include/numeric/../numeric/dyadic.h \
  ../../include/numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../utility/stlio.h \
  ../../include/linear_algebra/interval_matrix.h \
  ../../include/linear_algebra/../linear_algebra/matrix.h \
  ../../include/linear_algebra/../linear_algebra/../declarations.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/integer.h \
  ../../include/linear_algebra/../declarations.h \
  ../../include/linear_algebra/../numeric/interval.h \
  ../../include/linear_algebra/../numeric/../declarations.h \
  ../../include/linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/linear_algebra/interval_matrix.tpl \
  ../../include/declarations.h \
  ../../include/linear_algebra/interval_matrix.h \
  ../../include/linear_algebra/../numeric/arithmetic.h \
  ../../include/linear_algebra/../numeric/../numeric/numerical_types.h \
  ../../include/linear_algebra/../numeric/../numeric/approximation.h \
  ../../include/linear_algebra/../numeric/../numeric/../declarations.h \
  ../../include/linear_algebra/../numeric/../numeric/../numeric/numerical_types.h \
  ../../include/linear_algebra/../linear_algebra/vector.h \
  ../../include/linear_algebra/../linear_algebra/interval_vector.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/interval.h \
  ../../include/linear_algebra/../linear_algebra/../linear_algebra/vector.h \
  ../../include/linear_algebra/../linear_algebra/lu_matrix.h \
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
  ../../include/linear_algebra/../linear_algebra/../base/array.h \
  ../../include/linear_algebra/../linear_algebra/../linear_algebra/matrix.h

constraint_system.o: constraint_system.cc ../../include/real_typedef.h \
  ../../include/numeric/numerical_types.h \
  ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/integer.h \
  ../../include/numeric/../numeric/../declarations.h \
  ../../include/numeric/../numeric/../utility/stlio.h \
  ../../include/numeric/../numeric/../utility/../base/array.h \
  ../../include/numeric/../numeric/../utility/../base/../declarations.h \
  ../../include/numeric/../numeric/dyadic.h \
  ../../include/numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../utility/stlio.h \
  ../../include/linear_algebra/constraint_system.h \
  /ufs/collins/include/ppl.hh \
  ../../include/linear_algebra/../declarations.h \
  ../../include/linear_algebra/../numeric/numerical_types.h \
  ../../include/linear_algebra/../linear_algebra/vector.h \
  ../../include/linear_algebra/../linear_algebra/../declarations.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/integer.h \
  ../../include/linear_algebra/../linear_algebra/matrix.h \
  ../../include/linear_algebra/constraint_system.tpl \
  ../../include/linear_algebra/constraint_system.h

generator_system.o: generator_system.cc ../../include/real_typedef.h \
  ../../include/numeric/numerical_types.h \
  ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/integer.h \
  ../../include/numeric/../numeric/../declarations.h \
  ../../include/numeric/../numeric/../utility/stlio.h \
  ../../include/numeric/../numeric/../utility/../base/array.h \
  ../../include/numeric/../numeric/../utility/../base/../declarations.h \
  ../../include/numeric/../numeric/dyadic.h \
  ../../include/numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../utility/stlio.h \
  ../../include/linear_algebra/generator_system.h \
  /ufs/collins/include/ppl.hh \
  ../../include/linear_algebra/../declarations.h \
  ../../include/linear_algebra/../linear_algebra/vector.h \
  ../../include/linear_algebra/../linear_algebra/../declarations.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/integer.h \
  ../../include/linear_algebra/../linear_algebra/matrix.h \
  ../../include/linear_algebra/generator_system.tpl \
  ../../include/linear_algebra/../numeric/arithmetic.h \
  ../../include/linear_algebra/../numeric/../numeric/numerical_types.h \
  ../../include/linear_algebra/../numeric/../numeric/approximation.h \
  ../../include/linear_algebra/../numeric/../numeric/../declarations.h \
  ../../include/linear_algebra/../numeric/../numeric/../numeric/numerical_types.h \
  ../../include/linear_algebra/generator_system.h

transformation_system.o: transformation_system.cc \
  ../../include/real_typedef.h ../../include/numeric/numerical_types.h \
  ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/integer.h \
  ../../include/numeric/../numeric/../declarations.h \
  ../../include/numeric/../numeric/../utility/stlio.h \
  ../../include/numeric/../numeric/../utility/../base/array.h \
  ../../include/numeric/../numeric/../utility/../base/../declarations.h \
  ../../include/numeric/../numeric/dyadic.h \
  ../../include/numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../utility/stlio.h \
  ../../include/linear_algebra/transformation_system.h \
  ../../include/linear_algebra/../linear_algebra/vector.h \
  ../../include/linear_algebra/../linear_algebra/../declarations.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/integer.h \
  ../../include/linear_algebra/../linear_algebra/matrix.h \
  ../../include/linear_algebra/../linear_algebra/interval_vector.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/interval.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/../declarations.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/linear_algebra/../linear_algebra/../linear_algebra/vector.h \
  ../../include/linear_algebra/../linear_algebra/interval_matrix.h \
  ../../include/linear_algebra/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/linear_algebra/transformation_system.tpl \
  ../../include/linear_algebra/transformation_system.h

linear_program.o: linear_program.cc ../../include/real_typedef.h \
  ../../include/numeric/numerical_types.h \
  ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/integer.h \
  ../../include/numeric/../numeric/../declarations.h \
  ../../include/numeric/../numeric/../utility/stlio.h \
  ../../include/numeric/../numeric/../utility/../base/array.h \
  ../../include/numeric/../numeric/../utility/../base/../declarations.h \
  ../../include/numeric/../numeric/dyadic.h \
  ../../include/numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../utility/stlio.h \
  ../../include/linear_algebra/linear_program.h \
  ../../include/linear_algebra/../linear_algebra/matrix.h \
  ../../include/linear_algebra/../linear_algebra/../declarations.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/integer.h \
  ../../include/linear_algebra/linear_program.tpl \
  ../../include/linear_algebra/linear_program.h \
  ../../include/linear_algebra/../linear_algebra/vector.h \
  ../../include/utility/stlio.h

