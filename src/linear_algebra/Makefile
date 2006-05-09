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
  /ufs/collins/include/tblas/xerbla.hpp \
  /ufs/collins/include/tblas/trcpy.hpp \
  /ufs/collins/include/tblas/gemm.hpp /ufs/collins/include/tblas/trmm.hpp \
  /ufs/collins/include/tlapack/laswp.hpp \
  /ufs/collins/include/tlapack/tlapack.hpp \
  /ufs/collins/include/tblas/tblas.hpp \
  /ufs/collins/include/tlapack/getrf.hpp \
  /ufs/collins/include/tblas/iamax.hpp \
  /ufs/collins/include/tblas/scal.hpp /ufs/collins/include/tblas/swap.hpp \
  /ufs/collins/include/tblas/ger.hpp \
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
  ../../include/numeric/../numeric/../utility/stlio.h \
  ../../include/numeric/../numeric/approximation.h \
  ../../include/numeric/../numeric/../numeric/numerical_types.h \
  ../../include/real_typedef.h ../../include/numeric/numerical_types.h \
  ../../include/linear_algebra/interval_vector.h \
  ../../include/linear_algebra/../declarations.h \
  ../../include/linear_algebra/../numeric/interval.h \
  ../../include/linear_algebra/../numeric/../declarations.h \
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
  ../../include/numeric/../utility/stlio.h \
  ../../include/linear_algebra/interval_matrix.h \
  ../../include/linear_algebra/../linear_algebra/matrix.h \
  ../../include/linear_algebra/../linear_algebra/../declarations.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/integer.h \
  ../../include/linear_algebra/../declarations.h \
  ../../include/linear_algebra/../numeric/interval.h \
  ../../include/linear_algebra/../numeric/../declarations.h \
  ../../include/linear_algebra/interval_matrix.tpl \
  ../../include/linear_algebra/interval_matrix.h \
  ../../include/linear_algebra/../numeric/arithmetic.h \
  ../../include/linear_algebra/../numeric/../numeric/numerical_types.h \
  ../../include/linear_algebra/../numeric/../numeric/approximation.h \
  ../../include/linear_algebra/../numeric/../numeric/../declarations.h \
  ../../include/linear_algebra/../numeric/../numeric/../numeric/numerical_types.h \
  ../../include/linear_algebra/../linear_algebra/vector.h \
  ../../include/linear_algebra/../linear_algebra/interval_vector.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/interval.h \
  ../../include/linear_algebra/../linear_algebra/../linear_algebra/vector.h

constraint_system.o: constraint_system.cc ../../include/real_typedef.h \
  ../../include/numeric/numerical_types.h \
  ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/integer.h \
  ../../include/numeric/../numeric/../declarations.h \
  ../../include/numeric/../numeric/../utility/stlio.h \
  ../../include/numeric/../numeric/../utility/../base/array.h \
  ../../include/numeric/../numeric/../utility/../base/../declarations.h \
  ../../include/numeric/../numeric/dyadic.h \
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
  ../../include/numeric/../utility/stlio.h \
  ../../include/linear_algebra/transformation_system.h \
  ../../include/linear_algebra/../linear_algebra/vector.h \
  ../../include/linear_algebra/../linear_algebra/../declarations.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/integer.h \
  ../../include/linear_algebra/../linear_algebra/matrix.h \
  ../../include/linear_algebra/../linear_algebra/interval_vector.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/interval.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/../declarations.h \
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
  ../../include/numeric/../utility/stlio.h \
  ../../include/linear_algebra/linear_program.h \
  ../../include/linear_algebra/../linear_algebra/matrix.h \
  ../../include/linear_algebra/../linear_algebra/../declarations.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/integer.h \
  ../../include/linear_algebra/linear_program.tpl \
  ../../include/linear_algebra/linear_program.h \
  ../../include/linear_algebra/../linear_algebra/vector.h \
  ../../include/utility/stlio.h

