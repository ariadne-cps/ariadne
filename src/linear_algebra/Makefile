include ../../config.mk

OBJ =   vector.o matrix.o tensor.o \
	interval_vector.o interval_matrix.o \
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
vector.o: vector.cc ../../include/linear_algebra/vector.h \
  ../../include/linear_algebra/../declarations.h \
  ../../include/linear_algebra/../numeric/integer.h \
  ../../include/linear_algebra/../numeric/../declarations.h \
  ../../include/linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/linear_algebra/../numeric/arithmetic.h \
  ../../include/linear_algebra/../numeric/../numeric/integer.h \
  ../../include/linear_algebra/vector.tpl \
  ../../include/linear_algebra/vector.h \
  ../../include/linear_algebra/../base/stlio.h \
  ../../include/linear_algebra/../base/../base/array.h \
  ../../include/linear_algebra/../base/../base/../declarations.h \
  ../../include/numeric/float64.h ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/function.h \
  ../../include/numeric/../numeric/integer.h \
  ../../include/numeric/../numeric/integer.h \
  ../../include/numeric/mpfloat.h \
  ../../include/numeric/../numeric/rational.h \
  ../../include/numeric/../numeric/../declarations.h \
  ../../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/../numeric/arithmetic.h \
  ../../include/numeric/../numeric/../numeric/function.h \
  ../../include/numeric/../numeric/../numeric/integer.h \
  ../../include/numeric/rational.h

matrix.o: matrix.cc ../../include/linear_algebra/matrix.h \
  ../../include/linear_algebra/../declarations.h \
  ../../include/linear_algebra/../numeric/integer.h \
  ../../include/linear_algebra/../numeric/../declarations.h \
  ../../include/linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/linear_algebra/matrix.tpl \
  ../../include/linear_algebra/matrix.h \
  ../../include/linear_algebra/../base/array.h \
  ../../include/linear_algebra/../base/../declarations.h \
  ../../include/linear_algebra/../numeric/arithmetic.h \
  ../../include/linear_algebra/../numeric/../numeric/integer.h \
  ../../include/linear_algebra/../linear_algebra/vector.h \
  ../../include/linear_algebra/../linear_algebra/../declarations.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/integer.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/arithmetic.h \
  ../../include/linear_algebra/../linear_algebra/matrix.h \
  ../../include/linear_algebra/../linear_algebra/lu_matrix.h \
  ../../include/linear_algebra/../linear_algebra/../base/array.h \
  ../../include/linear_algebra/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/numeric/float64.h ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/function.h \
  ../../include/numeric/../numeric/integer.h \
  ../../include/numeric/../numeric/integer.h \
  ../../include/numeric/mpfloat.h \
  ../../include/numeric/../numeric/rational.h \
  ../../include/numeric/../numeric/../declarations.h \
  ../../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/../numeric/arithmetic.h \
  ../../include/numeric/../numeric/../numeric/function.h \
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
  ../../include/linear_algebra/../linear_algebra/../numeric/arithmetic.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/../numeric/integer.h \
  ../../include/linear_algebra/../linear_algebra/matrix.h \
  ../../include/numeric/float64.h ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/function.h \
  ../../include/numeric/../numeric/integer.h \
  ../../include/numeric/../numeric/integer.h \
  ../../include/numeric/mpfloat.h \
  ../../include/numeric/../numeric/rational.h \
  ../../include/numeric/../numeric/../declarations.h \
  ../../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/../numeric/arithmetic.h \
  ../../include/numeric/../numeric/../numeric/function.h \
  ../../include/numeric/../numeric/../numeric/integer.h \
  ../../include/numeric/rational.h

interval_vector.o: interval_vector.cc ../../include/numeric/arithmetic.h \
  ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/integer.h \
  ../../include/numeric/../numeric/../declarations.h \
  ../../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/linear_algebra/interval_vector.h \
  ../../include/linear_algebra/../declarations.h \
  ../../include/linear_algebra/../numeric/interval.h \
  ../../include/linear_algebra/../numeric/../declarations.h \
  ../../include/linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/linear_algebra/../numeric/../numeric/arithmetic.h \
  ../../include/linear_algebra/../numeric/../numeric/rounding.h \
  ../../include/linear_algebra/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/linear_algebra/../numeric/../numeric/../numeric/function.h \
  ../../include/linear_algebra/../numeric/../numeric/../numeric/integer.h \
  ../../include/linear_algebra/../linear_algebra/vector.h \
  ../../include/linear_algebra/../linear_algebra/../declarations.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/integer.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/arithmetic.h \
  ../../include/linear_algebra/interval_vector.tpl \
  ../../include/linear_algebra/interval_vector.h \
  ../../include/linear_algebra/vector.tpl \
  ../../include/linear_algebra/vector.h \
  ../../include/linear_algebra/../base/stlio.h \
  ../../include/linear_algebra/../base/../base/array.h \
  ../../include/linear_algebra/../base/../base/../declarations.h \
  ../../include/numeric/float64.h \
  ../../include/numeric/../numeric/function.h \
  ../../include/numeric/mpfloat.h \
  ../../include/numeric/../numeric/rational.h \
  ../../include/numeric/../numeric/../numeric/arithmetic.h \
  ../../include/numeric/../numeric/../numeric/function.h \
  ../../include/numeric/../numeric/../numeric/integer.h \
  ../../include/numeric/rational.h

interval_matrix.o: interval_matrix.cc \
  ../../include/linear_algebra/interval_matrix.h \
  ../../include/linear_algebra/../linear_algebra/matrix.h \
  ../../include/linear_algebra/../linear_algebra/../declarations.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/integer.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/../declarations.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/linear_algebra/../declarations.h \
  ../../include/linear_algebra/../numeric/interval.h \
  ../../include/linear_algebra/../numeric/../declarations.h \
  ../../include/linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/linear_algebra/../numeric/../numeric/arithmetic.h \
  ../../include/linear_algebra/../numeric/../numeric/../declarations.h \
  ../../include/linear_algebra/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/linear_algebra/../numeric/../numeric/../numeric/integer.h \
  ../../include/linear_algebra/../numeric/../numeric/rounding.h \
  ../../include/linear_algebra/../numeric/../numeric/../numeric/function.h \
  ../../include/linear_algebra/../numeric/../numeric/../numeric/integer.h \
  ../../include/linear_algebra/interval_matrix.tpl \
  ../../include/declarations.h \
  ../../include/linear_algebra/interval_matrix.h \
  ../../include/linear_algebra/../numeric/arithmetic.h \
  ../../include/linear_algebra/../numeric/approximation.h \
  ../../include/linear_algebra/../linear_algebra/vector.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/arithmetic.h \
  ../../include/linear_algebra/../linear_algebra/interval_vector.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/interval.h \
  ../../include/linear_algebra/../linear_algebra/../linear_algebra/vector.h \
  ../../include/linear_algebra/../linear_algebra/lu_matrix.h \
  ../../include/linear_algebra/../linear_algebra/../base/array.h \
  ../../include/linear_algebra/../linear_algebra/../base/../declarations.h \
  ../../include/linear_algebra/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/linear_algebra/matrix.tpl \
  ../../include/linear_algebra/matrix.h \
  ../../include/linear_algebra/../base/array.h \
  ../../include/linear_algebra/../numeric/integer.h \
  ../../include/numeric/float64.h ../../include/numeric/../declarations.h \
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
  ../../include/linear_algebra/linear_program.tpl \
  ../../include/linear_algebra/linear_program.h \
  ../../include/linear_algebra/../linear_algebra/vector.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/arithmetic.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/../numeric/integer.h \
  ../../include/linear_algebra/../base/stlio.h \
  ../../include/linear_algebra/../base/../base/array.h \
  ../../include/linear_algebra/../base/../base/../declarations.h \
  ../../include/numeric/rational.h \
  ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/arithmetic.h \
  ../../include/numeric/../numeric/function.h \
  ../../include/numeric/../numeric/integer.h \
  ../../include/numeric/../numeric/integer.h

matrix_function.o: matrix_function.cc \
  ../../include/linear_algebra/matrix_function.h \
  ../../include/linear_algebra/../declarations.h \
  ../../include/linear_algebra/matrix_function.tpl \
  ../../include/linear_algebra/matrix_function.h \
  ../../include/linear_algebra/../linear_algebra/matrix.h \
  ../../include/linear_algebra/../linear_algebra/../declarations.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/integer.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/../declarations.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../../include/linear_algebra/../linear_algebra/interval_matrix.h \
  ../../include/linear_algebra/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/interval.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/../numeric/arithmetic.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/../numeric/../declarations.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/../numeric/../numeric/integer.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/../numeric/rounding.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/../numeric/../numeric/function.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/../numeric/../numeric/integer.h \
  ../../include/numeric/float64.h ../../include/numeric/../declarations.h \
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

