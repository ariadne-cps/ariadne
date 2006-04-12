include ../../config.mk

OBJ =   vector.o interval_vector.o zonotopic_vector.o \
	matrix.o interval_matrix.o \
	tensor.o \
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
  ../../include/numeric/../numeric/dyadic.h \
  ../../include/numeric/../utility/stlio.h \
  ../../include/numeric/../utility/../base/array.h \
  ../../include/numeric/../utility/../base/../declarations.h \
  ../../include/linear_algebra/vector.h \
  ../../include/linear_algebra/../declarations.h \
  ../../include/linear_algebra/../numeric/integer.h \
  ../../include/linear_algebra/../numeric/../declarations.h \
  ../../include/linear_algebra/../numeric/../utility/stlio.h \
  ../../include/linear_algebra/vector.tpl \
  ../../include/linear_algebra/vector.h \
  ../../include/linear_algebra/../utility/stlio.h

interval_vector.o: interval_vector.cc ../../include/numeric/arithmetic.h \
  ../../include/numeric/../numeric/numerical_types.h \
  ../../include/numeric/../numeric/../declarations.h \
  ../../include/numeric/../numeric/../numeric/dyadic.h \
  ../../include/numeric/../numeric/../utility/stlio.h \
  ../../include/numeric/../numeric/../utility/../base/array.h \
  ../../include/numeric/../numeric/../utility/../base/../declarations.h \
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
  ../../include/linear_algebra/../linear_algebra/../numeric/../declarations.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/../utility/stlio.h \
  ../../include/linear_algebra/interval_vector.tpl \
  ../../include/linear_algebra/interval_vector.h

zonotopic_vector.o: zonotopic_vector.cc ../../include/real_typedef.h \
  ../../include/numeric/numerical_types.h \
  ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/dyadic.h \
  ../../include/numeric/../utility/stlio.h \
  ../../include/numeric/../utility/../base/array.h \
  ../../include/numeric/../utility/../base/../declarations.h \
  ../../include/linear_algebra/zonotopic_vector.h \
  ../../include/linear_algebra/../linear_algebra/vector.h \
  ../../include/linear_algebra/../linear_algebra/../declarations.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/integer.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/../declarations.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/../utility/stlio.h \
  ../../include/linear_algebra/../linear_algebra/matrix.h \
  ../../include/linear_algebra/../linear_algebra/interval_vector.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/interval.h \
  ../../include/linear_algebra/../linear_algebra/../linear_algebra/vector.h \
  ../../include/linear_algebra/../linear_algebra/interval_matrix.h \
  ../../include/linear_algebra/../linear_algebra/../linear_algebra/matrix.h \
  ../../include/linear_algebra/zonotopic_vector.tpl \
  ../../include/linear_algebra/zonotopic_vector.h

matrix.o: matrix.cc ../../include/real_typedef.h \
  ../../include/numeric/numerical_types.h \
  ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/dyadic.h \
  ../../include/numeric/../utility/stlio.h \
  ../../include/numeric/../utility/../base/array.h \
  ../../include/numeric/../utility/../base/../declarations.h \
  ../../include/linear_algebra/matrix.h \
  ../../include/linear_algebra/../declarations.h \
  ../../include/linear_algebra/../numeric/integer.h \
  ../../include/linear_algebra/../numeric/../declarations.h \
  ../../include/linear_algebra/../numeric/../utility/stlio.h \
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
  ../../include/linear_algebra/../linear_algebra/matrix.h

interval_matrix.o: interval_matrix.cc ../../include/real_typedef.h \
  ../../include/numeric/numerical_types.h \
  ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/dyadic.h \
  ../../include/numeric/../utility/stlio.h \
  ../../include/numeric/../utility/../base/array.h \
  ../../include/numeric/../utility/../base/../declarations.h \
  ../../include/linear_algebra/interval_matrix.h \
  ../../include/linear_algebra/../linear_algebra/matrix.h \
  ../../include/linear_algebra/../linear_algebra/../declarations.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/integer.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/../declarations.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/../utility/stlio.h \
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

tensor.o: tensor.cc ../../include/real_typedef.h \
  ../../include/numeric/numerical_types.h \
  ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/dyadic.h \
  ../../include/numeric/../utility/stlio.h \
  ../../include/numeric/../utility/../base/array.h \
  ../../include/numeric/../utility/../base/../declarations.h \
  ../../include/linear_algebra/tensor.h \
  ../../include/linear_algebra/../declarations.h \
  ../../include/linear_algebra/../base/array.h \
  ../../include/linear_algebra/tensor.tpl \
  ../../include/linear_algebra/tensor.h \
  ../../include/linear_algebra/../linear_algebra/vector.h \
  ../../include/linear_algebra/../linear_algebra/../declarations.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/integer.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/../declarations.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/../utility/stlio.h \
  ../../include/linear_algebra/../linear_algebra/matrix.h

linear_program.o: linear_program.cc ../../include/real_typedef.h \
  ../../include/numeric/numerical_types.h \
  ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/dyadic.h \
  ../../include/numeric/../utility/stlio.h \
  ../../include/numeric/../utility/../base/array.h \
  ../../include/numeric/../utility/../base/../declarations.h \
  ../../include/linear_algebra/linear_program.h \
  ../../include/linear_algebra/../linear_algebra/matrix.h \
  ../../include/linear_algebra/../linear_algebra/../declarations.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/integer.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/../declarations.h \
  ../../include/linear_algebra/../linear_algebra/../numeric/../utility/stlio.h \
  ../../include/linear_algebra/linear_program.tpl \
  ../../include/linear_algebra/linear_program.h \
  ../../include/linear_algebra/../linear_algebra/vector.h \
  ../../include/utility/stlio.h

