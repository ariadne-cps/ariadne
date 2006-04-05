include ../../config.mk

OBJ = \
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
matrix.o: matrix.cc ../../include/linear_algebra/matrix.h \
  ../../include/linear_algebra/../base/basic_type.h \
  ../../include/linear_algebra/../base/../base/array.h \
  ../../include/linear_algebra/../base/../base/sequence.h \
  ../../include/linear_algebra/../base/numerical_type.h \
  ../../include/linear_algebra/../base/../base/dyadic.h \
  ../../include/linear_algebra/../base/arithmetic.h \
  ../../include/linear_algebra/../base/../base/numerical_type.h \
  ../../include/linear_algebra/matrix.tpl \
  ../../include/linear_algebra/../linear_algebra/matrix.h

interval_matrix.o: interval_matrix.cc \
  ../../include/linear_algebra/interval_matrix.h \
  ../../include/linear_algebra/../base/basic_type.h \
  ../../include/linear_algebra/../base/../base/array.h \
  ../../include/linear_algebra/../base/../base/sequence.h \
  ../../include/linear_algebra/../base/numerical_type.h \
  ../../include/linear_algebra/../base/../base/dyadic.h \
  ../../include/linear_algebra/../base/interval.h \
  ../../include/linear_algebra/../linear_algebra/vector.h \
  ../../include/linear_algebra/../linear_algebra/../base/basic_type.h \
  ../../include/linear_algebra/../linear_algebra/../base/numerical_type.h \
  ../../include/linear_algebra/../linear_algebra/../base/interval.h \
  ../../include/linear_algebra/../linear_algebra/matrix.h \
  ../../include/linear_algebra/../linear_algebra/../base/arithmetic.h \
  ../../include/linear_algebra/../linear_algebra/../base/../base/numerical_type.h \
  ../../include/linear_algebra/../linear_algebra/interval_vector.h \
  ../../include/linear_algebra/interval_matrix.tpl \
  ../../include/linear_algebra/../linear_algebra/interval_matrix.h

tensor.o: tensor.cc ../../include/linear_algebra/tensor.h \
  ../../include/linear_algebra/../base/basic_type.h \
  ../../include/linear_algebra/../base/../base/array.h \
  ../../include/linear_algebra/../base/../base/sequence.h \
  ../../include/linear_algebra/../base/array.h \
  ../../include/linear_algebra/../linear_algebra/vector.h \
  ../../include/linear_algebra/../linear_algebra/../base/basic_type.h \
  ../../include/linear_algebra/../linear_algebra/../base/numerical_type.h \
  ../../include/linear_algebra/../linear_algebra/../base/../base/dyadic.h \
  ../../include/linear_algebra/../linear_algebra/../base/interval.h \
  ../../include/linear_algebra/../linear_algebra/matrix.h \
  ../../include/linear_algebra/../linear_algebra/../base/arithmetic.h \
  ../../include/linear_algebra/../linear_algebra/../base/../base/numerical_type.h \
  ../../include/linear_algebra/tensor.tpl \
  ../../include/linear_algebra/tensor.h \
  ../../include/linear_algebra/../base/numerical_type.h

linear_program.o: linear_program.cc \
  ../../include/linear_algebra/linear_program.h \
  ../../include/linear_algebra/../base/basic_type.h \
  ../../include/linear_algebra/../base/../base/array.h \
  ../../include/linear_algebra/../base/../base/sequence.h \
  ../../include/linear_algebra/../linear_algebra/matrix.h \
  ../../include/linear_algebra/../linear_algebra/../base/basic_type.h \
  ../../include/linear_algebra/../linear_algebra/../base/numerical_type.h \
  ../../include/linear_algebra/../linear_algebra/../base/../base/dyadic.h \
  ../../include/linear_algebra/../linear_algebra/../base/arithmetic.h \
  ../../include/linear_algebra/../linear_algebra/../base/../base/numerical_type.h \
  ../../include/linear_algebra/linear_program.tpl \
  ../../include/linear_algebra/linear_program.h \
  ../../include/base/utility.h ../../include/base/../base/array.h

