include ../../config.mk

OBJ = linear_program.o \

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
linear_program.o: linear_program.cc \
  ../../include/linear_algebra/linear_program.h \
  ../../include/linear_algebra/../base/basic_type.h \
  ../../include/linear_algebra/../base/../base/array.h \
  ../../include/linear_algebra/../base/../base/sequence.h \
  ../../include/linear_algebra/../linear_algebra/matrix.h \
  ../../include/linear_algebra/../linear_algebra/../base/basic_type.h \
  ../../include/linear_algebra/../linear_algebra/../base/numerical_type.h \
  ../../include/linear_algebra/../linear_algebra/../base/../base/dyadic.h \
  ../../include/linear_algebra/../linear_algebra/../base/interval.h \
  ../../include/linear_algebra/../linear_algebra/../base/utility.h \
  ../../include/linear_algebra/../linear_algebra/../base/../base/array.h \
  ../../include/linear_algebra/linear_program.tpl \
  ../../include/linear_algebra/linear_program.h \
  ../../include/base/utility.h

