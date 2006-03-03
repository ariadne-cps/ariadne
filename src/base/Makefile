include ../../config.mk

OBJ = binary_tree.o array_operations.o

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
binary_tree.o: binary_tree.cc ../../include/base/arithmetic.h \
  ../../include/base/../base/numerical_type.h \
  ../../include/base/../base/../base/dyadic.h \
  ../../include/base/utility.h ../../include/base/../base/array.h \
  ../../include/base/binary_tree.h \
  ../../include/base/../base/basic_type.h \
  ../../include/base/../base/../base/array.h \
  ../../include/base/../base/../base/sequence.h \
  ../../include/base/../base/binary_word.h \
  ../../include/base/../base/../base/basic_type.h

array_operations.o: array_operations.cc \
  ../../include/base/array_operations.h \
  ../../include/base/../base/array.h \
  ../../include/base/../base/basic_type.h \
  ../../include/base/../base/../base/array.h \
  ../../include/base/../base/../base/sequence.h

