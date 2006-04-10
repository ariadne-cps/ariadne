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
binary_tree.o: binary_tree.cc ../../include/base/binary_tree.h \
  ../../include/base/../declarations.h ../../include/base/../base/array.h \
  ../../include/base/../base/../declarations.h \
  ../../include/base/../base/binary_word.h ../../include/utility/stlio.h \
  ../../include/utility/../base/array.h ../../include/base/array.h \
  ../../include/numeric/arithmetic.h \
  ../../include/numeric/../numeric/numerical_types.h \
  ../../include/numeric/../numeric/../declarations.h \
  ../../include/numeric/../numeric/../base/utility.h \
  ../../include/numeric/../numeric/../base/../declarations.h \
  ../../include/numeric/../numeric/../numeric/dyadic.h \
  ../../include/numeric/../numeric/approximation.h \
  ../../include/numeric/../numeric/../numeric/numerical_types.h

array_operations.o: array_operations.cc \
  ../../include/base/array_operations.h \
  ../../include/base/../declarations.h ../../include/base/array.h

