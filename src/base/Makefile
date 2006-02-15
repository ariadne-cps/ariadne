include ../../config.mk

OBJ = binary_tree.o

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
  ../../include/base/basic_type.h ../../include/base/array.h \
  ../../include/base/sequence.h ../../include/base/binary_word.h \
  ../../include/base/utility.h ../../include/geometry/rectangle.h \
  ../../include/base/utility.h ../../include/base/interval.h \
  ../../include/base/binary_word.h ../../include/geometry/point.h \
  ../../include/linear_algebra/vector.h ../../include/base/basic_type.h \
  ../../include/base/numerical_type.h ../../include/base/dyadic.h \
  ../../include/base/interval.h \
  ../../include/geometry/geometry_declarations.h

