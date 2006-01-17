include ../config.mk

OBJ = binary_tree.o grid_operations.o partition_tree_operations.o

all: ${LIBS}

libariadne.so: ${OBJ}
	${CXX} ${LDXXFLAGS} -shared -o libariadne.so ${OBJ} ${LIBPPL} ${LIBGMPXX} ${RPATHFLAGS}
  
libariadne.a: ${OBJ}
	${CXX} ${LDXXFLAGS} -o libariadne.a ${OBJ} -lppl -lgmpxx -lgmp

# ${CXX} ${INCLUDE} -M 

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
binary_tree.o: binary_tree.cc ../include/binary_tree.h ../include/array.h \
  ../include/utility.h ../include/sequence.h ../include/utility.h \
  ../include/binary_word.h ../include/utility.h ../include/array.h \
  ../include/rectangle.h ../include/ariadne.h ../include/utility.h \
  ../include/interval.h ../include/state.h ../include/utility.h \
  ../include/linear_algebra.h ../include/numerical_type.h \
  ../include/binary_word.h

grid_operations.o: grid_operations.cc ../include/grid_operations.h \
  ../include/array.h ../include/utility.h ../include/rectangle.h \
  ../include/ariadne.h ../include/utility.h ../include/interval.h \
  ../include/state.h ../include/utility.h ../include/linear_algebra.h \
  ../include/numerical_type.h ../include/binary_word.h \
  ../include/utility.h ../include/array.h

partition_tree_operations.o: partition_tree_operations.cc \
  ../include/partition_tree_set.h ../include/array.h ../include/utility.h \
  ../include/sequence.h ../include/utility.h ../include/binary_word.h \
  ../include/utility.h ../include/array.h ../include/binary_tree.h \
  ../include/array.h ../include/sequence.h ../include/utility.h \
  ../include/binary_word.h ../include/rectangle.h ../include/ariadne.h \
  ../include/utility.h ../include/interval.h ../include/state.h \
  ../include/utility.h ../include/linear_algebra.h \
  ../include/numerical_type.h ../include/binary_word.h \
  ../include/rectangle.h ../include/list_set.h

