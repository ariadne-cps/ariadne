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
