include ../config.mk

OBJ = binary_tree.o grid_operations.o partition_tree_operations.o

all: ${LIBS}

libariadne.so: ${OBJ}
	${CXX} ${LDXXFLAGS} -shared -o libariadne.so ${OBJ}

libariadne.a: ${OBJ}
	${CXX} ${LDXXFLAGS} -o libariadne.a ${OBJ}

# ${CXX} ${INCLUDE} -M 

.cc.o: 
	${CXX} ${INCLUDE} ${CXXFLAGS} -c $<

clean:                                                                          
	rm -f *.bak *.o *~ *% core libariadne.*
