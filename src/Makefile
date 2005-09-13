include ../config.mk

all: ${LIBS}

libariadne.so: grid_operations.o
	${CXX} ${LDXXFLAGS} -shared -o libariadne.so *.o

libariadne.a: grid_operations.o
	${CXX} ${LDXXFLAGS} -o libariadne.a *.o

# ${CXX} ${INCLUDE} -M 

.cc.o: 
	${CXX} ${INCLUDE} ${CXXFLAGS} -c $<

clean:                                                                          
	rm -f *.bak *.o *~ *% core libariadne.*
