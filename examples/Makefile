include ../config.mk

all: model2

model2: model2.o
	${CXX} ${LDXXFLAGS} model2.o -o model2


.cc.o: model2.cc
	${CXX} ${CXXFLAGS} -c $<

clean:                                                                          
	rm -f *.bak *.o *~ *% core model2
