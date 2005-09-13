include ../config.mk

CXXFLAGS = -I../${INCLUDEDIR} -O2 -fpic -Wall
LDXXFLAGS += -lboost_python

OBJ=ariadne.o

all: python

python: ariadne.so

ariadne.so: ${OBJ}
	${CXX} -shared ${OBJ} ${LDXXFLAGS} -o ariadne.so

.cc.o: ariadne.cc
	${CXX} ${CXXFLAGS} -I${PYTHONINCLUDEDIR} -c $<

clean:
	rm -f *.o *.so
