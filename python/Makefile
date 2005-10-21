include ../config.mk

CXXFLAGS = -I../${INCLUDEDIR} -O2 -fpic -Wall
LDXXFLAGS += -lboost_python -lsuperlu -lblas -lmpfr -lgmpxx -lgmp

OBJ=ariadne.o 

all: python

python: ariadne.so numeric.so interval.so function.so geometry.so

ariadne.so: ${OBJ}
	${CXX} -shared ${OBJ} ${LDXXFLAGS} -o ariadne.so

numeric.so: numeric.o
	${CXX} -shared numeric.o ${LDXXFLAGS} -o numeric.so

interval.so: interval.o
	${CXX} -shared interval.o ${LDXXFLAGS} -o interval.so

function.so: function.o
	${CXX} -shared function.o ${LDXXFLAGS} -o function.so

geometry.so: geometry.o
	${CXX} -shared geometry.o ${LDXXFLAGS} -o geometry.so

.cc.o: ariadne.cc interval.cc function.cc
	${CXX} ${CXXFLAGS} -I${PYTHONINCLUDEDIR} -c $<

clean:
	rm -f *.o *.so
