include ../config.mk

all: test_array test_interval test_rectangle test_state
	@test_array;
	@test_interval;
	@test_rectangle;
	@test_state;

test_array: test_array.cc
	${CXX} ${CXXFLAGS} ${LDXXFLAGS} test_array.cc -o test_array

test_interval: test_interval.cc
	${CXX} ${CXXFLAGS} ${LDXXFLAGS} test_interval.cc -o test_interval

test_rectangle: test_rectangle.cc
	${CXX} ${CXXFLAGS} ${LDXXFLAGS} test_rectangle.cc -o test_rectangle

test_state: test_state.cc
	${CXX} ${CXXFLAGS} ${LDXXFLAGS} test_state.cc -o test_state

clean:                                                                          
	rm -f *.bak *.o *~ *% core test_array test_interval test_rectangle test_state
