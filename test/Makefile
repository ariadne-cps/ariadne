TESTS = test_array test_binary_word test_denotable_set test_grid test_interval test_polyhedron test_rectangle test_state

include ../config.mk

all: ${TESTS}
	@test_array 
	@test_binary_word 
	@test_denotable_set 
	@test_grid 
	@test_interval 
	@test_rectangle 
	@test_state
	@test_polyhedron 


test_array: test_array.cc ../include/array.h
	${CXX} ${CXXFLAGS} ${LDXXFLAGS} test_array.cc -o test_array

test_binary_word: test_binary_word.cc ../include/binary_word.h
	${CXX} ${CXXFLAGS} ${LDXXFLAGS} test_binary_word.cc -o test_binary_word

test_denotable_set: test_denotable_set.cc ../include/denotable_set.h
	${CXX} ${CXXFLAGS} ${LDXXFLAGS} test_denotable_set.cc -o test_denotable_set

test_grid: test_grid.cc ../include/grid.h ../include/grid_rectangle.h ../include/grid_denotable_set.h
	${CXX} ${CXXFLAGS} ${LDXXFLAGS} test_grid.cc -o test_grid

test_interval: test_interval.cc ../include/interval.h
	${CXX} ${CXXFLAGS} ${LDXXFLAGS} test_interval.cc -o test_interval

test_polyhedron: test_polyhedron.cc ../include/polyhedron.h
	${CXX} ${CXXFLAGS} ${LDXXFLAGS} test_polyhedron.cc -o test_polyhedron

test_rectangle: test_rectangle.cc ../include/rectangle.h
	${CXX} ${CXXFLAGS} ${LDXXFLAGS} test_rectangle.cc -o test_rectangle

test_state: test_state.cc ../include/state.h
	${CXX} ${CXXFLAGS} ${LDXXFLAGS} test_state.cc -o test_state

clean:                                                                          
	rm -f *.bak *.o *~ *% core ${TESTS}
