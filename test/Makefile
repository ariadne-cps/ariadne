include ../config.mk

TESTS = test_array test_arithmetic test_binary_word \
  test_grid test_interval test_list_set \
  test_partition_tree test_polyhedron test_rectangle test_state
#  test_epsfstream \

CXXFLAGS += -Wall
LDXXFLAGS += -L../src/ -lariadne

all: ${TESTS}
	@test_array
	@test_arithmetic
	@test_binary_word
#	@test_epsfstream
	@test_grid
	@test_interval
	@test_list_set
	@test_partition_tree
	test_rectangle
	@test_state
	@test_polyhedron


test_array: test_array.cc ../include/array.h
	${CXX} ${CXXFLAGS} -Wall test_array.cc -o test_array  ${LDXXFLAGS}

test_arithmetic: test_arithmetic.cc ../include/arithmetic.h
	${CXX} ${CXXFLAGS} -Wall test_arithmetic.cc -o test_arithmetic  ${LDXXFLAGS}

test_binary_word: test_binary_word.cc ../include/binary_word.h
	${CXX} ${CXXFLAGS} -Wall test_binary_word.cc -o test_binary_word  ${LDXXFLAGS}

test_denotable_set: test_denotable_set.cc ../include/denotable_set.h
	${CXX} ${CXXFLAGS} -Wall test_denotable_set.cc -o test_denotable_set  ${LDXXFLAGS}

test_epsfstream: test_epsfstream.cc ../include/epsfstream.h
	${CXX} ${CXXFLAGS} -Wall test_epsfstream.cc -o test_epsfstream  ${LDXXFLAGS}

test_grid: test_grid.cc ../include/grid_set.h ../include/grid_operations.h
	${CXX} ${CXXFLAGS} -Wall test_grid.cc -o test_grid  ${LDXXFLAGS}

test_interval: test_interval.cc ../include/interval.h
	${CXX} ${CXXFLAGS} -Wall test_interval.cc -o test_interval  ${LDXXFLAGS}

test_list_set: test_list_set.cc ../include/list_set.h
	${CXX} ${CXXFLAGS} -Wall test_list_set.cc -o test_list_set  ${LDXXFLAGS}

test_partition_tree: test_partition_tree.cc ../include/partition_tree.h
	${CXX} ${CXXFLAGS} -Wall test_partition_tree.cc -o test_partition_tree  ${LDXXFLAGS}

test_polyhedron: test_polyhedron.cc ../include/polyhedron.h
	${CXX} ${CXXFLAGS} -Wall test_polyhedron.cc -o test_polyhedron  ${LDXXFLAGS}

test_rectangle: test_rectangle.cc ../include/rectangle.h
	${CXX} ${CXXFLAGS} -Wall test_rectangle.cc -o test_rectangle  ${LDXXFLAGS}

test_state: test_state.cc ../include/state.h
	${CXX} ${CXXFLAGS} -Wall test_state.cc -o test_state 

clean:
	rm -f *.bak *.log *.o *~ *% core ${TESTS}
