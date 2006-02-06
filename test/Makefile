include ../config.mk

TESTS = test_array test_arithmetic test_binary_word \
  test_grid test_interval test_list_set \
  test_polyhedron test_rectangle test_point \
  test_linear_algebra
#  test_partition_tree \
#  test_epsfstream \

CXXFLAGS += -Wall
LDXXFLAGS += -L../src/ -lariadne -lgmp -lgmpxx -lppl 

dep:  
	@sed -e '/^### Do Not edit this line$$/,$$d' Makefile \
		> MMM.$$$$ && mv MMM.$$$$ Makefile
	@echo "### Do Not edit this line" >> Makefile
	@echo "### Everything below is auto-generated" >> Makefile
	@echo "" >> Makefile 
	@echo "all: $$""{TESTS}" >> Makefile
	@for ff in ${TESTS}; do echo "	@./$$ff"; done >> Makefile
	@echo "" >> Makefile 
	@for ff in ${TESTS}; do g++ -MM -MT$$ff ${CXXFLAGS} $$ff.cc; \
	echo "	$$""{CXX} $$""{CXXFLAGS} $$""{LDXXFLAGS} $$ff.cc -o $$ff"; \
	echo ""; done >> Makefile

clean:
	@rm -rf core ${TESTS}

### Do Not edit this line
### Everything below is auto-generated

all: ${TESTS}
	@./test_array
	@./test_arithmetic
	@./test_binary_word
	@./test_grid
	@./test_interval
	@./test_list_set
	@./test_polyhedron
	@./test_rectangle
	@./test_point
	@./test_linear_algebra

test_array: test_array.cc ../include/array.h ../include/utility.h \
  ../include/ariadne.h ../include/basic_type.h ../include/array.h \
  ../include/sequence.h ../include/numerical_type.h \
  ../include/synaps/dyadic.h test.h
	${CXX} ${CXXFLAGS} ${LDXXFLAGS} test_array.cc -o test_array

test_arithmetic: test_arithmetic.cc ../include/arithmetic.h \
  ../include/numerical_type.h ../include/synaps/dyadic.h test.h
	${CXX} ${CXXFLAGS} ${LDXXFLAGS} test_arithmetic.cc -o test_arithmetic

test_binary_word: test_binary_word.cc ../include/ariadne.h \
  ../include/basic_type.h ../include/array.h ../include/utility.h \
  ../include/sequence.h ../include/exception.h ../include/utility.h \
  ../include/numerical_type.h ../include/synaps/dyadic.h \
  ../include/point.h ../include/linear_algebra.h \
  ../include/numerical_type.h ../include/rectangle.h ../include/ariadne.h \
  ../include/geometry.h ../include/interval.h ../include/point.h \
  ../include/binary_word.h ../include/binary_word.h \
  ../include/binary_tree.h ../include/rectangle.h test.h
	${CXX} ${CXXFLAGS} ${LDXXFLAGS} test_binary_word.cc -o test_binary_word

test_grid: test_grid.cc ../include/ariadne.h ../include/basic_type.h \
  ../include/array.h ../include/utility.h ../include/sequence.h \
  ../include/exception.h ../include/utility.h ../include/numerical_type.h \
  ../include/synaps/dyadic.h ../include/point.h \
  ../include/linear_algebra.h ../include/numerical_type.h \
  ../include/rectangle.h ../include/ariadne.h ../include/geometry.h \
  ../include/interval.h ../include/point.h ../include/binary_word.h \
  ../include/binary_word.h ../include/grid_set.h ../include/rectangle.h \
  ../include/grid_operations.h ../include/list_set.h test.h
	${CXX} ${CXXFLAGS} ${LDXXFLAGS} test_grid.cc -o test_grid

test_interval: test_interval.cc ../include/ariadne.h \
  ../include/basic_type.h ../include/array.h ../include/utility.h \
  ../include/sequence.h ../include/numerical_type.h \
  ../include/synaps/dyadic.h ../include/interval.h test.h
	${CXX} ${CXXFLAGS} ${LDXXFLAGS} test_interval.cc -o test_interval

test_list_set: test_list_set.cc ../include/ariadne.h \
  ../include/basic_type.h ../include/array.h ../include/utility.h \
  ../include/sequence.h ../include/exception.h ../include/utility.h \
  ../include/numerical_type.h ../include/synaps/dyadic.h \
  ../include/point.h ../include/linear_algebra.h \
  ../include/numerical_type.h ../include/rectangle.h ../include/ariadne.h \
  ../include/geometry.h ../include/interval.h ../include/point.h \
  ../include/binary_word.h ../include/list_set.h test.h
	${CXX} ${CXXFLAGS} ${LDXXFLAGS} test_list_set.cc -o test_list_set

test_polyhedron: test_polyhedron.cc ../include/ariadne.h \
  ../include/basic_type.h ../include/array.h ../include/utility.h \
  ../include/sequence.h ../include/exception.h ../include/utility.h \
  ../include/numerical_type.h ../include/synaps/dyadic.h \
  ../include/point.h ../include/linear_algebra.h \
  ../include/numerical_type.h ../include/polyhedron.h \
  ../include/constraint.h ../include/ariadne.h ../include/generator.h \
  ../include/point.h ../include/rectangle.h ../include/geometry.h \
  ../include/interval.h ../include/binary_word.h
	${CXX} ${CXXFLAGS} ${LDXXFLAGS} test_polyhedron.cc -o test_polyhedron

test_rectangle: test_rectangle.cc ../include/ariadne.h \
  ../include/basic_type.h ../include/array.h ../include/utility.h \
  ../include/sequence.h ../include/exception.h ../include/utility.h \
  ../include/rectangle.h ../include/ariadne.h ../include/geometry.h \
  ../include/interval.h ../include/point.h ../include/linear_algebra.h \
  ../include/numerical_type.h ../include/synaps/dyadic.h \
  ../include/binary_word.h ../include/list_set.h test.h
	${CXX} ${CXXFLAGS} ${LDXXFLAGS} test_rectangle.cc -o test_rectangle

test_point: test_point.cc ../include/point.h ../include/utility.h \
  ../include/linear_algebra.h ../include/basic_type.h ../include/array.h \
  ../include/sequence.h ../include/numerical_type.h \
  ../include/synaps/dyadic.h test.h
	${CXX} ${CXXFLAGS} ${LDXXFLAGS} test_point.cc -o test_point

test_linear_algebra: test_linear_algebra.cc ../include/numerical_type.h \
  ../include/synaps/dyadic.h ../include/linear_algebra.h \
  ../include/basic_type.h ../include/array.h ../include/utility.h \
  ../include/sequence.h ../include/numerical_type.h test.h
	${CXX} ${CXXFLAGS} ${LDXXFLAGS} test_linear_algebra.cc -o test_linear_algebra

