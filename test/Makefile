include ../config.mk

TESTS = test_array test_arithmetic test_binary_word \
  test_polyhedron test_rectangle test_point \
  test_linear_algebra test_interval test_list_set test_grid \
  test_lattice
#  test_partition_tree \
#  test_epsfstream \

CXXFLAGS += -Wall -I../include/ -I../
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
	@./test_polyhedron
	@./test_rectangle
	@./test_point
	@./test_linear_algebra
	@./test_interval
	@./test_list_set
	@./test_grid

test_array: test_array.cc test.h
	${CXX} ${CXXFLAGS} ${LDXXFLAGS} test_array.cc -o test_array

test_arithmetic: test_arithmetic.cc test.h
	${CXX} ${CXXFLAGS} ${LDXXFLAGS} test_arithmetic.cc -o test_arithmetic

test_binary_word: test_binary_word.cc test.h
	${CXX} ${CXXFLAGS} ${LDXXFLAGS} test_binary_word.cc -o test_binary_word

	${CXX} ${CXXFLAGS} ${LDXXFLAGS} test_polyhedron.cc -o test_polyhedron

	${CXX} ${CXXFLAGS} ${LDXXFLAGS} test_rectangle.cc -o test_rectangle

	${CXX} ${CXXFLAGS} ${LDXXFLAGS} test_point.cc -o test_point

	${CXX} ${CXXFLAGS} ${LDXXFLAGS} test_linear_algebra.cc -o test_linear_algebra

	${CXX} ${CXXFLAGS} ${LDXXFLAGS} test_interval.cc -o test_interval

	${CXX} ${CXXFLAGS} ${LDXXFLAGS} test_list_set.cc -o test_list_set

	${CXX} ${CXXFLAGS} ${LDXXFLAGS} test_grid.cc -o test_grid

