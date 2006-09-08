include ../config.mk

TESTS = \
  test_array test_binary_word \
  test_lattice \
  test_mpfr test_arithmetic test_function test_interval \
  test_vector test_matrix test_linear_algebra \
  test_point \
  test_rectangle test_parallelotope test_zonotope \
  test_simplex test_polyhedron \
  test_list_set test_grid_set test_partition_tree_set \
  test_epsfstream \
  test_interval_newton \
  test_chainreach \
  test_integration_step 

OBJS=$(addsuffix .o,$(TESTS))

CXXFLAGS += -Wall -I../include/ 
LDLIBS = ../src/libariadne.so -lmpfr -lgmp -lgmpxx -lppl #-lblas

all: ${TESTS}

check: $(TESTS)
	@(for test in $(TESTS); do ./$$test; done);

$(TESTS): %: %.o ../src/libariadne.so
	$(CXX) $@.o $(LDXXFLAGS) $(LDLIBS) -o $@

.cc.o:
	${CXX} ${INCLUDE} ${CXXFLAGS} -c $<

clean:
	-rm -rf core *.log *.eps *.o ${TESTS} 

dep:  
	@sed -e '/^### Do Not edit this line$$/,$$d' Makefile \
		> MMM.$$$$ && mv MMM.$$$$ Makefile
	@echo "### Do Not edit this line" >> Makefile
	@echo "### Everything below is auto-generated" >> Makefile
	@echo "" >> Makefile 
	@for ff in ${TESTS}; do g++ -MM ${CXXFLAGS} $$ff.cc; \
	echo "	-$$""{CXX} $$""{CXXFLAGS} -c $$ff.cc "; \
	echo ""; done >> Makefile

depclean:
	@sed -e '/^### Do Not edit this line$$/,$$d' Makefile \
		> MMM.$$$$ && mv MMM.$$$$ Makefile
	@echo "### Do Not edit this line" >> Makefile

### Do Not edit this line
### Everything below is auto-generated

test_array.o: test_array.cc ../include/ariadne.h ../include/base/stlio.h \
  ../include/base/../base/array.h \
  ../include/base/../base/../declarations.h ../include/base/array.h \
  ../include/numeric/numerical_types.h ../include/numeric/mpfloat.h \
  ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/function.h \
  ../include/numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/rational.h \
  ../include/numeric/../numeric/../numeric/arithmetic.h \
  ../include/numeric/../numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../numeric/function.h \
  ../include/numeric/../numeric/../numeric/integer.h test.h
	-${CXX} ${CXXFLAGS} -c test_array.cc 

test_binary_word.o: test_binary_word.cc ../include/ariadne.h \
  ../include/base/stlio.h ../include/base/../base/array.h \
  ../include/base/../base/../declarations.h ../include/base/exception.h \
  ../include/base/binary_word.h ../include/base/../declarations.h \
  ../include/base/binary_tree.h ../include/base/../base/binary_word.h \
  ../include/numeric/numerical_types.h ../include/numeric/mpfloat.h \
  ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/function.h \
  ../include/numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/rational.h \
  ../include/numeric/../numeric/../numeric/arithmetic.h \
  ../include/numeric/../numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../numeric/function.h \
  ../include/numeric/../numeric/../numeric/integer.h \
  ../include/base/utility.h ../include/geometry/point.h \
  ../include/geometry/../declarations.h \
  ../include/geometry/../base/array.h \
  ../include/geometry/../linear_algebra/vector.h \
  ../include/geometry/../linear_algebra/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/integer.h \
  ../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../include/geometry/../linear_algebra/matrix.h \
  ../include/geometry/rectangle.h \
  ../include/geometry/../numeric/arithmetic.h \
  ../include/geometry/../numeric/interval.h \
  ../include/geometry/../numeric/../declarations.h \
  ../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../numeric/../numeric/arithmetic.h \
  ../include/geometry/../numeric/../numeric/rounding.h \
  ../include/geometry/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../numeric/../numeric/../numeric/function.h \
  ../include/geometry/../linear_algebra/interval_vector.h \
  ../include/geometry/../linear_algebra/../numeric/interval.h \
  ../include/geometry/../linear_algebra/../linear_algebra/vector.h \
  ../include/geometry/../geometry/ppl_polyhedron.h \
  ../include/geometry/../geometry/../declarations.h \
  ../include/geometry/../geometry/../numeric/rational.h \
  ../include/geometry/../geometry/../linear_algebra/vector.h \
  ../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../include/geometry/../geometry/../linear_algebra/interval_vector.h \
  ../include/geometry/../geometry/point.h test.h
	-${CXX} ${CXXFLAGS} -c test_binary_word.cc 

test_lattice.o: test_lattice.cc ../include/declarations.h \
  ../include/real_typedef.h ../include/numeric/numerical_types.h \
  ../include/numeric/mpfloat.h ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/function.h \
  ../include/numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/rational.h \
  ../include/numeric/../numeric/../numeric/arithmetic.h \
  ../include/numeric/../numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../numeric/function.h \
  ../include/numeric/../numeric/../numeric/integer.h \
  ../include/geometry/lattice_set.h ../include/geometry/../base/array.h \
  ../include/geometry/../base/../declarations.h \
  ../include/geometry/../base/iterator.h \
  ../include/geometry/../base/../base/array.h \
  ../include/geometry/../numeric/interval.h \
  ../include/geometry/../numeric/../declarations.h \
  ../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../numeric/../numeric/arithmetic.h \
  ../include/geometry/../numeric/../numeric/rounding.h \
  ../include/geometry/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../numeric/../numeric/../numeric/function.h \
  ../include/geometry/../base/array_operations.h \
  ../include/system/lattice_map.h \
  ../include/system/../geometry/lattice_set.h
	-${CXX} ${CXXFLAGS} -c test_lattice.cc 

test_mpfr.o: test_mpfr.cc
	-${CXX} ${CXXFLAGS} -c test_mpfr.cc 

test_arithmetic.o: test_arithmetic.cc ../include/numeric/arithmetic.h \
  ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/interval.h \
  ../include/numeric/../numeric/arithmetic.h \
  ../include/numeric/../numeric/rounding.h \
  ../include/numeric/../numeric/../numeric/function.h \
  ../include/numeric/../numeric/../numeric/integer.h \
  ../include/numeric/float64.h ../include/numeric/../numeric/function.h \
  ../include/numeric/mpfloat.h ../include/numeric/../numeric/rational.h \
  ../include/numeric/../numeric/../numeric/arithmetic.h \
  ../include/numeric/../numeric/../numeric/integer.h \
  ../include/numeric/rational.h test.h
	-${CXX} ${CXXFLAGS} -c test_arithmetic.cc 

test_function.o: test_function.cc ../include/numeric/function.h \
  ../include/numeric/integer.h ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/float64.h ../include/numeric/../numeric/function.h \
  ../include/numeric/../numeric/integer.h ../include/numeric/mpfloat.h \
  ../include/numeric/../numeric/rational.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../numeric/arithmetic.h \
  ../include/numeric/../numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../numeric/function.h \
  ../include/numeric/../numeric/../numeric/integer.h test.h
	-${CXX} ${CXXFLAGS} -c test_function.cc 

test_interval.o: test_interval.cc ../include/numeric/interval.h \
  ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/arithmetic.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/rounding.h \
  ../include/numeric/../numeric/../numeric/function.h \
  ../include/numeric/../numeric/../numeric/integer.h \
  ../include/numeric/float64.h ../include/numeric/../numeric/function.h \
  ../include/numeric/../numeric/integer.h ../include/numeric/mpfloat.h \
  ../include/numeric/../numeric/rational.h \
  ../include/numeric/../numeric/../numeric/arithmetic.h \
  ../include/numeric/rational.h test.h
	-${CXX} ${CXXFLAGS} -c test_interval.cc 

test_vector.o: test_vector.cc ../include/real_typedef.h \
  ../include/numeric/numerical_types.h ../include/numeric/mpfloat.h \
  ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/function.h \
  ../include/numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/rational.h \
  ../include/numeric/../numeric/../numeric/arithmetic.h \
  ../include/numeric/../numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../numeric/function.h \
  ../include/numeric/../numeric/../numeric/integer.h \
  ../include/numeric/numerical_types.h ../include/linear_algebra/vector.h \
  ../include/linear_algebra/../declarations.h \
  ../include/linear_algebra/../numeric/integer.h \
  ../include/linear_algebra/../numeric/arithmetic.h \
  ../include/linear_algebra/interval_vector.h \
  ../include/linear_algebra/../numeric/interval.h \
  ../include/linear_algebra/../numeric/../declarations.h \
  ../include/linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../include/linear_algebra/../numeric/../numeric/arithmetic.h \
  ../include/linear_algebra/../numeric/../numeric/rounding.h \
  ../include/linear_algebra/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/linear_algebra/../numeric/../numeric/../numeric/function.h \
  ../include/linear_algebra/../linear_algebra/vector.h \
  ../include/linear_algebra/vector.tpl ../include/linear_algebra/vector.h \
  ../include/linear_algebra/../base/stlio.h \
  ../include/linear_algebra/../base/../base/array.h \
  ../include/linear_algebra/../base/../base/../declarations.h
	-${CXX} ${CXXFLAGS} -c test_vector.cc 

test_matrix.o: test_matrix.cc ../include/declarations.h \
  ../include/real_typedef.h ../include/numeric/numerical_types.h \
  ../include/numeric/mpfloat.h ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/function.h \
  ../include/numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/rational.h \
  ../include/numeric/../numeric/../numeric/arithmetic.h \
  ../include/numeric/../numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../numeric/function.h \
  ../include/numeric/../numeric/../numeric/integer.h \
  ../include/numeric/numerical_types.h ../include/linear_algebra/vector.h \
  ../include/linear_algebra/../declarations.h \
  ../include/linear_algebra/../numeric/integer.h \
  ../include/linear_algebra/../numeric/arithmetic.h \
  ../include/linear_algebra/vector.tpl ../include/linear_algebra/vector.h \
  ../include/linear_algebra/../base/stlio.h \
  ../include/linear_algebra/../base/../base/array.h \
  ../include/linear_algebra/../base/../base/../declarations.h \
  ../include/linear_algebra/matrix.h ../include/linear_algebra/matrix.tpl \
  ../include/linear_algebra/matrix.h \
  ../include/linear_algebra/../base/array.h \
  ../include/linear_algebra/../linear_algebra/vector.h \
  ../include/linear_algebra/../linear_algebra/matrix.h \
  ../include/linear_algebra/../linear_algebra/lu_matrix.h \
  ../include/linear_algebra/../linear_algebra/../declarations.h \
  ../include/linear_algebra/../linear_algebra/../base/array.h \
  ../include/linear_algebra/../linear_algebra/../linear_algebra/matrix.h \
  ../include/linear_algebra/interval_matrix.h \
  ../include/linear_algebra/../numeric/interval.h \
  ../include/linear_algebra/../numeric/../declarations.h \
  ../include/linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../include/linear_algebra/../numeric/../numeric/arithmetic.h \
  ../include/linear_algebra/../numeric/../numeric/rounding.h \
  ../include/linear_algebra/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/linear_algebra/../numeric/../numeric/../numeric/function.h \
  ../include/linear_algebra/interval_matrix.tpl ../include/declarations.h \
  ../include/linear_algebra/interval_matrix.h \
  ../include/linear_algebra/../numeric/approximation.h \
  ../include/linear_algebra/../linear_algebra/interval_vector.h \
  ../include/linear_algebra/../linear_algebra/../numeric/interval.h \
  ../include/linear_algebra/../linear_algebra/../linear_algebra/vector.h
	-${CXX} ${CXXFLAGS} -c test_matrix.cc 

test_linear_algebra.o: test_linear_algebra.cc ../include/real_typedef.h \
  ../include/numeric/numerical_types.h ../include/numeric/mpfloat.h \
  ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/function.h \
  ../include/numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/rational.h \
  ../include/numeric/../numeric/../numeric/arithmetic.h \
  ../include/numeric/../numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../numeric/function.h \
  ../include/numeric/../numeric/../numeric/integer.h \
  ../include/numeric/numerical_types.h ../include/linear_algebra/vector.h \
  ../include/linear_algebra/../declarations.h \
  ../include/linear_algebra/../numeric/integer.h \
  ../include/linear_algebra/../numeric/arithmetic.h \
  ../include/linear_algebra/matrix.h \
  ../include/linear_algebra/interval_matrix.h \
  ../include/linear_algebra/../linear_algebra/matrix.h \
  ../include/linear_algebra/../numeric/interval.h \
  ../include/linear_algebra/../numeric/../declarations.h \
  ../include/linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../include/linear_algebra/../numeric/../numeric/arithmetic.h \
  ../include/linear_algebra/../numeric/../numeric/rounding.h \
  ../include/linear_algebra/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/linear_algebra/../numeric/../numeric/../numeric/function.h \
  ../include/linear_algebra/lu_matrix.h \
  ../include/linear_algebra/../base/array.h \
  ../include/linear_algebra/../base/../declarations.h \
  ../include/linear_algebra/qr_matrix.h \
  ../include/linear_algebra/svd_matrix.h \
  ../include/linear_algebra/../linear_algebra/vector.h test.h
	-${CXX} ${CXXFLAGS} -c test_linear_algebra.cc 

test_point.o: test_point.cc ../include/real_typedef.h \
  ../include/numeric/numerical_types.h ../include/numeric/mpfloat.h \
  ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/function.h \
  ../include/numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/rational.h \
  ../include/numeric/../numeric/../numeric/arithmetic.h \
  ../include/numeric/../numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../numeric/function.h \
  ../include/numeric/../numeric/../numeric/integer.h \
  ../include/geometry/point.h ../include/geometry/../declarations.h \
  ../include/geometry/../base/array.h \
  ../include/geometry/../base/../declarations.h \
  ../include/geometry/../linear_algebra/vector.h \
  ../include/geometry/../linear_algebra/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/integer.h \
  ../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../include/geometry/../linear_algebra/matrix.h test.h
	-${CXX} ${CXXFLAGS} -c test_point.cc 

test_rectangle.o: test_rectangle.cc ../include/ariadne.h \
  ../include/real_typedef.h ../include/numeric/numerical_types.h \
  ../include/numeric/mpfloat.h ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/function.h \
  ../include/numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/rational.h \
  ../include/numeric/../numeric/../numeric/arithmetic.h \
  ../include/numeric/../numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../numeric/function.h \
  ../include/numeric/../numeric/../numeric/integer.h \
  ../include/base/exception.h ../include/base/utility.h \
  ../include/base/../declarations.h ../include/geometry/rectangle.h \
  ../include/geometry/../declarations.h \
  ../include/geometry/../numeric/arithmetic.h \
  ../include/geometry/../numeric/interval.h \
  ../include/geometry/../numeric/../declarations.h \
  ../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../numeric/../numeric/arithmetic.h \
  ../include/geometry/../numeric/../numeric/rounding.h \
  ../include/geometry/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../numeric/../numeric/../numeric/function.h \
  ../include/geometry/../linear_algebra/interval_vector.h \
  ../include/geometry/../linear_algebra/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/interval.h \
  ../include/geometry/../linear_algebra/../linear_algebra/vector.h \
  ../include/geometry/../linear_algebra/../linear_algebra/../declarations.h \
  ../include/geometry/../linear_algebra/../linear_algebra/../numeric/integer.h \
  ../include/geometry/../linear_algebra/../linear_algebra/../numeric/arithmetic.h \
  ../include/geometry/../geometry/ppl_polyhedron.h \
  ../include/geometry/../geometry/../declarations.h \
  ../include/geometry/../geometry/../numeric/rational.h \
  ../include/geometry/../geometry/../linear_algebra/vector.h \
  ../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../include/geometry/../geometry/../linear_algebra/../declarations.h \
  ../include/geometry/../geometry/../linear_algebra/../numeric/integer.h \
  ../include/geometry/../geometry/../linear_algebra/interval_vector.h \
  ../include/geometry/../geometry/point.h \
  ../include/geometry/../geometry/../base/array.h \
  ../include/geometry/../geometry/../base/../declarations.h \
  ../include/geometry/list_set.h ../include/geometry/../base/utility.h \
  test.h
	-${CXX} ${CXXFLAGS} -c test_rectangle.cc 

test_parallelotope.o: test_parallelotope.cc ../include/ariadne.h \
  ../include/real_typedef.h ../include/numeric/numerical_types.h \
  ../include/numeric/mpfloat.h ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/function.h \
  ../include/numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/rational.h \
  ../include/numeric/../numeric/../numeric/arithmetic.h \
  ../include/numeric/../numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../numeric/function.h \
  ../include/numeric/../numeric/../numeric/integer.h \
  ../include/base/exception.h ../include/base/utility.h \
  ../include/base/../declarations.h ../include/numeric/numerical_types.h \
  ../include/geometry/point.h ../include/geometry/../declarations.h \
  ../include/geometry/../base/array.h \
  ../include/geometry/../base/../declarations.h \
  ../include/geometry/../linear_algebra/vector.h \
  ../include/geometry/../linear_algebra/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/integer.h \
  ../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../include/geometry/../linear_algebra/matrix.h \
  ../include/geometry/rectangle.h \
  ../include/geometry/../numeric/arithmetic.h \
  ../include/geometry/../numeric/interval.h \
  ../include/geometry/../numeric/../declarations.h \
  ../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../numeric/../numeric/arithmetic.h \
  ../include/geometry/../numeric/../numeric/rounding.h \
  ../include/geometry/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../numeric/../numeric/../numeric/function.h \
  ../include/geometry/../linear_algebra/interval_vector.h \
  ../include/geometry/../linear_algebra/../numeric/interval.h \
  ../include/geometry/../linear_algebra/../linear_algebra/vector.h \
  ../include/geometry/../geometry/ppl_polyhedron.h \
  ../include/geometry/../geometry/../declarations.h \
  ../include/geometry/../geometry/../numeric/rational.h \
  ../include/geometry/../geometry/../linear_algebra/vector.h \
  ../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../include/geometry/../geometry/../linear_algebra/interval_vector.h \
  ../include/geometry/../geometry/point.h \
  ../include/geometry/parallelotope.h \
  ../include/geometry/../geometry/zonotope.h \
  ../include/geometry/../geometry/../numeric/interval.h \
  ../include/geometry/../geometry/../geometry/ppl_polyhedron.h \
  ../include/geometry/../geometry/../geometry/point.h
	-${CXX} ${CXXFLAGS} -c test_parallelotope.cc 

test_zonotope.o: test_zonotope.cc ../include/ariadne.h \
  ../include/real_typedef.h ../include/numeric/numerical_types.h \
  ../include/numeric/mpfloat.h ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/function.h \
  ../include/numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/rational.h \
  ../include/numeric/../numeric/../numeric/arithmetic.h \
  ../include/numeric/../numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../numeric/function.h \
  ../include/numeric/../numeric/../numeric/integer.h \
  ../include/base/exception.h ../include/base/utility.h \
  ../include/base/../declarations.h ../include/numeric/numerical_types.h \
  ../include/geometry/point.h ../include/geometry/../declarations.h \
  ../include/geometry/../base/array.h \
  ../include/geometry/../base/../declarations.h \
  ../include/geometry/../linear_algebra/vector.h \
  ../include/geometry/../linear_algebra/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/integer.h \
  ../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../include/geometry/../linear_algebra/matrix.h \
  ../include/geometry/rectangle.h \
  ../include/geometry/../numeric/arithmetic.h \
  ../include/geometry/../numeric/interval.h \
  ../include/geometry/../numeric/../declarations.h \
  ../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../numeric/../numeric/arithmetic.h \
  ../include/geometry/../numeric/../numeric/rounding.h \
  ../include/geometry/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../numeric/../numeric/../numeric/function.h \
  ../include/geometry/../linear_algebra/interval_vector.h \
  ../include/geometry/../linear_algebra/../numeric/interval.h \
  ../include/geometry/../linear_algebra/../linear_algebra/vector.h \
  ../include/geometry/../geometry/ppl_polyhedron.h \
  ../include/geometry/../geometry/../declarations.h \
  ../include/geometry/../geometry/../numeric/rational.h \
  ../include/geometry/../geometry/../linear_algebra/vector.h \
  ../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../include/geometry/../geometry/../linear_algebra/interval_vector.h \
  ../include/geometry/../geometry/point.h \
  ../include/geometry/parallelotope.h \
  ../include/geometry/../geometry/zonotope.h \
  ../include/geometry/../geometry/../numeric/interval.h \
  ../include/geometry/../geometry/../geometry/ppl_polyhedron.h \
  ../include/geometry/../geometry/../geometry/point.h \
  ../include/geometry/zonotope.h ../include/output/epsfstream.h \
  ../include/output/../numeric/approximation.h \
  ../include/output/../numeric/../declarations.h \
  ../include/output/../numeric/../numeric/numerical_traits.h \
  ../include/output/../geometry/point.h \
  ../include/output/../geometry/rectangle.h \
  ../include/output/../geometry/list_set.h \
  ../include/output/../geometry/../declarations.h \
  ../include/output/../geometry/../base/utility.h \
  ../include/output/../geometry/grid_set.h \
  ../include/output/../geometry/../base/array.h \
  ../include/output/../geometry/../geometry/rectangle.h \
  ../include/output/../geometry/../geometry/grid.h \
  ../include/output/../geometry/../geometry/../declarations.h \
  ../include/output/../geometry/../geometry/../base/array.h \
  ../include/output/../geometry/../geometry/../numeric/arithmetic.h \
  ../include/output/../geometry/../geometry/../geometry/lattice_set.h \
  ../include/output/../geometry/../geometry/../geometry/../base/array.h \
  ../include/output/../geometry/../geometry/../geometry/../base/iterator.h \
  ../include/output/../geometry/../geometry/../geometry/../base/../base/array.h \
  ../include/output/../geometry/../geometry/../geometry/../numeric/interval.h \
  ../include/output/../geometry/../geometry/../geometry/../base/array_operations.h \
  ../include/output/../geometry/../geometry/../geometry/../base/../declarations.h \
  ../include/output/../geometry/../geometry/lattice_set.h \
  ../include/output/../geometry/parallelotope.h \
  ../include/output/../geometry/zonotope.h \
  ../include/output/../geometry/polyhedron.h \
  ../include/output/../geometry/../linear_algebra/vector.h \
  ../include/output/../geometry/../linear_algebra/matrix.h \
  ../include/output/../geometry/../linear_algebra/interval_vector.h \
  ../include/output/../geometry/../geometry/ppl_polyhedron.h \
  ../include/output/../geometry/../geometry/point.h \
  ../include/output/../geometry/../geometry/point_list.h \
  ../include/output/../geometry/../geometry/../linear_algebra/vector.h \
  ../include/output/../geometry/../geometry/../linear_algebra/matrix.h \
  ../include/output/../geometry/../geometry/../geometry/point.h \
  ../include/output/../geometry/partition_tree_set.h \
  ../include/output/../geometry/../base/binary_word.h \
  ../include/output/../geometry/../base/../declarations.h \
  ../include/output/../geometry/../base/binary_tree.h \
  ../include/output/../geometry/../base/../base/array.h \
  ../include/output/../geometry/../base/../base/binary_word.h \
  ../include/output/../geometry/../base/iterator.h \
  ../include/output/../geometry/../geometry/subdivision_tree_set.h \
  ../include/output/../geometry/../geometry/../base/sequence.h \
  ../include/output/../geometry/../geometry/../base/binary_word.h \
  ../include/output/../geometry/../geometry/../base/binary_tree.h \
  ../include/output/../geometry/../geometry/../base/iterator.h \
  ../include/output/../geometry/../geometry/../geometry/rectangle.h \
  ../include/output/../linear_algebra/matrix.h \
  ../include/output/../system/affine_map.h \
  ../include/output/../system/../declarations.h \
  ../include/output/../system/../linear_algebra/vector.h \
  ../include/output/../system/../linear_algebra/matrix.h \
  ../include/output/../system/../system/map.h \
  ../include/output/../system/../system/../declarations.h test.h
	-${CXX} ${CXXFLAGS} -c test_zonotope.cc 

test_simplex.o: test_simplex.cc ../include/ariadne.h \
  ../include/real_typedef.h ../include/numeric/numerical_types.h \
  ../include/numeric/mpfloat.h ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/function.h \
  ../include/numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/rational.h \
  ../include/numeric/../numeric/../numeric/arithmetic.h \
  ../include/numeric/../numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../numeric/function.h \
  ../include/numeric/../numeric/../numeric/integer.h \
  ../include/base/exception.h ../include/base/utility.h \
  ../include/base/../declarations.h ../include/numeric/numerical_types.h \
  ../include/geometry/point.h ../include/geometry/../declarations.h \
  ../include/geometry/../base/array.h \
  ../include/geometry/../base/../declarations.h \
  ../include/geometry/../linear_algebra/vector.h \
  ../include/geometry/../linear_algebra/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/integer.h \
  ../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../include/geometry/../linear_algebra/matrix.h \
  ../include/geometry/simplex.h ../include/geometry/../geometry/point.h \
  ../include/geometry/../geometry/point_list.h \
  ../include/geometry/../geometry/../declarations.h \
  ../include/geometry/../geometry/../linear_algebra/vector.h \
  ../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../include/geometry/../geometry/../geometry/point.h \
  ../include/geometry/../geometry/polyhedron.h \
  ../include/geometry/../geometry/../linear_algebra/interval_vector.h \
  ../include/geometry/../geometry/../linear_algebra/../declarations.h \
  ../include/geometry/../geometry/../linear_algebra/../numeric/interval.h \
  ../include/geometry/../geometry/../linear_algebra/../numeric/../declarations.h \
  ../include/geometry/../geometry/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../geometry/../linear_algebra/../numeric/../numeric/arithmetic.h \
  ../include/geometry/../geometry/../linear_algebra/../numeric/../numeric/rounding.h \
  ../include/geometry/../geometry/../linear_algebra/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../geometry/../linear_algebra/../numeric/../numeric/../numeric/function.h \
  ../include/geometry/../geometry/../linear_algebra/../linear_algebra/vector.h \
  ../include/geometry/../geometry/../geometry/ppl_polyhedron.h \
  ../include/geometry/../geometry/../geometry/../declarations.h \
  ../include/geometry/../geometry/../geometry/../numeric/rational.h \
  ../include/geometry/../geometry/../geometry/../linear_algebra/vector.h \
  ../include/geometry/../geometry/../geometry/../linear_algebra/matrix.h \
  ../include/geometry/../geometry/../geometry/../linear_algebra/interval_vector.h \
  ../include/geometry/../geometry/../geometry/point_list.h \
  ../include/geometry/../geometry/../geometry/rectangle.h \
  ../include/geometry/../geometry/../geometry/../numeric/arithmetic.h \
  ../include/geometry/../geometry/../geometry/../numeric/interval.h \
  ../include/geometry/../geometry/../geometry/../geometry/ppl_polyhedron.h \
  ../include/geometry/../geometry/../geometry/../geometry/point.h \
  ../include/geometry/polyhedron.h test.h
	-${CXX} ${CXXFLAGS} -c test_simplex.cc 

test_polyhedron.o: test_polyhedron.cc ../include/ariadne.h \
  ../include/base/exception.h ../include/base/utility.h \
  ../include/base/../declarations.h ../include/numeric/float64.h \
  ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/function.h \
  ../include/numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/integer.h ../include/numeric/mpfloat.h \
  ../include/numeric/../numeric/rational.h \
  ../include/numeric/../numeric/../numeric/arithmetic.h \
  ../include/numeric/../numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../numeric/function.h \
  ../include/numeric/../numeric/../numeric/integer.h \
  ../include/numeric/rational.h ../include/geometry/point.h \
  ../include/geometry/../declarations.h \
  ../include/geometry/../base/array.h \
  ../include/geometry/../base/../declarations.h \
  ../include/geometry/../linear_algebra/vector.h \
  ../include/geometry/../linear_algebra/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/integer.h \
  ../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../include/geometry/../linear_algebra/matrix.h \
  ../include/geometry/point_list.h \
  ../include/geometry/../geometry/point.h \
  ../include/geometry/polyhedron.h \
  ../include/geometry/../linear_algebra/interval_vector.h \
  ../include/geometry/../linear_algebra/../numeric/interval.h \
  ../include/geometry/../linear_algebra/../numeric/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/arithmetic.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/rounding.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/../numeric/function.h \
  ../include/geometry/../linear_algebra/../linear_algebra/vector.h \
  ../include/geometry/../geometry/ppl_polyhedron.h \
  ../include/geometry/../geometry/../declarations.h \
  ../include/geometry/../geometry/../numeric/rational.h \
  ../include/geometry/../geometry/../linear_algebra/vector.h \
  ../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../include/geometry/../geometry/../linear_algebra/interval_vector.h \
  ../include/geometry/../geometry/point_list.h \
  ../include/geometry/../geometry/rectangle.h \
  ../include/geometry/../geometry/../numeric/arithmetic.h \
  ../include/geometry/../geometry/../numeric/interval.h \
  ../include/geometry/../geometry/../geometry/ppl_polyhedron.h \
  ../include/geometry/../geometry/../geometry/point.h \
  ../include/geometry/polyhedron.tpl ../include/geometry/polyhedron.h \
  ../include/geometry/../base/stlio.h \
  ../include/geometry/../base/../base/array.h \
  ../include/geometry/../numeric/interval.h \
  ../include/geometry/../numeric/arithmetic.h test.h
	-${CXX} ${CXXFLAGS} -c test_polyhedron.cc 

test_list_set.o: test_list_set.cc ../include/ariadne.h \
  ../include/real_typedef.h ../include/numeric/numerical_types.h \
  ../include/numeric/mpfloat.h ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/function.h \
  ../include/numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/rational.h \
  ../include/numeric/../numeric/../numeric/arithmetic.h \
  ../include/numeric/../numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../numeric/function.h \
  ../include/numeric/../numeric/../numeric/integer.h \
  ../include/base/exception.h ../include/base/utility.h \
  ../include/base/../declarations.h ../include/numeric/numerical_types.h \
  ../include/geometry/point.h ../include/geometry/../declarations.h \
  ../include/geometry/../base/array.h \
  ../include/geometry/../base/../declarations.h \
  ../include/geometry/../linear_algebra/vector.h \
  ../include/geometry/../linear_algebra/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/integer.h \
  ../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../include/geometry/../linear_algebra/matrix.h \
  ../include/geometry/rectangle.h \
  ../include/geometry/../numeric/arithmetic.h \
  ../include/geometry/../numeric/interval.h \
  ../include/geometry/../numeric/../declarations.h \
  ../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../numeric/../numeric/arithmetic.h \
  ../include/geometry/../numeric/../numeric/rounding.h \
  ../include/geometry/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../numeric/../numeric/../numeric/function.h \
  ../include/geometry/../linear_algebra/interval_vector.h \
  ../include/geometry/../linear_algebra/../numeric/interval.h \
  ../include/geometry/../linear_algebra/../linear_algebra/vector.h \
  ../include/geometry/../geometry/ppl_polyhedron.h \
  ../include/geometry/../geometry/../declarations.h \
  ../include/geometry/../geometry/../numeric/rational.h \
  ../include/geometry/../geometry/../linear_algebra/vector.h \
  ../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../include/geometry/../geometry/../linear_algebra/interval_vector.h \
  ../include/geometry/../geometry/point.h ../include/geometry/list_set.h \
  ../include/geometry/../base/utility.h test.h
	-${CXX} ${CXXFLAGS} -c test_list_set.cc 

test_grid_set.o: test_grid_set.cc ../include/ariadne.h \
  ../include/real_typedef.h ../include/numeric/numerical_types.h \
  ../include/numeric/mpfloat.h ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/function.h \
  ../include/numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/rational.h \
  ../include/numeric/../numeric/../numeric/arithmetic.h \
  ../include/numeric/../numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../numeric/function.h \
  ../include/numeric/../numeric/../numeric/integer.h \
  ../include/base/exception.h ../include/base/utility.h \
  ../include/base/../declarations.h ../include/base/binary_word.h \
  ../include/numeric/numerical_types.h ../include/geometry/point.h \
  ../include/geometry/../declarations.h \
  ../include/geometry/../base/array.h \
  ../include/geometry/../base/../declarations.h \
  ../include/geometry/../linear_algebra/vector.h \
  ../include/geometry/../linear_algebra/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/integer.h \
  ../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../include/geometry/../linear_algebra/matrix.h \
  ../include/geometry/rectangle.h \
  ../include/geometry/../numeric/arithmetic.h \
  ../include/geometry/../numeric/interval.h \
  ../include/geometry/../numeric/../declarations.h \
  ../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../numeric/../numeric/arithmetic.h \
  ../include/geometry/../numeric/../numeric/rounding.h \
  ../include/geometry/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../numeric/../numeric/../numeric/function.h \
  ../include/geometry/../linear_algebra/interval_vector.h \
  ../include/geometry/../linear_algebra/../numeric/interval.h \
  ../include/geometry/../linear_algebra/../linear_algebra/vector.h \
  ../include/geometry/../geometry/ppl_polyhedron.h \
  ../include/geometry/../geometry/../declarations.h \
  ../include/geometry/../geometry/../numeric/rational.h \
  ../include/geometry/../geometry/../linear_algebra/vector.h \
  ../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../include/geometry/../geometry/../linear_algebra/interval_vector.h \
  ../include/geometry/../geometry/point.h ../include/geometry/grid_set.h \
  ../include/geometry/../geometry/rectangle.h \
  ../include/geometry/../geometry/grid.h \
  ../include/geometry/../geometry/../base/array.h \
  ../include/geometry/../geometry/../numeric/arithmetic.h \
  ../include/geometry/../geometry/../geometry/lattice_set.h \
  ../include/geometry/../geometry/../geometry/../base/array.h \
  ../include/geometry/../geometry/../geometry/../base/iterator.h \
  ../include/geometry/../geometry/../geometry/../base/../base/array.h \
  ../include/geometry/../geometry/../geometry/../numeric/interval.h \
  ../include/geometry/../geometry/../geometry/../base/array_operations.h \
  ../include/geometry/../geometry/../geometry/../base/../declarations.h \
  ../include/geometry/../geometry/lattice_set.h \
  ../include/geometry/list_set.h ../include/geometry/../base/utility.h \
  test.h
	-${CXX} ${CXXFLAGS} -c test_grid_set.cc 

test_partition_tree_set.o: test_partition_tree_set.cc \
  ../include/ariadne.h ../include/real_typedef.h \
  ../include/numeric/numerical_types.h ../include/numeric/mpfloat.h \
  ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/function.h \
  ../include/numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/rational.h \
  ../include/numeric/../numeric/../numeric/arithmetic.h \
  ../include/numeric/../numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../numeric/function.h \
  ../include/numeric/../numeric/../numeric/integer.h \
  ../include/base/exception.h ../include/base/stlio.h \
  ../include/base/../base/array.h \
  ../include/base/../base/../declarations.h ../include/base/utility.h \
  ../include/base/../declarations.h ../include/base/binary_word.h \
  ../include/numeric/numerical_types.h ../include/linear_algebra/matrix.h \
  ../include/linear_algebra/../declarations.h \
  ../include/linear_algebra/../numeric/integer.h \
  ../include/linear_algebra/vector.h \
  ../include/linear_algebra/../numeric/arithmetic.h \
  ../include/geometry/point.h ../include/geometry/../declarations.h \
  ../include/geometry/../base/array.h \
  ../include/geometry/../linear_algebra/vector.h \
  ../include/geometry/../linear_algebra/matrix.h \
  ../include/geometry/rectangle.h \
  ../include/geometry/../numeric/arithmetic.h \
  ../include/geometry/../numeric/interval.h \
  ../include/geometry/../numeric/../declarations.h \
  ../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../numeric/../numeric/arithmetic.h \
  ../include/geometry/../numeric/../numeric/rounding.h \
  ../include/geometry/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../numeric/../numeric/../numeric/function.h \
  ../include/geometry/../linear_algebra/interval_vector.h \
  ../include/geometry/../linear_algebra/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/interval.h \
  ../include/geometry/../linear_algebra/../linear_algebra/vector.h \
  ../include/geometry/../geometry/ppl_polyhedron.h \
  ../include/geometry/../geometry/../declarations.h \
  ../include/geometry/../geometry/../numeric/rational.h \
  ../include/geometry/../geometry/../linear_algebra/vector.h \
  ../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../include/geometry/../geometry/../linear_algebra/interval_vector.h \
  ../include/geometry/../geometry/point.h \
  ../include/geometry/parallelotope.h \
  ../include/geometry/../geometry/zonotope.h \
  ../include/geometry/../geometry/../numeric/interval.h \
  ../include/geometry/../geometry/../geometry/ppl_polyhedron.h \
  ../include/geometry/../geometry/../geometry/point.h \
  ../include/geometry/list_set.h ../include/geometry/../base/utility.h \
  ../include/geometry/partition_tree_set.h \
  ../include/geometry/../base/binary_word.h \
  ../include/geometry/../base/binary_tree.h \
  ../include/geometry/../base/../declarations.h \
  ../include/geometry/../base/../base/array.h \
  ../include/geometry/../base/../base/binary_word.h \
  ../include/geometry/../base/iterator.h \
  ../include/geometry/../geometry/subdivision_tree_set.h \
  ../include/geometry/../geometry/../base/array.h \
  ../include/geometry/../geometry/../base/sequence.h \
  ../include/geometry/../geometry/../base/binary_word.h \
  ../include/geometry/../geometry/../base/binary_tree.h \
  ../include/geometry/../geometry/../base/iterator.h \
  ../include/geometry/../geometry/../geometry/rectangle.h \
  ../include/output/epsfstream.h \
  ../include/output/../numeric/approximation.h \
  ../include/output/../numeric/../declarations.h \
  ../include/output/../numeric/../numeric/numerical_traits.h \
  ../include/output/../geometry/point.h \
  ../include/output/../geometry/rectangle.h \
  ../include/output/../geometry/list_set.h \
  ../include/output/../geometry/grid_set.h \
  ../include/output/../geometry/../declarations.h \
  ../include/output/../geometry/../base/array.h \
  ../include/output/../geometry/../geometry/rectangle.h \
  ../include/output/../geometry/../geometry/grid.h \
  ../include/output/../geometry/../geometry/../declarations.h \
  ../include/output/../geometry/../geometry/../base/array.h \
  ../include/output/../geometry/../geometry/../numeric/arithmetic.h \
  ../include/output/../geometry/../geometry/../geometry/lattice_set.h \
  ../include/output/../geometry/../geometry/../geometry/../base/array.h \
  ../include/output/../geometry/../geometry/../geometry/../base/iterator.h \
  ../include/output/../geometry/../geometry/../geometry/../numeric/interval.h \
  ../include/output/../geometry/../geometry/../geometry/../base/array_operations.h \
  ../include/output/../geometry/../geometry/../geometry/../base/../declarations.h \
  ../include/output/../geometry/../geometry/lattice_set.h \
  ../include/output/../geometry/parallelotope.h \
  ../include/output/../geometry/zonotope.h \
  ../include/output/../geometry/polyhedron.h \
  ../include/output/../geometry/../linear_algebra/vector.h \
  ../include/output/../geometry/../linear_algebra/matrix.h \
  ../include/output/../geometry/../linear_algebra/interval_vector.h \
  ../include/output/../geometry/../geometry/ppl_polyhedron.h \
  ../include/output/../geometry/../geometry/point.h \
  ../include/output/../geometry/../geometry/point_list.h \
  ../include/output/../geometry/../geometry/../linear_algebra/vector.h \
  ../include/output/../geometry/../geometry/../linear_algebra/matrix.h \
  ../include/output/../geometry/../geometry/../geometry/point.h \
  ../include/output/../geometry/partition_tree_set.h \
  ../include/output/../linear_algebra/matrix.h \
  ../include/output/../system/affine_map.h \
  ../include/output/../system/../declarations.h \
  ../include/output/../system/../linear_algebra/vector.h \
  ../include/output/../system/../linear_algebra/matrix.h \
  ../include/output/../system/../system/map.h \
  ../include/output/../system/../system/../declarations.h test.h
	-${CXX} ${CXXFLAGS} -c test_partition_tree_set.cc 

test_epsfstream.o: test_epsfstream.cc ../include/ariadne.h \
  ../include/real_typedef.h ../include/numeric/numerical_types.h \
  ../include/numeric/mpfloat.h ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/function.h \
  ../include/numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/rational.h \
  ../include/numeric/../numeric/../numeric/arithmetic.h \
  ../include/numeric/../numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../numeric/function.h \
  ../include/numeric/../numeric/../numeric/integer.h \
  ../include/base/exception.h ../include/base/utility.h \
  ../include/base/../declarations.h ../include/numeric/numerical_types.h \
  ../include/geometry/point.h ../include/geometry/../declarations.h \
  ../include/geometry/../base/array.h \
  ../include/geometry/../base/../declarations.h \
  ../include/geometry/../linear_algebra/vector.h \
  ../include/geometry/../linear_algebra/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/integer.h \
  ../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../include/geometry/../linear_algebra/matrix.h \
  ../include/geometry/rectangle.h \
  ../include/geometry/../numeric/arithmetic.h \
  ../include/geometry/../numeric/interval.h \
  ../include/geometry/../numeric/../declarations.h \
  ../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../numeric/../numeric/arithmetic.h \
  ../include/geometry/../numeric/../numeric/rounding.h \
  ../include/geometry/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../numeric/../numeric/../numeric/function.h \
  ../include/geometry/../linear_algebra/interval_vector.h \
  ../include/geometry/../linear_algebra/../numeric/interval.h \
  ../include/geometry/../linear_algebra/../linear_algebra/vector.h \
  ../include/geometry/../geometry/ppl_polyhedron.h \
  ../include/geometry/../geometry/../declarations.h \
  ../include/geometry/../geometry/../numeric/rational.h \
  ../include/geometry/../geometry/../linear_algebra/vector.h \
  ../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../include/geometry/../geometry/../linear_algebra/interval_vector.h \
  ../include/geometry/../geometry/point.h ../include/geometry/list_set.h \
  ../include/geometry/../base/utility.h ../include/output/epsfstream.h \
  ../include/output/../numeric/approximation.h \
  ../include/output/../numeric/../declarations.h \
  ../include/output/../numeric/../numeric/numerical_traits.h \
  ../include/output/../geometry/point.h \
  ../include/output/../geometry/rectangle.h \
  ../include/output/../geometry/list_set.h \
  ../include/output/../geometry/grid_set.h \
  ../include/output/../geometry/../declarations.h \
  ../include/output/../geometry/../base/array.h \
  ../include/output/../geometry/../geometry/rectangle.h \
  ../include/output/../geometry/../geometry/grid.h \
  ../include/output/../geometry/../geometry/../declarations.h \
  ../include/output/../geometry/../geometry/../base/array.h \
  ../include/output/../geometry/../geometry/../numeric/arithmetic.h \
  ../include/output/../geometry/../geometry/../geometry/lattice_set.h \
  ../include/output/../geometry/../geometry/../geometry/../base/array.h \
  ../include/output/../geometry/../geometry/../geometry/../base/iterator.h \
  ../include/output/../geometry/../geometry/../geometry/../base/../base/array.h \
  ../include/output/../geometry/../geometry/../geometry/../numeric/interval.h \
  ../include/output/../geometry/../geometry/../geometry/../base/array_operations.h \
  ../include/output/../geometry/../geometry/../geometry/../base/../declarations.h \
  ../include/output/../geometry/../geometry/lattice_set.h \
  ../include/output/../geometry/parallelotope.h \
  ../include/output/../geometry/../numeric/interval.h \
  ../include/output/../geometry/../linear_algebra/vector.h \
  ../include/output/../geometry/../linear_algebra/matrix.h \
  ../include/output/../geometry/../geometry/point.h \
  ../include/output/../geometry/../geometry/zonotope.h \
  ../include/output/../geometry/../geometry/../linear_algebra/vector.h \
  ../include/output/../geometry/../geometry/../linear_algebra/matrix.h \
  ../include/output/../geometry/../geometry/../numeric/interval.h \
  ../include/output/../geometry/../geometry/../geometry/ppl_polyhedron.h \
  ../include/output/../geometry/../geometry/../geometry/point.h \
  ../include/output/../geometry/zonotope.h \
  ../include/output/../geometry/polyhedron.h \
  ../include/output/../geometry/../linear_algebra/interval_vector.h \
  ../include/output/../geometry/../geometry/ppl_polyhedron.h \
  ../include/output/../geometry/../geometry/point_list.h \
  ../include/output/../geometry/partition_tree_set.h \
  ../include/output/../geometry/../base/binary_word.h \
  ../include/output/../geometry/../base/../declarations.h \
  ../include/output/../geometry/../base/binary_tree.h \
  ../include/output/../geometry/../base/../base/array.h \
  ../include/output/../geometry/../base/../base/binary_word.h \
  ../include/output/../geometry/../base/iterator.h \
  ../include/output/../geometry/../geometry/subdivision_tree_set.h \
  ../include/output/../geometry/../geometry/../base/sequence.h \
  ../include/output/../geometry/../geometry/../base/binary_word.h \
  ../include/output/../geometry/../geometry/../base/binary_tree.h \
  ../include/output/../geometry/../geometry/../base/iterator.h \
  ../include/output/../geometry/../geometry/../geometry/rectangle.h \
  ../include/output/../linear_algebra/matrix.h \
  ../include/output/../system/affine_map.h \
  ../include/output/../system/../declarations.h \
  ../include/output/../system/../linear_algebra/vector.h \
  ../include/output/../system/../linear_algebra/matrix.h \
  ../include/output/../system/../system/map.h \
  ../include/output/../system/../system/../declarations.h test.h
	-${CXX} ${CXXFLAGS} -c test_epsfstream.cc 

test_interval_newton.o: test_interval_newton.cc \
  ../include/geometry/point.h ../include/geometry/../declarations.h \
  ../include/geometry/../base/array.h \
  ../include/geometry/../base/../declarations.h \
  ../include/geometry/../linear_algebra/vector.h \
  ../include/geometry/../linear_algebra/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/integer.h \
  ../include/geometry/../linear_algebra/../numeric/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../include/geometry/../linear_algebra/../numeric/../numeric/integer.h \
  ../include/geometry/../linear_algebra/matrix.h \
  ../include/geometry/rectangle.h \
  ../include/geometry/../numeric/arithmetic.h \
  ../include/geometry/../numeric/interval.h \
  ../include/geometry/../numeric/../declarations.h \
  ../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../numeric/../numeric/arithmetic.h \
  ../include/geometry/../numeric/../numeric/rounding.h \
  ../include/geometry/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../numeric/../numeric/../numeric/function.h \
  ../include/geometry/../numeric/../numeric/../numeric/integer.h \
  ../include/geometry/../linear_algebra/interval_vector.h \
  ../include/geometry/../linear_algebra/../numeric/interval.h \
  ../include/geometry/../linear_algebra/../linear_algebra/vector.h \
  ../include/geometry/../geometry/ppl_polyhedron.h \
  ../include/geometry/../geometry/../declarations.h \
  ../include/geometry/../geometry/../numeric/rational.h \
  ../include/geometry/../geometry/../numeric/../declarations.h \
  ../include/geometry/../geometry/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../geometry/../numeric/../numeric/arithmetic.h \
  ../include/geometry/../geometry/../numeric/../numeric/function.h \
  ../include/geometry/../geometry/../numeric/../numeric/integer.h \
  ../include/geometry/../geometry/../linear_algebra/vector.h \
  ../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../include/geometry/../geometry/../linear_algebra/interval_vector.h \
  ../include/geometry/../geometry/point.h ../include/system/henon_map.h \
  ../include/system/../linear_algebra/matrix.h \
  ../include/system/../linear_algebra/interval_matrix.h \
  ../include/system/../linear_algebra/../linear_algebra/matrix.h \
  ../include/system/../linear_algebra/../declarations.h \
  ../include/system/../linear_algebra/../numeric/interval.h \
  ../include/system/../geometry/point.h \
  ../include/system/../geometry/rectangle.h \
  ../include/system/../geometry/parallelotope.h \
  ../include/system/../geometry/../declarations.h \
  ../include/system/../geometry/../numeric/interval.h \
  ../include/system/../geometry/../linear_algebra/vector.h \
  ../include/system/../geometry/../linear_algebra/matrix.h \
  ../include/system/../geometry/../geometry/point.h \
  ../include/system/../geometry/../geometry/zonotope.h \
  ../include/system/../geometry/../geometry/../linear_algebra/vector.h \
  ../include/system/../geometry/../geometry/../linear_algebra/matrix.h \
  ../include/system/../geometry/../geometry/../numeric/interval.h \
  ../include/system/../geometry/../geometry/../geometry/ppl_polyhedron.h \
  ../include/system/../geometry/../geometry/../geometry/point.h \
  ../include/system/../system/map.h \
  ../include/system/../system/../declarations.h \
  ../include/system/../evaluation/apply.h \
  ../include/system/../evaluation/../declarations.h \
  ../include/evaluation/newton.h ../include/evaluation/../declarations.h \
  ../include/linear_algebra/interval_vector.h \
  ../include/linear_algebra/interval_matrix.h \
  ../include/geometry/rectangle.h ../include/system/map.h \
  ../include/system/vector_field.h ../include/system/../declarations.h \
  test.h ../include/real_typedef.h ../include/numeric/numerical_types.h \
  ../include/numeric/mpfloat.h ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/function.h \
  ../include/numeric/../numeric/rational.h
	-${CXX} ${CXXFLAGS} -c test_interval_newton.cc 

test_chainreach.o: test_chainreach.cc ../include/real_typedef.h \
  ../include/numeric/numerical_types.h ../include/numeric/mpfloat.h \
  ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/function.h \
  ../include/numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/rational.h \
  ../include/numeric/../numeric/../numeric/arithmetic.h \
  ../include/numeric/../numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../numeric/function.h \
  ../include/numeric/../numeric/../numeric/integer.h \
  ../include/geometry/point.h ../include/geometry/../declarations.h \
  ../include/geometry/../base/array.h \
  ../include/geometry/../base/../declarations.h \
  ../include/geometry/../linear_algebra/vector.h \
  ../include/geometry/../linear_algebra/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/integer.h \
  ../include/geometry/../linear_algebra/../numeric/arithmetic.h \
  ../include/geometry/../linear_algebra/matrix.h \
  ../include/geometry/rectangle.h \
  ../include/geometry/../numeric/arithmetic.h \
  ../include/geometry/../numeric/interval.h \
  ../include/geometry/../numeric/../declarations.h \
  ../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../numeric/../numeric/arithmetic.h \
  ../include/geometry/../numeric/../numeric/rounding.h \
  ../include/geometry/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../numeric/../numeric/../numeric/function.h \
  ../include/geometry/../linear_algebra/interval_vector.h \
  ../include/geometry/../linear_algebra/../numeric/interval.h \
  ../include/geometry/../linear_algebra/../linear_algebra/vector.h \
  ../include/geometry/../geometry/ppl_polyhedron.h \
  ../include/geometry/../geometry/../declarations.h \
  ../include/geometry/../geometry/../numeric/rational.h \
  ../include/geometry/../geometry/../linear_algebra/vector.h \
  ../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../include/geometry/../geometry/../linear_algebra/interval_vector.h \
  ../include/geometry/../geometry/point.h ../include/geometry/grid.h \
  ../include/geometry/../geometry/lattice_set.h \
  ../include/geometry/../geometry/../base/array.h \
  ../include/geometry/../geometry/../base/iterator.h \
  ../include/geometry/../geometry/../base/../base/array.h \
  ../include/geometry/../geometry/../numeric/interval.h \
  ../include/geometry/../geometry/../base/array_operations.h \
  ../include/geometry/../geometry/../base/../declarations.h \
  ../include/geometry/grid_set.h \
  ../include/geometry/../geometry/rectangle.h \
  ../include/geometry/../geometry/grid.h \
  ../include/geometry/partition_tree_set.h \
  ../include/geometry/../base/binary_word.h \
  ../include/geometry/../base/binary_tree.h \
  ../include/geometry/../base/../base/array.h \
  ../include/geometry/../base/../base/binary_word.h \
  ../include/geometry/../base/iterator.h \
  ../include/geometry/../geometry/subdivision_tree_set.h \
  ../include/geometry/../geometry/../base/sequence.h \
  ../include/geometry/../geometry/../base/binary_word.h \
  ../include/geometry/../geometry/../base/binary_tree.h \
  ../include/geometry/../geometry/../geometry/rectangle.h \
  ../include/system/henon_map.h \
  ../include/system/../linear_algebra/matrix.h \
  ../include/system/../linear_algebra/interval_matrix.h \
  ../include/system/../linear_algebra/../linear_algebra/matrix.h \
  ../include/system/../linear_algebra/../declarations.h \
  ../include/system/../linear_algebra/../numeric/interval.h \
  ../include/system/../geometry/point.h \
  ../include/system/../geometry/rectangle.h \
  ../include/system/../geometry/parallelotope.h \
  ../include/system/../geometry/../declarations.h \
  ../include/system/../geometry/../numeric/interval.h \
  ../include/system/../geometry/../linear_algebra/vector.h \
  ../include/system/../geometry/../linear_algebra/matrix.h \
  ../include/system/../geometry/../geometry/point.h \
  ../include/system/../geometry/../geometry/zonotope.h \
  ../include/system/../geometry/../geometry/../linear_algebra/vector.h \
  ../include/system/../geometry/../geometry/../linear_algebra/matrix.h \
  ../include/system/../geometry/../geometry/../numeric/interval.h \
  ../include/system/../geometry/../geometry/../geometry/ppl_polyhedron.h \
  ../include/system/../geometry/../geometry/../geometry/point.h \
  ../include/system/../system/map.h \
  ../include/system/../system/../declarations.h \
  ../include/system/../evaluation/apply.h \
  ../include/system/../evaluation/../declarations.h \
  ../include/evaluation/apply.h ../include/output/epsfstream.h \
  ../include/output/../numeric/approximation.h \
  ../include/output/../numeric/../declarations.h \
  ../include/output/../numeric/../numeric/numerical_traits.h \
  ../include/output/../geometry/point.h \
  ../include/output/../geometry/rectangle.h \
  ../include/output/../geometry/list_set.h \
  ../include/output/../geometry/../declarations.h \
  ../include/output/../geometry/../base/utility.h \
  ../include/output/../geometry/../base/../declarations.h \
  ../include/output/../geometry/grid_set.h \
  ../include/output/../geometry/parallelotope.h \
  ../include/output/../geometry/zonotope.h \
  ../include/output/../geometry/polyhedron.h \
  ../include/output/../geometry/../linear_algebra/vector.h \
  ../include/output/../geometry/../linear_algebra/matrix.h \
  ../include/output/../geometry/../linear_algebra/interval_vector.h \
  ../include/output/../geometry/../geometry/ppl_polyhedron.h \
  ../include/output/../geometry/../geometry/point.h \
  ../include/output/../geometry/../geometry/point_list.h \
  ../include/output/../geometry/../geometry/../declarations.h \
  ../include/output/../geometry/../geometry/../linear_algebra/vector.h \
  ../include/output/../geometry/../geometry/../linear_algebra/matrix.h \
  ../include/output/../geometry/../geometry/../geometry/point.h \
  ../include/output/../geometry/../geometry/rectangle.h \
  ../include/output/../geometry/partition_tree_set.h \
  ../include/output/../linear_algebra/matrix.h \
  ../include/output/../system/affine_map.h \
  ../include/output/../system/../declarations.h \
  ../include/output/../system/../linear_algebra/vector.h \
  ../include/output/../system/../linear_algebra/matrix.h \
  ../include/output/../system/../system/map.h test.h
	-${CXX} ${CXXFLAGS} -c test_chainreach.cc 

test_integration_step.o: test_integration_step.cc ../include/ariadne.h \
  ../include/real_typedef.h ../include/numeric/numerical_types.h \
  ../include/numeric/mpfloat.h ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/function.h \
  ../include/numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/rational.h \
  ../include/numeric/../numeric/../numeric/arithmetic.h \
  ../include/numeric/../numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/numeric/../numeric/../numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../numeric/function.h \
  ../include/numeric/../numeric/../numeric/integer.h \
  ../include/numeric/numerical_types.h ../include/linear_algebra/vector.h \
  ../include/linear_algebra/../declarations.h \
  ../include/linear_algebra/../numeric/integer.h \
  ../include/linear_algebra/../numeric/arithmetic.h \
  ../include/linear_algebra/matrix.h ../include/geometry/rectangle.h \
  ../include/geometry/../declarations.h \
  ../include/geometry/../numeric/arithmetic.h \
  ../include/geometry/../numeric/interval.h \
  ../include/geometry/../numeric/../declarations.h \
  ../include/geometry/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../numeric/../numeric/arithmetic.h \
  ../include/geometry/../numeric/../numeric/rounding.h \
  ../include/geometry/../numeric/../numeric/../numeric/numerical_traits.h \
  ../include/geometry/../numeric/../numeric/../numeric/function.h \
  ../include/geometry/../linear_algebra/interval_vector.h \
  ../include/geometry/../linear_algebra/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/interval.h \
  ../include/geometry/../linear_algebra/../linear_algebra/vector.h \
  ../include/geometry/../geometry/ppl_polyhedron.h \
  ../include/geometry/../geometry/../declarations.h \
  ../include/geometry/../geometry/../numeric/rational.h \
  ../include/geometry/../geometry/../linear_algebra/vector.h \
  ../include/geometry/../geometry/../linear_algebra/matrix.h \
  ../include/geometry/../geometry/../linear_algebra/interval_vector.h \
  ../include/geometry/../geometry/point.h \
  ../include/geometry/../geometry/../base/array.h \
  ../include/geometry/../geometry/../base/../declarations.h \
  ../include/geometry/parallelotope.h \
  ../include/geometry/../linear_algebra/vector.h \
  ../include/geometry/../linear_algebra/matrix.h \
  ../include/geometry/../geometry/zonotope.h \
  ../include/geometry/../geometry/../numeric/interval.h \
  ../include/geometry/../geometry/../geometry/ppl_polyhedron.h \
  ../include/geometry/../geometry/../geometry/point.h \
  ../include/system/affine_vector_field.h \
  ../include/system/../linear_algebra/vector.h \
  ../include/system/../linear_algebra/matrix.h \
  ../include/system/../system/vector_field.h \
  ../include/system/../system/../declarations.h \
  ../include/evaluation/lohner_integrator.h \
  ../include/evaluation/../declarations.h \
  ../include/evaluation/../evaluation/integrator.h \
  ../include/evaluation/../evaluation/../declarations.h test.h
	-${CXX} ${CXXFLAGS} -c test_integration_step.cc 

