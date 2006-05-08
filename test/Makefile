include ../config.mk

TESTS = test_array test_arithmetic test_binary_word \
  test_matrix \
  test_polyhedron test_rectangle test_point \
  test_linear_algebra test_interval test_list_set test_grid \
  test_lattice
#  test_partition_tree \
#  test_epsfstream \

CXXFLAGS += -Wall -I../include/ -I../
LDXXFLAGS += -L../src/ -lariadne -lgmp -lgmpxx -lppl -lcblas -lblas

all: ${TESTS}

clean:
	-rm -rf core ${TESTS}

dep:  
	@sed -e '/^### Do Not edit this line$$/,$$d' Makefile \
		> MMM.$$$$ && mv MMM.$$$$ Makefile
	@echo "### Do Not edit this line" >> Makefile
	@echo "### Everything below is auto-generated" >> Makefile
	@echo "" >> Makefile 
	@echo "all: $$""{TESTS}" >> Makefile
	@for ff in ${TESTS}; do echo "	-./$$ff"; done >> Makefile
	@echo "" >> Makefile 
	@for ff in ${TESTS}; do g++ -MM -MT$$ff ${CXXFLAGS} $$ff.cc; \
	echo "	-$$""{CXX} $$""{CXXFLAGS} $$ff.cc  $$""{LDXXFLAGS} -o $$ff"; \
	echo ""; done >> Makefile


### Do Not edit this line
### Everything below is auto-generated

all: ${TESTS}
	-./test_array
	-./test_arithmetic
	-./test_binary_word
	-./test_matrix
	-./test_polyhedron
	-./test_rectangle
	-./test_point
	-./test_linear_algebra
	-./test_interval
	-./test_list_set
	-./test_grid
	-./test_lattice

test_array: test_array.cc ../include/ariadne.h ../include/base/array.h \
  ../include/base/../declarations.h ../include/numeric/numerical_types.h \
  ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../utility/stlio.h \
  ../include/numeric/../numeric/../utility/../base/array.h \
  ../include/numeric/../numeric/dyadic.h \
  ../include/numeric/../utility/stlio.h test.h
	-${CXX} ${CXXFLAGS} test_array.cc  ${LDXXFLAGS} -o test_array

test_arithmetic: test_arithmetic.cc ../include/numeric/arithmetic.h \
  ../include/numeric/../numeric/numerical_types.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../numeric/../utility/stlio.h \
  ../include/numeric/../numeric/../numeric/../utility/../base/array.h \
  ../include/numeric/../numeric/../numeric/../utility/../base/../declarations.h \
  ../include/numeric/../numeric/../numeric/dyadic.h \
  ../include/numeric/../numeric/../utility/stlio.h \
  ../include/numeric/../numeric/approximation.h \
  ../include/numeric/../numeric/../numeric/numerical_types.h test.h
	-${CXX} ${CXXFLAGS} test_arithmetic.cc  ${LDXXFLAGS} -o test_arithmetic

test_binary_word: test_binary_word.cc ../include/ariadne.h \
  ../include/base/exception.h ../include/base/binary_word.h \
  ../include/base/../declarations.h ../include/base/binary_tree.h \
  ../include/base/../base/array.h \
  ../include/base/../base/../declarations.h \
  ../include/base/../base/binary_word.h \
  ../include/numeric/numerical_types.h \
  ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../utility/stlio.h \
  ../include/numeric/../numeric/../utility/../base/array.h \
  ../include/numeric/../numeric/dyadic.h \
  ../include/numeric/../utility/stlio.h ../include/base/utility.h \
  ../include/geometry/point.h ../include/geometry/../declarations.h \
  ../include/geometry/../base/array.h \
  ../include/geometry/../linear_algebra/vector.h \
  ../include/geometry/../linear_algebra/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/integer.h \
  ../include/geometry/rectangle.h \
  ../include/geometry/../numeric/arithmetic.h \
  ../include/geometry/../numeric/../numeric/numerical_types.h \
  ../include/geometry/../numeric/../numeric/approximation.h \
  ../include/geometry/../numeric/../numeric/../declarations.h \
  ../include/geometry/../numeric/../numeric/../numeric/numerical_types.h \
  ../include/geometry/../numeric/interval.h \
  ../include/geometry/../numeric/../declarations.h \
  ../include/geometry/../geometry/point.h test.h
	-${CXX} ${CXXFLAGS} test_binary_word.cc  ${LDXXFLAGS} -o test_binary_word

test_matrix: test_matrix.cc ../include/blas/blas.hpp \
  ../include/blas/fblas.hpp ../include/blas/blas.hpp \
  /ufs/collins/include/cblas.h ../include/blas/xerbla.hpp \
  ../include/blas/blas_output.hpp ../include/lapack/lapack.hpp \
  ../include/numeric/numerical_types.h \
  ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../utility/stlio.h \
  ../include/numeric/../numeric/../utility/../base/array.h \
  ../include/numeric/../numeric/../utility/../base/../declarations.h \
  ../include/numeric/../numeric/dyadic.h \
  ../include/numeric/../utility/stlio.h \
  ../include/linear_algebra/vector.h \
  ../include/linear_algebra/../declarations.h \
  ../include/linear_algebra/../numeric/integer.h \
  ../include/linear_algebra/vector.tpl ../include/linear_algebra/vector.h \
  ../include/linear_algebra/../utility/stlio.h \
  ../include/linear_algebra/matrix.h ../include/linear_algebra/matrix.tpl \
  ../include/linear_algebra/matrix.h \
  ../include/linear_algebra/../base/array.h \
  ../include/linear_algebra/../numeric/arithmetic.h \
  ../include/linear_algebra/../numeric/../numeric/numerical_types.h \
  ../include/linear_algebra/../numeric/../numeric/approximation.h \
  ../include/linear_algebra/../numeric/../numeric/../declarations.h \
  ../include/linear_algebra/../numeric/../numeric/../numeric/numerical_types.h \
  ../include/linear_algebra/../numeric/function.h \
  /ufs/collins/include/mpfr.h \
  ../include/linear_algebra/../linear_algebra/vector.h \
  ../include/linear_algebra/../linear_algebra/matrix.h \
  ../include/linear_algebra/../linear_algebra/lu_matrix.h \
  ../include/blas/geset.hpp ../include/blas/trcpy.hpp \
  ../include/blas/gemm.hpp ../include/blas/trmm.hpp \
  ../include/lapack/laswp.hpp ../include/lapack/lapack.hpp \
  ../include/lapack/getrf.hpp ../include/blas/iamax.hpp \
  ../include/blas/scal.hpp ../include/blas/swap.hpp \
  ../include/blas/ger.hpp ../include/lapack/getrs.hpp \
  ../include/linear_algebra/../linear_algebra/../declarations.h \
  ../include/linear_algebra/../linear_algebra/../base/array.h \
  ../include/linear_algebra/../linear_algebra/../linear_algebra/matrix.h \
  ../include/linear_algebra/lu_matrix.h \
  ../include/linear_algebra/qr_matrix.h ../include/blas/gecpy.hpp \
  ../include/blas/trset.hpp ../include/lapack/larf.hpp \
  ../include/lapack/larfg.hpp ../include/blas/nrm2.hpp \
  ../include/lapack/geqrf.hpp ../include/lapack/larf.hpp \
  ../include/lapack/larfg.hpp ../include/lapack/orgqr.hpp \
  ../include/linear_algebra/svd_matrix.h ../include/blas/copy.hpp \
  ../include/lapack/gesvd.hpp
	-${CXX} ${CXXFLAGS} test_matrix.cc  ${LDXXFLAGS} -o test_matrix

test_polyhedron: test_polyhedron.cc ../include/ariadne.h \
  ../include/base/exception.h ../include/base/utility.h \
  ../include/base/../declarations.h ../include/numeric/numerical_types.h \
  ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../utility/stlio.h \
  ../include/numeric/../numeric/../utility/../base/array.h \
  ../include/numeric/../numeric/../utility/../base/../declarations.h \
  ../include/numeric/../numeric/dyadic.h \
  ../include/numeric/../utility/stlio.h ../include/geometry/point.h \
  ../include/geometry/../declarations.h \
  ../include/geometry/../base/array.h \
  ../include/geometry/../linear_algebra/vector.h \
  ../include/geometry/../linear_algebra/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/integer.h \
  ../include/geometry/polyhedron.h /ufs/collins/include/ppl.hh \
  ../include/geometry/../linear_algebra/constraint_system.h \
  ../include/geometry/../linear_algebra/../numeric/numerical_types.h \
  ../include/geometry/../linear_algebra/../linear_algebra/vector.h \
  ../include/geometry/../linear_algebra/../linear_algebra/matrix.h \
  ../include/geometry/../linear_algebra/../linear_algebra/../declarations.h \
  ../include/geometry/../linear_algebra/../linear_algebra/../numeric/integer.h \
  ../include/geometry/../linear_algebra/matrix.h \
  ../include/geometry/../linear_algebra/generator_system.h \
  ../include/geometry/../geometry/point.h \
  ../include/geometry/../geometry/rectangle.h \
  ../include/geometry/../geometry/../declarations.h \
  ../include/geometry/../geometry/../numeric/arithmetic.h \
  ../include/geometry/../geometry/../numeric/../numeric/numerical_types.h \
  ../include/geometry/../geometry/../numeric/../numeric/approximation.h \
  ../include/geometry/../geometry/../numeric/../numeric/../declarations.h \
  ../include/geometry/../geometry/../numeric/../numeric/../numeric/numerical_types.h \
  ../include/geometry/../geometry/../numeric/interval.h \
  ../include/geometry/../geometry/../numeric/../declarations.h \
  ../include/geometry/../geometry/../geometry/point.h
	-${CXX} ${CXXFLAGS} test_polyhedron.cc  ${LDXXFLAGS} -o test_polyhedron

test_rectangle: test_rectangle.cc ../include/ariadne.h \
  ../include/base/exception.h ../include/base/utility.h \
  ../include/base/../declarations.h ../include/geometry/rectangle.h \
  ../include/geometry/../declarations.h \
  ../include/geometry/../numeric/arithmetic.h \
  ../include/geometry/../numeric/../numeric/numerical_types.h \
  ../include/geometry/../numeric/../numeric/../declarations.h \
  ../include/geometry/../numeric/../numeric/../numeric/integer.h \
  ../include/geometry/../numeric/../numeric/../numeric/../declarations.h \
  ../include/geometry/../numeric/../numeric/../numeric/../utility/stlio.h \
  ../include/geometry/../numeric/../numeric/../numeric/../utility/../base/array.h \
  ../include/geometry/../numeric/../numeric/../numeric/../utility/../base/../declarations.h \
  ../include/geometry/../numeric/../numeric/../numeric/dyadic.h \
  ../include/geometry/../numeric/../numeric/../utility/stlio.h \
  ../include/geometry/../numeric/../numeric/approximation.h \
  ../include/geometry/../numeric/../numeric/../numeric/numerical_types.h \
  ../include/geometry/../numeric/interval.h \
  ../include/geometry/../numeric/../declarations.h \
  ../include/geometry/../geometry/point.h \
  ../include/geometry/../geometry/../declarations.h \
  ../include/geometry/../geometry/../base/array.h \
  ../include/geometry/../geometry/../linear_algebra/vector.h \
  ../include/geometry/../geometry/../linear_algebra/../declarations.h \
  ../include/geometry/../geometry/../linear_algebra/../numeric/integer.h \
  ../include/geometry/list_set.h ../include/geometry/../base/utility.h \
  test.h
	-${CXX} ${CXXFLAGS} test_rectangle.cc  ${LDXXFLAGS} -o test_rectangle

test_point: test_point.cc ../include/geometry/point.h \
  ../include/geometry/../declarations.h \
  ../include/geometry/../base/array.h \
  ../include/geometry/../base/../declarations.h \
  ../include/geometry/../linear_algebra/vector.h \
  ../include/geometry/../linear_algebra/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/integer.h \
  ../include/geometry/../linear_algebra/../numeric/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/../utility/stlio.h \
  ../include/geometry/../linear_algebra/../numeric/../utility/../base/array.h \
  test.h
	-${CXX} ${CXXFLAGS} test_point.cc  ${LDXXFLAGS} -o test_point

test_linear_algebra: test_linear_algebra.cc \
  ../include/numeric/numerical_types.h \
  ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../utility/stlio.h \
  ../include/numeric/../numeric/../utility/../base/array.h \
  ../include/numeric/../numeric/../utility/../base/../declarations.h \
  ../include/numeric/../numeric/dyadic.h \
  ../include/numeric/../utility/stlio.h \
  ../include/linear_algebra/vector.h \
  ../include/linear_algebra/../declarations.h \
  ../include/linear_algebra/../numeric/integer.h \
  ../include/linear_algebra/matrix.h test.h
	-${CXX} ${CXXFLAGS} test_linear_algebra.cc  ${LDXXFLAGS} -o test_linear_algebra

test_interval: test_interval.cc ../include/ariadne.h \
  ../include/numeric/numerical_types.h \
  ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../utility/stlio.h \
  ../include/numeric/../numeric/../utility/../base/array.h \
  ../include/numeric/../numeric/../utility/../base/../declarations.h \
  ../include/numeric/../numeric/dyadic.h \
  ../include/numeric/../utility/stlio.h ../include/numeric/interval.h \
  test.h
	-${CXX} ${CXXFLAGS} test_interval.cc  ${LDXXFLAGS} -o test_interval

test_list_set: test_list_set.cc ../include/ariadne.h \
  ../include/base/exception.h ../include/base/utility.h \
  ../include/base/../declarations.h ../include/numeric/numerical_types.h \
  ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../utility/stlio.h \
  ../include/numeric/../numeric/../utility/../base/array.h \
  ../include/numeric/../numeric/../utility/../base/../declarations.h \
  ../include/numeric/../numeric/dyadic.h \
  ../include/numeric/../utility/stlio.h ../include/geometry/point.h \
  ../include/geometry/../declarations.h \
  ../include/geometry/../base/array.h \
  ../include/geometry/../linear_algebra/vector.h \
  ../include/geometry/../linear_algebra/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/integer.h \
  ../include/geometry/rectangle.h \
  ../include/geometry/../numeric/arithmetic.h \
  ../include/geometry/../numeric/../numeric/numerical_types.h \
  ../include/geometry/../numeric/../numeric/approximation.h \
  ../include/geometry/../numeric/../numeric/../declarations.h \
  ../include/geometry/../numeric/../numeric/../numeric/numerical_types.h \
  ../include/geometry/../numeric/interval.h \
  ../include/geometry/../numeric/../declarations.h \
  ../include/geometry/../geometry/point.h ../include/geometry/list_set.h \
  ../include/geometry/../base/utility.h test.h
	-${CXX} ${CXXFLAGS} test_list_set.cc  ${LDXXFLAGS} -o test_list_set

test_grid: test_grid.cc ../include/ariadne.h ../include/base/exception.h \
  ../include/base/utility.h ../include/base/../declarations.h \
  ../include/numeric/numerical_types.h \
  ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../utility/stlio.h \
  ../include/numeric/../numeric/../utility/../base/array.h \
  ../include/numeric/../numeric/../utility/../base/../declarations.h \
  ../include/numeric/../numeric/dyadic.h \
  ../include/numeric/../utility/stlio.h ../include/geometry/point.h \
  ../include/geometry/../declarations.h \
  ../include/geometry/../base/array.h \
  ../include/geometry/../linear_algebra/vector.h \
  ../include/geometry/../linear_algebra/../declarations.h \
  ../include/geometry/../linear_algebra/../numeric/integer.h \
  ../include/geometry/rectangle.h \
  ../include/geometry/../numeric/arithmetic.h \
  ../include/geometry/../numeric/../numeric/numerical_types.h \
  ../include/geometry/../numeric/../numeric/approximation.h \
  ../include/geometry/../numeric/../numeric/../declarations.h \
  ../include/geometry/../numeric/../numeric/../numeric/numerical_types.h \
  ../include/geometry/../numeric/interval.h \
  ../include/geometry/../numeric/../declarations.h \
  ../include/geometry/../geometry/point.h ../include/base/binary_word.h \
  ../include/geometry/grid_set.h \
  ../include/geometry/../geometry/rectangle.h \
  ../include/geometry/../geometry/grid.h \
  ../include/geometry/../geometry/../declarations.h \
  ../include/geometry/../geometry/../base/array.h \
  ../include/geometry/../geometry/../numeric/arithmetic.h \
  ../include/geometry/../geometry/../geometry/lattice_set.h \
  ../include/geometry/../geometry/../geometry/../base/array.h \
  ../include/geometry/../geometry/../geometry/../base/iterator.h \
  ../include/geometry/../geometry/../geometry/../base/../base/array.h \
  ../include/geometry/../geometry/../geometry/../numeric/interval.h \
  ../include/geometry/../geometry/../geometry/../base/array_operations.h \
  ../include/geometry/../geometry/../geometry/../base/../declarations.h \
  ../include/geometry/../geometry/lattice_set.h test.h
	-${CXX} ${CXXFLAGS} test_grid.cc  ${LDXXFLAGS} -o test_grid

test_lattice: test_lattice.cc ../include/declarations.h \
  ../include/real_typedef.h ../include/numeric/numerical_types.h \
  ../include/numeric/../declarations.h \
  ../include/numeric/../numeric/integer.h \
  ../include/numeric/../numeric/../declarations.h \
  ../include/numeric/../numeric/../utility/stlio.h \
  ../include/numeric/../numeric/../utility/../base/array.h \
  ../include/numeric/../numeric/../utility/../base/../declarations.h \
  ../include/numeric/../numeric/dyadic.h \
  ../include/numeric/../utility/stlio.h ../include/geometry/lattice_set.h \
  ../include/geometry/../base/array.h \
  ../include/geometry/../base/iterator.h \
  ../include/geometry/../base/../base/array.h \
  ../include/geometry/../numeric/interval.h \
  ../include/geometry/../numeric/../declarations.h \
  ../include/geometry/../base/array_operations.h \
  ../include/geometry/../base/../declarations.h \
  ../include/evaluation/lattice_map.h \
  ../include/evaluation/../geometry/lattice_set.h
	-${CXX} ${CXXFLAGS} test_lattice.cc  ${LDXXFLAGS} -o test_lattice

