#pragma once
// C++ program to demostrate working of Guassian Elimination
// method
#include <bits/stdc++.h>
#include <eigen3/Eigen/Dense>

#include <glpk.h>

<<<<<<< HEAD
namespace Ariadne {
class IndefiniteMatrixException : public std::runtime_error {
public:
  IndefiniteMatrixException(const StringType &what)
      : std::runtime_error(what) {}
};

//  perform simplex algorithm to find minimum with only lower constriants
//    using glpk library
template <class X>
Vector<X> lp_min(const Vector<X> &C, const Matrix<X> &A, const Vector<X> &b,
                 const Vector<X> &lb, int &errnum);
=======
namespace Ariadne
{
  class IndefiniteMatrixException : public std::runtime_error
  {
  public:
    IndefiniteMatrixException(const StringType &what) :
      std::runtime_error(what){}
  };
>>>>>>> Small fixes. Implemented temporary __feasible__ function to test barrier method.

// Compute the orthogonal decomposition A=QR with or without column pivoting.
// The matrix Q is built up as a composition of elementary Householder
// transformations H=I-vv' with |v|=1. Note that inv(H)=H'=H. The vector v is
// chosen to be a multiple of the first working column of A.
template <class X>
Tuple<Matrix<X>, Matrix<X>, PivotMatrix>
orthogonal_decomposition(const Matrix<X> &A, Bool allow_pivoting);

// Compute the pseudo-inverse matrix operation
template <class X> Matrix<X> pinv(const Matrix<X> &G);

// compute the null basis of a matrix A
template <class X> Tuple<Matrix<X>, unsigned> null(const Matrix<X> &A);

// compute the Cholesky factorization
template <class X> Matrix<X> chol(const Matrix<X> &A);

// Compute the norm_2 on a vector
template <class X> X norm2(const Vector<X> &A);

// Compute the norm_1 on a vector
template <class X> X norm1(const Vector<X> &A);

// Jacobi computation of eigenvalues
template <class X> Tuple<Matrix<X>, Vector<X>> jacobi_eigs(const Matrix<X> &A);

// Power method computation of eigenvalues
template <class X>
Tuple<Vector<X>, X> power_eigs(const Matrix<X> &A, const unsigned itmax = 100);

// Inverse power method computation of eigenvalues:
//    compute power method with (A-muI)^-1
//    taken from https://en.wikiversity.org/wiki/Shifted_inverse_iteration
template <class X>
Tuple<Vector<X>, X> inverse_power_eigs(const Matrix<X> &A, const X &mu,
                                       const unsigned itmax = 100);

// Eigen computation of eigenvalues
template <class X> Tuple<Vector<X>, X> eigen_eigs(const Matrix<X> &A);

// Eigen computation of null space!
// @Warning if kernel has dimension 0, return a vector of only 0
template <class X> Tuple<Matrix<X>, unsigned> eigen_null(const Matrix<X> &G);

// Compute the pseudo-inverse matrix operation
template <class X> Matrix<X> eigen_pinv(const Matrix<X> &G);

<<<<<<< HEAD
// compute the Cholesky factorization using Eigen library
template <class X> Matrix<X> eigen_chol(const Matrix<X> &A);

} // namespace Ariadne
=======
  // Eigen computation of null space!
  // @Warning if kernel has dimension 0, return a vector of only 0
  template<class X> Tuple<Matrix<X>, unsigned>
  eigen_null(const Matrix<X> &G);

  // Compute the pseudo-inverse matrix operation
  template<class X> Matrix<X>
  eigen_pinv(const Matrix<X> &G);

  // compute the Cholesky factorization using Eigen library
  template<class X> Matrix<X>
  eigen_chol(const Matrix<X> &A);

} //namespace Ariadne
>>>>>>> Small fixes. Implemented temporary __feasible__ function to test barrier method.

#include "impl/AlgebraAddOn.i.hpp"
