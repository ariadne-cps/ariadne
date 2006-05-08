#include <blas/blas_output.hpp>

#include "lapack.hpp"
#include "gesvd.hpp"

using std::cout;
using std::cerr;
using std::endl;

using BLAS::matrix;
using BLAS::vector;

int main() {
  int m=3;
  int n=4;
  
  const double A[12]={2.0,2.0,1.0,4.0,  2.0,1.0,1.0,0.0,  1.0,0.5,4.0,1.0};
  //const double A[4]={4.0,1.0,  3.0,2.0};
  //const double A[6]={4.0,1.0,-1.0,  3.0,2.0,1.0};
  double S[BLAS::max(m,n)];
  double U[m*m];
  double V[n*n];
  const double I[9]={1.0,0.0,0.0,  0.0,1.0,0.0, 0.0,0.0,1.0};

  double D[m*n];
  double B[m*n];
  BLAS::copy(m*n,A,1,B,1);
  
  BLAS::set(BLAS::max(m,n),0.0,S,1);
  BLAS::set(m*m,0.0,U,1);
  BLAS::set(n*n,0.0,V,1);
  cout << matrix(m,n,A) << endl;
  LAPACK::gesvd(BLAS::RowMajor,m,n,B,n,S,U,m,V,n);
  cout << "B=\n" << matrix(m,n,B) << endl << endl;
  cout << "S=\n" << vector(m,S) << endl << endl;
  cout << "U=\n" << matrix(m,m,U) << endl;
  cout << "V=\n" << matrix(n,n,V) << endl;
  BLAS::set(m*n,0.0,D,1);
  BLAS::copy(m,S,1,D,n+1);
  cout << "D=\n" << matrix(m,n,D) << endl;
  double T[m*n];
  double C[m*n];
  gemm(BLAS::RowMajor,BLAS::NoTrans,BLAS::NoTrans,m,n,m,1.0,U,m,D,n,0.0,T,n);
  gemm(BLAS::RowMajor,BLAS::NoTrans,BLAS::Trans,m,n,n,1.0,T,n,V,n,0.0,C,n);
  cout << "C=\n" << matrix(m,n,C) << endl;
  cout << "A=\n" << matrix(m,n,A) << endl;
  double E[m*n];
  BLAS::copy(m*n,A,1,E,1);
  BLAS::axpy(m*n,-1.0,C,1,E,1);
  cout << "Error=" << BLAS::amax(m*n,E,1) << endl;

  return 0;
}
