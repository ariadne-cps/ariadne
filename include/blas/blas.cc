#include "blas.hpp"
#include "blas_output.hpp"
#include "fblas.hpp"

#include "amax.hpp"
#include "asum.hpp"
#include "axpy.hpp"
#include "copy.hpp"
#include "dot.hpp"
#include "gbmv.hpp"
#include "gemm.hpp"
#include "gemv.hpp"
#include "ger.hpp"
#include "iamax.hpp"
#include "nrm2.hpp"
#include "rotg.hpp"
#include "rot.hpp"
#include "rotmg.hpp"
#include "rotm.hpp"
#include "sbmv.hpp"
#include "scal.hpp"
#include "set.hpp"
#include "spmv.hpp"
#include "spr2.hpp"
#include "spr.hpp"
#include "swap.hpp"
#include "symm.hpp"
#include "symv.hpp"
#include "syr2.hpp"
#include "syr2k.hpp"
#include "syr.hpp"
#include "syrk.hpp"
#include "tbmv.hpp"
#include "tbsv.hpp"
#include "tpmv.hpp"
#include "tpsv.hpp"
#include "trmm.hpp"
#include "trmv.hpp"
#include "trsm.hpp"
#include "trsv.hpp"
#include "zdotc.hpp"
#include "zdotu.hpp"
#include "zgerc.hpp"
#include "zgeru.hpp"
#include "zhbmv.hpp"
#include "zhemm.hpp"
#include "zhemv.hpp"
#include "zher2.hpp"
#include "zher2k.hpp"
#include "zher.hpp"
#include "zherk.hpp"
#include "zhpmv.hpp"
#include "zhpr2.hpp"
#include "zhpr.hpp"

namespace BLAS {

template void set(const int N, const double alpha, 
                  double *X, const int incX);
template void set(const int N, const complex<double> alpha, 
                  complex<double> *X, const int incX);

template double amax(const int N, const double *X, const int incX);
template double amax(const int N, const complex<double> *X, const int incX);

/*
 * ===========================================================================
 * Instantiation of level 1 BLAS functions
 * ===========================================================================
 */

template double nrm2(const int N, const double *X, const int incX);
template double asum(const int N, const double *X, const int incX);
template int iamax(const int N, const double *X, const int incX);
template double dot(const int N, const double *X, const int incX,
                    const double *Y, const int incY);

template double nrm2(const int N, const complex<double> *X, const int incX);
template double asum(const int N, const complex<double> *X, const int incX);
template int iamax(const int N, const complex<double> *X, const int incX);
template complex<double> dotu(const int N, 
                              const complex<double> *X, const int incX,
                              const complex<double> *Y, const int incY);
template complex<double> dotc(const int N, 
                              const complex<double> *X, const int incX,
                              const complex<double> *Y, const int incY);


/*
 * ===========================================================================
 * Instantiation of level 1 BLAS routines
 * ===========================================================================
 */

template void swap(const int N, double *X, const int incX, 
                   double *Y, const int incY);
template void copy(const int N, const double *X, const int incX, 
                   double *Y, const int incY);
template void axpy(const int N, const double alpha, const double *X,
                   const int incX, double *Y, const int incY);
template void scal(const int N, const double alpha, double *X, const int incX);

template void rotg(double *a, double *b, double *c, double *s);
template void rotmg(double *d1, double *d2, double *b1, 
                    const double b2, double *P);
template void rot(const int N, double *X, const int incX,
                  double *Y, const int incY, const double c, const double  s);
template void rotm(const int N, double *X, const int incX,
                   double *Y, const int incY, const double *P);


template void swap(const int N, complex<double> *X, const int incX, 
                   complex<double> *Y, const int incY);
template void copy(const int N, const complex<double> *X, const int incX, 
                   complex<double> *Y, const int incY);
template void axpy(const int N, const complex<double> alpha, 
                   const complex<double> *X, const int incX, 
                   complex<double> *Y, const int incY);
template void scal(const int N, const complex<double> alpha, 
                   complex<double> *X, const int incX);


/*
 * ===========================================================================
 * Prototypes for level 2 BLAS
 * ===========================================================================
 */

template void gemv(const enum ORDER order,
                 const enum TRANSPOSE TransA, const int M, const int N,
                 const double alpha, const double *A, const int lda,
                 const double *X, const int incX, const double beta,
                 double *Y, const int incY);
template void gbmv(const enum ORDER order,
                 const enum TRANSPOSE TransA, const int M, const int N,
                 const int KL, const int KU, const double alpha,
                 const double *A, const int lda, const double *X,
                 const int incX, const double beta, double *Y, const int incY);
template void trmv(const enum ORDER order, const enum UPLO Uplo,
                 const enum TRANSPOSE TransA, const enum DIAG Diag,
                 const int N, const double *A, const int lda, 
                 double *X, const int incX);
template void tbmv(const enum ORDER order, const enum UPLO Uplo,
                 const enum TRANSPOSE TransA, const enum DIAG Diag,
                 const int N, const int K, const double *A, const int lda, 
                 double *X, const int incX);
template void tpmv(const enum ORDER order, const enum UPLO Uplo,
                 const enum TRANSPOSE TransA, const enum DIAG Diag,
                 const int N, const double *Ap, double *X, const int incX);
template void trsv(const enum ORDER order, const enum UPLO Uplo,
                 const enum TRANSPOSE TransA, const enum DIAG Diag,
                 const int N, const double *A, const int lda, double *X,
                 const int incX);
template void tbsv(const enum ORDER order, const enum UPLO Uplo,
                 const enum TRANSPOSE TransA, const enum DIAG Diag,
                 const int N, const int K, const double *A, const int lda,
                 double *X, const int incX);
template void tpsv(const enum ORDER order, const enum UPLO Uplo,
                 const enum TRANSPOSE TransA, const enum DIAG Diag,
                 const int N, const double *Ap, double *X, const int incX);

template void symv(const enum ORDER order, const enum UPLO Uplo,
                 const int N, const double alpha, const double *A,
                 const int lda, const double *X, const int incX,
                 const double beta, double *Y, const int incY);
template void sbmv(const enum ORDER order, const enum UPLO Uplo,
                 const int N, const int K, const double alpha, const double *A,
                 const int lda, const double *X, const int incX,
                 const double beta, double *Y, const int incY);
template void spmv(const enum ORDER order, const enum UPLO Uplo,
                 const int N, const double alpha, const double *Ap,
                 const double *X, const int incX,
                 const double beta, double *Y, const int incY);
template void ger(const enum ORDER order, const int M, const int N,
                const double alpha, const double *X, const int incX,
                const double *Y, const int incY, double *A, const int lda);
template void syr(const enum ORDER order, const enum UPLO Uplo,
                const int N, const double alpha, const double *X,
                const int incX, double *A, const int lda);
template void spr(const enum ORDER order, const enum UPLO Uplo,
                const int N, const double alpha, const double *X,
                const int incX, double *Ap);
template void syr2(const enum ORDER order, const enum UPLO Uplo,
                const int N, const double alpha, const double *X,
                const int incX, const double *Y, const int incY, double *A,
                const int lda);
template void spr2(const enum ORDER order, const enum UPLO Uplo,
                const int N, const double alpha, const double *X,
                const int incX, const double *Y, const int incY, double *A);



template void gemv(const enum ORDER order, const enum TRANSPOSE TransA, 
                   const int M, const int N, const complex<double> alpha, 
                   const complex<double> *A, const int lda,
                   const complex<double> *X, const int incX, 
                   const complex<double> beta,
                   complex<double> *Y, const int incY);
template void gbmv(const enum ORDER order, const enum TRANSPOSE TransA, 
                   const int M, const int N, const int KL, const int KU, 
                   const complex<double> alpha,
                   const complex<double> *A, const int lda, 
                   const complex<double> *X, const int incX, 
                   const complex<double> beta, 
                   complex<double> *Y, const int incY);
template void trmv(const enum ORDER order, const enum UPLO Uplo,
                   const enum TRANSPOSE TransA, const enum DIAG Diag,
                   const int N, const complex<double> *A, const int lda, 
                   complex<double> *X, const int incX);
template void tbmv(const enum ORDER order, const enum UPLO Uplo,
                   const enum TRANSPOSE TransA, const enum DIAG Diag,
                   const int N, const int K, 
                   const complex<double> *A, const int lda, 
                   complex<double> *X, const int incX);
template void tpmv(const enum ORDER order, const enum UPLO Uplo,
                   const enum TRANSPOSE TransA, const enum DIAG Diag,
                   const int N, const complex<double> *Ap,
                   complex<double> *X, const int incX);
template void trsv(const enum ORDER order, const enum UPLO Uplo,
                   const enum TRANSPOSE TransA, const enum DIAG Diag,
                   const int N, const complex<double> *A, const int lda, 
                   complex<double> *X, const int incX);
template void tbsv(const enum ORDER order, const enum UPLO Uplo,
                   const enum TRANSPOSE TransA, const enum DIAG Diag,
                   const int N, const int K, 
                   const complex<double> *A, const int lda,
                   complex<double> *X, const int incX);
template void tpsv(const enum ORDER order, const enum UPLO Uplo,
                   const enum TRANSPOSE TransA, const enum DIAG Diag,
                   const int N, const complex<double> *Ap, 
                   complex<double> *X, const int incX);


template void hemv(const enum ORDER order, const enum UPLO Uplo,
                 const int N, const complex<double> alpha, const complex<double> *A,
                 const int lda, const complex<double> *X, const int incX,
                 const complex<double> beta, complex<double> *Y, const int incY);
template void hbmv(const enum ORDER order, const enum UPLO Uplo,
                 const int N, const int K, const complex<double> alpha, const complex<double> *A,
                 const int lda, const complex<double> *X, const int incX,
                 const complex<double> beta, complex<double> *Y, const int incY);
template void hpmv(const enum ORDER order, const enum UPLO Uplo,
                 const int N, const complex<double> alpha, const complex<double> *Ap,
                 const complex<double> *X, const int incX,
                 const complex<double> beta, complex<double> *Y, const int incY);
template void geru(const enum ORDER order, const int M, const int N,
                 const complex<double> alpha, const complex<double> *X, const int incX,
                 const complex<double> *Y, const int incY, complex<double> *A, const int lda);
template void gerc(const enum ORDER order, const int M, const int N,
                 const complex<double> alpha, const complex<double> *X, const int incX,
                 const complex<double> *Y, const int incY, complex<double> *A, const int lda);
template void her(const enum ORDER order, const enum UPLO Uplo,
                const int N, const double alpha, const complex<double> *X, const int incX,
                complex<double> *A, const int lda);
template void hpr(const enum ORDER order, const enum UPLO Uplo,
                const int N, const double alpha, const complex<double> *X,
                const int incX, complex<double> *A);
template void her2(const enum ORDER order, const enum UPLO Uplo, const int N,
                const complex<double> alpha, const complex<double> *X, const int incX,
                const complex<double> *Y, const int incY, complex<double> *A, const int lda);
template void hpr2(const enum ORDER order, const enum UPLO Uplo, const int N,
                const complex<double> alpha, const complex<double> *X, const int incX,
                const complex<double> *Y, const int incY, complex<double> *Ap);

/*
 * ===========================================================================
 * Prototypes for level 3 BLAS
 * ===========================================================================
 */

/* 
 * Routines with standard 4 prefixes (S, D, C, Z)
 */
template void gemm(const enum ORDER Order, const enum TRANSPOSE TransA,
                 const enum TRANSPOSE TransB, const int M, const int N,
                 const int K, const double alpha, const double *A,
                 const int lda, const double *B, const int ldb,
                 const double beta, double *C, const int ldc);
template void symm(const enum ORDER Order, const enum SIDE Side,
                 const enum UPLO Uplo, const int M, const int N,
                 const double alpha, const double *A, const int lda,
                 const double *B, const int ldb, const double beta,
                 double *C, const int ldc);
template void syrk(const enum ORDER Order, const enum UPLO Uplo,
                 const enum TRANSPOSE Trans, const int N, const int K,
                 const double alpha, const double *A, const int lda,
                 const double beta, double *C, const int ldc);
template void syr2k(const enum ORDER Order, const enum UPLO Uplo,
                  const enum TRANSPOSE Trans, const int N, const int K,
                  const double alpha, const double *A, const int lda,
                  const double *B, const int ldb, const double beta,
                  double *C, const int ldc);
template void trmm(const enum ORDER Order, const enum SIDE Side,
                 const enum UPLO Uplo, const enum TRANSPOSE TransA,
                 const enum DIAG Diag, const int M, const int N,
                 const double alpha, const double *A, const int lda,
                 double *B, const int ldb);
template void trsm(const enum ORDER Order, const enum SIDE Side,
                 const enum UPLO Uplo, const enum TRANSPOSE TransA,
                 const enum DIAG Diag, const int M, const int N,
                 const double alpha, const double *A, const int lda,
                 double *B, const int ldb);



/* 
 * Routines with prefixes C and Z only
 */
template void hemm(const enum ORDER Order, const enum SIDE Side,
                 const enum UPLO Uplo, const int M, const int N,
                 const complex<double> alpha, const complex<double> *A, const int lda,
                 const complex<double> *B, const int ldb, const complex<double> beta,
                 complex<double> *C, const int ldc);
template void herk(const enum ORDER Order, const enum UPLO Uplo,
                 const enum TRANSPOSE Trans, const int N, const int K,
                 const double alpha, const complex<double> *A, const int lda,
                 const double beta, complex<double> *C, const int ldc);
template void her2k(const enum ORDER Order, const enum UPLO Uplo,
                  const enum TRANSPOSE Trans, const int N, const int K,
                  const complex<double> alpha, const complex<double> *A, const int lda,
                  const complex<double> *B, const int ldb, const double beta,
                  complex<double> *C, const int ldc);

}

