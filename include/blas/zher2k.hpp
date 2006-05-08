#ifndef __BLAS_ZHER2K_HPP__
#define __BLAS_ZHER2K_HPP__

#include "blas.hpp"

template<typename real>
void
BLAS::her2k (const enum ORDER Order, const enum UPLO Uplo,
              const enum TRANSPOSE Trans, const int N, const int K,
              const complex<real> alpha, const complex<real> *A, const int lda, const complex<real> *B,
              const int ldb, const real beta, complex<real> *C, const int ldc)
{
{
  int i, j, k;
  int uplo, trans;
  real alpha_real = alpha.real(); 
  real alpha_imag = alpha.imag(); 
  if (beta == 1.0 && ((alpha_real == 0.0 && alpha_imag == 0.0) || K == 0))
    return;
  if (Order == RowMajor) {
    uplo = Uplo;
    trans = Trans;
  } else {
    uplo = (Uplo == Upper) ? Lower : Upper;
    trans = (Trans == NoTrans) ? ConjTrans : NoTrans;
    alpha_imag *= -1;
  }
  if (beta == 0.0) {
    if (uplo == Upper) {
      for (i = 0; i < N; i++) {
        for (j = i; j < N; j++) {
          (((real *) C)[2*(ldc * i + j)]) = 0.0;
          (((real *) C)[2*(ldc * i + j)+1]) = 0.0;
        }
      }
    } else {
      for (i = 0; i < N; i++) {
        for (j = 0; j <= i; j++) {
          (((real *) C)[2*(ldc * i + j)]) = 0.0;
          (((real *) C)[2*(ldc * i + j)+1]) = 0.0;
        }
      }
    }
  } else if (beta != 1.0) {
    if (uplo == Upper) {
      for (i = 0; i < N; i++) {
        (((real *) C)[2*(ldc * i + i)]) *= beta;
        (((real *) C)[2*(ldc * i + i)+1]) = 0.0;
        for (j = i + 1; j < N; j++) {
          (((real *) C)[2*(ldc * i + j)]) *= beta;
          (((real *) C)[2*(ldc * i + j)+1]) *= beta;
        }
      }
    } else {
      for (i = 0; i < N; i++) {
        for (j = 0; j < i; j++) {
          (((real *) C)[2*(ldc * i + j)]) *= beta;
          (((real *) C)[2*(ldc * i + j)+1]) *= beta;
        }
        (((real *) C)[2*(ldc * i + i)]) *= beta;
        (((real *) C)[2*(ldc * i + i)+1]) = 0.0;
      }
    }
  } else {
    for (i = 0; i < N; i++) {
      (((real *) C)[2*(ldc * i + i)+1]) = 0.0;
    }
  }
  if (alpha_real == 0.0 && alpha_imag == 0.0)
    return;
  if (uplo == Upper && trans == NoTrans) {
    for (i = 0; i < N; i++) {
      {
        real temp_real = 0.0;
        for (k = 0; k < K; k++) {
          const real Aik_real = (((const real *) A)[2*(i * lda + k)]);
          const real Aik_imag = (((const real *) A)[2*(i * lda + k)+1]);
          const real temp1_real = alpha_real * Aik_real - alpha_imag * Aik_imag;
          const real temp1_imag = alpha_real * Aik_imag + alpha_imag * Aik_real;
          const real Bik_real = (((const real *) B)[2*(i * ldb + k)]);
          const real Bik_imag = (((const real *) B)[2*(i * ldb + k)+1]);
          temp_real += temp1_real * Bik_real + temp1_imag * Bik_imag;
        }
        (((real *) C)[2*(i * ldc + i)]) += 2 * temp_real;
        (((real *) C)[2*(i * ldc + i)+1]) = 0.0;
      }
      for (j = i + 1; j < N; j++) {
        real temp_real = 0.0;
        real temp_imag = 0.0;
        for (k = 0; k < K; k++) {
          const real Aik_real = (((const real *) A)[2*(i * lda + k)]);
          const real Aik_imag = (((const real *) A)[2*(i * lda + k)+1]);
          const real temp1_real = alpha_real * Aik_real - alpha_imag * Aik_imag;
          const real temp1_imag = alpha_real * Aik_imag + alpha_imag * Aik_real;
          const real Bik_real = (((const real *) B)[2*(i * ldb + k)]);
          const real Bik_imag = (((const real *) B)[2*(i * ldb + k)+1]);
          const real Ajk_real = (((const real *) A)[2*(j * lda + k)]);
          const real Ajk_imag = (((const real *) A)[2*(j * lda + k)+1]);
          const real temp2_real = alpha_real * Ajk_real - alpha_imag * Ajk_imag;
          const real temp2_imag = alpha_real * Ajk_imag + alpha_imag * Ajk_real;
          const real Bjk_real = (((const real *) B)[2*(j * ldb + k)]);
          const real Bjk_imag = (((const real *) B)[2*(j * ldb + k)+1]);
          temp_real += ((temp1_real * Bjk_real + temp1_imag * Bjk_imag)
                        + (Bik_real * temp2_real + Bik_imag * temp2_imag));
          temp_imag += ((temp1_real * (-Bjk_imag) + temp1_imag * Bjk_real)
                        + (Bik_real * (-temp2_imag) + Bik_imag * temp2_real));
        }
        (((real *) C)[2*(i * ldc + j)]) += temp_real;
        (((real *) C)[2*(i * ldc + j)+1]) += temp_imag;
      }
    }
  } else if (uplo == Upper && trans == ConjTrans) {
    for (k = 0; k < K; k++) {
      for (i = 0; i < N; i++) {
        real Aki_real = (((const real *) A)[2*(k * lda + i)]);
        real Aki_imag = (((const real *) A)[2*(k * lda + i)+1]);
        real Bki_real = (((const real *) B)[2*(k * ldb + i)]);
        real Bki_imag = (((const real *) B)[2*(k * ldb + i)+1]);
        real temp1_real = alpha_real * Aki_real - alpha_imag * (-Aki_imag);
        real temp1_imag = alpha_real * (-Aki_imag) + alpha_imag * Aki_real;
        real temp2_real = alpha_real * Bki_real - alpha_imag * Bki_imag;
        real temp2_imag = -(alpha_real * Bki_imag + alpha_imag * Bki_real);
        {
          (((real *) C)[2*(i * lda + i)]) += 2 * (temp1_real * Bki_real - temp1_imag * Bki_imag);
          (((real *) C)[2*(i * lda + i)+1]) = 0.0;
        }
        for (j = i + 1; j < N; j++) {
          real Akj_real = (((const real *) A)[2*(k * lda + j)]);
          real Akj_imag = (((const real *) A)[2*(k * lda + j)+1]);
          real Bkj_real = (((const real *) B)[2*(k * ldb + j)]);
          real Bkj_imag = (((const real *) B)[2*(k * ldb + j)+1]);
          (((real *) C)[2*(i * lda + j)]) += (temp1_real * Bkj_real - temp1_imag * Bkj_imag)
              + (temp2_real * Akj_real - temp2_imag * Akj_imag);
          (((real *) C)[2*(i * lda + j)+1]) += (temp1_real * Bkj_imag + temp1_imag * Bkj_real)
              + (temp2_real * Akj_imag + temp2_imag * Akj_real);
        }
      }
    }
  } else if (uplo == Lower && trans == NoTrans) {
    for (i = 0; i < N; i++) {
      for (j = 0; j < i; j++) {
        real temp_real = 0.0;
        real temp_imag = 0.0;
        for (k = 0; k < K; k++) {
          const real Aik_real = (((const real *) A)[2*(i * lda + k)]);
          const real Aik_imag = (((const real *) A)[2*(i * lda + k)+1]);
          const real temp1_real = alpha_real * Aik_real - alpha_imag * Aik_imag;
          const real temp1_imag = alpha_real * Aik_imag + alpha_imag * Aik_real;
          const real Bik_real = (((const real *) B)[2*(i * ldb + k)]);
          const real Bik_imag = (((const real *) B)[2*(i * ldb + k)+1]);
          const real Ajk_real = (((const real *) A)[2*(j * lda + k)]);
          const real Ajk_imag = (((const real *) A)[2*(j * lda + k)+1]);
          const real temp2_real = alpha_real * Ajk_real - alpha_imag * Ajk_imag;
          const real temp2_imag = alpha_real * Ajk_imag + alpha_imag * Ajk_real;
          const real Bjk_real = (((const real *) B)[2*(j * ldb + k)]);
          const real Bjk_imag = (((const real *) B)[2*(j * ldb + k)+1]);
          temp_real += ((temp1_real * Bjk_real + temp1_imag * Bjk_imag)
                        + (Bik_real * temp2_real + Bik_imag * temp2_imag));
          temp_imag += ((temp1_real * (-Bjk_imag) + temp1_imag * Bjk_real)
                        + (Bik_real * (-temp2_imag) + Bik_imag * temp2_real));
        }
        (((real *) C)[2*(i * ldc + j)]) += temp_real;
        (((real *) C)[2*(i * ldc + j)+1]) += temp_imag;
      }
      {
        real temp_real = 0.0;
        for (k = 0; k < K; k++) {
          const real Aik_real = (((const real *) A)[2*(i * lda + k)]);
          const real Aik_imag = (((const real *) A)[2*(i * lda + k)+1]);
          const real temp1_real = alpha_real * Aik_real - alpha_imag * Aik_imag;
          const real temp1_imag = alpha_real * Aik_imag + alpha_imag * Aik_real;
          const real Bik_real = (((const real *) B)[2*(i * ldb + k)]);
          const real Bik_imag = (((const real *) B)[2*(i * ldb + k)+1]);
          temp_real += temp1_real * Bik_real + temp1_imag * Bik_imag;
        }
        (((real *) C)[2*(i * ldc + i)]) += 2 * temp_real;
        (((real *) C)[2*(i * ldc + i)+1]) = 0.0;
      }
    }
  } else if (uplo == Lower && trans == ConjTrans) {
    for (k = 0; k < K; k++) {
      for (i = 0; i < N; i++) {
        real Aki_real = (((const real *) A)[2*(k * lda + i)]);
        real Aki_imag = (((const real *) A)[2*(k * lda + i)+1]);
        real Bki_real = (((const real *) B)[2*(k * ldb + i)]);
        real Bki_imag = (((const real *) B)[2*(k * ldb + i)+1]);
        real temp1_real = alpha_real * Aki_real - alpha_imag * (-Aki_imag);
        real temp1_imag = alpha_real * (-Aki_imag) + alpha_imag * Aki_real;
        real temp2_real = alpha_real * Bki_real - alpha_imag * Bki_imag;
        real temp2_imag = -(alpha_real * Bki_imag + alpha_imag * Bki_real);
        for (j = 0; j < i; j++) {
          real Akj_real = (((const real *) A)[2*(k * lda + j)]);
          real Akj_imag = (((const real *) A)[2*(k * lda + j)+1]);
          real Bkj_real = (((const real *) B)[2*(k * ldb + j)]);
          real Bkj_imag = (((const real *) B)[2*(k * ldb + j)+1]);
          (((real *) C)[2*(i * lda + j)]) += (temp1_real * Bkj_real - temp1_imag * Bkj_imag)
              + (temp2_real * Akj_real - temp2_imag * Akj_imag);
          (((real *) C)[2*(i * lda + j)+1]) += (temp1_real * Bkj_imag + temp1_imag * Bkj_real)
              + (temp2_real * Akj_imag + temp2_imag * Akj_real);
        }
        {
          (((real *) C)[2*(i * lda + i)]) += 2 * (temp1_real * Bki_real - temp1_imag * Bki_imag);
          (((real *) C)[2*(i * lda + i)+1]) = 0.0;
        }
      }
    }
  } else {
    xerbla(0, "source_her2k.h", "unrecognized operation");;
  }
}
}

#endif // __BLAS_ZHER2K_HPP__
