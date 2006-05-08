/*
 * Copyright (C) 2006 Pieter Collins <Pieter.Collins@cwi.nl>
 *
 * Based on the routine in LAPACK (version 3.0)
 *   Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
 *   Courant Institute, Argonne National Lab, and Rice University
 */

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef __LAPACK_GESVD_HPP__
#define __LAPACK_GESVD_HPP__

#include "lapack.hpp"

template<typename real>
void
LAPACK::gesvd(BLAS::ORDER order, int m, int n, real *A, int ldA,
              real *S, real *U, int ldU, real *V, int ldV)
{
  std::cerr << "LAPACK::gesvd\n";
  assert(order==BLAS::RowMajor);

  real zero=0;
  real one=1;
  
  real *Ain=new real[m*ldA];
  BLAS::copy(m*ldA,A,1,Ain,1);
  
  int nu = min(m,n);
  real *E=new real[n];
  real *work=new real[m];

  bool wantu = true;
  bool wantv = true;
  
  // Reduce A to bidiagonal form, storing the diagonal elements
  // in S and the super-diagonal elements in E.
  
  int nct = min(m-1,n);
  int nrt = max(0,min(n-2,m));
  for (int k = 0; k < max(nct,nrt); k++) {
    if (k < nct) {
      
      // Compute the transformation for the k-th column and
      // place the k-th diagonal in S[k].
      // Compute 2-norm of k-th column without under/overflow.
      S[k] = BLAS::nrm2(m-k,&A[k*ldA+k],ldA);
      if (S[k] != zero) {
        if (A[k*ldA+k] < zero) {
          S[k] = -S[k];
        }
        for (int i = k; i < m; i++) {
          A[i*ldA+k] /= S[k];
        }
        A[k*ldA+k] += one;
      }
      S[k] = -S[k];
      //std::cerr << "k=" << k << "  A[k:m,k]=" << vector(m-k,&A[k*ldA+k],ldA) << "  S[k]=" << S[k] << "\n";
    }
    for (int j = k+1; j < n; j++) {
      if ((k < nct) && (S[k] != zero))  {
        // Apply the transformation.
        real t = zero;
        for (int i = k; i < m; i++) {
          t += A[i*ldA+k]*A[i*ldA+j];
        }
        t = -t/A[k*ldA+k];
        for (int i = k; i < m; i++) {
          A[i*ldA+j] += t*A[i*ldA+k];
        }
      }
    }
    if (wantu & (k < nct)) {
      // Place the transformation in U for subsequent back
      // multiplication.
      for (int i = k; i < m; i++) {
        U[i*ldU+k] = A[i*ldA+k];
      }
    }
 
    if (k < nrt) {      
      // Compute the k-th row transformation and place the
      // k-th super-diagonal in E[k].
      // Compute 2-norm without under/overflow.
      E[k] = BLAS::nrm2(n-k-1,&A[k*ldA+k+1],1);
      if (E[k] != zero) {
        if (A[k*ldA+k+1] < 0) {
          E[k] = -E[k];
        }
        for (int i = k+1; i < n; i++) {
          A[k*ldA+i] /= E[k];
        }
        A[k*ldA+k+1] += one;
      }
      E[k] = -E[k];
      //std::cerr << "k=" << k << "  A[k,k+1,n]=" << vector(n-k-1,&A[k*ldA+k+1],1) << "  E[k]=" << E[k] << "\n";
      if ((k+1 < m) & (E[k] != zero)) {
        // Apply the transformation.
        for (int i = k+1; i < m; i++) {
          work[i] = zero;
        }
        for (int j = k+1; j < n; j++) {
          for (int i = k+1; i < m; i++) {
            work[i] += A[k*ldA+j]*A[i*ldA+j];
          }
        }
        for (int j = k+1; j < n; j++) {
          real t = -A[k*ldA+j]/A[k*ldA+k+1];
          for (int i = k+1; i < m; i++) {
            A[i*ldA+j] += t*work[i];
          }
        }
      }
      if (wantv) {
        // Place the transformation in V for subsequent
        // back multiplication.
        for (int i = k+1; i < n; i++) {
          V[i*ldV+k+1] = A[k*ldA+i];
        }
      }
    }
    
  }
  
  // Set up the final bidiagonal matrix of order p.
  
  int p = min(n,m+1);
  if (nct < n) {
    S[nct] = A[nct*ldA+nct];
  }
  if (m < p) {
    S[p-1] = zero;
  }
  if (nrt+1 < p) {
    E[nrt] = A[nrt*ldA+p-1];
  }
  E[p-1] = zero;
  
  //std::cerr << "Working values\n";
  //std::cerr << "A=\n" << matrix(m,n,A,ldA);
  //std::cerr << "U=\n" << matrix(m,m,U,ldU);
  //std::cerr << "V=\n" << matrix(n,n,V,ldV);
  //std::cerr << "S=" << vector(n,S,1) << "\n";
  //std::cerr << "E=" << vector(n,E,1) << "\n";
  
  // If required, generate U.
  if (wantu) {
    //std::cerr << "U=\n" << matrix(m,m,U,ldU);
    for (int j = nct; j < nu; j++) {
      for (int i = 0; i < m; i++) {
        U[i*ldU+j] = zero;
      }
      U[j*ldU+j] = one;
    }
    //std::cerr << "U=\n" << matrix(m,m,U,ldU);
    for (int k = nct-1; k >= 0; k--) {
      //std::cerr << "nct=" << nct << "  k=" << k << "S[k]=" << S[k] << "\n";
      if (S[k] != zero) {
        for (int j = k+1; j < nu; j++) {
          real t = zero;
          for (int i = k; i < m; i++) {
            t += U[i*ldU+k]*U[i*ldU+j];
          }
          t = -t/U[k*ldU+k];
          for (int i = k; i < m; i++) {
            U[i*ldU+j] += t*U[i*ldU+k];
          }
        }
        for (int i = k; i < m; i++ ) {
          U[i*ldU+k] = -U[i*ldU+k];
        }
        U[k*ldU+k] = one + U[k*ldU+k];
        for (int i = 0; i < k; i++) {
          U[i*ldU+k] = zero;
        }
      } else {
        for (int i = 0; i < m; i++) {
          U[i*ldU+k] = zero;
        }
        U[k*ldU+k] = one;
      }
      //std::cerr << "U=\n" << matrix(m,m,U,ldU);
    }
  }
  
  // If required, generate V.
  
  if (wantv) {
    //std::cerr << "V=\n" << matrix(n,n,V,ldV);
    for (int k = n-1; k >= 0; k--) {
      //std::cerr << "nrt=" << nrt << "  k=" << k << "  E[k]=" << E[k] << "\n";
      if ((k!=0) && (k <= nrt) && (E[k-1] != zero)) {
        for (int j = k+1; j < n; j++) {
          real t = zero;
          for (int i = k; i < n; i++) {
            t += V[i*ldV+k]*V[i*ldV+j];
          }
          t = -t/V[k*ldV+k];
          for (int i = k; i < n; i++) {
            V[i*ldV+j] += t*V[i*ldV+k];
          }
        }
        for (int i = k; i < n; i++ ) {
          V[i*ldV+k] = -V[i*ldV+k];
        }
        V[k*ldV+k] = one + V[k*ldV+k];
        for (int i = 0; i < k; i++) {
          V[i*ldV+k] = zero;
        }
      } else {
        for (int i = 0; i < n; i++) {
          V[i*ldV+k] = zero;
        }
        V[k*ldV+k] = one;
      } 
      //std::cerr << "V=\n" << matrix(n,n,V,ldV);
    }
  }
  
  /*
  std::cerr << "Intermediate values:\n";
  std::cerr << "U=\n" << BLAS::matrix(m,m,U,ldU);
  std::cerr << "V=\n" << BLAS::matrix(n,n,V,ldV);
  real *BD=new real(m*n);
  int ldBD=n;
  for(int i=0; i!=m; ++i) {
    for(int j=0; j!=n; ++j) {
      BD[i*ldBD+j]=0;
    }
    if(i<n) { BD[i*ldBD+i]=S[i]; }
    if(i+1<n) { BD[i*ldBD+i+1]=E[i]; }
  }
  std::cerr << "BD=\n" << BLAS::matrix(m,n,BD,ldBD) << std::endl << std::endl;
  real *I=new real(m*m);
  BLAS::gemm(BLAS::RowMajor,BLAS::Trans,BLAS::NoTrans,m,m,m,1,U,ldU,U,ldU,0,I,m);
  std::cerr << "I=\n" << BLAS::matrix(m,m,I) << std::endl;
  real *J=new real(n*n);
  BLAS::gemm(BLAS::RowMajor,BLAS::Trans,BLAS::NoTrans,n,n,n,1,V,ldV,V,ldV,0,J,n);
  std::cerr << "J=\n" << BLAS::matrix(n,n,J) << std::endl;
  BLAS::axpy(n,real(-1),&one,0,J,n+1);
  std::cerr << "error=" << BLAS::amax(n*n,J,1) << std::endl;
  delete[] BD; delete[] I; delete[] J;
  */
  
   
  // Main iteration loop for the singular values.
  int pp = p-1;
  int iter = 0;
  real eps = pow(2.0,-52);
  //std::cerr << "eps=" << eps << "\n";
  while (p > 0) {
    int k=0;
    int kase=0;
    
    // Here is where a test for too many iterations would go.
    
    // This section of the program inspects for
    // negligible elements in the S and E arrays.  On
    // completion the variables kase and k are set as follows.
    
    // kase = 1     if S(p) and E[k-1] are negligible and k<p
    // kase = 2     if S(k) is negligible and k<p
    // kase = 3     if E[k-1] is negligible, k<p, and
    //              S(k), ..., S(p) are not negligible (qr step).
    // kase = 4     if E(p-1) is negligible (convergence).
    
    for (k = p-2; k >= -1; k--) {
      if (k == -1) {
        break;
      }
      if (abs(E[k]) <= eps*(abs(S[k]) + abs(S[k+1]))) {
        E[k] = zero;
        break;
      }
    }
    if (k == p-2) {
      kase = 4;
    } else {
      int ks;
      for (ks = p-1; ks >= k; ks--) {
        if (ks == k) {
          break;
        }
        real t = (ks != p ? abs(E[ks]) : 0.) + 
          (ks != k+1 ? abs(E[ks-1]) : 0.);
        if (abs(S[ks]) <= eps*t)  {
          S[ks] = zero;
          break;
        }
      }
      if (ks == k) {
        kase = 3;
      } else if (ks == p-1) {
        kase = 1;
      } else {
        kase = 2;
        k = ks;
      }
    }
    k++;
    
    // Perform the task indicated by kase.
  
    switch (kase) {
      
      // Deflate negligible S(p).
      case 1: {
        real f = E[p-2];
        E[p-2] = zero;
        for (int j = p-2; j >= k; j--) {
          real t = hypot(S[j],f);
          real cs = S[j]/t;
          real sn = f/t;
          S[j] = t;
          if (j != k) {
            f = -sn*E[j-1];
            E[j-1] = cs*E[j-1];
          }
          if (wantv) {
            for (int i = 0; i < n; i++) {
              t = cs*V[i*ldV+j] + sn*V[i*ldV+p-1];
              V[i*ldV+p-1] = -sn*V[i*ldV+j] + cs*V[i*ldV+p-1];
              V[i*ldV+j] = t;
            }
          }
        }
      }
      break;
      
      // Split at negligible S(k).
      case 2: {
        real f = E[k-1];
        E[k-1] = zero;
        for (int j = k; j < p; j++) {
          real t = hypot(S[j],f);
          real cs = S[j]/t;
          real sn = f/t;
          S[j] = t;
          f = -sn*E[j];
          E[j] = cs*E[j];
          if (wantu) {
            for (int i = 0; i < m; i++) {
              t = cs*U[i*ldU+j] + sn*U[i*ldU+k-1];
              U[i*ldU+k-1] = -sn*U[i*ldU+j] + cs*U[i*ldU+k-1];
              U[i*ldU+j] = t;
            }
          }
        }
      }
      break;
      
      // Perform one qr step.      
      case 3: {
        
        // Calculate the shift.
        real scale = max(max(max(max(
                                       abs(S[p-1]),abs(S[p-2])),abs(E[p-2])), 
                               abs(S[k])),abs(E[k]));
        real sp = S[p-1]/scale;
        real spm1 = S[p-2]/scale;
        real epm1 = E[p-2]/scale;
        real sk = S[k]/scale;
        real ek = E[k]/scale;
        real b = ((spm1 + sp)*(spm1 - sp) + epm1*epm1)/2.0;
        real c = (sp*epm1)*(sp*epm1);
        real shift = zero;
        if ((b != zero) || (c != zero)) {
          shift = sqrt(b*b + c);
          if (b < 0) {
            shift = -shift;
          }
          shift = c/(b + shift);
        }
        real f = (sk + sp)*(sk - sp) + shift;
        real g = sk*ek;
        
        // Chase zeros.
        for (int j = k; j < p-1; j++) {
          real t = hypot(f,g);
          real cs = f/t;
          real sn = g/t;
          if (j != k) {
            E[j-1] = t;
          }
          f = cs*S[j] + sn*E[j];
          E[j] = cs*E[j] - sn*S[j];
          g = sn*S[j+1];
          S[j+1] = cs*S[j+1];
          if (wantv) {
            for (int i = 0; i < n; i++) {
              t = cs*V[i*ldV+j] + sn*V[i*ldV+j+1];
              V[i*ldV+j+1] = -sn*V[i*ldV+j] + cs*V[i*ldV+j+1];
              V[i*ldV+j] = t;
            }
          }
          t = hypot(f,g);
          cs = f/t;
          sn = g/t;
          S[j] = t;
          f = cs*E[j] + sn*S[j+1];
          S[j+1] = -sn*E[j] + cs*S[j+1];
          g = sn*E[j+1];
          E[j+1] = cs*E[j+1];
          if (wantu && (j < m-1)) {
            for (int i = 0; i < m; i++) {
              t = cs*U[i*ldU+j] + sn*U[i*ldU+j+1];
              U[i*ldU+j+1] = -sn*U[i*ldU+j] + cs*U[i*ldU+j+1];
              U[i*ldU+j] = t;
            }
          }
        }
        E[p-2] = f;
        iter = iter + 1;
      }
      break;

      // Convergence.      
      case 4: {
        
        // Make the singular values positive.
        if (S[k] <= zero) {
          S[k] = (S[k] < 0 ? -S[k] : real(0));
          if (wantv) {
            for (int i = 0; i <= pp; i++) {
              V[i*ldV+k] = -V[i*ldV+k];
            }
          }
        }
        
        // Order the singular values.
        while (k < pp) {
          if (S[k] >= S[k+1]) {
            break;
          }
          real t = S[k];
          S[k] = S[k+1];
          S[k+1] = t;
          if (wantv && (k < n-1)) {
            for (int i = 0; i < n; i++) {
              t = V[i*ldV+k+1]; V[i*ldV+k+1] = V[i*ldV+k]; V[i*ldV+k] = t;
            }
          }
          if (wantu && (k < m-1)) {
            for (int i = 0; i < m; i++) {
              t = U[i*ldU+k+1]; U[i*ldU+k+1] = U[i*ldU+k]; U[i*ldU+k] = t;
            }
          }
          k++;
        }
        iter = 0;
        p--;
      }
      break;
    }
  
  }
  
  delete[] E;
  delete[] work;
}

#endif // __LAPACK_GESVD_HPP__
