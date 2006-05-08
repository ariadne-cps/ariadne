/*
 * Copyright (C) 2006 Pieter Collins <Pieter.Collins@cwi.nl>
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

#include "lapack.hpp"
#include <blas/blas_output.hpp>

#include "gesv.hpp"
#include "getrf.hpp"
#include "getrs.hpp"
#include "laswp.hpp"

using std::cout;
using std::cerr;
using std::endl;

using BLAS::matrix;
using BLAS::vector;

int main() {
  double A[9]={2.0,1.0,0.1,  1.0,1.0,2.0,  5.0,0.5,4.0};
  double LU[9];
  const double I[9]={1.0,0.0,0.0,  0.0,1.0,0.0,  0.0,0.0,1.0};
  double J[9]={1.0,0.0,0.0,  0.0,1.0,0.0,  0.0,0.0,1.0};
  double Z[9];
  double E[9];
  int piv[3]={0,0,0};
  int info=0;
  
  BLAS::copy(9,A,1,LU,1);
  
  cout << matrix(3,3,A) << endl;
  cout << matrix(3,3,J) << endl;
  LAPACK::gesv(BLAS::RowMajor,3,3,LU,3,piv,J,3);
  cout << matrix(3,3,LU) << endl;
  cout << matrix(3,3,J) << endl;
  cout << vector(3,piv) << endl;
  
  BLAS::gemm(BLAS::RowMajor, BLAS::NoTrans, BLAS::NoTrans, 3, 3, 3, 1.0, A, 3, J, 3, 0.0, Z, 3);
  cout << matrix(3,3,Z) << endl;
  
  BLAS::copy(9,I,1,E,1);
  BLAS::axpy(9,-1.0,Z,1,E,1);

  cout << BLAS::amax(9,E,1) << endl;

  return 0;
}
