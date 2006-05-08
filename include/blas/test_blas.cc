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

#include "blas.hpp"
#include "blas_output.hpp"

using namespace BLAS;
using std::cout;

int main() {
  complex<double> X[3] = {complex<double>(1,1),complex<double>(0,2),3};
  complex<double> Y[3] = {4,5,6};
  complex<double> A[12] = {1,2,3,4,5,6,7,8,9,10,11,12};
  int m=3;
  int n=3;

  //axpy(3,0.1,X,1,Y,1);

  cout << vector(3,X) << "\n\n";
  cout << vector(3,Y) << "\n\n";
  cout << matrix(3,3,A) << "\n";
  cout << vector(3,X) << "\n";
  cout << vector(3,Y) << "\n";

  gemv(RowMajor, NoTrans, 3,3,complex<double>(0.01),A,3,X,1,complex<double>(2.0),Y,1);
  cout << vector(3,Y) << "\n\n";

  cout << "Exiting\n";

  return 0;
}
  
