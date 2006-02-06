/***************************************************************************
 *            test_linear_algebra.cc
 *
 *  31 Jan 2006
 *  Copyright  2006  Pieter Collins, Alberto Casagrande
 *  Email Pieter.Collins@cwi.nl, casagrande@dimi.uniud.it
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include <iostream>
#include <string>

#include "numerical_type.h"
#include "linear_algebra.h"

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::LinearAlgebra;
using namespace std;

int main() {
    cout << "test_linear_algebra: " << flush;

    Ariadne::LinearAlgebra::matrix< Rational > A(3,3), LU(3,3);
    Ariadne::LinearAlgebra::matrix< Dyadic > Aq(3,3), QT, R;
    Ariadne::LinearAlgebra::vector< Rational > b(3), x(3);
    Ariadne::LinearAlgebra::vector< dimension_type > pivot(3);
    
    A(0,0)=Rational(7,3);A(0,1)=4;A(0,2)=-1;
    A(1,0)=-13;A(1,1)=Rational(10,9);A(1,2)=6;
    A(2,0)=0;A(2,1)=Rational(1,13);A(2,2)=9;
   
    b(0)=Rational(16,3);
    b(1)=-Rational(53,9);
    b(2)=Rational(118,13);
   
    test_assert(1==13*A(2,1),"mutiplication");
   
    LU=lu_decompose(A,pivot);
    
    x=lu_solve(LU,pivot,b);
	   
    test_assert(x(0)==x(1)==x(2)==1,"lu_solve");

 /*   for (unsigned int i=0; i< A.size1(); i++) {
       for (unsigned int j=0; j< A.size2(); j++) {
	   Aq(i,j)=(Dyadic)A(i,j);
       }
    }
   
    QT=hermitian(Householder_QR(Aq));
    R=prod(QT,Aq); */
    
    cout << "PASS\n";

    return 0;
}
