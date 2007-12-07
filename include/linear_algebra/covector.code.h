/***************************************************************************
 *            covector.code.h
 *
 *  Copyright  2007  Pieter Collins
 *  Pieter.Collins@cwi.nl
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
 
#include "vector.h"
#include "covector.h"
#include "matrix.h"

#include <iostream>
#include <string>
#include <sstream>
#include <exception>

#include "../debug.h"
#include "../base/stlio.h"

namespace Ariadne {

template<class R>
void
LinearAlgebra::Covector<R>::instantiate()
{
  typedef typename Numeric::traits<R>::arithmetic_type I;
  std::ostream* os=0;
  Vector<R>* v=0;
  Vector<I>* iv=0;
  Matrix<R>* A=0;
  Matrix<I>* iA=0;
  Covector<R>* cv=0;
  Covector<I>* icv=0;
  *cv * *v;
  *cv * *iv;
  *icv * *v;
  *icv * *iv;
  *v * *cv;
  *v * *icv;
  *iv * *cv;
  *iv * *icv;
  *cv * *A;
  *cv * *iA;
  *icv * *A;
  *icv * *iA;
  *os << *cv;
  *os << *icv;
}

template<class R1, class R2>
typename Numeric::traits<R1,R2>::arithmetic_type
LinearAlgebra::operator*(const Covector<R1>& cv1, const Vector<R2>& v2)
{
  ARIADNE_CHECK_EQUAL_SIZES(cv1,v2,"Scalar operator*(Covector,Vector)");
  typedef typename Numeric::traits<R1,R2>::arithmetic_type R0;
  R0 s(0);
  for(size_type i=0; i!=cv1.size(); ++i) {
    s+=cv1(i)*v2(i);
  }
  return s;
}

template<class R1, class R2>
LinearAlgebra::Covector<typename Numeric::traits<R1,R2>::arithmetic_type>
LinearAlgebra::operator*(const Covector<R1>& cv1, const Matrix<R2>& A2)
{
  //ARIADNE_ASSERT(cv1.size()==A2.number_of_rows());
  assert(cv1.size()==A2.number_of_rows());
  typedef typename Numeric::traits<R1,R2>::arithmetic_type R0;
  LinearAlgebra::Covector<R0> cv0(A2.number_of_columns());
  for(size_type i=0; i!=A2.number_of_rows(); ++i) {
    for(size_type j=0; j!=A2.number_of_columns(); ++j) {
      cv0(j)+=cv1(i)*A2(i,j);
    }
  }
  return cv0;
}

template<class R1, class R2>
LinearAlgebra::Matrix<typename Numeric::traits<R1,R2>::arithmetic_type>
LinearAlgebra::operator*(const Vector<R1>& v1, const Covector<R2>& cv2)
{
  typedef typename Numeric::traits<R1,R2>::arithmetic_type R0;
  Matrix<R0> A(v1.size(),cv2.size());
  for(size_type i=0; i!=v1.size(); ++i) {
    for(size_type j=0; j!=cv2.size(); ++j) {
      A(i,j)=v1(i)*cv2(j);
    }
  }
  return A;
}


template<class R>
std::ostream&
LinearAlgebra::operator<<(std::ostream& os, const Covector<R>& cv) 
{  
  os << "[";
  if(cv.size()>0) {
    os << cv(0);
    for(uint i=1; i!=cv.size(); ++i) {
      os << "," << cv(i);
    }
  }
  os << "]";
  return os;
}      



}
