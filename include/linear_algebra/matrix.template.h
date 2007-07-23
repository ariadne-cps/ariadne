/***************************************************************************
 *            matrix.template.h
 *
 *  Copyright  2004-7  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it Pieter.Collins@cwi.nl
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
 

namespace Ariadne {


template<class R> 
LinearAlgebra::Matrix<R>
LinearAlgebra::over_approximation(const Matrix< Numeric::Interval<R> >& A)
{
  if(A.number_of_rows()!=A.number_of_columns()) {
    ARIADNE_THROW(NotImplemented,"Matrix<Real> over_approximation(Matrix<Interval>)","A="<<A<<" (only implemented for square matrices)"); 
  }
  dimension_type n=A.number_of_rows();
  
  Matrix<R> Amid(n,n);
  for(size_type i=0; i!=n; ++i) {
    for(size_type j=0; j!=n; ++j) {
      Amid(i,j)=(A(i,j).upper()+A(i,j).lower())/2;
    }
  }
  Matrix< Numeric::Interval<R> > I=LinearAlgebra::inverse(Matrix< Numeric::Interval<R> >(Amid))*A;
  
  R excess=LinearAlgebra::norm(I).upper();
  
  // FIXME: Outer bound on multiplication
  return excess*Amid;
}


} // namespace Ariadne
