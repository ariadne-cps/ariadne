/***************************************************************************
 *            affine_function.code.h
 *
 *  Copyright  2005-7  Alberto Casagrande, Pieter Collins
 *
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
 
#include<iostream>

#include "affine_function.h"

#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "differentiation/taylor_derivative.h"
#include "differentiation/sparse_differential.h"


namespace Ariadne {

template<class R>
Vector<typename AffineFunction<R>::F> 
AffineFunction<R>::evaluate(const Vector<F>& x) const
{ 
  return this->_a*x+this->_b; 
}

template<class R>
Matrix<typename AffineFunction<R>::F> 
AffineFunction<R>::jacobian(const Vector<F>& x) const 
{ 
  return this->_a; 
}


template<class R>
TaylorDerivative<typename AffineFunction<R>::F> 
AffineFunction<R>::derivative(const Vector<F>& x, const smoothness_type& s) const
{
  TaylorDerivative<F> result(this->result_size(),this->argument_size(),s);
  Vector<F> value=this->evaluate(x);
  for(size_type i=0; i!=this->result_size(); ++i) {
    array<F>& data=result[i].data();
    data[0]=value[i];
    if(s>0) {
      for(size_type j=0; j!=this->argument_size(); ++j) {
        data[j+1u]=this->_a(i,j);
      }
    }
  }
  return result;
}

template<class R>
SparseDifferentialVector<typename AffineFunction<R>::AA> 
AffineFunction<R>::expansion(const Vector<AA>& x, const smoothness_type& s) const
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}

template<class R>
std::ostream& 
AffineFunction<R>::write(std::ostream& os) const
{
  return os << "AffineFunction( A=" << this->A()
            << ", b=" << this->b() << " )";
}



} // namespace Ariadne

