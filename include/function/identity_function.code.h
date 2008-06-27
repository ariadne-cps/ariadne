/***************************************************************************
 *            identity_function.code.h
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

#include "identity_function.h"

#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "differentiation/taylor_derivative.h"
#include "differentiation/sparse_differential.h"


namespace Ariadne {

template<class R>
Vector<typename IdentityFunction<R>::F> 
IdentityFunction<R>::evaluate(const Vector<F>& x) const
{ 
  return x; 
}

template<class R>
Matrix<typename IdentityFunction<R>::F> 
IdentityFunction<R>::jacobian(const Vector<F>& x) const 
{ 
  return Matrix<F>::identity(this->_dimension);
}


template<class R>
TaylorDerivative<typename IdentityFunction<R>::F> 
IdentityFunction<R>::derivative(const Vector<F>& x, const smoothness_type& s) const
{
  TaylorDerivative<F> result(this->result_size(),this->argument_size(),s);
  for(size_type i=0; i!=this->result_size(); ++i) {
    array<F>& data=result[i].data();
    data[0]=x[i];
    if(s>0) {
      data[i+1u]=1;
    }
  }
  return result;
}

template<class R>
SparseDifferentialVector<typename IdentityFunction<R>::A> 
IdentityFunction<R>::expansion(const Vector<A>& x, const smoothness_type& s) const
{
  SparseDifferentialVector<A> result(this->result_size(),this->argument_size(),s);
  for(size_type i=0; i!=this->result_size(); ++i) {
    result[i]=x[i];
    result[i][i]=A(1);
  }
  return result;
}

template<class R>
std::ostream& 
IdentityFunction<R>::write(std::ostream& os) const
{
  return os << "IdentityFunction( dimension=" << this->_dimension << " )";
}



} // namespace Ariadne
