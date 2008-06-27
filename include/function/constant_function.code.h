/***************************************************************************
 *            constant_function.code.h
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

#include "constant_function.h"
#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "differentiation/taylor_derivative.h"
#include "differentiation/sparse_differential.h"


namespace Ariadne {

template<class R>
Vector<typename ConstantFunction<R>::F> 
ConstantFunction<R>::evaluate(const Vector<F>& x) const 
{ 
  return this->_c; 
}


template<class R>
Matrix<typename ConstantFunction<R>::F> 
ConstantFunction<R>::jacobian(const Vector<F>& x) const 
{ 
  return Matrix<F>(this->result_size(),this->argument_size()); 
}


template<class R>
TaylorDerivative<typename ConstantFunction<R>::F> 
ConstantFunction<R>::derivative(const Vector<F>& x, const smoothness_type& s) const 
{
  return TaylorDerivative<F>::constant(this->result_size(),this->argument_size(),s,this->c());
}


template<class R>
SparseDifferentialVector<typename ConstantFunction<R>::A> 
ConstantFunction<R>::expansion(const Vector<A>& x, const smoothness_type& s) const 
{
  SparseDifferentialVector<A> result(this->result_size(),this->argument_size(),s);
  for(uint i=0; i!=x.size(); ++i) { result[i]=A(this->_c[i]); }
  return result;
}


template<class R>
std::ostream& 
ConstantFunction<R>::write(std::ostream& os) const
{
  return os << "ConstantFunction( c=" << this->_c << " )";
}



} // namespace Ariadne

