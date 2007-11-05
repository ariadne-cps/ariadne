/***************************************************************************
 *            multivariable_derivative.inline.h
 *
 *  Copyright 2007  Alberto Casagrande, Pieter Collins
 *  Email casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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
 
#include "multi_index.h"
#include "scalar_derivative.h"

namespace Ariadne {

namespace {

inline uint compute_polynomial_data_size(uint rs, uint as, uint d) {
  return rs*Numeric::choose(d+as,as);
}

}



template<class X> inline
Function::MultivariableDerivative<X>::MultivariableDerivative()
  : _result_size(0), _argument_size(0), _degree(0), _data() 
{
}

template<class X> inline
Function::MultivariableDerivative<X>::MultivariableDerivative(uint r, uint a, uint d)
  : _result_size(r), _argument_size(a), _degree(d), _data(compute_polynomial_data_size(r,a,d))
{
}

template<class X> template<class XX> inline
Function::MultivariableDerivative<X>::MultivariableDerivative(uint r, uint a, uint d, const XX* ptr)
  : _result_size(r), _argument_size(a), _degree(d), _data(compute_polynomial_data_size(r,a,d)) 
{
  for(uint i=0; i!=this->_data.size(); ++i) {
    this->_data[i]=ptr[i];
  }
}


template<class X> template<class XX> inline
Function::MultivariableDerivative<X>::MultivariableDerivative(const MultivariableDerivative<XX>& other) 
  : _result_size(other._result_size), _argument_size(other._argument_size), _degree(other._degree), _data(other._data) 
{
}
  
template<class X> template<class XX> inline
Function::MultivariableDerivative<X>& 
Function::MultivariableDerivative<X>::operator=(const MultivariableDerivative<XX>& other) 
{
  this->_result_size=other._result_size;
  this->_argument_size=other._argument_size;
  this->_degree=other._degree;
  this->_data=other._data;
  return *this;
}


template<class X> template<class XX> inline
bool 
Function::MultivariableDerivative<X>::operator==(const MultivariableDerivative<XX>& other) 
{
  return this->_result_size==other._result_size
    && this->_argument_size==other->_argument_size
    && this->_degree==other._degree
    && this->_data==other._data; 
}

template<class X> template<class XX> inline
bool 
Function::MultivariableDerivative<X>::operator!=(const MultivariableDerivative<XX>& other) 
{
  return !(*this==other); 
}



template<class X> inline
uint 
Function::MultivariableDerivative<X>::result_size() const 
{ 
  return this->_result_size;
}

template<class X> inline
uint 
Function::MultivariableDerivative<X>::argument_size() const 
{ 
  return this->_argument_size;
}

template<class X> inline
uint 
Function::MultivariableDerivative<X>::degree() const 
{ 
  return this->_degree;
}

template<class X> inline
const LinearAlgebra::Vector<X>
Function::MultivariableDerivative<X>::value() const 
{ 
  return LinearAlgebra::Vector<X>(this->_result_size,this->_data.begin());
}

template<class X> inline
array<X>& 
Function::MultivariableDerivative<X>::data()
{
  return this->_data; 
}

template<class X> inline
const array<X>& 
Function::MultivariableDerivative<X>::data() const 
{
  return this->_data; 
}




template<class X> inline
const X& 
Function::MultivariableDerivative<X>::get(const uint& i, const MultiIndex& j) const
{
  return this->_data[i+this->_result_size*j.position()];
}

template<class X> template<class XX> inline
void 
Function::MultivariableDerivative<X>::set(const uint& i, const MultiIndex& j, const XX& x) 
{
  this->_data[i+this->_result_size*j.position()]=x;
}

template<class X> template<class XX> inline
void 
Function::MultivariableDerivative<X>::set(const uint& i, const TaylorDerivative<XX>& x) 
{
  for(MultiIndex j(this->argument_size()); j.degree()<=this->degree(); ++j) {
    this->set(i,j,x[j]);
  }
}







template<class X> inline
Function::TaylorDerivative<X>& 
Function::MultivariableDerivative<X>::operator[](const uint& i) 
{ 
  throw NotImplemented(__PRETTY_FUNCTION__);
}

template<class X> inline
const Function::TaylorDerivative<X> 
Function::MultivariableDerivative<X>::operator[](const uint& i) const 
{ 
  uint m=this->_data.size()/this->_result_size;
  LinearAlgebra::Vector<X> v(m,this->_data.begin()+i,this->_result_size);
  return TaylorDerivative<X>(this->_argument_size,this->_degree,v.data().begin());
}

template<class X> inline
X& 
Function::MultivariableDerivative<X>::operator()(const uint& i, const MultiIndex& j) 
{ 
  return this->_data[i+this->_result_size*j.position()];
}

template<class X> inline
const X&
Function::MultivariableDerivative<X>::operator()(const uint& i, const MultiIndex& j) const 
{ 
  return this->_data[i+this->_result_size*j.position()];
}










template<class X> inline
Function::MultivariableDerivative<X> 
Function::neg(const MultivariableDerivative<X>& x1)
{
  MultivariableDerivative<X> x0(x1.result_size(),x1.argument_size(),x1.degree());
  LinearAlgebra::Vector<X>& v0=reinterpret_cast<LinearAlgebra::Vector<X>&>(x0.data());
  const LinearAlgebra::Vector<X>& v1=reinterpret_cast<const LinearAlgebra::Vector<X>&>(x1.data());
  v0=-v1;
  return x0;
}

template<class X> inline
Function::MultivariableDerivative<X> 
Function::add(const MultivariableDerivative<X>& x1, const MultivariableDerivative<X>& x2)
{
  assert(x1.result_size()==x2.result_size());
  assert(x1.argument_size()==x2.argument_size());
  assert(x1.degree()==x2.degree());
  MultivariableDerivative<X> x0(x1.result_size(),x1.argument_size(),std::min(x1.degree(),x2.degree()));
  LinearAlgebra::Vector<X>& v0=reinterpret_cast<LinearAlgebra::Vector<X>&>(x0.data());
  const LinearAlgebra::Vector<X>& v1=reinterpret_cast<const LinearAlgebra::Vector<X>&>(x1.data());
  const LinearAlgebra::Vector<X>& v2=reinterpret_cast<const LinearAlgebra::Vector<X>&>(x2.data());
  v0=v1+v2;
  return x0;
}

template<class X> inline
Function::MultivariableDerivative<X> 
Function::sub(const MultivariableDerivative<X>& x1, const MultivariableDerivative<X>& x2)
{
  assert(x1.result_size()==x2.result_size());
  assert(x1.argument_size()==x2.argument_size());
  assert(x1.degree()==x2.degree());
  MultivariableDerivative<X> x0(x1.result_size(),x1.argument_size(),std::min(x1.degree(),x2.degree()));
  LinearAlgebra::Vector<X>& v0=reinterpret_cast<LinearAlgebra::Vector<X>&>(x0.data());
  const LinearAlgebra::Vector<X>& v1=reinterpret_cast<const LinearAlgebra::Vector<X>&>(x1.data());
  const LinearAlgebra::Vector<X>& v2=reinterpret_cast<const LinearAlgebra::Vector<X>&>(x2.data());
  v0=v1-v2;
  return x0;
}






template<class X> inline
Function::MultivariableDerivative<X> 
Function::compose(const MultivariableDerivative<X>& y, const MultivariableDerivative<X>& x)
{
  assert(y.argument_size()==x.result_size());
  MultivariableDerivative<X> z(y.result_size(),x.argument_size(),std::min(x.degree(),y.degree()));
  compute_composition(z,y,x);
  return z;
}



template<class X> inline
Function::MultivariableDerivative<X> 
Function::operator-(const MultivariableDerivative<X>& x)
{
  return neg(x);
}

template<class X> inline
Function::MultivariableDerivative<X> 
Function::operator+(const MultivariableDerivative<X>& x, const MultivariableDerivative<X>& y)
{
  return add(x,y);
}

template<class X> inline
Function::MultivariableDerivative<X> 
Function::operator-(const MultivariableDerivative<X>& x, const MultivariableDerivative<X>& y)
{
  return sub(x,y);
}

template<class X, class R> inline
Function::MultivariableDerivative<X> 
Function::operator+(const MultivariableDerivative<X>& x, const R& c)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}

template<class X, class R> inline
Function::MultivariableDerivative<X> 
Function::operator+(const R& c, const MultivariableDerivative<X>& x)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}

template<class X, class R> inline
Function::MultivariableDerivative<X> 
Function::operator-(const MultivariableDerivative<X>& x, const R& c)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}

template<class X, class R> inline
Function::MultivariableDerivative<X> 
Function::operator-(const R& c, const MultivariableDerivative<X>& x)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}

template<class X, class R> inline
Function::MultivariableDerivative<X> 
Function::operator*(const MultivariableDerivative<X>& x, const R& c)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}

template<class X, class R> inline
Function::MultivariableDerivative<X> 
Function::operator*(const R& c, const MultivariableDerivative<X>& x)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}

template<class X, class R> inline
Function::MultivariableDerivative<X> 
Function::operator/(const MultivariableDerivative<X>& x, const R& c)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}

template<class X, class R> inline
Function::MultivariableDerivative<X> 
Function::operator/(const R& c, const MultivariableDerivative<X>& x)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}


}
