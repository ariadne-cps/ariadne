/***************************************************************************
 *            taylor_derivative.inline.h
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

inline uint compute_polynomial_data_size(uint as, uint d) {
  return Numeric::choose(d+as,as);
}

}



template<class X> inline
Function::TaylorDerivative<X>::TaylorDerivative()
  : _argument_size(0), _degree(0), _data(1u) 
{
}

template<class X> inline
Function::TaylorDerivative<X>::TaylorDerivative(uint a, uint d)
  : _argument_size(a), _degree(d), _data(compute_polynomial_data_size(a,d))
{
}

template<class X> inline
Function::TaylorDerivative<X>::TaylorDerivative(uint a, uint d, uint i, const X& x)
  : _argument_size(a), _degree(d), _data(compute_polynomial_data_size(a,d)) 
{
  this->_data[0]=x; this->_data[i+1u]=1;
}

template<class X> inline
Function::TaylorDerivative<X>::TaylorDerivative(uint a, uint d, const X& c)
  : _argument_size(a), _degree(d), _data(compute_polynomial_data_size(a,d)) 
{
  this->_data[0]=c; 
}

template<class X> template<class X0> inline
Function::TaylorDerivative<X>::TaylorDerivative(uint a, uint d, const X0* ptr)
  : _argument_size(a), _degree(d), _data(compute_polynomial_data_size(a,d)) 
{
  for(size_type i=0; i!=this->_data.size(); ++i) {
    this->_data[i]=ptr[i];
  }
}


template<class X> template<class XX> inline
Function::TaylorDerivative<X>::TaylorDerivative(const TaylorDerivative<XX>& other) 
  : _argument_size(other._argument_size), _degree(other._degree), _data(other._data) 
{
}
  
template<class X> template<class XX> inline
Function::TaylorDerivative<X>& 
Function::TaylorDerivative<X>::operator=(const TaylorDerivative<XX>& other) 
{
  this->_argument_size=other._argument_size;
  this->_degree=other._degree;
  this->_data=other._data;
  return *this;
}


template<class X> template<class XX> inline
bool 
Function::TaylorDerivative<X>::operator==(const TaylorDerivative<XX>& other) 
{
  return this->_argument_size==other->_argument_size
    && this->_degree==other._degree
    && this->_data==other._data; 
}

template<class X> template<class XX> inline
bool 
Function::TaylorDerivative<X>::operator!=(const TaylorDerivative<XX>& other) 
{
  return !(*this==other); 
}


template<class X> inline
Function::TaylorDerivative<X> 
Function::TaylorDerivative<X>::constant(uint a, uint d, const X& c) 
{
  return TaylorDerivative<X>(a,d,c);
}

template<class X> inline
Function::TaylorDerivative<X> 
Function::TaylorDerivative<X>::variable(uint a, uint d, uint i, const X& x) 
{
  return TaylorDerivative<X>(a,d,i,x);
}



template<class X> inline
uint 
Function::TaylorDerivative<X>::argument_size() const 
{ 
  return this->_argument_size;
}

template<class X> inline
uint 
Function::TaylorDerivative<X>::degree() const 
{ 
  return this->_degree;
}

template<class X> inline
const X&
Function::TaylorDerivative<X>::value() const 
{ 
  return this->_data[0];
}

template<class X> inline
X&
Function::TaylorDerivative<X>::value()  
{ 
  return this->_data[0];
}

template<class X> inline
array<X>& 
Function::TaylorDerivative<X>::data()
{
  return this->_data; 
}

template<class X> inline
const array<X>& 
Function::TaylorDerivative<X>::data() const 
{
  return this->_data; 
}


template<class X> inline
X& 
Function::TaylorDerivative<X>::operator[](const MultiIndex& a) 
{ 
  return this->_data[a.position()]; 
}

template<class X> inline
const X& 
Function::TaylorDerivative<X>::operator[](const MultiIndex& a) const 
{ 
  return this->_data[a.position()]; 
}















template<class X> inline
Function::TaylorDerivative<X> 
Function::compose(const ScalarDerivative<X>& y, const TaylorDerivative<X>& x)
{
  TaylorDerivative<X> z(x.argument_size(),std::min(x.degree(),y.degree()));
  compute_composition(z,y,x);
  return z;
}


template<class X> inline 
Function::TaylorDerivative<X> 
Function::min(const TaylorDerivative<X>& x1, const TaylorDerivative<X>& x2) 
{
  if(x1.value()==x2.value()) {
    ARIADNE_THROW(std::runtime_error,"min(TaylorDerivative x1, TaylorDerivative x2)","x1[0]==x2[0]");
  }
  return x1.value()<x2.value() ? x1 : x2;
}

template<class X> inline 
Function::TaylorDerivative<X> 
Function::max(const TaylorDerivative<X>& x1,const TaylorDerivative<X>& x2) 
{
  if(x1.value()==x2.value()) { 
    ARIADNE_THROW(std::runtime_error,"max(TaylorDerivative x1, TaylorDerivative x2)","x1[0]==x2[0]"); 
  }
  return x1.value()>x2.value() ? x1 : x2;
}

template<class X> inline
Function::TaylorDerivative<X> 
Function::pos(const TaylorDerivative<X>& x)
{
  return x;
}

template<class X> inline
Function::TaylorDerivative<X> 
Function::neg(const TaylorDerivative<X>& x)
{
  TaylorDerivative<X> y(x.argument_size(),x.degree());
  for(uint n=0; n<=y.data().size(); ++n) {
    y.data()[n] = -x.data()[n];
  }
  return y;
}

template<class X>  
Function::TaylorDerivative<X> 
Function::abs(const TaylorDerivative<X>& x) 
{
  if(x.value()==0) { 
    ARIADNE_THROW(std::runtime_error,"abs(TaylorDerivative x)","x[0]==0"); 
  }
  return x.value()>0 ? pos(x) : neg(x); 
}

template<class X> inline
Function::TaylorDerivative<X> 
Function::inv(const TaylorDerivative<X>& x)
{
  ScalarDerivative<X> y(x.degree());
  X mr = X(-1)/x.value();
  for(uint i=0; i<=y.degree(); ++i) {
    y[i]=(-Numeric::factorial<int>(i))*pow(mr,i+1);
  }
  //std::cerr << y << std::endl;
  return compose(y,x);
}

template<class X> inline
Function::TaylorDerivative<X> 
Function::add(const TaylorDerivative<X>& x, const TaylorDerivative<X>& y)
{
  assert(x.argument_size()==y.argument_size());
  TaylorDerivative<X> z(x.argument_size(),std::min(x.degree(),y.degree()));
  for(uint n=0; n<z.data().size(); ++n) {
    z.data()[n] = x.data()[n]+y.data()[n];
  }
  return z;
}

template<class X> inline
Function::TaylorDerivative<X> 
Function::sub(const TaylorDerivative<X>& x, const TaylorDerivative<X>& y)
{
  assert(x.argument_size()==y.argument_size());
  TaylorDerivative<X> z(x.argument_size(),std::min(x.degree(),y.degree()));
  for(uint n=0; n<z.data().size(); ++n) {
    z.data()[n] = x.data()[n]-y.data()[n];
  }
  return z;
}

template<class X> inline
Function::TaylorDerivative<X> 
Function::mul(const TaylorDerivative<X>& x, const TaylorDerivative<X>& y)
{
  TaylorDerivative<X> z(x.argument_size(),std::min(x.degree(),y.degree()));
  compute_product(z,x,y);
  return z;
}

template<class X> inline
Function::TaylorDerivative<X> 
Function::div(const TaylorDerivative<X>& x, const TaylorDerivative<X>& y)
{
  return mul(x,inv(y));
}

template<class X, class N> inline
Function::TaylorDerivative<X> 
Function::pow(const TaylorDerivative<X>& x, N k)
{
  uint n=k;
  ScalarDerivative<X> y(x.degree());
  for(uint i=0; i<=std::min(y.degree(),n); ++i) {
    int j=n-i;
    y[i]=(Numeric::factorial<int>(n)/Numeric::factorial<int>(j))*Numeric::pow(x.value(),j);
  }
  return compose(y,x);
}

template<class X>  
Function::TaylorDerivative<X> 
Function::sqrt(const TaylorDerivative<X>& x) 
{
  ScalarDerivative<X> y(x.degree());
  y[0]=sqrt(x.value());
  X mhr=(-0.5)/x.value();
  for(uint i=1; i<=y.degree(); ++i) {
    y[i]=(2*int(i)-3)*mhr*y[i-1];
  }
  return compose(y,x);
}

template<class X>  
Function::TaylorDerivative<X> 
Function::exp(const TaylorDerivative<X>& x) 
{
  ScalarDerivative<X> y(x.degree());
  y[0]=exp(x.value());
  for(uint i=1; i<=y.degree(); ++i) {
    y[i]=y[0];
  }
  return compose(y,x);
}

template<class X>  
Function::TaylorDerivative<X> 
Function::log(const TaylorDerivative<X>& x) 
{
  ScalarDerivative<X> y(x.degree());
  y[0]=log(x.value());
  X mr=(-1)/x.value();
  for(uint i=1; i!=y.degree();++i) {
    y[i]=(-Numeric::factorial<int>(i-1))*pow(x.value(),i);
  }
  return compose(y,x);
}

template<class X> 
Function::TaylorDerivative<X> 
Function::sin(const TaylorDerivative<X>& x)
{
  uint d=x.degree();
  ScalarDerivative<X> y(d);
  y[0]=sin(x.value());
  y[1]=cos(x.value());
  for(uint i=2; i!=d; ++i) {
    y[i]=-y[i-2];
  }
  return compose(y,x);
}

template<class X>  
Function::TaylorDerivative<X> 
Function::cos(const TaylorDerivative<X>& x) 
{
  uint d=x.degree();
  ScalarDerivative<X> y(d);
  y[0]=cos(x.value());
  y[1]=-sin(x.value());
  for(uint i=2; i!=d; ++i) {
    y[i]=-y[i-2];
  }
  return compose(y,x);
}

template<class X>  
Function::TaylorDerivative<X> 
Function::tan(const TaylorDerivative<X>& x) 
{
  return sin(x)/cos(x);
}

template<class X>  
Function::TaylorDerivative<X> 
Function::asin(const TaylorDerivative<X>& x) 
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}

template<class X>  
Function::TaylorDerivative<X> 
Function::acos(const TaylorDerivative<X>& x) 
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}

template<class X>  
Function::TaylorDerivative<X> 
Function::atan(const TaylorDerivative<X>& x) 
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}


template<class X> inline
Function::TaylorDerivative<X> 
Function::operator-(const TaylorDerivative<X>& x)
{
  return neg(x);
}

template<class X> inline
Function::TaylorDerivative<X> 
Function::operator+(const TaylorDerivative<X>& x, const TaylorDerivative<X>& y)
{
  return add(x,y);
}

template<class X> inline
Function::TaylorDerivative<X> 
Function::operator-(const TaylorDerivative<X>& x, const TaylorDerivative<X>& y)
{
  return sub(x,y);
}

template<class X> inline
Function::TaylorDerivative<X> 
Function::operator*(const TaylorDerivative<X>& x, const TaylorDerivative<X>& y)
{
  return mul(x,y);
}

template<class X> inline
Function::TaylorDerivative<X> 
Function::operator/(const TaylorDerivative<X>& x, const TaylorDerivative<X>& y)
{
  return div(x,y);
}

template<class X, class R> inline
Function::TaylorDerivative<X> 
Function::operator+(const TaylorDerivative<X>& x, const R& c)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}

template<class X, class R> inline
Function::TaylorDerivative<X> 
Function::operator+(const R& c, const TaylorDerivative<X>& x)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}

template<class X, class R> inline
Function::TaylorDerivative<X> 
Function::operator-(const TaylorDerivative<X>& x, const R& c)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}

template<class X, class R> inline
Function::TaylorDerivative<X> 
Function::operator-(const R& c, const TaylorDerivative<X>& x)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}

template<class X, class R> inline
Function::TaylorDerivative<X> 
Function::operator*(const TaylorDerivative<X>& x, const R& c)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}

template<class X, class R> inline
Function::TaylorDerivative<X> 
Function::operator*(const R& c, const TaylorDerivative<X>& x)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}

template<class X, class R> inline
Function::TaylorDerivative<X> 
Function::operator/(const TaylorDerivative<X>& x, const R& c)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}

template<class X, class R> inline
Function::TaylorDerivative<X> 
Function::operator/(const R& c, const TaylorDerivative<X>& x)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}


}
