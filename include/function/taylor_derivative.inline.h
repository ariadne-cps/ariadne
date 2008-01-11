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
 
#include <cassert>

#include "multi_index.h"
#include "taylor_series.h"
#include "taylor_variable.h"

namespace Ariadne {

namespace {

inline size_type compute_polynomial_data_size(size_type rs, size_type as,smoothness_type d) {
  return rs*Numeric::bin(d+as,as);
}

}

namespace Function { 
template<class X> class TaylorVariableReference { 
 public:
  TaylorVariableReference(TaylorDerivative<X>& td, const size_type& i) : _td(td), _i(i) { };
  void operator=(const TaylorVariable<X>& tv) { _td.set(_i,tv); }
  template<class XX> void operator=(const XX& x) { _td.data()[_i*_increment()]=x; }
 private:
  size_type _increment() { return compute_polynomial_data_size(1u,_td.argument_size(),_td.degree()); }
  TaylorDerivative<X>& _td; const size_type _i;
};

}

template<class X> inline
Function::TaylorDerivative<X>::TaylorDerivative()
  : _result_size(0), _argument_size(0), _degree(0), _variables() 
{
}

template<class X> inline
Function::TaylorDerivative<X>::TaylorDerivative(size_type r, size_type a, smoothness_type d)
  : _result_size(r), _argument_size(a), _degree(d), _variables(r,TaylorVariable<X>(a,d))
{
}

template<class X> template<class XX> inline
Function::TaylorDerivative<X>::TaylorDerivative(size_type r, size_type a, smoothness_type d, const XX* ptr)
  : _result_size(r), _argument_size(a), _degree(d), _variables(r,TaylorVariable<X>(a,d))
{
  for(size_type i=0; i!=r; ++i) {
    array<X>& tvd=this->_variables[i].data();
    for(size_type j=0; j!=compute_polynomial_data_size(a,d); ++j) {
      tvd[j]=*ptr;
      ++ptr;
    }
  }
}


template<class X> template<class XX> 
Function::TaylorDerivative<X>::TaylorDerivative(const TaylorSeries<XX>& ts) 
  : _result_size(1u), _argument_size(1u), _degree(ts.degree()), _variables(1u,TaylorVariable<X>(ts))
{
}
  
template<class X> template<class XX> 
Function::TaylorDerivative<X>::TaylorDerivative(const TaylorVariable<XX>& tv) 
  : _result_size(1u), _argument_size(tv.argument_size()), _degree(tv.degree()), _variables(tv.data())
{
}
  
template<class X> template<class XX> inline
Function::TaylorDerivative<X>::TaylorDerivative(const TaylorDerivative<XX>& other) 
  :  _result_size(other._result_size), _argument_size(other._argument_size), _degree(other._degree), _variables(other._variables) 
{
}
  
template<class X> template<class XX> inline
Function::TaylorDerivative<X>& 
Function::TaylorDerivative<X>::operator=(const TaylorDerivative<XX>& other) 
{
  this->_argument_size=other._argument_size;
  this->_degree=other._degree;
  this->_variables=other._variables;
  return *this;
}


template<class X> template<class XX> inline
bool 
Function::TaylorDerivative<X>::operator==(const TaylorDerivative<XX>& other) const
{
  return this->_argument_size==other._argument_size
    && this->_degree==other._degree
    && this->_variables==other._variables; 
}

template<class X> template<class XX> inline
bool 
Function::TaylorDerivative<X>::operator!=(const TaylorDerivative<XX>& other) const
{
  return !(*this==other); 
}




template<class X> inline
size_type 
Function::TaylorDerivative<X>::_increment() const 
{ 
  return compute_polynomial_data_size(1u,this->_argument_size,this->_degree); 
}

template<class X> inline
size_type 
Function::TaylorDerivative<X>::result_size() const 
{ 
  return this->_result_size;
}

template<class X> inline
size_type 
Function::TaylorDerivative<X>::argument_size() const 
{ 
  return this->_argument_size;
}

template<class X> inline
smoothness_type 
Function::TaylorDerivative<X>::degree() const 
{ 
  return this->_degree;
}

template<class X> inline
const X&
Function::TaylorDerivative<X>::get(const size_type& i, const MultiIndex& j) const 
{ 
  return this->_variables[i]._data[j.position()];
}

template<class X> inline
const Function::TaylorVariable<X>&
Function::TaylorDerivative<X>::get(const size_type& i) const 
{ 
  return this->_variables[i];
}

template<class X> template<class XX> inline
void
Function::TaylorDerivative<X>::set(const size_type& i, const MultiIndex& j, const XX& x)  
{ 
  this->_variables[i*this->_increment()+j.position()]=x;
}

template<class X> template<class XX> inline
void
Function::TaylorDerivative<X>::set(const size_type& i, const TaylorVariable<XX>& tv)
{ 
  assert(tv.argument_size())==this->_argument_size;
  assert(tv.degree()>=this->_degree);
  size_type n=compute_polynomial_data_size(this->_argument_size,this->_degree);
  for(size_type j=0; j!=n; ++j) { this->_variables[i].data()[j]=tv.data()[j]; } 
}

/*
template<class X> inline
array<X>& 
Function::TaylorDerivative<X>::data()
{
  return this->_variables; 
}

template<class X> inline
const array<X>& 
Function::TaylorDerivative<X>::data() const 
{
  return this->_variables; 
}
*/


template<class X> inline
Function::TaylorVariable<X>&
Function::TaylorDerivative<X>::operator[](const size_type& i) 
{ 
  return this->_variables[i];
}

template<class X> inline
Function::TaylorVariable<X> const&
Function::TaylorDerivative<X>::operator[](const size_type& i) const 
{ 
  return this->_variables[i];
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
  TaylorDerivative<X> y(x.result_size(),x.argument_size(),x.degree());
  for(size_type n=0; n<=y.result_size(); ++n) {
    y[n] = -x[n];
  }
  return y;
}


template<class X> inline
Function::TaylorDerivative<X> 
Function::add(const TaylorDerivative<X>& x, const TaylorDerivative<X>& y)
{
  assert(x.result_size()==y.result_size());
  assert(x.argument_size()==y.argument_size());
  TaylorDerivative<X> z(x.result_size(),x.argument_size(),std::min(x.degree(),y.degree()));
  for(size_type n=0; n<z.result_size(); ++n) {
    z[n] = x[n]+y[n];
  }
  return z;
}

template<class X> inline
Function::TaylorDerivative<X> 
Function::sub(const TaylorDerivative<X>& x, const TaylorDerivative<X>& y)
{
  assert(x.result_size()==y.result_size());
  assert(x.argument_size()==y.argument_size());
  TaylorDerivative<X> z(x.result_size(),x.argument_size(),std::min(x.degree(),y.degree()));
  for(size_type n=0; n<z.result_size(); ++n) {
    z[n] = x[n]-y[n];
  }
  return z;
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
  TaylorDerivative<X> y=x;
  X m(c);
  for(uint i=0; i!=y.result_size(); ++i) {
    y[i]*=m;
  }
  return y;
}

template<class X, class R> inline
Function::TaylorDerivative<X> 
Function::operator*(const R& c, const TaylorDerivative<X>& x)
{
  TaylorDerivative<X> y=x;
  X m(c);
  for(uint i=0; i!=y.result_size(); ++i) {
    y[i]*=m;
  }
  return y;
}

template<class X, class R> inline
Function::TaylorDerivative<X> 
Function::operator/(const TaylorDerivative<X>& x, const R& c)
{
  TaylorDerivative<X> y=x;
  for(uint i=0; i!=y.result_size(); ++i) {
    y[i]/=c;
  }
  return y;
}


}
