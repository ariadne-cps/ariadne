/***************************************************************************
 *            differential_vector.inline.h
 *
 *  Copyright 2007  Pieter Collins
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
 
#include <cassert>

#include "multi_index.h"
#include "power_series.h"
#include "differential.h"

namespace Ariadne {

namespace {

inline size_type compute_polynomial_data_size(size_type rs, size_type as,smoothness_type d) {
  return rs*bin(d+as,as);
}

} // namespace


template<class X> class DifferentialReference { 
 public:
  DifferentialReference(DifferentialVector<X>& td, const size_type& i) : _td(td), _i(i) { };
  void operator=(const Differential<X>& tv) { _td.set(_i,tv); }
  template<class XX> void operator=(const XX& x) { _td.data()[_i*_increment()]=x; }
 private:
  size_type _increment() { return compute_polynomial_data_size(1u,_td.argument_size(),_td.degree()); }
  DifferentialVector<X>& _td; const size_type _i;
};


template<class X> inline
DifferentialVector<X>::DifferentialVector()
  : _result_size(0), _argument_size(0), _degree(0), _variables() 
{
}

template<class X> inline
DifferentialVector<X>::DifferentialVector(size_type r, size_type a, smoothness_type d)
  : _result_size(r), _argument_size(a), _degree(d), _variables(r,Differential<X>(a,d))
{
}

template<class X> template<class XX> inline
DifferentialVector<X>::DifferentialVector(size_type r, size_type a, smoothness_type d, const XX* ptr)
  : _result_size(r), _argument_size(a), _degree(d), _variables(r,Differential<X>(a,d))
{
  for(size_type i=0; i!=r; ++i) {
    array<X>& tvd=this->_variables[i].data();
    assert(tvd.size()==compute_polynomial_data_size(a,d));
    for(size_type j=0; j!=tvd.size(); ++j) {
      tvd[j]=*ptr;
      ++ptr;
    }
  }
}


template<class X> template<class XX> 
DifferentialVector<X>::DifferentialVector(const PowerSeries<XX>& ts) 
  : _result_size(1u), _argument_size(1u), _degree(ts.degree()), _variables(1u,Differential<X>(ts))
{
}
  
template<class X> template<class XX> 
DifferentialVector<X>::DifferentialVector(const Differential<XX>& tv) 
  : _result_size(1u), _argument_size(tv.argument_size()), _degree(tv.degree()), _variables(tv.data())
{
}
  
template<class X> template<class XX> inline
DifferentialVector<X>::DifferentialVector(const DifferentialVector<XX>& other) 
  :  _result_size(other.result_size()), _argument_size(other.argument_size()), _degree(other.degree()), _variables(other.variables()) 
{
}
  
template<class X> template<class XX> inline
DifferentialVector<X>& 
DifferentialVector<X>::operator=(const DifferentialVector<XX>& other) 
{
  this->_argument_size=other._argument_size;
  this->_degree=other._degree;
  this->_variables=other._variables;
  return *this;
}


template<class X> template<class XX> inline
bool 
DifferentialVector<X>::operator==(const DifferentialVector<XX>& other) const
{
  return this->_argument_size==other._argument_size
    && this->_degree==other._degree
    && this->_variables==other._variables; 
}

template<class X> template<class XX> inline
bool 
DifferentialVector<X>::operator!=(const DifferentialVector<XX>& other) const
{
  return !(*this==other); 
}




template<class X> inline
size_type 
DifferentialVector<X>::_increment() const 
{ 
  return compute_polynomial_data_size(1u,this->_argument_size,this->_degree); 
}

// Synonym for "result_size" for templated evaluate function
template<class X> inline
size_type 
DifferentialVector<X>::size() const 
{ 
  return this->_result_size;
}

template<class X> inline
size_type 
DifferentialVector<X>::result_size() const 
{ 
  return this->_result_size;
}

template<class X> inline
size_type 
DifferentialVector<X>::argument_size() const 
{ 
  return this->_argument_size;
}

template<class X> inline
smoothness_type 
DifferentialVector<X>::degree() const 
{ 
  return this->_degree;
}

template<class X> inline
const X&
DifferentialVector<X>::get(const size_type& i, const MultiIndex& j) const 
{ 
  return this->_variables[i]._data[j.position()];
}

template<class X> inline
const Differential<X>&
DifferentialVector<X>::get(const size_type& i) const 
{ 
  return this->_variables[i];
}

template<class X> template<class XX> inline
void
DifferentialVector<X>::set(const size_type& i, const MultiIndex& j, const XX& x)  
{ 
  this->_variables[i*this->_increment()+j.position()]=x;
}

template<class X> template<class XX> inline
void
DifferentialVector<X>::set(const size_type& i, const Differential<XX>& tv)
{ 
  assert(tv.argument_size())==this->_argument_size;
  assert(tv.degree()>=this->_degree);
  size_type n=compute_polynomial_data_size(this->_argument_size,this->_degree);
  for(size_type j=0; j!=n; ++j) { this->_variables[i].data()[j]=tv.data()[j]; } 
}

/*
template<class X> inline
array<X>& 
DifferentialVector<X>::data()
{
  return this->_variables; 
}

template<class X> inline
const array<X>& 
DifferentialVector<X>::data() const 
{
  return this->_variables; 
}
*/


template<class X> inline
Differential<X>&
DifferentialVector<X>::operator[](const size_type& i) 
{ 
  return this->_variables[i];
}

template<class X> inline
Differential<X> const&
DifferentialVector<X>::operator[](const size_type& i) const 
{ 
  return this->_variables[i];
}
















template<class X> inline 
DifferentialVector<X> 
min(const DifferentialVector<X>& x1, const DifferentialVector<X>& x2) 
{
  if(x1.value()==x2.value()) {
    ARIADNE_THROW(std::runtime_error,"min(DifferentialVector x1, DifferentialVector x2)","x1[0]==x2[0]");
  }
  return x1.value()<x2.value() ? x1 : x2;
}

template<class X> inline 
DifferentialVector<X> 
max(const DifferentialVector<X>& x1,const DifferentialVector<X>& x2) 
{
  if(x1.value()==x2.value()) { 
    ARIADNE_THROW(std::runtime_error,"max(DifferentialVector x1, DifferentialVector x2)","x1[0]==x2[0]"); 
  }
  return x1.value()>x2.value() ? x1 : x2;
}

template<class X> inline
DifferentialVector<X> 
pos(const DifferentialVector<X>& x)
{
  return x;
}

template<class X> inline
DifferentialVector<X> 
neg(const DifferentialVector<X>& x)
{
  DifferentialVector<X> y(x.result_size(),x.argument_size(),x.degree());
  for(size_type n=0; n<=y.result_size(); ++n) {
    y[n] = -x[n];
  }
  return y;
}


template<class X> inline
DifferentialVector<X> 
add(const DifferentialVector<X>& x, const DifferentialVector<X>& y)
{
  assert(x.result_size()==y.result_size());
  assert(x.argument_size()==y.argument_size());
  DifferentialVector<X> z(x.result_size(),x.argument_size(),std::min(x.degree(),y.degree()));
  for(size_type n=0; n<z.result_size(); ++n) {
    z[n] = x[n]+y[n];
  }
  return z;
}

template<class X> inline
DifferentialVector<X> 
sub(const DifferentialVector<X>& x, const DifferentialVector<X>& y)
{
  assert(x.result_size()==y.result_size());
  assert(x.argument_size()==y.argument_size());
  DifferentialVector<X> z(x.result_size(),x.argument_size(),std::min(x.degree(),y.degree()));
  for(size_type n=0; n<z.result_size(); ++n) {
    z[n] = x[n]-y[n];
  }
  return z;
}


template<class X> inline
DifferentialVector<X> 
operator-(const DifferentialVector<X>& x)
{
  return neg(x);
}

template<class X> inline
DifferentialVector<X> 
operator+(const DifferentialVector<X>& x, const DifferentialVector<X>& y)
{
  return add(x,y);
}

template<class X> inline
DifferentialVector<X> 
operator-(const DifferentialVector<X>& x, const DifferentialVector<X>& y)
{
  return sub(x,y);
}


template<class X, class R> inline
DifferentialVector<X> 
operator+(const DifferentialVector<X>& x, const R& c)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}

template<class X, class R> inline
DifferentialVector<X> 
operator+(const R& c, const DifferentialVector<X>& x)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}

template<class X, class R> inline
DifferentialVector<X> 
operator-(const DifferentialVector<X>& x, const R& c)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}

template<class X, class R> inline
DifferentialVector<X> 
operator-(const R& c, const DifferentialVector<X>& x)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}

template<class X, class R> inline
DifferentialVector<X> 
operator*(const DifferentialVector<X>& x, const R& c)
{
  DifferentialVector<X> y=x;
  X m(c);
  for(uint i=0; i!=y.result_size(); ++i) {
    y[i]*=m;
  }
  return y;
}

template<class X, class R> inline
DifferentialVector<X> 
operator*(const R& c, const DifferentialVector<X>& x)
{
  DifferentialVector<X> y=x;
  X m(c);
  for(uint i=0; i!=y.result_size(); ++i) {
    y[i]*=m;
  }
  return y;
}

template<class X, class R> inline
DifferentialVector<X> 
operator/(const DifferentialVector<X>& x, const R& c)
{
  DifferentialVector<X> y=x;
  for(uint i=0; i!=y.result_size(); ++i) {
    y[i]/=c;
  }
  return y;
}


} // namespace Ariadne
