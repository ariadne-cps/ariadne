/***************************************************************************
 *            power_series.inline.h
 *
 *  Copyright 2007  A Pieter Collins
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
 
#include "linear_algebra/vector.h"
#include "differentiation/function_series.h"

namespace Ariadne {

template<class X> inline
PowerSeries<X>::PowerSeries() 
  : _data(1u)
{
}

template<class X> inline
PowerSeries<X>::PowerSeries(smoothness_type d) 
  : _data(d+1u)
{
}

template<class X> template<class XX> inline
PowerSeries<X>::PowerSeries(smoothness_type d, const XX* ptr) 
  : _data(ptr,ptr+d+1)
{
}

template<class X> template<class XX> inline
PowerSeries<X>::PowerSeries(const PowerSeries<XX>& ts) 
  : _data(ts.data())
{
}

template<class X> template<class XX> inline
PowerSeries<X>&
PowerSeries<X>::operator=(const PowerSeries<XX>& ts) 
{
  this->_data=ts.data();
  return *this;
}

template<class X> template<class XX> inline
PowerSeries<X>&
PowerSeries<X>::operator=(const XX& c) 
{
  this->_data[0]=c;
  return *this;
}


template<class X> inline
PowerSeries<X>
PowerSeries<X>::zero(smoothness_type d)
{
  return PowerSeries<X>(d);
}

template<class X> inline
PowerSeries<X>
PowerSeries<X>::constant(smoothness_type d, const X& c)
{
  PowerSeries<X> result(d);
  result[0]=c;
  return result;
}

template<class X> inline
PowerSeries<X>
PowerSeries<X>::variable(smoothness_type d, const X& c)
{
  PowerSeries<X> result(d);
  result[0]=c;
  if(d>=1) {
    result[1]=1;
  }
  return result;
}

template<class X> inline
smoothness_type 
PowerSeries<X>::degree() const 
{
  return this->_data.size()-1;
}

template<class X> inline 
const X& 
PowerSeries<X>::value() const
{
  return this->_data[0];
}

template<class X> inline 
X& 
PowerSeries<X>::value() 
{
  return this->_data[0];
}

template<class X> inline
const array<X>& 
PowerSeries<X>::data() const
{
  return this->_data;
}

template<class X> inline
array<X>& 
PowerSeries<X>::data() 
{
  return this->_data;
}

template<class X> inline 
const X& 
PowerSeries<X>::operator[](const smoothness_type& j) const
{
  return this->_data[j];
}

template<class X> inline
X& 
PowerSeries<X>::operator[](const smoothness_type& j)
{
  return this->_data[j];
}




template<class X> inline 
PowerSeries<X> 
min(const PowerSeries<X>& x1, const PowerSeries<X>& x2) 
{
  if(x1[0]==x2[0]) {
    ARIADNE_THROW(std::runtime_error,"min(PowerSeries x1, PowerSeries x2)","x1[0]==x2[0]");
  }
  return x1[0]<x2[0] ? x1 : x2;
}

template<class X> inline 
PowerSeries<X> 
max(const PowerSeries<X>& x1,const PowerSeries<X>& x2) 
{
  if(x1[0]==x2[0]) { 
    ARIADNE_THROW(std::runtime_error,"max(PowerSeries x1, PowerSeries x2)","x1[0]==x2[0]"); 
  }
  return x1[0]>x2[0] ? x1 : x2;
}

template<class X> inline
PowerSeries<X> 
pos(const PowerSeries<X>& x)
{
  return x;
}

template<class X> inline
PowerSeries<X> 
neg(const PowerSeries<X>& x)
{
  PowerSeries<X> result(x.degree());
  for(size_type n=0; n<=result.degree(); ++n) {
    result[n] = -x[n];
  }
  return result;
}

template<class X> inline 
PowerSeries<X> 
abs(const PowerSeries<X>& x) 
{
  if(x[0]==0) { 
    ARIADNE_THROW(std::runtime_error,"abs(PowerSeries x)","x[0]==0"); 
  }
  return x[0]>0 ? pos(x) : neg(x); 
}

template<class X> inline
PowerSeries<X> 
rec(const PowerSeries<X>& x)
{
  return compose(FunctionSeries<X>::rec(x.degree(),x.value()),x);
}



template<class X> inline
PowerSeries<X> 
add(const PowerSeries<X>& x, const PowerSeries<X>& y)
{
  PowerSeries<X> result(std::min(x.degree(),y.degree()));
  for(size_type n=0; n<=result.degree(); ++n) {
    result[n] = x[n]+y[n];
  }
  return result;
}

template<class X> inline
PowerSeries<X> 
sub(const PowerSeries<X>& x, const PowerSeries<X>& y)
{
  PowerSeries<X> result(std::min(x.degree(),y.degree()));
  for(size_type n=0; n<=result.degree(); ++n) {
    result[n] = x[n]-y[n];
  }
  return result;
}


template<class X> inline
PowerSeries<X> 
div(const PowerSeries<X>& x, const PowerSeries<X>& y)
{
  return x*rec(y);
}

template<class X> inline
PowerSeries<X> 
pow(const PowerSeries<X>& x, const uint& k)
{
  return compose(FunctionSeries<X>::pow(x.degree(),x.value(),k),x);
}

template<class X> inline
PowerSeries<X> 
pow(const PowerSeries<X>& x, const int& k)
{
  return compose(FunctionSeries<X>::pow(x.degree(),x.value(),uint(k)),x);
}

template<class X>  
PowerSeries<X> 
sqrt(const PowerSeries<X>& x) 
{
  return compose(FunctionSeries<X>::sqrt(x.degree(),x.value()),x);
}

template<class X>  
PowerSeries<X> 
exp(const PowerSeries<X>& x) 
{
  return compose(FunctionSeries<X>::exp(x.degree(),x.value()),x);
}

template<class X>  
PowerSeries<X> 
log(const PowerSeries<X>& x) 
{
  return compose(FunctionSeries<X>::log(x.degree(),x.value()),x);
}

template<class X>  
PowerSeries<X> 
sin(const PowerSeries<X>& x) 
{
  return compose(FunctionSeries<X>::sin(x.degree(),x.value()),x);
}

template<class X>  
PowerSeries<X> 
cos(const PowerSeries<X>& x) 
{
  return compose(FunctionSeries<X>::cos(x.degree(),x.value()),x);
}

template<class X>  
PowerSeries<X> 
tan(const PowerSeries<X>& x) 
{
  return compose(FunctionSeries<X>::tan(x.degree(),x.value()),x);
}

template<class X>  
PowerSeries<X> 
asin(const PowerSeries<X>& x) 
{
  return compose(FunctionSeries<X>::asin(x.degree(),x.value()),x);
}

template<class X>  
PowerSeries<X> 
acos(const PowerSeries<X>& x) 
{
  return compose(FunctionSeries<X>::acos(x.degree(),x.value()),x);
}

template<class X>  
PowerSeries<X> 
atan(const PowerSeries<X>& x) 
{
  return compose(FunctionSeries<X>::atan(x.degree(),x.value()),x);
}


template<class X> inline
PowerSeries<X> 
operator-(const PowerSeries<X>& x)
{
  return neg(x);
}

template<class X> inline
PowerSeries<X> 
operator+(const PowerSeries<X>& x, const PowerSeries<X>& y)
{
  return add(x,y);
}

template<class X> inline
PowerSeries<X> 
operator-(const PowerSeries<X>& x, const PowerSeries<X>& y)
{
  return sub(x,y);
}

template<class X> inline
PowerSeries<X> 
operator*(const PowerSeries<X>& x, const PowerSeries<X>& y)
{
  return mul(x,y);
}

template<class X> inline
PowerSeries<X> 
operator/(const PowerSeries<X>& x, const PowerSeries<X>& y)
{
  return div(x,y);
}


template<class X> inline
PowerSeries<X>& 
operator+=(PowerSeries<X>& x, const PowerSeries<X>& y)
{
  for(size_type n=0; n<=std::min(x.degree(),y.degree()); ++n) {
    x[n] += y[n];
  }
  return x;
}



template<class X> inline
PowerSeries<X> 
operator+(const PowerSeries<X>& x, const X& c)
{
  return add(x,PowerSeries<X>::constant(x.degree(),c));
}

template<class X> inline
PowerSeries<X> 
operator+(const X& c, const PowerSeries<X>& x)
{
  return add(PowerSeries<X>::constant(x.degree(),c),x);
}

template<class X> inline
PowerSeries<X> 
operator-(const PowerSeries<X>& x, const X& c)
{
  return sub(x,PowerSeries<X>::constant(x.degree(),c));
}

template<class X> inline
PowerSeries<X> 
operator-(const X& c, const PowerSeries<X>& x)
{
  return sub(PowerSeries<X>::constant(x.degree(),c),x);
}

template<class X> inline
PowerSeries<X> 
operator*(const PowerSeries<X>& x, const X& c)
{
  return mul(x,PowerSeries<X>::constant(x.degree(),c));
}

template<class X> inline
PowerSeries<X> 
operator*(const X& c, const PowerSeries<X>& x)
{
  return mul(PowerSeries<X>::constant(x.degree(),c),x);
}

template<class X> inline
PowerSeries<X> 
operator/(const PowerSeries<X>& x, const X& c)
{
  return div(x,PowerSeries<X>::constant(x.degree(),c));
}

template<class X> inline
PowerSeries<X> 
operator/(const X& c, const PowerSeries<X>& x)
{
  return div(PowerSeries<X>::constant(x.degree(),c),x);
}

template<class X> inline
PowerSeries<X> 
operator/(const double& c, const PowerSeries<X>& x)
{
  return X(c)/x;
}

template<class X>  
PowerSeries<X>&
operator+=(PowerSeries<X>& x, const X& c)
{
  x[0]+=c;
  return x;
}

template<class X>  
PowerSeries<X>&
operator-=(PowerSeries<X>& x, const X& c)
{
  x[0]-=c;
  return x;
}

template<class X>  
PowerSeries<X>&
operator*=(PowerSeries<X>& x, const X& c)
{
  reinterpret_cast< Vector<X>& >(x.data())*=c;
  return x;
}

template<class X>  
PowerSeries<X>&
operator/=(PowerSeries<X>& x, const X& c)
{
  reinterpret_cast< Vector<X>& >(x.data())/=c;
  return x;
}


template<class X>  
PowerSeries<X>&
operator+=(PowerSeries<X>& x, const double& c)
{
  X& v=x[0];
  v+=c;
  return x;
}

template<class X>  
PowerSeries<X>&
operator-=(PowerSeries<X>& x, const double& c)
{
  x[0]-=c;
  return x;
}

template<class X>  
PowerSeries<X>&
operator*=(PowerSeries<X>& x, const double& c)
{
  reinterpret_cast< Vector<X>& >(x.data())*=c;
  return x;
}

template<class X>  
PowerSeries<X>&
operator/=(PowerSeries<X>& x, const double& c)
{
  reinterpret_cast< Vector<X>& >(x.data())/=c;
  return x;
}


template<class X>  
PowerSeries<X>
operator+(const PowerSeries<X>& x, const double& c)
{
  PowerSeries<X> r(x); r+=c; return r;
}

template<class X>  
PowerSeries<X>
operator-(const PowerSeries<X>& x, const double& c)
{
  PowerSeries<X> r(x); r-=c; return r;
}

template<class X, class R>  
PowerSeries<X>&
operator*=(PowerSeries<X>& x, const R& c)
{
  reinterpret_cast< Vector<X>& >(x.data())*=c;
  return x;
}







} // namespace Ariadne
