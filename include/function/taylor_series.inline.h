/***************************************************************************
 *            taylor_series.inline.h
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
 

namespace Ariadne {

template<class X> inline
Function::TaylorSeries<X>::TaylorSeries() 
  : _data(1u)
{
}

template<class X> inline
Function::TaylorSeries<X>::TaylorSeries(smoothness_type d) 
  : _data(d+1u)
{
}

template<class X> template<class XX> inline
Function::TaylorSeries<X>::TaylorSeries(smoothness_type d, const XX* ptr) 
  : _data(ptr,ptr+d+1)
{
}

template<class X> template<class XX> inline
Function::TaylorSeries<X>::TaylorSeries(const TaylorSeries<XX>& ts) 
  : _data(ts.data())
{
}

template<class X> template<class XX> inline
Function::TaylorSeries<X>&
Function::TaylorSeries<X>::operator=(const TaylorSeries<XX>& ts) 
{
  this->_data=ts.data();
  return *this;
}

template<class X> template<class XX> inline
Function::TaylorSeries<X>&
Function::TaylorSeries<X>::operator=(const XX& c) 
{
  this->_data[0]=c;
  return *this;
}

template<class X> inline
smoothness_type 
Function::TaylorSeries<X>::degree() const 
{
  return this->_data.size()-1;
}

template<class X> inline 
const X& 
Function::TaylorSeries<X>::value() const
{
  return this->_data[0];
}

template<class X> inline 
X& 
Function::TaylorSeries<X>::value() 
{
  return this->_data[0];
}

template<class X> inline
const array<X>& 
Function::TaylorSeries<X>::data() const
{
  return this->_data;
}

template<class X> inline
array<X>& 
Function::TaylorSeries<X>::data() 
{
  return this->_data;
}

template<class X> inline 
const X& 
Function::TaylorSeries<X>::operator[](const smoothness_type& j) const
{
  return this->_data[j];
}

template<class X> inline
X& 
Function::TaylorSeries<X>::operator[](const smoothness_type& j)
{
  return this->_data[j];
}




template<class X> inline 
Function::TaylorSeries<X> 
Function::min(const TaylorSeries<X>& x1, const TaylorSeries<X>& x2) 
{
  if(x1[0]==x2[0]) {
    ARIADNE_THROW(std::runtime_error,"min(TaylorSeries x1, TaylorSeries x2)","x1[0]==x2[0]");
  }
  return x1[0]<x2[0] ? x1 : x2;
}

template<class X> inline 
Function::TaylorSeries<X> 
Function::max(const TaylorSeries<X>& x1,const TaylorSeries<X>& x2) 
{
  if(x1[0]==x2[0]) { 
    ARIADNE_THROW(std::runtime_error,"max(TaylorSeries x1, TaylorSeries x2)","x1[0]==x2[0]"); 
  }
  return x1[0]>x2[0] ? x1 : x2;
}

template<class X> inline
Function::TaylorSeries<X> 
Function::pos(const TaylorSeries<X>& x)
{
  return x;
}

template<class X> inline
Function::TaylorSeries<X> 
Function::neg(const TaylorSeries<X>& x)
{
  TaylorSeries<X> result(x.degree());
  for(size_type n=0; n<=result.degree(); ++n) {
    result[n] = -x[n];
  }
  return result;
}

template<class X> inline 
Function::TaylorSeries<X> 
Function::abs(const TaylorSeries<X>& x) 
{
  if(x[0]==0) { 
    ARIADNE_THROW(std::runtime_error,"abs(TaylorSeries x)","x[0]==0"); 
  }
  return x[0]>0 ? pos(x) : neg(x); 
}

template<class X> inline
Function::TaylorSeries<X> 
Function::rec(const TaylorSeries<X>& x)
{
  return compose(TaylorSeries<X>::rec(x.degree(),x.value()),x);
}



template<class X> inline
Function::TaylorSeries<X> 
Function::add(const TaylorSeries<X>& x, const TaylorSeries<X>& y)
{
  TaylorSeries<X> result(std::min(x.degree(),y.degree()));
  for(size_type n=0; n<=result.degree(); ++n) {
    result[n] = x[n]+y[n];
  }
  return result;
}

template<class X> inline
Function::TaylorSeries<X> 
Function::sub(const TaylorSeries<X>& x, const TaylorSeries<X>& y)
{
  TaylorSeries<X> result(std::min(x.degree(),y.degree()));
  for(size_type n=0; n<=result.degree(); ++n) {
    result[n] = x[n]-y[n];
  }
  return result;
}

template<class X> inline
Function::TaylorSeries<X> 
Function::mul(const TaylorSeries<X>& x, const TaylorSeries<X>& y)
{
  TaylorSeries<X> result(std::min(x.degree(),y.degree()));
  for(size_type n=0; n<=result.degree(); ++n) {
    for(size_type i=0; i<=n; ++i) {
      result[n] += Numeric::bin<int>(n,i)*x[i]*y[n-i];
    }
  }
  return result;
}

template<class X> inline
Function::TaylorSeries<X> 
Function::div(const TaylorSeries<X>& x, const TaylorSeries<X>& y)
{
  return x*rec(y);
}

template<class X> inline
Function::TaylorSeries<X> 
Function::pow(const TaylorSeries<X>& x, const uint& k)
{
  return compose(x,TaylorSeries<X>::pow(x.degree(),x.value(),k));
}

template<class X> inline
Function::TaylorSeries<X> 
Function::pow(const TaylorSeries<X>& x, const int& k)
{
  return compose(x,TaylorSeries<X>::pow(x.degree(),x.value(),uint(k)));
}

template<class X>  
Function::TaylorSeries<X> 
Function::sqrt(const TaylorSeries<X>& x) 
{
  return compose(x,TaylorSeries<X>::sqrt(x.degree(),x.value()));
}

template<class X>  
Function::TaylorSeries<X> 
Function::exp(const TaylorSeries<X>& x) 
{
  return compose(x,TaylorSeries<X>::exp(x.degree(),x.value()));
}

template<class X>  
Function::TaylorSeries<X> 
Function::log(const TaylorSeries<X>& x) 
{
  return compose(x,TaylorSeries<X>::log(x.degree(),x.value()));
}

template<class X>  
Function::TaylorSeries<X> 
Function::sin(const TaylorSeries<X>& x) 
{
  return compose(x,TaylorSeries<X>::sin(x.degree(),x.value()));
}

template<class X>  
Function::TaylorSeries<X> 
Function::cos(const TaylorSeries<X>& x) 
{
  return compose(x,TaylorSeries<X>::cos(x.degree(),x.value()));
}

template<class X>  
Function::TaylorSeries<X> 
Function::tan(const TaylorSeries<X>& x) 
{
  return compose(x,TaylorSeries<X>::tan(x.degree(),x.value()));
}

template<class X>  
Function::TaylorSeries<X> 
Function::asin(const TaylorSeries<X>& x) 
{
  return compose(x,TaylorSeries<X>::asin(x.degree(),x.value()));
}

template<class X>  
Function::TaylorSeries<X> 
Function::acos(const TaylorSeries<X>& x) 
{
  return compose(x,TaylorSeries<X>::acos(x.degree(),x.value()));
}

template<class X>  
Function::TaylorSeries<X> 
Function::atan(const TaylorSeries<X>& x) 
{
  return compose(x,TaylorSeries<X>::atan(x.degree(),x.value()));
}


template<class X> inline
Function::TaylorSeries<X> 
Function::operator-(const TaylorSeries<X>& x)
{
  return neg(x);
}

template<class X> inline
Function::TaylorSeries<X> 
Function::operator+(const TaylorSeries<X>& x, const TaylorSeries<X>& y)
{
  return add(x,y);
}

template<class X> inline
Function::TaylorSeries<X> 
Function::operator-(const TaylorSeries<X>& x, const TaylorSeries<X>& y)
{
  return sub(x,y);
}

template<class X> inline
Function::TaylorSeries<X> 
Function::operator*(const TaylorSeries<X>& x, const TaylorSeries<X>& y)
{
  return mul(x,y);
}

template<class X> inline
Function::TaylorSeries<X> 
Function::operator/(const TaylorSeries<X>& x, const TaylorSeries<X>& y)
{
  return div(x,y);
}

template<class X, class R> inline
Function::TaylorSeries<X> 
Function::operator+(const TaylorSeries<X>& x, const R& c)
{
  return add(x,TaylorSeries<X>::constant(x.degree(),c));
}

template<class X, class R> inline
Function::TaylorSeries<X> 
Function::operator+(const R& c, const TaylorSeries<X>& x)
{
  return add(TaylorSeries<X>::constant(x.degree(),c),x);
}

template<class X, class R> inline
Function::TaylorSeries<X> 
Function::operator-(const TaylorSeries<X>& x, const R& c)
{
  return sub(x,TaylorSeries<X>::constant(x.degree(),c));
}

template<class X, class R> inline
Function::TaylorSeries<X> 
Function::operator-(const R& c, const TaylorSeries<X>& x)
{
  return sub(TaylorSeries<X>::constant(x.degree(),c),x);
}

template<class X, class R> inline
Function::TaylorSeries<X> 
Function::operator*(const TaylorSeries<X>& x, const R& c)
{
  return mul(x,TaylorSeries<X>::constant(x.degree(),c));
}

template<class X, class R> inline
Function::TaylorSeries<X> 
Function::operator*(const R& c, const TaylorSeries<X>& x)
{
  return mul(TaylorSeries<X>::constant(x.degree(),c),x);
}

template<class X, class R> inline
Function::TaylorSeries<X> 
Function::operator/(const TaylorSeries<X>& x, const R& c)
{
  return div(x,TaylorSeries<X>::constant(x.degree(),c));
}

template<class X, class R> inline
Function::TaylorSeries<X> 
Function::operator/(const R& c, const TaylorSeries<X>& x)
{
  return div(TaylorSeries<X>::constant(x.degree(),c),x);
}






inline
Function::TaylorSeries<Numeric::Rational>::TaylorSeries() 
  : _data(1u)
{
}

inline
Function::TaylorSeries<Numeric::Rational>::TaylorSeries(smoothness_type d) 
  : _data(d+1u)
{
}

template<class XX> inline
Function::TaylorSeries<Numeric::Rational>::TaylorSeries(smoothness_type d, const XX* ptr) 
  : _data(ptr,ptr+d+1)
{
}

template<class XX> inline
Function::TaylorSeries<Numeric::Rational>::TaylorSeries(const TaylorSeries<XX>& ts) 
  : _data(ts.data())
{
}

template<class XX> inline
Function::TaylorSeries<Numeric::Rational>&
Function::TaylorSeries<Numeric::Rational>::operator=(const TaylorSeries<XX>& ts) 
{
  this->_data=ts.data();
  return *this;
}

template<class XX> inline
Function::TaylorSeries<Numeric::Rational>&
Function::TaylorSeries<Numeric::Rational>::operator=(const XX& c) 
{
  this->_data[0]=c;
  return *this;
}

inline
smoothness_type 
Function::TaylorSeries<Numeric::Rational>::degree() const 
{
  return this->_data.size()-1;
}

inline 
const Numeric::Rational& 
Function::TaylorSeries<Numeric::Rational>::value() const
{
  return this->_data[0];
}

inline 
Numeric::Rational& 
Function::TaylorSeries<Numeric::Rational>::value() 
{
  return this->_data[0];
}

inline
const array<Numeric::Rational>& 
Function::TaylorSeries<Numeric::Rational>::data() const
{
  return this->_data;
}

inline
array<Numeric::Rational>& 
Function::TaylorSeries<Numeric::Rational>::data() 
{
  return this->_data;
}

inline 
const Numeric::Rational& 
Function::TaylorSeries<Numeric::Rational>::operator[](const smoothness_type& j) const
{
  return this->_data[j];
}

inline
Numeric::Rational& 
Function::TaylorSeries<Numeric::Rational>::operator[](const smoothness_type& j)
{
  return this->_data[j];
}

inline
Function::TaylorSeries<Numeric::Rational>
Function::TaylorSeries<Numeric::Rational>::constant(smoothness_type d, const Numeric::Rational& c)
{
  TaylorSeries<Numeric::Rational> result(d);
  result[0]=c;
  return result;
}

inline
Function::TaylorSeries<Numeric::Rational>
Function::TaylorSeries<Numeric::Rational>::variable(smoothness_type d, const Numeric::Rational& c)
{
  TaylorSeries<Numeric::Rational> result(d);
  result[0]=c;
  result[1]=1;
  return result;
}





} // namespace Ariadne
