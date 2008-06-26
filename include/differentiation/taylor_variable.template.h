/***************************************************************************
 *            taylor_variable.template.h
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
 
#include "linear_algebra/vector.h"
#include "differentiation/multi_index.h"
#include "differentiation/taylor_series.h"

namespace Ariadne {

template<class X> template<class XX> 
TaylorVariable<X>::TaylorVariable(size_type a, smoothness_type d, const XX* ptr)
  : _argument_size(a), _degree(d), _data(compute_polynomial_data_size(a,d)) 
{
  for(size_type i=0; i!=this->_data.size(); ++i) {
    this->_data[i]=ptr[i];
  }
}


template<class X> template<class XX> 
TaylorVariable<X>::TaylorVariable(const TaylorSeries<XX>& ts) 
  : _argument_size(1u), _degree(ts.degree()), _data(ts.data())
{
}
  

template<class X> template<class XX> 
TaylorVariable<X>::TaylorVariable(const TaylorVariable<XX>& other) 
  : _argument_size(other._argument_size), _degree(other._degree), _data(other._data) 
{
}

  
template<class X> template<class XX> 
TaylorVariable<X>& 
TaylorVariable<X>::operator=(const TaylorVariable<XX>& other) 
{
  this->_argument_size=other.argument_size();
  this->_degree=other.degree();
  this->_data=other.data();
  return *this;
}


template<class X> template<class XX> 
TaylorVariable<X>& 
TaylorVariable<X>::operator=(const XX& c) 
{
  this->_data[0]=c;
  for(size_type i=1; i!=this->_data.size(); ++i) {
    this->_data[i]=0;
  }
  return *this;
}


template<class X> template<class XX> 
bool 
TaylorVariable<X>::operator==(const TaylorVariable<XX>& other) const
{
  return this->_argument_size==other.argument_size()
    && this->_degree==other.degree()
    && this->_data==other.data(); 
}


template<class X> template<class XX> 
bool 
TaylorVariable<X>::operator!=(const TaylorVariable<XX>& other) const
{
  return !(*this==other); 
}


template<class X> template<class XX> 
TaylorVariable<X> 
TaylorVariable<X>::constant(size_type a, smoothness_type d, const XX& c) 
{
  TaylorVariable<X> result(a,d);
  result._data[0]=c; 
  return result;
}


template<class X> template<class XX> 
TaylorVariable<X> 
TaylorVariable<X>::variable(size_type a, smoothness_type d, const XX& x, size_type i) 
{
  ARIADNE_ASSERT(d>=1);
  TaylorVariable<X> result(a,d);
  result._data[0]=x; 
  result._data[i+1u]=1; 
  return result;
}


template<class X, class XX> 
XX
evaluate(const TaylorVariable<X>& y, const array<XX>& x)
{
  using namespace std;
  ARIADNE_ASSERT(y.argument_size()==x.size());
  //std::cerr << "y=" << y << std::endl;
  //std::cerr << "x=" << x << std::endl;
  size_type d=y.degree();
  size_type ms=x.size();
  assert(d>=1);

  XX zero = x[0]; zero*=0;
  XX one = zero; one+=1;

  // Use inefficient brute-force approach with lots of storage...
  array< array< XX > > val(ms, array< XX >(d+1));
  for(uint j=0; j!=ms; ++j) {
    val[j][0]=one;
    val[j][1]=x[j];
    for(uint k=2; k<=d; ++k) {
      val[j][k]=val[j][k-1]*x[j];
    }
  }

  XX r(zero);
  for(MultiIndex j(ms); j.degree()<=d; ++j) {
    X sf=fac(j);
    XX t=one;
    for(uint k=0; k!=ms; ++k) {
      t=t*val[k][j[k]];
    }
    t*=X(y[j]/sf);
    r+=t;
  }
  return r;
}





template<class X>  
TaylorVariable<X> 
min(const TaylorVariable<X>& x1, const TaylorVariable<X>& x2) 
{
  if(x1.value()==x2.value()) {
    ARIADNE_THROW(std::runtime_error,"min(TaylorVariable x1, TaylorVariable x2)","x1[0]==x2[0]");
  }
  return x1.value()<x2.value() ? x1 : x2;
}

template<class X>  
TaylorVariable<X> 
max(const TaylorVariable<X>& x1,const TaylorVariable<X>& x2) 
{
  if(x1.value()==x2.value()) { 
    ARIADNE_THROW(std::runtime_error,"max(TaylorVariable x1, TaylorVariable x2)","x1[0]==x2[0]"); 
  }
  return x1.value()>x2.value() ? x1 : x2;
}

template<class X> 
TaylorVariable<X> 
pos(const TaylorVariable<X>& x)
{
  return x;
}

template<class X> 
TaylorVariable<X> 
neg(const TaylorVariable<X>& x)
{
  TaylorVariable<X> y(x.argument_size(),x.degree());
  for(size_type n=0; n<y.data().size(); ++n) {
    y.data()[n] = -x.data()[n];
  }
  return y;
}

template<class X>  
TaylorVariable<X> 
abs(const TaylorVariable<X>& x) 
{
  if(x.value()==0) { 
    ARIADNE_THROW(std::runtime_error,"abs(TaylorVariable x)","x[0]==0"); 
  }
  return x.value()>0 ? pos(x) : neg(x); 
}

template<class X> 
TaylorVariable<X> 
add(const TaylorVariable<X>& x, const TaylorVariable<X>& y)
{
  assert(x.argument_size()==y.argument_size());
  TaylorVariable<X> z(x.argument_size(),std::min(x.degree(),y.degree()));
  for(size_type n=0; n<z.data().size(); ++n) {
    z.data()[n] = x.data()[n]+y.data()[n];
  }
  return z;
}

template<class X> 
TaylorVariable<X> 
sub(const TaylorVariable<X>& x, const TaylorVariable<X>& y)
{
  assert(x.argument_size()==y.argument_size());
  TaylorVariable<X> z(x.argument_size(),std::min(x.degree(),y.degree()));
  for(size_type n=0; n<z.data().size(); ++n) {
    z.data()[n] = x.data()[n]-y.data()[n];
  }
  return z;
}

template<class X> 
TaylorVariable<X> 
mul(const TaylorVariable<X>& x, const TaylorVariable<X>& y)
{
  assert(x.argument_size()==y.argument_size());
  TaylorVariable<X> z(x.argument_size(),std::min(x.degree(),y.degree()));
  acc(z,x,y);
  return z;
}

template<class X> 
TaylorVariable<X> 
div(const TaylorVariable<X>& x, const TaylorVariable<X>& y)
{
  return mul(x,rec(y));
}

template<class X, class N> 
TaylorVariable<X> 
pow(const TaylorVariable<X>& x, N k)
{
  return compose(FunctionSeries<X>::pow(x.degree(),x.value(),k),x);
}

template<class X> 
TaylorVariable<X> 
rec(const TaylorVariable<X>& x)
{
  return compose(FunctionSeries<X>::rec(x.degree(),x.value()),x);
}

template<class X>  
TaylorVariable<X> 
sqrt(const TaylorVariable<X>& x) 
{
  return compose(FunctionSeries<X>::sqrt(x.degree(),x.value()),x);
}

template<class X>  
TaylorVariable<X> 
exp(const TaylorVariable<X>& x) 
{
  return compose(FunctionSeries<X>::exp(x.degree(),x.value()),x);
}

template<class X>  
TaylorVariable<X> 
log(const TaylorVariable<X>& x) 
{
  return compose(FunctionSeries<X>::log(x.degree(),x.value()),x);
}

template<class X> 
TaylorVariable<X> 
sin(const TaylorVariable<X>& x)
{
  return compose(FunctionSeries<X>::sin(x.degree(),x.value()),x);
}

template<class X>  
TaylorVariable<X> 
cos(const TaylorVariable<X>& x) 
{
  return compose(FunctionSeries<X>::cos(x.degree(),x.value()),x);
}

template<class X>  
TaylorVariable<X> 
tan(const TaylorVariable<X>& x) 
{
  return compose(FunctionSeries<X>::tan(x.degree(),x.value()),x);
}

template<class X>  
TaylorVariable<X> 
asin(const TaylorVariable<X>& x) 
{
  return compose(FunctionSeries<X>::asin(x.degree(),x.value()),x);
}

template<class X>  
TaylorVariable<X> 
acos(const TaylorVariable<X>& x) 
{
  return compose(FunctionSeries<X>::acos(x.degree(),x.value()),x);
}

template<class X>  
TaylorVariable<X> 
atan(const TaylorVariable<X>& x) 
{
  return compose(FunctionSeries<X>::atan(x.degree(),x.value()),x);
}


template<class X> 
TaylorVariable<X> 
operator+(const TaylorVariable<X>& x)
{
  return pos(x);
}

template<class X> 
TaylorVariable<X> 
operator-(const TaylorVariable<X>& x)
{
  return neg(x);
}

template<class X> 
TaylorVariable<X> 
operator+(const TaylorVariable<X>& x, const TaylorVariable<X>& y)
{
  return add(x,y);
}

template<class X> 
TaylorVariable<X> 
operator-(const TaylorVariable<X>& x, const TaylorVariable<X>& y)
{
  return sub(x,y);
}

template<class X> 
TaylorVariable<X> 
operator*(const TaylorVariable<X>& x, const TaylorVariable<X>& y)
{
  return mul(x,y);
}

template<class X> 
TaylorVariable<X> 
operator/(const TaylorVariable<X>& x, const TaylorVariable<X>& y)
{
  return div(x,y);
}





template<class X, class R> 
TaylorVariable<X> 
operator+(const TaylorVariable<X>& x, const R& c)
{
  TaylorVariable<X> r=x; r.data()[0]+=c; return r;
}

template<class X, class R> 
TaylorVariable<X> 
operator+(const R& c, const TaylorVariable<X>& x)
{
  TaylorVariable<X> r=x; r.data()[0]+=c; return r;
}

template<class X, class R> 
TaylorVariable<X> 
operator-(const TaylorVariable<X>& x, const R& c)
{
  TaylorVariable<X> r=x; r.data()[0]-=c; return r;
}

template<class X, class R> 
TaylorVariable<X> 
operator-(const R& c, const TaylorVariable<X>& x)
{

  TaylorVariable<X> r=-x; 
  r.data()[0]+=c; 
  return r;
}

template<class X, class R> 
TaylorVariable<X> 
operator*(const TaylorVariable<X>& x, const R& c)
{
  TaylorVariable<X> r(x.argument_size(),x.degree()); 
  reinterpret_cast<Vector<X>&>(r.data())=c*reinterpret_cast<const Vector<X>&>(x.data());
  return r;
}

template<class X, class R> 
TaylorVariable<X> 
operator*(const R& c, const TaylorVariable<X>& x)
{
  return x*c;
}

template<class X, class R> 
TaylorVariable<X> 
operator/(const TaylorVariable<X>& x, const R& c)
{
  return x*X(X(1)/c); // Careful! If R=int then 1/c gives integer division
}

template<class X, class R> 
TaylorVariable<X> 
operator/(const R& c, const TaylorVariable<X>& x)
{
  return c*rec(x);
}



template<class X,class R>  
TaylorVariable<X>&
operator+=(TaylorVariable<X>& x, const R& c)
{
  x.value()+=c; 
  return x;
}

template<class X,class R>  
TaylorVariable<X>&
operator-=(TaylorVariable<X>& x, const R& c)
{
  x.value()-=c; 
  return x;
}

template<class X,class R>  
TaylorVariable<X>&
operator*=(TaylorVariable<X>& x, const R& c)
{
  reinterpret_cast< Vector<X>& >(x.data())*=X(c);
  return x;
}

template<class X,class R>  
TaylorVariable<X>&
operator/=(TaylorVariable<X>& x, const R& c)
{
  reinterpret_cast< Vector<X>& >(x.data())/=X(c);
  return x;
}


} // namespace Ariadne
