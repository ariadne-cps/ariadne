/***************************************************************************
 *            taylor_variable.code.h
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
#include "taylor_series.h"

namespace Ariadne {

namespace {

template<class X0, class X1, class X2> 
void
compute_product(Function::TaylorVariable<X0>& x0, const Function::TaylorVariable<X1>& x1, const Function::TaylorVariable<X2>& x2)
{
  using namespace Function;
  assert(x0.argument_size()==x1.argument_size());
  assert(x0.argument_size()==x2.argument_size());
  for(MultiIndex i1(x1.argument_size()); i1.degree() <= x1.degree(); ++i1) {
    for(MultiIndex i2(x2.argument_size()); i2.degree() <= std::min(x2.degree(),smoothness_type(x0.degree()-i1.degree())); ++i2) {
      MultiIndex i0=i1+i2;
      //std::cout << "i0=" << i0 << ", i1=" << i1 << ", i2=" << i2 << std::endl;
      // FIXME: Use Integer
      //Numeric::Integer c=i0.factorial()/(i1.factorial()*i2.factorial());
      unsigned int c=choose(i0,i1);
      x0[i0]+=X0(c)*x1[i1]*x2[i2];
    }
  }
}

template<class X0, class X1, class X2>
void 
compute_composition(Function::TaylorVariable<X0>& z, const Function::TaylorSeries<X1>& y, const Function::TaylorVariable<X2>& x)
{
  using namespace Function;
  //std::cerr << "y=" << y << std::endl;
  //std::cerr << "z=" << z << std::endl;
  assert(z.degree()==x.degree());
  assert(z.degree()==y.degree());
  size_type d=z.degree();
  TaylorVariable<X2> w=x;
  w.value()=0;
  //std::cerr << "w=" << w << std::endl;
  TaylorVariable<X0> t=TaylorVariable<X0>::constant(x.argument_size(),0,y[d]);
  //std::cerr << "t[0]=" << t << std::endl;
  for(uint n=1; n<=d; ++n) {
    TaylorVariable<X0> u(x.argument_size(),n);
    compute_product(u,t,w);
    u.value()=y[d-n];
    t=u;
    //std::cerr << "t[" << n << "]=" << t << std::endl;
  }
  z=t;
}


template<class X> 
void
instantiate_taylor_variable()
{
  Function::TaylorSeries<X>* ts = 0;
  Function::TaylorVariable<X>* tv = 0;
  - *tv;
  *tv + *tv;
  *tv - *tv;
  *tv * *tv;
  *tv / *tv;
  min(*tv,*tv);
  max(*tv,*tv);
  abs(*tv);  
  pow(*tv,0);
  pow(*tv,0u);
  sqrt(*tv);
  exp(*tv);
  log(*tv);
  sin(*tv);
  cos(*tv);
  tan(*tv);
  asin(*tv);
  acos(*tv);
  atan(*tv);
  
  compose(*ts,*tv);
}

template<>
void
instantiate_taylor_variable<Numeric::Rational>()
{
  typedef Numeric::Rational X;
  Function::TaylorSeries<X>* ts = 0;
  Function::TaylorVariable<X>* tv = 0;
  - *tv;
  *tv + *tv;
  *tv - *tv;
  *tv * *tv;
  *tv / *tv;
  min(*tv,*tv);
  max(*tv,*tv);
  abs(*tv);
  pow(*tv,0u);
  pow(*tv,0);
  
  compose(*ts,*tv);
}


} // namespace




template<class X> 
Function::TaylorVariable<X>::TaylorVariable()
  : _argument_size(0), _degree(0), _data(1u) 
{
}

template<class X> 
Function::TaylorVariable<X>::TaylorVariable(size_type a, smoothness_type d)
  : _argument_size(a), _degree(d), _data(compute_polynomial_data_size(a,d))
{
}



template<class X> 
size_type 
Function::TaylorVariable<X>::argument_size() const 
{ 
  return this->_argument_size;
}

template<class X> 
smoothness_type 
Function::TaylorVariable<X>::degree() const 
{ 
  return this->_degree;
}

template<class X> 
const X&
Function::TaylorVariable<X>::value() const 
{ 
  return this->_data[0];
}

template<class X> 
X&
Function::TaylorVariable<X>::value()  
{ 
  return this->_data[0];
}

template<class X> 
array<X>& 
Function::TaylorVariable<X>::data()
{
  return this->_data; 
}

template<class X> 
const array<X>& 
Function::TaylorVariable<X>::data() const 
{
  return this->_data; 
}


template<class X> 
X& 
Function::TaylorVariable<X>::operator[](const MultiIndex& a) 
{ 
  return this->_data[a.position()]; 
}

template<class X> 
const X& 
Function::TaylorVariable<X>::operator[](const MultiIndex& a) const 
{ 
  return this->_data[a.position()]; 
}















template<class X> 
Function::TaylorVariable<X> 
Function::compose(const TaylorSeries<X>& y, const TaylorVariable<X>& x)
{
  TaylorVariable<X> z(x.argument_size(),std::min(x.degree(),y.degree()));
  compute_composition(z,y,x);
  return z;
}


template<class X>  
Function::TaylorVariable<X> 
Function::min(const TaylorVariable<X>& x1, const TaylorVariable<X>& x2) 
{
  if(x1.value()==x2.value()) {
    ARIADNE_THROW(std::runtime_error,"min(TaylorVariable x1, TaylorVariable x2)","x1[0]==x2[0]");
  }
  return x1.value()<x2.value() ? x1 : x2;
}

template<class X>  
Function::TaylorVariable<X> 
Function::max(const TaylorVariable<X>& x1,const TaylorVariable<X>& x2) 
{
  if(x1.value()==x2.value()) { 
    ARIADNE_THROW(std::runtime_error,"max(TaylorVariable x1, TaylorVariable x2)","x1[0]==x2[0]"); 
  }
  return x1.value()>x2.value() ? x1 : x2;
}

template<class X> 
Function::TaylorVariable<X> 
Function::pos(const TaylorVariable<X>& x)
{
  return x;
}

template<class X> 
Function::TaylorVariable<X> 
Function::neg(const TaylorVariable<X>& x)
{
  TaylorVariable<X> y(x.argument_size(),x.degree());
  for(size_type n=0; n<=y.data().size(); ++n) {
    y.data()[n] = -x.data()[n];
  }
  return y;
}

template<class X>  
Function::TaylorVariable<X> 
Function::abs(const TaylorVariable<X>& x) 
{
  if(x.value()==0) { 
    ARIADNE_THROW(std::runtime_error,"abs(TaylorVariable x)","x[0]==0"); 
  }
  return x.value()>0 ? pos(x) : neg(x); 
}

template<class X> 
Function::TaylorVariable<X> 
Function::inv(const TaylorVariable<X>& x)
{
  TaylorSeries<X> y(x.degree());
  X mr = X(-1)/x.value();
  for(size_type i=0; i<=y.degree(); ++i) {
    y[i]=(-Numeric::factorial<int>(i))*pow(mr,i+1);
  }
  //std::cerr << y << std::endl;
  return compose(y,x);
}

template<class X> 
Function::TaylorVariable<X> 
Function::add(const TaylorVariable<X>& x, const TaylorVariable<X>& y)
{
  assert(x.argument_size()==y.argument_size());
  TaylorVariable<X> z(x.argument_size(),std::min(x.degree(),y.degree()));
  for(size_type n=0; n<z.data().size(); ++n) {
    z.data()[n] = x.data()[n]+y.data()[n];
  }
  return z;
}

template<class X> 
Function::TaylorVariable<X> 
Function::sub(const TaylorVariable<X>& x, const TaylorVariable<X>& y)
{
  assert(x.argument_size()==y.argument_size());
  TaylorVariable<X> z(x.argument_size(),std::min(x.degree(),y.degree()));
  for(size_type n=0; n<z.data().size(); ++n) {
    z.data()[n] = x.data()[n]-y.data()[n];
  }
  return z;
}

template<class X> 
Function::TaylorVariable<X> 
Function::mul(const TaylorVariable<X>& x, const TaylorVariable<X>& y)
{
  TaylorVariable<X> z(x.argument_size(),std::min(x.degree(),y.degree()));
  compute_product(z,x,y);
  return z;
}

template<class X> 
Function::TaylorVariable<X> 
Function::div(const TaylorVariable<X>& x, const TaylorVariable<X>& y)
{
  return mul(x,inv(y));
}

template<class X, class N> 
Function::TaylorVariable<X> 
Function::pow(const TaylorVariable<X>& x, N k)
{
  return compose(TaylorSeries<X>::pow(x.degree(),x.value(),k),x);
}

template<class X>  
Function::TaylorVariable<X> 
Function::sqrt(const TaylorVariable<X>& x) 
{
  return compose(TaylorSeries<X>::sqrt(x.degree(),x.value()),x);
}

template<class X>  
Function::TaylorVariable<X> 
Function::exp(const TaylorVariable<X>& x) 
{
  return compose(TaylorSeries<X>::exp(x.degree(),x.value()),x);
}

template<class X>  
Function::TaylorVariable<X> 
Function::log(const TaylorVariable<X>& x) 
{
  return compose(TaylorSeries<X>::log(x.degree(),x.value()),x);
}

template<class X> 
Function::TaylorVariable<X> 
Function::sin(const TaylorVariable<X>& x)
{
  return compose(TaylorSeries<X>::sin(x.degree(),x.value()),x);
}

template<class X>  
Function::TaylorVariable<X> 
Function::cos(const TaylorVariable<X>& x) 
{
  return compose(TaylorSeries<X>::cos(x.degree(),x.value()),x);
}

template<class X>  
Function::TaylorVariable<X> 
Function::tan(const TaylorVariable<X>& x) 
{
  return compose(TaylorSeries<X>::tan(x.degree(),x.value()),x);
}

template<class X>  
Function::TaylorVariable<X> 
Function::asin(const TaylorVariable<X>& x) 
{
  return compose(TaylorSeries<X>::asin(x.degree(),x.value()),x);
}

template<class X>  
Function::TaylorVariable<X> 
Function::acos(const TaylorVariable<X>& x) 
{
  return compose(TaylorSeries<X>::acos(x.degree(),x.value()),x);
}

template<class X>  
Function::TaylorVariable<X> 
Function::atan(const TaylorVariable<X>& x) 
{
  return compose(TaylorSeries<X>::atan(x.degree(),x.value()),x);
}


template<class X> 
Function::TaylorVariable<X> 
Function::operator-(const TaylorVariable<X>& x)
{
  return neg(x);
}

template<class X> 
Function::TaylorVariable<X> 
Function::operator+(const TaylorVariable<X>& x, const TaylorVariable<X>& y)
{
  return add(x,y);
}

template<class X> 
Function::TaylorVariable<X> 
Function::operator-(const TaylorVariable<X>& x, const TaylorVariable<X>& y)
{
  return sub(x,y);
}

template<class X> 
Function::TaylorVariable<X> 
Function::operator*(const TaylorVariable<X>& x, const TaylorVariable<X>& y)
{
  return mul(x,y);
}

template<class X> 
Function::TaylorVariable<X> 
Function::operator/(const TaylorVariable<X>& x, const TaylorVariable<X>& y)
{
  return div(x,y);
}



template<class X> 
void
Function::TaylorVariable<X>::instantiate()
{ 
  instantiate_taylor_variable<X>();
}

}
