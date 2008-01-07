/***************************************************************************
 *            polynomial_variable.template.h
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
#include "function/multi_index.h"
#include "function/polynomial_variable.h"

namespace Ariadne {

namespace {


template<class X0, class X1, class X2> 
void
compute_product(Function::PolynomialVariable<X0>& x0, const Function::PolynomialVariable<X1>& x1, const Function::PolynomialVariable<X2>& x2)
{
  using namespace Function;
  assert(x0.argument_size()==x1.argument_size());
  assert(x0.argument_size()==x2.argument_size());
  for(MultiIndex i1(x1.argument_size()); i1.degree() <= x1.degree(); ++i1) {
    for(MultiIndex i2(x2.argument_size()); i2.degree() <= std::min(x2.degree(),smoothness_type(x0.degree()-i1.degree())); ++i2) {
      MultiIndex i0=i1+i2;
      x0[i0]+=x1[i1]*x2[i2];
    }
  }
}

template<class X0, class X1, class X2>
void 
compute_composition(Function::PolynomialVariable<X0>& z, 
                    const Function::PolynomialVariable<X1>& y, 
                    const Function::PolynomialVariable<X2>& x)
{
  using namespace Function;
  assert(y.argument_size()==1);
  size_type as=x.argument_size();
  size_type d=z.degree();

  using namespace std;
  PolynomialVariable<X2> w=x;
  w.value()=0;
  cerr<<"x="<<x<<"\ny="<<y<<"\nw="<<w<<"\n"<<endl;
  PolynomialVariable<X0> t(as,d);
  PolynomialVariable<X0> u(as,d);
  t.value()=y.data()[d];
  for(uint n=1; n<=d; ++n) {
    cerr<<"t="<<t<<"\n";
    u=t*w;
    cerr<<"u="<<u<<"\n";
    t=u+y.data()[d-n];
  };
  cerr<<"t="<<t<<"\n\n"<<endl;
  z=t;
  return;
}

} // namespace




template<class X> 
Function::PolynomialVariable<X>::PolynomialVariable()
  : _argument_size(0), _degree(0), _data(1u) 
{
}

template<class X> 
Function::PolynomialVariable<X>::PolynomialVariable(size_type a, smoothness_type d)
  : _argument_size(a), _degree(d), _data(compute_polynomial_data_size(a,d))
{
}

template<class X>
Function::PolynomialVariable<X>::PolynomialVariable(const PolynomialVariable<X>& other) 
  : _argument_size(other._argument_size),
    _degree(other._degree),
    _data(other._data)
{
}

template<class X>
Function::PolynomialVariable<X>&
Function::PolynomialVariable<X>::operator=(const PolynomialVariable<X>& other) 
{ 
  this->_argument_size=other._argument_size;
  this->_degree=other._degree;
  this->_data=other._data;
  return *this;
}

template<class X>
bool 
Function::PolynomialVariable<X>::operator==(const PolynomialVariable<X>& other) 
{
  return this->_argument_size==other._argument_size
    && this->_degree==other._degree
    && this->_data==other._data; 
}

template<class X> 
bool 
Function::PolynomialVariable<X>::operator!=(const PolynomialVariable<X>& other) 
{
  return !(*this==other); 
}


template<class X> 
size_type 
Function::PolynomialVariable<X>::argument_size() const 
{ 
  return this->_argument_size;
}

template<class X> 
smoothness_type 
Function::PolynomialVariable<X>::degree() const 
{ 
  return this->_degree;
}

template<class X> 
const X&
Function::PolynomialVariable<X>::value() const 
{ 
  return this->_data[0];
}

template<class X> 
X&
Function::PolynomialVariable<X>::value()  
{ 
  return this->_data[0];
}

template<class X> 
array<X>& 
Function::PolynomialVariable<X>::data()
{
  return this->_data; 
}

template<class X> 
const array<X>& 
Function::PolynomialVariable<X>::data() const 
{
  return this->_data; 
}


template<class X> 
X& 
Function::PolynomialVariable<X>::operator[](const MultiIndex& a) 
{ 
  return this->_data[a.position()]; 
}

template<class X> 
const X& 
Function::PolynomialVariable<X>::operator[](const MultiIndex& a) const 
{ 
  return this->_data[a.position()]; 
}






template<class X> 
Function::PolynomialVariable<X> 
Function::PolynomialVariable<X>::rec(smoothness_type d, const X& c) 
{
  PolynomialVariable<X> y(1u,d);
  X mr = -1/c;
  for(size_type i=0; i<=y.degree(); ++i) {
    y.data()[i]=-Numeric::pow(mr,i+1u);
  }
  
  return y;
}


template<class X> 
Function::PolynomialVariable<X> 
Function::PolynomialVariable<X>::pow(smoothness_type d, const X& c, const uint& k)
{
  size_type n=k;
  PolynomialVariable<X> y(1u,d);
  for(size_type i=0; i<=std::min(size_type(d),n); ++i) {
    int j=n-i;
    y.data()[i]=Numeric::bin<int>(n,i)*Numeric::pow(c,j);
  }
  //std::cout << "pow("<<d<<","<<c<<","<<k<<")="<<y<<std::endl;
  return y;
}

template<class X> 
Function::PolynomialVariable<X> 
Function::PolynomialVariable<X>::sqrt(smoothness_type d, const X& c)
{
  PolynomialVariable<X> y(1u,d);
  y.data()[0]=Numeric::sqrt(c);
  X mhr=(-0.5)/c;
  for(size_type i=1; i<=y.degree(); ++i) {
    y.data()[i]=(2*int(i)-3)*mhr*y.data()[i-1]/Numeric::fac< Numeric::Integer>(i);
  }
  return y;
}

template<class X> 
Function::PolynomialVariable<X> 
Function::PolynomialVariable<X>::exp(smoothness_type d, const X& c)
{
  PolynomialVariable<X> y(1u,d);
  y.data()[0]=Numeric::exp(c);
  for(size_type i=1; i<=y.degree(); ++i) {
    y.data()[i]=y.data()[0]/Numeric::fac< Numeric::Integer>(i);
  }
  return y;
}

template<class X>  
Function::PolynomialVariable<X> 
Function::PolynomialVariable<X>::log(smoothness_type d, const X& c)
{
  PolynomialVariable<X> y(1u,d);
  y.data()[0]=Numeric::log(c);
  X mr=(-1)/c;
  for(size_type i=1; i<=y.degree();++i) {
    y.data()[i]=-Numeric::pow(mr,i)/i;
  }
  return y;
}

template<class X> 
Function::PolynomialVariable<X> 
Function::PolynomialVariable<X>::sin(smoothness_type d, const X& c)
{
  PolynomialVariable<X> y(1u,d);
  y.data()[0]=Numeric::sin(c);
  y.data()[1]=Numeric::cos(c);
  for(size_type i=2; i!=d; ++i) {
    y.data()[i]=-y.data()[i-2]/Numeric::fac< Numeric::Integer>(i);
  }
  return y;
}

template<class X> 
Function::PolynomialVariable<X> 
Function::PolynomialVariable<X>::cos(smoothness_type d, const X& c)
{
  PolynomialVariable<X> y(1u,d);
  y.data()[0]=Numeric::cos(c);
  y.data()[1]=-Numeric::sin(c);
  for(size_type i=2; i!=d; ++i) {
    y.data()[i]=-y.data()[i-2]/Numeric::fac< Numeric::Integer>(i);
  }
  return y;
}





template<class X> 
Function::PolynomialVariable<X> 
Function::PolynomialVariable<X>::tan(smoothness_type d, const X& c)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}

template<class X> 
Function::PolynomialVariable<X> 
Function::PolynomialVariable<X>::asin(smoothness_type d, const X& c)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}


template<class X> 
Function::PolynomialVariable<X> 
Function::PolynomialVariable<X>::acos(smoothness_type d, const X& c)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}

template<class X> 
Function::PolynomialVariable<X> 
Function::PolynomialVariable<X>::atan(smoothness_type d, const X& c)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}













template<class X> 
Function::PolynomialVariable<X> 
Function::compose(const PolynomialVariable<X>& y, const PolynomialVariable<X>& x)
{
  assert(y.argument_size()==1);
  PolynomialVariable<X> z(x.argument_size(),std::min(x.degree(),y.degree()));
  compute_composition(z,y,x);
  return z;
}


template<class X>  
Function::PolynomialVariable<X> 
Function::min(const PolynomialVariable<X>& x1, const PolynomialVariable<X>& x2) 
{
  if(x1.value()==x2.value()) {
    ARIADNE_THROW(std::runtime_error,"min(PolynomialVariable x1, PolynomialVariable x2)","x1[0]==x2[0]");
  }
  return x1.value()<x2.value() ? x1 : x2;
}

template<class X>  
Function::PolynomialVariable<X> 
Function::max(const PolynomialVariable<X>& x1,const PolynomialVariable<X>& x2) 
{
  if(x1.value()==x2.value()) { 
    ARIADNE_THROW(std::runtime_error,"max(PolynomialVariable x1, PolynomialVariable x2)","x1[0]==x2[0]"); 
  }
  return x1.value()>x2.value() ? x1 : x2;
}

template<class X> 
Function::PolynomialVariable<X> 
Function::pos(const PolynomialVariable<X>& x)
{
  return x;
}

template<class X> 
Function::PolynomialVariable<X> 
Function::neg(const PolynomialVariable<X>& x)
{
  PolynomialVariable<X> y(x.argument_size(),x.degree());
  for(size_type n=0; n<=y.data().size(); ++n) {
    y.data()[n] = -x.data()[n];
  }
  return y;
}

template<class X>  
Function::PolynomialVariable<X> 
Function::abs(const PolynomialVariable<X>& x) 
{
  if(x.value()==0) { 
    ARIADNE_THROW(std::runtime_error,"abs(PolynomialVariable x)","x[0]==0"); 
  }
  return x.value()>0 ? pos(x) : neg(x); 
}

template<class X> 
Function::PolynomialVariable<X> 
Function::rec(const PolynomialVariable<X>& x)
{
  return compose(PolynomialVariable<X>::rec(x.degree(),x.value()),x);
}

template<class X> 
Function::PolynomialVariable<X> 
Function::add(const PolynomialVariable<X>& x, const PolynomialVariable<X>& y)
{
  assert(x.argument_size()==y.argument_size());
  PolynomialVariable<X> z(x.argument_size(),std::min(x.degree(),y.degree()));
  for(size_type n=0; n<z.data().size(); ++n) {
    z.data()[n] = x.data()[n]+y.data()[n];
  }
  return z;
}

template<class X> 
Function::PolynomialVariable<X> 
Function::sub(const PolynomialVariable<X>& x, const PolynomialVariable<X>& y)
{
  assert(x.argument_size()==y.argument_size());
  PolynomialVariable<X> z(x.argument_size(),std::min(x.degree(),y.degree()));
  for(size_type n=0; n<z.data().size(); ++n) {
    z.data()[n] = x.data()[n]-y.data()[n];
  }
  return z;
}

template<class X> 
Function::PolynomialVariable<X> 
Function::mul(const PolynomialVariable<X>& x, const PolynomialVariable<X>& y)
{
  PolynomialVariable<X> z(x.argument_size(),std::min(x.degree(),y.degree()));
  compute_product(z,x,y);
  return z;
}

template<class X> 
Function::PolynomialVariable<X> 
Function::div(const PolynomialVariable<X>& x, const PolynomialVariable<X>& y)
{
  return mul(x,rec(y));
}

template<class X, class N> 
Function::PolynomialVariable<X> 
Function::pow(const PolynomialVariable<X>& x, N k)
{
  return compose(PolynomialVariable<X>::pow(x.degree(),x.value(),k),x);
}

template<class X>  
Function::PolynomialVariable<X> 
Function::sqrt(const PolynomialVariable<X>& x) 
{
  return compose(PolynomialVariable<X>::sqrt(x.degree(),x.value()),x);
}

template<class X>  
Function::PolynomialVariable<X> 
Function::exp(const PolynomialVariable<X>& x) 
{
  return compose(PolynomialVariable<X>::exp(x.degree(),x.value()),x);
}

template<class X>  
Function::PolynomialVariable<X> 
Function::log(const PolynomialVariable<X>& x) 
{
  return compose(PolynomialVariable<X>::log(x.degree(),x.value()),x);
}

template<class X> 
Function::PolynomialVariable<X> 
Function::sin(const PolynomialVariable<X>& x)
{
  return compose(PolynomialVariable<X>::sin(x.degree(),x.value()),x);
}

template<class X>  
Function::PolynomialVariable<X> 
Function::cos(const PolynomialVariable<X>& x) 
{
  return compose(PolynomialVariable<X>::cos(x.degree(),x.value()),x);
}

template<class X>  
Function::PolynomialVariable<X> 
Function::tan(const PolynomialVariable<X>& x) 
{
  return compose(PolynomialVariable<X>::tan(x.degree(),x.value()),x);
}

template<class X>  
Function::PolynomialVariable<X> 
Function::asin(const PolynomialVariable<X>& x) 
{
  return compose(PolynomialVariable<X>::asin(x.degree(),x.value()),x);
}

template<class X>  
Function::PolynomialVariable<X> 
Function::acos(const PolynomialVariable<X>& x) 
{
  return compose(PolynomialVariable<X>::acos(x.degree(),x.value()),x);
}

template<class X>  
Function::PolynomialVariable<X> 
Function::atan(const PolynomialVariable<X>& x) 
{
  return compose(PolynomialVariable<X>::atan(x.degree(),x.value()),x);
}


template<class X> 
Function::PolynomialVariable<X> 
Function::operator-(const PolynomialVariable<X>& x)
{
  return neg(x);
}

template<class X> 
Function::PolynomialVariable<X> 
Function::operator+(const PolynomialVariable<X>& x, const PolynomialVariable<X>& y)
{
  return add(x,y);
}

template<class X> 
Function::PolynomialVariable<X> 
Function::operator-(const PolynomialVariable<X>& x, const PolynomialVariable<X>& y)
{
  return sub(x,y);
}

template<class X> 
Function::PolynomialVariable<X> 
Function::operator*(const PolynomialVariable<X>& x, const PolynomialVariable<X>& y)
{
  return mul(x,y);
}

template<class X> 
Function::PolynomialVariable<X> 
Function::operator/(const PolynomialVariable<X>& x, const PolynomialVariable<X>& y)
{
  return div(x,y);
}









template<class X> 
std::ostream& 
Function::operator<<(std::ostream& os, const PolynomialVariable<X>& x) {
  //  return os << "PolynomialVariable( argument_size=" << x.argument_size() << ", degree=" << x.degree() << ", data=" << x.data() << ")";
  os << "PolynomialVariable(";
  size_type degree=0;
  for(MultiIndex i(x.argument_size()); i.degree()<=x.degree(); ++i) {
    if(i.degree()==0) {
      os << '[';
    } else if(i.degree()==degree) {
      os << ',';
    } else {
      degree=i.degree();
      os << ';';
    }
    os << x[i];
  }
  os << ']';
  os << ")";
  return os;
}


}
