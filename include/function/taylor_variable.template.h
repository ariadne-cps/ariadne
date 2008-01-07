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

namespace Ariadne {

namespace {

template<class X0, class X1, class X2> 
void
add_product(Function::TaylorVariable<X0>& x0, const Function::TaylorVariable<X1>& x1, const Function::TaylorVariable<X2>& x2)
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
      uint c=Function::bin(i0,i1);
      x0[i0]+=X0(c)*x1[i1]*x2[i2];
    }
  }
}



template<class X0, class X1, class X2>
void 
compute_composition(Function::TaylorVariable<X0>& z, 
                    const Function::TaylorSeries<X1>& y, 
                    const Function::TaylorVariable<X2>& x)
{
  using namespace Function;
  size_type as=x.argument_size();
  size_type d=z.degree();

  TaylorVariable<X2> w=x;
  w.value()=0;
  TaylorVariable<X0> t(as,d);
  t.value()=y.data()[d]/Numeric::fac<Numeric::Integer>(d);
  for(uint n=1; n<=d; ++n) {
    TaylorVariable<X0> u(as,d);
    add_product(u,t,w);
    t=u+y.data()[d-n]/Numeric::fac<Numeric::Integer>(d-n);
  };
  z=t;
  return;
}

} // namespace




template<class X> 
Function::TaylorVariable<X>::TaylorVariable()
  : _argument_size(1), _degree(0), _data(1u,X(0)) 
{
}

template<class X> 
Function::TaylorVariable<X>::TaylorVariable(size_type a, smoothness_type d)
  : _argument_size(a), _degree(d), _data(compute_polynomial_data_size(a,d),X(0))
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
Function::TaylorVariable<X>&
Function::TaylorVariable<X>::operator+=(const TaylorVariable<X>& x)
{
  ARIADNE_ASSERT(this->argument_size()==x.argument_size());
  ARIADNE_ASSERT(this->degree()==x.degree());
  for(uint i=0; i!=this->_data.size(); ++i) {
    this->_data[i]+=x._data[i];
  }
  //reinterpret_cast<LinearAlgebra::Vector<X>&>(this->_data)
  //  += reinterpret_cast<LinearAlgebra::Vector<X>const&>(x._data);
  return *this;
}

template<class X> 
Function::TaylorVariable<X>&
Function::TaylorVariable<X>::operator*=(const X& x)
{
  LinearAlgebra::Vector<X>& v=reinterpret_cast<LinearAlgebra::Vector<X>&>(this->_data);
  v *= x;
  return *this;
}

template<class X> 
Function::TaylorVariable<X>&
Function::TaylorVariable<X>::operator/=(const X& x)
{
  LinearAlgebra::Vector<X>& v=reinterpret_cast<LinearAlgebra::Vector<X>&>(this->_data);
  v /= x;
  return *this;
}














template<class X> 
Function::TaylorVariable<X> 
Function::compose(const TaylorVariable<X>& y, const TaylorVariable<X>& x)
{
  assert(y.argument_size()==1);
  TaylorSeries<X> t(y.degree(),y.data().begin());
  TaylorVariable<X> z(x.argument_size(),std::min(x.degree(),t.degree()));
  compute_composition(z,t,x);
  return z;
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
  for(size_type n=0; n<y.data().size(); ++n) {
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
  assert(x.argument_size()==y.argument_size());
  TaylorVariable<X> z(x.argument_size(),std::min(x.degree(),y.degree()));
  add_product(z,x,y);
  return z;
}

template<class X> 
Function::TaylorVariable<X> 
Function::div(const TaylorVariable<X>& x, const TaylorVariable<X>& y)
{
  return mul(x,rec(y));
}

template<class X> 
Function::TaylorVariable<X> 
Function::rec(const TaylorVariable<X>& x)
{
  return compose(FunctionSeries<X>::rec(x.degree(),x.value()),x);
}

template<class X, class N> 
Function::TaylorVariable<X> 
Function::pow(const TaylorVariable<X>& x, N k)
{
  return compose(FunctionSeries<X>::pow(x.degree(),x.value(),k),x);
}

template<class X>  
Function::TaylorVariable<X> 
Function::sqrt(const TaylorVariable<X>& x) 
{
  return compose(FunctionSeries<X>::sqrt(x.degree(),x.value()),x);
}

template<class X>  
Function::TaylorVariable<X> 
Function::exp(const TaylorVariable<X>& x) 
{
  return compose(FunctionSeries<X>::exp(x.degree(),x.value()),x);
}

template<class X>  
Function::TaylorVariable<X> 
Function::log(const TaylorVariable<X>& x) 
{
  return compose(FunctionSeries<X>::log(x.degree(),x.value()),x);
}

template<class X> 
Function::TaylorVariable<X> 
Function::sin(const TaylorVariable<X>& x)
{
  return compose(FunctionSeries<X>::sin(x.degree(),x.value()),x);
}

template<class X>  
Function::TaylorVariable<X> 
Function::cos(const TaylorVariable<X>& x) 
{
  return compose(FunctionSeries<X>::cos(x.degree(),x.value()),x);
}

template<class X>  
Function::TaylorVariable<X> 
Function::tan(const TaylorVariable<X>& x) 
{
  return compose(FunctionSeries<X>::tan(x.degree(),x.value()),x);
}

template<class X>  
Function::TaylorVariable<X> 
Function::asin(const TaylorVariable<X>& x) 
{
  return compose(FunctionSeries<X>::asin(x.degree(),x.value()),x);
}

template<class X>  
Function::TaylorVariable<X> 
Function::acos(const TaylorVariable<X>& x) 
{
  return compose(FunctionSeries<X>::acos(x.degree(),x.value()),x);
}

template<class X>  
Function::TaylorVariable<X> 
Function::atan(const TaylorVariable<X>& x) 
{
  return compose(FunctionSeries<X>::atan(x.degree(),x.value()),x);
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



template<class X> template<class XX> 
Function::TaylorVariable<X>::TaylorVariable(size_type a, smoothness_type d, const XX* ptr)
  : _argument_size(a), _degree(d), _data(compute_polynomial_data_size(a,d)) 
{
  for(size_type i=0; i!=this->_data.size(); ++i) {
    this->_data[i]=ptr[i];
  }
}



template<class X> template<class XX> 
Function::TaylorVariable<X>::TaylorVariable(const TaylorSeries<XX>& ts) 
  : _argument_size(1u), _degree(ts.degree()), _data(ts.data())
{
}
  
template<class X> template<class XX> 
Function::TaylorVariable<X>::TaylorVariable(const TaylorVariable<XX>& other) 
  : _argument_size(other._argument_size), _degree(other._degree), _data(other._data) 
{
}
  
template<class X> template<class XX> 
Function::TaylorVariable<X>& 
Function::TaylorVariable<X>::operator=(const TaylorVariable<XX>& other) 
{
  this->_argument_size=other._argument_size;
  this->_degree=other._degree;
  this->_data=other._data;
  return *this;
}


template<class X> template<class XX> 
Function::TaylorVariable<X>& 
Function::TaylorVariable<X>::operator=(const XX& c) 
{
  this->_data[0]=c;
  for(size_type i=1; i!=this->_data.size(); ++i) {
    this->_data[i]=0;
  }
  return *this;
}


template<class X> template<class XX> 
bool 
Function::TaylorVariable<X>::operator==(const TaylorVariable<XX>& other) 
{
  return this->_argument_size==other._argument_size
    && this->_degree==other._degree
    && this->_data==other._data; 
}

template<class X> template<class XX> 
bool 
Function::TaylorVariable<X>::operator!=(const TaylorVariable<XX>& other) 
{
  return !(*this==other); 
}


template<class X> template<class XX> 
Function::TaylorVariable<X> 
Function::TaylorVariable<X>::constant(size_type a, smoothness_type d, const XX& c) 
{
  TaylorVariable<X> result(a,d);
  result._data[0]=c; 
  return result;
}

template<class X> template<class XX> 
Function::TaylorVariable<X> 
Function::TaylorVariable<X>::variable(size_type a, smoothness_type d, const XX& x, size_type i) 
{
  TaylorVariable<X> result(a,d);
  result._data[0]=x; 
  result._data[i+1u]=1; 
  return result;
}





template<class X> 
Function::TaylorVariable<X> 
Function::reduce(const TaylorVariable<X>& x, const size_type& d)
{
  assert(x.degree()>=d);
  Function::TaylorVariable<X> r(x.argument_size(),d);
  for(MultiIndex i(x.argument_size()); i.degree() <= x.degree(); ++i) {
    r[i]=x[i];
  }
}


template<class X> 
Function::TaylorVariable<X> 
Function::derivative(TaylorVariable<X>& x, const size_type& k)
{
  if(x.degree()==0) {
    return TaylorVariable<X>(x.argument_size(),0);
  } 
  TaylorVariable<X> r(x.argument_size(),x.degree()-1);
  MultiIndex e(x.argument_size());
  e.set(k,1);
  for(MultiIndex j(r.argument_size()); j.degree() <= r.degree(); ++j) {
    r[j]=x[j+e];
  }
}







template<class X> 
std::ostream& 
Function::operator<<(std::ostream& os, const TaylorVariable<X>& x) {
  //  return os << "TaylorVariable( argument_size=" << x.argument_size() << ", degree=" << x.degree() << ", data=" << x.data() << ")";
  os << "TaylorVariable(";
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

//  return os << "TaylorVariable( argument_size=" << x.argument_size() << ", degree=" << x.degree() << ", data=" << x.data() << ")";
}




template<class X, class R> 
Function::TaylorVariable<X> 
Function::operator+(const TaylorVariable<X>& x, const R& c)
{
  TaylorVariable<X> r=x; r.data()[0]+=c; return r;
}

template<class X, class R> 
Function::TaylorVariable<X> 
Function::operator+(const R& c, const TaylorVariable<X>& x)
{
  TaylorVariable<X> r=x; r.data()[0]+=c; return r;
}

template<class X, class R> 
Function::TaylorVariable<X> 
Function::operator-(const TaylorVariable<X>& x, const R& c)
{
  TaylorVariable<X> r=x; r.data()[0]-=c; return r;
}

template<class X, class R> 
Function::TaylorVariable<X> 
Function::operator-(const R& c, const TaylorVariable<X>& x)
{

  TaylorVariable<X> r=-x; 
  r.data()[0]+=c; 
  return r;
}

template<class X, class R> 
Function::TaylorVariable<X> 
Function::operator*(const TaylorVariable<X>& x, const R& c)
{
  using LinearAlgebra::Vector;
  TaylorVariable<X> r(x.argument_size(),x.degree()); 
  reinterpret_cast<Vector<X>&>(r.data())=c*reinterpret_cast<const Vector<X>&>(x.data());
  return r;
}

template<class X, class R> 
Function::TaylorVariable<X> 
Function::operator*(const R& c, const TaylorVariable<X>& x)
{
  return x*c;
}

template<class X, class R> 
Function::TaylorVariable<X> 
Function::operator/(const TaylorVariable<X>& x, const R& c)
{
  return x*X(X(1)/c); // Careful! If R=int then 1/c gives integer division
}

template<class X, class R> 
Function::TaylorVariable<X> 
Function::operator/(const R& c, const TaylorVariable<X>& x)
{
  return c*rec(x);
}

template<class X,class R>  
Function::TaylorVariable<X>&
Function::operator/=(const TaylorVariable<X>& x, const R& c)
{
  reinterpret_cast< LinearAlgebra::Vector<X>& >(x.data())/=X(c);
  return x;
}




} //namespace Ariadne
