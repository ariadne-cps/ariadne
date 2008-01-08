/***************************************************************************
 *            taylor_variable.code.h
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
#include "function/taylor_series.h"

namespace Ariadne {

namespace {

template<class X> 
void
add_product(Function::TaylorVariable<X>& x0, const Function::TaylorVariable<X>& x1, const Function::TaylorVariable<X>& x2)
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
      x0[i0]+=X(c)*x1[i1]*x2[i2];
    }
  }
}



template<class X>
void 
compute_composition(Function::TaylorVariable<X>& z, 
                    const Function::TaylorSeries<X>& y, 
                    const Function::TaylorVariable<X>& x)
{
  using namespace Function;
  size_type as=x.argument_size();
  size_type d=z.degree();

  TaylorVariable<X> w=x;
  w.value()=0;
  TaylorVariable<X> t(as,d);
  t.value()=y.data()[d]/Numeric::fac<Numeric::Integer>(d);
  for(uint n=1; n<=d; ++n) {
    TaylorVariable<X> u(as,d);
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
Function::TaylorVariable<X> 
Function::mul(const TaylorVariable<X>& x, const TaylorVariable<X>& y)
{
  assert(x.argument_size()==y.argument_size());
  TaylorVariable<X> z(x.argument_size(),std::min(x.degree(),y.degree()));
  add_product(z,x,y);
  return z;
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
Function::derivative(const TaylorVariable<X>& x, const size_type& k)
{
  TaylorVariable<X> r(x.argument_size(),x.degree()-1);
  MultiIndex e(x.argument_size());
  e.set(k,1);
  for(MultiIndex j(r.argument_size()); j.degree() <= r.degree(); ++j) {
    r[j]=x[j+e];
  }
  return r;
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


template<class X> 
void
Function::TaylorVariable<X>::instantiate() 
{
  TaylorSeries<X>* ts=0;
  TaylorVariable<X>* tv=0;
  std::ostream* os = 0;

  mul(*tv,*tv);
  compose(*ts,*tv);
  derivative(*tv,0u);
  + *tv;
  - *tv;
  *tv + *tv;
  *tv - *tv;
  *tv / *tv;
  *tv * *tv;
  
  operator<<(*os,*tv);
}


} //namespace Ariadne
