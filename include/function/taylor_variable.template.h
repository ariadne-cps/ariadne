/***************************************************************************
 *            taylor_derivative.template.h
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
 
#include "linear_algebra/vector.h"

namespace Ariadne {




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
  TaylorVariable<X> r=-x; r.data()[0]+=c; return r;
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
