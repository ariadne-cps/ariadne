/***************************************************************************
 *            affine_variable.inline.h
 *
 *  Copyright  2007  Pieter Collins
 *  Pieter.Collins@cwi.nl
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
 
#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

namespace Ariadne {

template<class X> inline
Function::AffineVariable<X>::~AffineVariable() 
{
  ARIADNE_LOG(3,"AffineVariable<X>::~AffineVariable()\n");
}


template<class X> inline
Function::AffineVariable<X>::AffineVariable() 
  : _x(0), _dx()
{
}

template<class X> inline
Function::AffineVariable<X>::AffineVariable(const dimension_type& ad) 
  :  _x(0), _dx(ad)
{
  ARIADNE_LOG(3,"AffineVariable<X>::AffineVariable(dimension_type d)\n");
}

template<class X> template<class XX> inline
Function::AffineVariable<X>::AffineVariable(const dimension_type& ad, const XX* ptr) 
  :  _x(*ptr), _dx(ad,ptr+1u)
{
  ARIADNE_LOG(3,"AffineVariable<X>::AffineVariable(dimension_type d)\n");
}

template<class X> template<class XX> inline
Function::AffineVariable<X>::AffineVariable(const uint& ad, const XX* ptr) 
  :  _x(*ptr), _dx(ad,ptr+1u)
{
  ARIADNE_LOG(3,"AffineVariable<X>::AffineVariable(dimension_type d)\n");
}


template<class X> inline
Function::AffineVariable<X>::AffineVariable(const X& x, const LinearAlgebra::Covector<X>& dx) 
  : _x(x), _dx(dx)
{
  ARIADNE_LOG(3,"AffineVariable<X>::AffineVariable(XTT)\n");
}

template<class X> inline
Function::AffineVariable<X>::AffineVariable(const AffineVariable<X>& av) 
  : _x(av._x), _dx(av._dx)
{
  ARIADNE_LOG(3,"AffineVariable<X>::AffineVariable(AffineVariable<X> av)\n");
}


template<class X> inline
Function::AffineVariable<X>& 
Function::AffineVariable<X>::operator=(const AffineVariable<X>& av) 
{
  ARIADNE_LOG(3,"AffineVariable<X>::operator=(AffineVariable<X> av)\n");
  if(this!=&av) {
    this->_x=av._x;
    this->_dx=av._dx;
  }
  return *this;
}

template<class X> template<class XX> inline
Function::AffineVariable<X>& 
Function::AffineVariable<X>::operator=(const XX& c) 
{
  ARIADNE_LOG(3,"AffineVariable<X>::operator=(XX c)\n");
  this->_x=c;
  for(size_type j=0; j!=this->_dx.size(); ++j) {
    this->_dx[j]=0;
  }
  return *this;
}


template<class X> inline
size_type
Function::AffineVariable<X>::argument_size() const
{
  return this->_dx.size();
}


template<class X> inline
Function::AffineVariable<X>
Function::AffineVariable<X>::constant(uint as, const X& c) 
{
  AffineVariable<X> result;
  result._x=c;
  result._dx=LinearAlgebra::Covector<X>(as);
  return result;
}

template<class X> inline
Function::AffineVariable<X>
Function::AffineVariable<X>::variable(uint as, uint i, const X& c) 
{
  AffineVariable<X> result;
  result._x=c;
  result._dx=LinearAlgebra::Covector<X>(as);
  result._dx[i]=1;
  return result;
}








template<class X> inline
smoothness_type 
Function::AffineVariable<X>::degree() const
{
  return 1u;
}

template<class X> template<class XT> inline
void 
Function::AffineVariable<X>::set(const XT& x) 
{
  this->_x=x;
}


template<class X> template<class XT> inline
void 
Function::AffineVariable<X>::set(dimension_type j, const XT& x) 
{
  assert(j<this->_dx.size());
  this->_dx[j]=x;
}




template<class X11, class X22> inline
bool
Function::operator==(const AffineVariable<X11>& x1, const AffineVariable<X22>& x2) 
{
  return x1.value()==x2.value() && x1.derivative()==x2.derivative();
}

template<class X11, class X22> inline
bool
Function::operator!=(const AffineVariable<X11>& x1, const AffineVariable<X22>& x2) 
{
  return !(x1==x2);
}



template<class C, class X> inline 
Function::AffineVariable<X>
Function::operator+(const C& c, const AffineVariable<X>& x) 
{
  return AffineVariable<X>(c+x.value(),x.derivative());
}

template<class C, class X> inline 
Function::AffineVariable<X>
Function::operator+(const AffineVariable<X>& x, const C& c) 
{
  return AffineVariable<X>(x.value()+c,x.derivative());
}

template<class C, class X> inline 
Function::AffineVariable<X>
Function::operator-(const C& c, const AffineVariable<X>& x) 
{
  return AffineVariable<X>(c-x.value(),-x.derivative());
}

template<class C, class X> inline 
Function::AffineVariable<X>
Function::operator-(const AffineVariable<X>& x, const C& c) 
{
  return AffineVariable<X>(x.value()-c,x.derivative());
}

template<class C, class X> inline 
Function::AffineVariable<X>
Function::operator*(const C& c, const AffineVariable<X>& x) 
{
  return AffineVariable<X>(c*x.value(),c*x.derivative());
}

template<class C, class X> inline 
Function::AffineVariable<X>
Function::operator*(const AffineVariable<X>& x, const C& c) 
{
  return AffineVariable<X>(c*x.value(),c*x.derivative());
}

template<class C, class X> inline 
Function::AffineVariable<X>
Function::operator/(const C& c, const AffineVariable<X>& x) 
{
  // Use this form to get right dimension of constant.
  return (c+0*x)/x;
}

template<class C, class X> inline 
Function::AffineVariable<X>
Function::operator/(const AffineVariable<X>& x, const C& c) 
{
  return AffineVariable<X>(x.value()/c,x.derivative()/c);
}



template<class X> inline 
Function::AffineVariable<X> 
Function::min(const AffineVariable<X>& x1,const AffineVariable<X>& x2) 
{
  if(x1.value()==x2.value()) {
    ARIADNE_THROW(std::runtime_error,"min(Differntial x1, AffineVariable x2)","x1.value()==x2.value()");
  }
  return x1.value()<x2.value() ? x1 : x2;
}

template<class X> inline 
Function::AffineVariable<X> 
Function::max(const AffineVariable<X>& x1,const AffineVariable<X>& x2) 
{
  if(x1.value()==x2.value()) { 
    ARIADNE_THROW(std::runtime_error,"max(Differntial x1, AffineVariable x2)","x1.value()==x2.value()"); 
  }
  return x1.value()>x2.value() ? x1 : x2;
}

template<class X> inline 
Function::AffineVariable<X> 
Function::abs(const AffineVariable<X>& x) 
{
  if(x.value()==0) { ARIADNE_THROW(std::runtime_error,"abs(Differntial x)","x.value()==0"); }
  if(x.value()>0) { return x; } else { return -x; }
}

template<class X> inline 
Function::AffineVariable<X> 
Function::pos(const AffineVariable<X>& x) 
{
  return AffineVariable<X>(x.value(),x.derivative());
}

template<class X> inline 
Function::AffineVariable<X> 
Function::neg(const AffineVariable<X>& x) 
{
  return AffineVariable<X>(-x.value(),-x.derivative());
}

template<class X> inline 
Function::AffineVariable<X> 
Function::add(const AffineVariable<X>& x1, const AffineVariable<X>& x2) 
{
  return AffineVariable<X>(x1.value()+x2.value(),x1.derivative()+x2.derivative());
}

template<class X> inline 
Function::AffineVariable<X> 
Function::sub(const AffineVariable<X>& x1, const AffineVariable<X>& x2) 
{
  return AffineVariable<X>(x1.value()-x2.value(),x1.derivative()-x2.derivative());
}

template<class X> inline 
Function::AffineVariable<X> 
Function::mul(const AffineVariable<X>& x1, const AffineVariable<X>& x2) 
{
  return AffineVariable<X>(x1.value()*x2.value(),x1.derivative()*x2.value()+x1.value()*x2.derivative());
}

template<class X> inline 
Function::AffineVariable<X> 
Function::div(const AffineVariable<X>& x1, const AffineVariable<X>& x2) 
{
  X y=x1.value()/x2.value();
  LinearAlgebra::Covector<X> v=x1.derivative()-y*x2.derivative();
  return AffineVariable<X>(y,v/x2.value());
}

template<class X> inline 
Function::AffineVariable<X> 
Function::pow(const AffineVariable<X>& x, int n) 
{
  //std::cerr << "pow(AffineVariable x, int n)"<<std::endl;
  if(n==0) {
    return AffineVariable<X>(Numeric::pow(x.value(),0),static_cast<X>(0)*x.derivative());
  } else {
    X y=Numeric::pow(x.value(),n-1);
    return AffineVariable<X>((x.value()*y),static_cast<X>(n*y)*x.derivative());
  }
}

template<class X> inline 
Function::AffineVariable<X> 
Function::sqrt(const AffineVariable<X>& x) 
{
  X y=Numeric::sqrt(x.value());
  return AffineVariable<X>(y,x.derivative()/(2*y));
}

template<class X> inline 
Function::AffineVariable<X> 
Function::exp(const AffineVariable<X>& x) 
{
  X y=exp(x.value());
  return AffineVariable<X>(y,y*x.derivative());
}

template<class X> inline 
Function::AffineVariable<X> 
Function::log(const AffineVariable<X>& x) 
{
  return AffineVariable<X>(log(x.value()),x.derivative()/x.value());
}

template<class X> inline 
Function::AffineVariable<X> 
Function::sin(const AffineVariable<X>& x) 
{
  return AffineVariable<X>(sin(x.value()),cos(x.value())*x.derivative());
}

template<class X> inline 
Function::AffineVariable<X> 
Function::cos(const AffineVariable<X>& x) 
{
  return AffineVariable<X>(cos(x.value()),-sin(x.value())*x.derivative());
}

template<class X> inline 
Function::AffineVariable<X> 
Function::tan(const AffineVariable<X>& x) 
{
  X y=tan(x.value());
  return AffineVariable<X>(y,(X(1)-y*y)*x.derivative());
}

template<class X> inline 
Function::AffineVariable<X> 
Function::asin(const AffineVariable<X>& x) 
{
  return AffineVariable<X>(asin(x.value()),x.derivative()/sqrt(X(1)-x.value()*x.value()));
}

template<class X> inline 
Function::AffineVariable<X> 
Function::acos(const AffineVariable<X>& x) 
{
  return AffineVariable<X>(acos(x.value()),-x.derivative()/sqrt(X(1)-x.value()*x.value()));
}

template<class X> inline 
Function::AffineVariable<X> 
Function::atan(const AffineVariable<X>& x) 
{
  return AffineVariable<X>(atan(x.value()),x.derivative()/(X(1)+x.value()*x.value()));
}



template<class X> inline 
Function::AffineVariable<X>
Function::operator+(const AffineVariable<X>& x) 
{
  return pos(x);
}

template<class X> inline 
Function::AffineVariable<X>
Function::operator-(const AffineVariable<X>& x) 
{
  return neg(x);
}

template<class X> inline 
Function::AffineVariable<X>
Function::operator+(const AffineVariable<X>& x1, const AffineVariable<X>& x2) 
{
  return add(x1,x2);
}

template<class X> inline 
Function::AffineVariable<X>
Function::operator-(const AffineVariable<X>& x1, const AffineVariable<X>& x2) 
{
  return sub(x1,x2);
}

template<class X> inline 
Function::AffineVariable<X>
Function::operator*(const AffineVariable<X>& x1, const AffineVariable<X>& x2) 
{
  return mul(x1,x2);
}

template<class X> inline 
Function::AffineVariable<X>
Function::operator/(const AffineVariable<X>& x1, const AffineVariable<X>& x2) 
{
  return div(x1,x2);
}



}

