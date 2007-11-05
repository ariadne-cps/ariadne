/***************************************************************************
 *            first_derivative.inline.h
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


namespace Ariadne {
    
template<class X1, class V1, class X2, class V2> inline
bool
Function::operator==(const FirstDerivative<X1,V1>& x1, const FirstDerivative<X2,V2>& x2) 
{
  return x1.value()==x2.value() && x1.derivative()==x2.derivative();
}

template<class X1, class V1, class X2, class V2> inline
bool
Function::operator!=(const FirstDerivative<X1,V1>& x1, const FirstDerivative<X2,V2>& x2) 
{
  return !(x1==x2);
}

template<class X, class V> inline
std::ostream& 
Function::operator<<(std::ostream& os, const FirstDerivative<X,V>& x) 
{
  return os << "("<<x.value()<<","<<x.derivative()<<")";
}



template<class C, class X, class V> inline 
Function::FirstDerivative<X,V>
Function::operator+(const C& c, const FirstDerivative<X,V>& x) 
{
  return FirstDerivative<X,V>(c+x.value(),x.derivative());
}

template<class C, class X, class V> inline 
Function::FirstDerivative<X,V>
Function::operator+(const FirstDerivative<X,V>& x, const C& c) 
{
  return FirstDerivative<X,V>(x.value()+c,x.derivative());
}

template<class C, class X, class V> inline 
Function::FirstDerivative<X,V>
Function::operator-(const C& c, const FirstDerivative<X,V>& x) 
{
  return FirstDerivative<X,V>(c-x.value(),-x.derivative());
}

template<class C, class X, class V> inline 
Function::FirstDerivative<X,V>
Function::operator-(const FirstDerivative<X,V>& x, const C& c) 
{
  return FirstDerivative<X,V>(x.value()-c,x.derivative());
}

template<class C, class X, class V> inline 
Function::FirstDerivative<X,V>
Function::operator*(const C& c, const FirstDerivative<X,V>& x) 
{
  return FirstDerivative<X,V>(c*x.value(),c*x.derivative());
}

template<class C, class X, class V> inline 
Function::FirstDerivative<X,V>
Function::operator*(const FirstDerivative<X,V>& x, const C& c) 
{
  return FirstDerivative<X,V>(c*x.value(),c*x.derivative());
}

template<class C, class X, class V> inline 
Function::FirstDerivative<X,V>
Function::operator/(const C& c, const FirstDerivative<X,V>& x) 
{
  // Use this form to get right dimension of constant.
  return (c+0*x)/x;
}

template<class C, class X, class V> inline 
Function::FirstDerivative<X,V>
Function::operator/(const FirstDerivative<X,V>& x, const C& c) 
{
  return FirstDerivative<X,V>(x.value()/c,x.derivative()/c);
}



template<class X, class V> inline 
Function::FirstDerivative<X,V> 
Function::min(const FirstDerivative<X,V>& x1,const FirstDerivative<X,V>& x2) 
{
  if(x1.value()==x2.value()) {
    ARIADNE_THROW(std::runtime_error,"min(Differntial x1, FirstDerivative x2)","x1.value()==x2.value()");
  }
  return x1.value()<x2.value() ? x1 : x2;
}

template<class X, class V> inline 
Function::FirstDerivative<X,V> 
Function::max(const FirstDerivative<X,V>& x1,const FirstDerivative<X,V>& x2) 
{
  if(x1.value()==x2.value()) { 
    ARIADNE_THROW(std::runtime_error,"max(Differntial x1, FirstDerivative x2)","x1.value()==x2.value()"); 
  }
  return x1.value()>x2.value() ? x1 : x2;
}

template<class X, class V> inline 
Function::FirstDerivative<X,V> 
Function::abs(const FirstDerivative<X,V>& x) 
{
  if(x.value()==0) { ARIADNE_THROW(std::runtime_error,"abs(Differntial x)","x.value()==0"); }
  if(x.value()>0) { return x; } else { return -x; }
}

template<class X, class V> inline 
Function::FirstDerivative<X,V> 
Function::pos(const FirstDerivative<X,V>& x) 
{
  return FirstDerivative<X,V>(x.value(),x.derivative());
}

template<class X, class V> inline 
Function::FirstDerivative<X,V> 
Function::neg(const FirstDerivative<X,V>& x) 
{
  return FirstDerivative<X,V>(-x.value(),-x.derivative());
}

template<class X, class V> inline 
Function::FirstDerivative<X,V> 
Function::add(const FirstDerivative<X,V>& x1, const FirstDerivative<X,V>& x2) 
{
  return FirstDerivative<X,V>(x1.value()+x2.value(),x1.derivative()+x2.derivative());
}

template<class X, class V> inline 
Function::FirstDerivative<X,V> 
Function::sub(const FirstDerivative<X,V>& x1, const FirstDerivative<X,V>& x2) 
{
  return FirstDerivative<X,V>(x1.value()-x2.value(),x1.derivative()-x2.derivative());
}

template<class X, class V> inline 
Function::FirstDerivative<X,V> 
Function::mul(const FirstDerivative<X,V>& x1, const FirstDerivative<X,V>& x2) 
{
  return FirstDerivative<X,V>(x1.value()*x2.value(),x1.derivative()*x2.value()+x1.value()*x2.derivative());
}

template<class X, class V> inline 
Function::FirstDerivative<X,V> 
Function::div(const FirstDerivative<X,V>& x1, const FirstDerivative<X,V>& x2) 
{
  X y=x1.value()/x2.value();
  V v=x1.derivative()-y*x2.derivative();
  return FirstDerivative<X,V>(y,v/x2.value());
}

template<class X, class V> inline 
Function::FirstDerivative<X,V> 
Function::pow(const FirstDerivative<X,V>& x, int n) 
{
  //std::cerr << "pow(FirstDerivative x, int n)"<<std::endl;
  if(n==0) {
    return FirstDerivative<X,V>(Numeric::pow(x.value(),0),static_cast<X>(0)*x.derivative());
  } else {
    X y=Numeric::pow(x.value(),n-1);
    return FirstDerivative<X,V>((x.value()*y),static_cast<X>(n*y)*x.derivative());
  }
}

template<class X, class V> inline 
Function::FirstDerivative<X,V> 
Function::sqrt(const FirstDerivative<X,V>& x) 
{
  X y=Numeric::sqrt(x.value());
  return FirstDerivative<X,V>(y,x.derivative()/(2*y));
}

template<class X, class V> inline 
Function::FirstDerivative<X,V> 
Function::exp(const FirstDerivative<X,V>& x) 
{
  X y=exp(x.value());
  return FirstDerivative<X,V>(y,y*x.derivative());
}

template<class X, class V> inline 
Function::FirstDerivative<X,V> 
Function::log(const FirstDerivative<X,V>& x) 
{
  return FirstDerivative<X,V>(log(x.value()),x.derivative()/x.value());
}

template<class X, class V> inline 
Function::FirstDerivative<X,V> 
Function::sin(const FirstDerivative<X,V>& x) 
{
  return FirstDerivative<X,V>(sin(x.value()),cos(x.value())*x.derivative());
}

template<class X, class V> inline 
Function::FirstDerivative<X,V> 
Function::cos(const FirstDerivative<X,V>& x) 
{
  return FirstDerivative<X,V>(cos(x.value()),-sin(x.value())*x.derivative());
}

template<class X, class V> inline 
Function::FirstDerivative<X,V> 
Function::tan(const FirstDerivative<X,V>& x) 
{
  X y=tan(x.value());
  return FirstDerivative<X,V>(y,(X(1)-y*y)*x.derivative());
}

template<class X, class V> inline 
Function::FirstDerivative<X,V> 
Function::asin(const FirstDerivative<X,V>& x) 
{
  return FirstDerivative<X,V>(asin(x.value()),x.derivative()/sqrt(X(1)-x.value()*x.value()));
}

template<class X, class V> inline 
Function::FirstDerivative<X,V> 
Function::acos(const FirstDerivative<X,V>& x) 
{
  return FirstDerivative<X,V>(acos(x.value()),-x.derivative()/sqrt(X(1)-x.value()*x.value()));
}

template<class X, class V> inline 
Function::FirstDerivative<X,V> 
Function::atan(const FirstDerivative<X,V>& x) 
{
  return FirstDerivative<X,V>(atan(x.value()),x.derivative()/(X(1)+x.value()*x.value()));
}



template<class X, class V> inline 
Function::FirstDerivative<X,V>
Function::operator+(const FirstDerivative<X,V>& x) 
{
  return pos(x);
}

template<class X, class V> inline 
Function::FirstDerivative<X,V>
Function::operator-(const FirstDerivative<X,V>& x) 
{
  return neg(x);
}

template<class X, class V> inline 
Function::FirstDerivative<X,V>
Function::operator+(const FirstDerivative<X,V>& x1, const FirstDerivative<X,V>& x2) 
{
  return add(x1,x2);
}

template<class X, class V> inline 
Function::FirstDerivative<X,V>
Function::operator-(const FirstDerivative<X,V>& x1, const FirstDerivative<X,V>& x2) 
{
  return sub(x1,x2);
}

template<class X, class V> inline 
Function::FirstDerivative<X,V>
Function::operator*(const FirstDerivative<X,V>& x1, const FirstDerivative<X,V>& x2) 
{
  return mul(x1,x2);
}

template<class X, class V> inline 
Function::FirstDerivative<X,V>
Function::operator/(const FirstDerivative<X,V>& x1, const FirstDerivative<X,V>& x2) 
{
  return div(x1,x2);
}


}
