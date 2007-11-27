/***************************************************************************
 *            scalar_derivative.inline.h
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


template<class X> inline
Function::ScalarDerivative<X> 
Function::compose(const ScalarDerivative<X>& y, const ScalarDerivative<X>& x)
{
  ScalarDerivative<X> result(std::min(x.degree(),y.degree()));
  for(size_type n=0; n<=result.degree(); ++n) { result[n]=y[n]; }
  compute_composition(result,x);
  return result;
}


template<class X> inline 
Function::ScalarDerivative<X> 
Function::min(const ScalarDerivative<X>& x1, const ScalarDerivative<X>& x2) 
{
  if(x1[0]==x2[0]) {
    ARIADNE_THROW(std::runtime_error,"min(ScalarDerivative x1, ScalarDerivative x2)","x1[0]==x2[0]");
  }
  return x1[0]<x2[0] ? x1 : x2;
}

template<class X> inline 
Function::ScalarDerivative<X> 
Function::max(const ScalarDerivative<X>& x1,const ScalarDerivative<X>& x2) 
{
  if(x1[0]==x2[0]) { 
    ARIADNE_THROW(std::runtime_error,"max(ScalarDerivative x1, ScalarDerivative x2)","x1[0]==x2[0]"); 
  }
  return x1[0]>x2[0] ? x1 : x2;
}

template<class X> inline
Function::ScalarDerivative<X> 
Function::pos(const ScalarDerivative<X>& x)
{
  return x;
}

template<class X> inline
Function::ScalarDerivative<X> 
Function::neg(const ScalarDerivative<X>& x)
{
  ScalarDerivative<X> result(x.degree());
  for(size_type n=0; n<=result.degree(); ++n) {
    result[n] = -x[n];
  }
  return result;
}

template<class X>  
Function::ScalarDerivative<X> 
Function::abs(const ScalarDerivative<X>& x) 
{
  if(x[0]==0) { 
    ARIADNE_THROW(std::runtime_error,"abs(ScalarDerivative x)","x[0]==0"); 
  }
  return x[0]>0 ? pos(x) : neg(x); 
}

template<class X> inline
Function::ScalarDerivative<X> 
Function::inv(const ScalarDerivative<X>& x)
{
  ScalarDerivative<X> y(x.degree());
  X mr = X(-1)/x[0];
  for(size_type i=0; i<=y.degree(); ++i) {
    y[i]=(-Numeric::factorial<int>(i))*pow(mr,i+1);
  }
  //std::cerr << y << std::endl;
  compute_composition(y,x);
  //std::cerr << y << std::endl;
  return y;
}

template<class X> inline
Function::ScalarDerivative<X> 
Function::add(const ScalarDerivative<X>& x, const ScalarDerivative<X>& y)
{
  ScalarDerivative<X> result(std::min(x.degree(),y.degree()));
  for(size_type n=0; n<=result.degree(); ++n) {
    result[n] = x[n]+y[n];
  }
  return result;
}

template<class X> inline
Function::ScalarDerivative<X> 
Function::sub(const ScalarDerivative<X>& x, const ScalarDerivative<X>& y)
{
  ScalarDerivative<X> result(std::min(x.degree(),y.degree()));
  for(size_type n=0; n<=result.degree(); ++n) {
    result[n] = x[n]-y[n];
  }
  return result;
}

template<class X> inline
Function::ScalarDerivative<X> 
Function::mul(const ScalarDerivative<X>& x, const ScalarDerivative<X>& y)
{
  ScalarDerivative<X> result(std::min(x.degree(),y.degree()));
  for(size_type n=0; n<=result.degree(); ++n) {
    for(size_type i=0; i<=n; ++i) {
      result[n] += Numeric::choose<int>(n,i)*x[i]*y[n-i];
    }
  }
  return result;
}

template<class X> inline
Function::ScalarDerivative<X> 
Function::div(const ScalarDerivative<X>& x, const ScalarDerivative<X>& y)
{
  return x*inv(y);
}

template<class X, class N> inline
Function::ScalarDerivative<X> 
Function::pow(const ScalarDerivative<X>& x, N k)
{
  size_type n=k;
  ScalarDerivative<X> result(x.degree());
  for(size_type i=0; i<=std::min(size_type(result.degree()),n); ++i) {
    int j=n-i;
    result[i]=(Numeric::factorial<int>(n)/Numeric::factorial<int>(j))*pow(x[0],j);
  }
  compute_composition(result,x);
  return result;
}

template<class X>  
Function::ScalarDerivative<X> 
Function::sqrt(const ScalarDerivative<X>& x) 
{
  ScalarDerivative<X> y(x.degree());
  y[0]=sqrt(x[0]);
  X mhr=(-0.5)/x[0];
  for(size_type i=1; i<=y.degree(); ++i) {
    y[i]=(2*int(i)-3)*mhr*y[i-1];
  }
  compute_composition(y,x);
  return y;
}

template<class X>  
Function::ScalarDerivative<X> 
Function::exp(const ScalarDerivative<X>& x) 
{
  ScalarDerivative<X> y(x.degree());
  y[0]=exp(x[0]);
  for(size_type i=1; i<=y.degree(); ++i) {
    y[i]=y[0];
  }
  compute_composition(y,x);
  return y;
}

template<class X>  
Function::ScalarDerivative<X> 
Function::log(const ScalarDerivative<X>& x) 
{
  ScalarDerivative<X> y(x.degree());
  y[0]=log(x[0]);
  X mr=(-1)/x[0];
  for(size_type i=1; i!=y.degree();++i) {
    y[i]=(-Numeric::factorial<int>(i-1))*pow(x[0],i);
  }
  compute_composition(y,x);
  return y;
}

template<class X>  
Function::ScalarDerivative<X> 
Function::sin(const ScalarDerivative<X>& x) 
{
  size_type d=x.degree();
  ScalarDerivative<X> y(d);
  y[0]=sin(x[0]);
  y[1]=cos(x[0]);
  for(size_type i=2; i!=d; ++i) {
    y[i]=-y[i-2];
  }
  compute_composition(y,x);
  return y;
}

template<class X>  
Function::ScalarDerivative<X> 
Function::cos(const ScalarDerivative<X>& x) 
{
  size_type d=x.degree();
  ScalarDerivative<X> y(d);
  y[0]=cos(x[0]);
  y[1]=-sin(x[0]);
  for(size_type i=2; i!=d; ++i) {
    y[i]=-y[i-2];
  }
  compute_composition(y,x);
  return y;
}

template<class X>  
Function::ScalarDerivative<X> 
Function::tan(const ScalarDerivative<X>& x) 
{
  return sin(x)/cos(x);
}

template<class X>  
Function::ScalarDerivative<X> 
Function::asin(const ScalarDerivative<X>& x) 
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}

template<class X>  
Function::ScalarDerivative<X> 
Function::acos(const ScalarDerivative<X>& x) 
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}

template<class X>  
Function::ScalarDerivative<X> 
Function::atan(const ScalarDerivative<X>& x) 
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}


template<class X> inline
Function::ScalarDerivative<X> 
Function::operator-(const ScalarDerivative<X>& x)
{
  return neg(x);
}

template<class X> inline
Function::ScalarDerivative<X> 
Function::operator+(const ScalarDerivative<X>& x, const ScalarDerivative<X>& y)
{
  return add(x,y);
}

template<class X> inline
Function::ScalarDerivative<X> 
Function::operator-(const ScalarDerivative<X>& x, const ScalarDerivative<X>& y)
{
  return sub(x,y);
}

template<class X> inline
Function::ScalarDerivative<X> 
Function::operator*(const ScalarDerivative<X>& x, const ScalarDerivative<X>& y)
{
  return mul(x,y);
}

template<class X> inline
Function::ScalarDerivative<X> 
Function::operator/(const ScalarDerivative<X>& x, const ScalarDerivative<X>& y)
{
  return div(x,y);
}

template<class X, class R> inline
Function::ScalarDerivative<X> 
Function::operator+(const ScalarDerivative<X>& x, const R& c)
{
  return add(x,ScalarDerivative<X>::constant(x.degree(),c));
}

template<class X, class R> inline
Function::ScalarDerivative<X> 
Function::operator+(const R& c, const ScalarDerivative<X>& x)
{
  return add(ScalarDerivative<X>::constant(x.degree(),c),x);
}

template<class X, class R> inline
Function::ScalarDerivative<X> 
Function::operator-(const ScalarDerivative<X>& x, const R& c)
{
  return sub(x,ScalarDerivative<X>::constant(x.degree(),c));
}

template<class X, class R> inline
Function::ScalarDerivative<X> 
Function::operator-(const R& c, const ScalarDerivative<X>& x)
{
  return sub(ScalarDerivative<X>::constant(x.degree(),c),x);
}

template<class X, class R> inline
Function::ScalarDerivative<X> 
Function::operator*(const ScalarDerivative<X>& x, const R& c)
{
  return mul(x,ScalarDerivative<X>::constant(x.degree(),c));
}

template<class X, class R> inline
Function::ScalarDerivative<X> 
Function::operator*(const R& c, const ScalarDerivative<X>& x)
{
  return mul(ScalarDerivative<X>::constant(x.degree(),c),x);
}

template<class X, class R> inline
Function::ScalarDerivative<X> 
Function::operator/(const ScalarDerivative<X>& x, const R& c)
{
  return div(x,ScalarDerivative<X>::constant(x.degree(),c));
}

template<class X, class R> inline
Function::ScalarDerivative<X> 
Function::operator/(const R& c, const ScalarDerivative<X>& x)
{
  return div(ScalarDerivative<X>::constant(x.degree(),c),x);
}

template<class X> inline
std::ostream& 
Function::operator<<(std::ostream& os, const ScalarDerivative<X>& x) {
  os << "ScalarDerivative";
  for(size_type i=0; i<=x.degree(); ++i) {
    os << (i==0 ? '(' : ',') << x[i]; 
  }
  os << ")";
  return os;
}


}
