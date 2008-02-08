/***************************************************************************
 *            function_series.template.h
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
 
#include "function/taylor_series.h"
#include "function/function_series.h"

namespace Ariadne {

template<class X> 
Function::TaylorSeries<X> 
Function::ArithmeticSeries<X>::rec(smoothness_type d, const X& c) 
{
  TaylorSeries<X> y(d);
  X mr = -1/c;
  for(uint i=0; i<=y.degree(); ++i) {
    y[i]=-Numeric::pow(mr,i+1u);
  }
  return y;
}


template<class X> 
Function::TaylorSeries<X> 
Function::ArithmeticSeries<X>::pow(smoothness_type d, const X& c, const uint& k)
{
  uint n=k;
  TaylorSeries<X> y(d);
  for(uint i=0; i<=std::min(uint(d),n); ++i) {
    uint j=n-i;
    y[i]=X(Numeric::bin<Numeric::Integer>(n,j))*Numeric::pow(c,j);
  }
  //std::cout << "pow("<<d<<","<<c<<","<<k<<")="<<y<<std::endl;
  return y;
}

template<class X> 
Function::TaylorSeries<X> 
Function::TranscendentalSeries<X>::sqrt(smoothness_type d, const X& c)
{
TaylorSeries<X> y(d);
  y[0]=Numeric::sqrt(c);
  X mhr=-0.5/c;
  for(uint i=1; i<=y.degree(); ++i) {
    y[i]=((2*i-3)*mhr)/i*y[i-1];
  }
  return y;
}

template<class X> 
Function::TaylorSeries<X> 
Function::TranscendentalSeries<X>::exp(smoothness_type d, const X& c)
{
TaylorSeries<X> y(d);
  y[0]=Numeric::exp(c);
  for(uint i=1; i<=y.degree(); ++i) {
    y[i]=y[i-1]/i;
  }
  return y;
}

template<class X>  
Function::TaylorSeries<X> 
Function::TranscendentalSeries<X>::log(smoothness_type d, const X& c)
{
TaylorSeries<X> y(d);
  y[0]=Numeric::log(c);
  X mr=(-1)/c;
  for(uint i=1; i<=y.degree();++i) {
    y[i]=-Numeric::pow(mr,i)/i;
  }
  return y;
}

template<class X> 
Function::TaylorSeries<X> 
Function::TranscendentalSeries<X>::sin(smoothness_type d, const X& c)
{
TaylorSeries<X> y(d);
  y[0]=Numeric::sin(c);
  y[1]=Numeric::cos(c);
  for(uint i=2; i!=d; ++i) {
    y[i]=-y[i-2]/(i*(i-1));
  }
  return y;
}

template<class X> 
Function::TaylorSeries<X> 
Function::TranscendentalSeries<X>::cos(smoothness_type d, const X& c)
{
  TaylorSeries<X> y(d);
  y[0]=Numeric::cos(c);
  y[1]=-Numeric::sin(c);
  for(uint i=2; i!=d; ++i) {
    y[i]=-y[i-2]/(i*(i-1));
  }
  return y;
}

template<class X> 
Function::TaylorSeries<X> 
Function::TranscendentalSeries<X>::tan(smoothness_type d, const X& c)
{
  return sin(d,c)/cos(d,c);
}

template<class X>  
Function::TaylorSeries<X> 
Function::TranscendentalSeries<X>::asin(smoothness_type d, const X& c)
{
  if(d==0) { return TaylorSeries<X>::constant(d,Numeric::atan(c)); }
  TaylorSeries<X> y = X(1)/Function::sqrt(X(1)-Function::pow(TaylorSeries<X>::variable(d-1,c),2));
  return antiderivative(y,Numeric::asin(c));
}

template<class X>  
Function::TaylorSeries<X> 
Function::TranscendentalSeries<X>::acos(smoothness_type d, const X& c)
{
  if(d==0) { return TaylorSeries<X>::constant(d,Numeric::atan(c)); }
  TaylorSeries<X> y = X(-1)/Function::sqrt(X(1)-Function::pow(TaylorSeries<X>::variable(d-1,c),2));
  return antiderivative(y,Numeric::acos(c));
}

template<class X>  
Function::TaylorSeries<X> 
Function::TranscendentalSeries<X>::atan(smoothness_type d, const X& c)
{
  if(d==0) { return TaylorSeries<X>::constant(d,Numeric::atan(c)); } 
  TaylorSeries<X> y = X(1)/(X(1)+Function::pow(TaylorSeries<X>::variable(d-1,c),2));
  return antiderivative(y,Numeric::atan(c));
}


}
