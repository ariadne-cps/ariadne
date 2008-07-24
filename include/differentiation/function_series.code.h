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
 
#include "differentiation/power_series.h"
#include "differentiation/function_series.h"

namespace Ariadne {

template<class X> 
PowerSeries<X> 
ArithmeticSeries<X>::rec(smoothness_type d, const X& c) 
{
  PowerSeries<X> y(d);
  X mr = (-1)/c;
  for(uint i=0; i<=y.degree(); ++i) {
    y[i]=-Ariadne::pow(mr,i+1u);
  }
  return y;
}


template<class X> 
PowerSeries<X> 
ArithmeticSeries<X>::pow(smoothness_type d, const X& c, const uint& k)
{
  uint n=k;
  PowerSeries<X> y(d);
  for(uint i=0; i<=std::min(uint(d),n); ++i) {
    uint j=n-i;
    y[i]=X(Ariadne::bin<Integer>(n,j))*Ariadne::pow(c,j);
  }
  return y;
}

template<class X> 
PowerSeries<X> 
TranscendentalSeries<X>::sqrt(smoothness_type d, const X& c)
{
  PowerSeries<X> y(d);
  y[0]=Ariadne::sqrt(c);
  X mhr=-0.5/c;
  for(uint i=1; i<=y.degree(); ++i) {
    // Need to convert uint to int to prevent wraparound for 2*1u-3
    y[i]=((2*int(i)-3)*mhr)/i*y[i-1];
  }
  return y;
}

template<class X> 
PowerSeries<X> 
TranscendentalSeries<X>::exp(smoothness_type d, const X& c)
{
  PowerSeries<X> y(d);
  y[0]=Ariadne::exp(c);
  for(uint i=1; i<=y.degree(); ++i) {
    y[i]=y[i-1]/i;
  }
  return y;
}

template<class X>  
PowerSeries<X> 
TranscendentalSeries<X>::log(smoothness_type d, const X& c)
{
  PowerSeries<X> y(d);
  y[0]=Ariadne::log(c);
  X mr=(-1)/c;
  for(uint i=1; i<=y.degree();++i) {
    y[i]=-Ariadne::pow(mr,i)/i;
  }
  return y;
}

template<class X> 
PowerSeries<X> 
TranscendentalSeries<X>::sin(smoothness_type d, const X& c)
{
  PowerSeries<X> y(d);
  y[0]=Ariadne::sin(c);
  if(d>=1) {
    y[1]=Ariadne::cos(c);
    for(uint i=2; i<=d; ++i) {
      y[i]=-y[i-2]/(i*(i-1));
    }
  }
  return y;
}

template<class X> 
PowerSeries<X> 
TranscendentalSeries<X>::cos(smoothness_type d, const X& c)
{
  PowerSeries<X> y(d);
  y[0]=Ariadne::cos(c);
  if(d>=1) {
    y[1]=-Ariadne::sin(c);
    for(uint i=2; i<=d; ++i) {
      y[i]=-y[i-2]/(i*(i-1));
    }
  }
  return y;
}

template<class X> 
PowerSeries<X> 
TranscendentalSeries<X>::tan(smoothness_type d, const X& c)
{
  return sin(d,c)/cos(d,c);
}

template<class X>  
PowerSeries<X> 
TranscendentalSeries<X>::asin(smoothness_type d, const X& c)
{
  if(d==0) { return PowerSeries<X>::constant(d,Ariadne::atan(c)); }
  PowerSeries<X> y = X(1)/Ariadne::sqrt(X(1)-Ariadne::pow(PowerSeries<X>::variable(d-1,c),2u));
  return antiderivative(y,Ariadne::asin(c));
}

template<class X>  
PowerSeries<X> 
TranscendentalSeries<X>::acos(smoothness_type d, const X& c)
{
  if(d==0) { return PowerSeries<X>::constant(d,Ariadne::atan(c)); }
  PowerSeries<X> y = X(-1)/Ariadne::sqrt(X(1)-pow(PowerSeries<X>::variable(d-1,c),2u));
  return antiderivative(y,Ariadne::acos(c));
}

template<class X>  
PowerSeries<X> 
TranscendentalSeries<X>::atan(smoothness_type d, const X& c)
{
  if(d==0) { return PowerSeries<X>::constant(d,Ariadne::atan(c)); } 
  PowerSeries<X> y = X(1)/(X(1)+pow(PowerSeries<X>::variable(d-1,c),2u));
  return antiderivative(y,Ariadne::atan(c));
}


}
