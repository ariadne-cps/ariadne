/***************************************************************************
 *            numeric/approximate_float.h
 *
 *  Copyright  2007 Pieter Collins
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
 
/*! \file numeric/approximate_float.h
 *  \brief Floating-point numbers with approximate arithmetic. */

#ifndef ARIADNE_NUMERIC_APPROXIMATE_FLOAT_H
#define ARIADNE_NUMERIC_APPROXIMATE_FLOAT_H

#include "numeric/float.h"
#include "numeric/interval.h"

namespace Ariadne {
  

    template<class T>
    class ApproximateFloat
    {
     public:
      Float<T> _value;
     public:
      ApproximateFloat() : _value(0) { }
      ApproximateFloat(const int& n) : _value(n) { }
      ApproximateFloat(const uint& n) : _value(n) { }
      ApproximateFloat(const double& x) : _value(x) { }
      explicit ApproximateFloat(const Integer& z) : _value(Float<T>(z)._value) { }
      explicit ApproximateFloat(const Float<T>& x) : _value(x) { }
      explicit ApproximateFloat(const Interval< Float<T> >& x) : _value(x.midpoint()) { }
      ApproximateFloat(const ApproximateFloat<T>& x) : _value(x._value) { }
  
      ApproximateFloat<T>& operator=(const int& n) { 
        _value=n; return *this; }
      ApproximateFloat<T>& operator=(const double& x) { 
        _value=x; return *this; }
      ApproximateFloat<T>& operator=(const ApproximateFloat<T>& x) { 
        _value=x._value; return *this; }
    };

    inline Float64::Float(const ApproximateFloat64& x) : _value(x._value._value) { }

    template<class T> inline
    bool
    operator==(const ApproximateFloat<T>& x, const int& y) {
      return x._value==y; }

    template<class T> inline
    bool
    operator==(const ApproximateFloat<T>& x, const double& y) {
      return x._value==y; }

    template<class T> inline
    bool
    operator==(const ApproximateFloat<T>& x, const ApproximateFloat<T>& y) {
      return x._value==y._value; }

    template<class T> inline
    bool
    operator!=(const ApproximateFloat<T>& x, const int& y) {
      return x._value!=y; }

    template<class T> inline
    bool
    operator!=(const ApproximateFloat<T>& x, const double& y) {
      return x._value!=y; }

    template<class T> inline
    bool
    operator!=(const ApproximateFloat<T>& x, const ApproximateFloat<T>& y) {
      return x._value!=y._value; }

    template<class T> inline
    bool
    operator>=(const ApproximateFloat<T>& x, const ApproximateFloat<T>& y) {
      return x._value>=y._value; }

    template<class T> inline
    bool
    operator<=(const ApproximateFloat<T>& x, const ApproximateFloat<T>& y) {
      return x._value<=y._value; }

    template<class T> inline
    bool
    operator> (const ApproximateFloat<T>& x, const ApproximateFloat<T>& y) {
      return x._value> y._value; }

    template<class T> inline
    bool
    operator< (const ApproximateFloat<T>& x, const ApproximateFloat<T>& y) {
      return x._value< y._value; }


    template<class T> inline
    ApproximateFloat<T> 
    abs(const ApproximateFloat<T>& x) {
      return ApproximateFloat<T>(Float<T>(abs(x._value))); }

    template<class T> inline
    ApproximateFloat<T> 
    operator+(const ApproximateFloat<T>& x) {
      return ApproximateFloat<T>(x._value); }

    template<class T> inline
    ApproximateFloat<T> 
    operator-(const ApproximateFloat<T>& x) {
      return ApproximateFloat<T>(Float<T>(-x._value)); }

    template<class T> inline
    ApproximateFloat<T> 
    operator+(const ApproximateFloat<T>& x, const ApproximateFloat<T>& y) {
      return ApproximateFloat<T>(add_approx(x._value,y._value)); }
  
    template<class T> inline
    ApproximateFloat<T> 
    operator-(const ApproximateFloat<T>& x, const ApproximateFloat<T>& y) {
      return ApproximateFloat<T>(sub_approx(x._value,y._value)); }

    template<class T> inline
    ApproximateFloat<T> 
    operator*(const ApproximateFloat<T>& x, const ApproximateFloat<T>& y) {
      return ApproximateFloat<T>(mul_approx(x._value,y._value)); }

    template<class T> inline
    ApproximateFloat<T> 
    operator/(const ApproximateFloat<T>& x, const ApproximateFloat<T>& y) {
      return ApproximateFloat<T>(div_approx(x._value,y._value)); }

    template<class T> inline
    ApproximateFloat<T> 
    pow(const ApproximateFloat<T>& x, const int& y) {
      return ApproximateFloat<T>(pow_approx(x._value,y)); }

    template<class T> inline
    ApproximateFloat<T> 
    pow(const ApproximateFloat<T>& x, const uint& y) {
      return ApproximateFloat<T>(pow_approx(x._value,y)); }



    template<class T> inline
    ApproximateFloat<T>&
    operator+=(ApproximateFloat<T>& x, const int& y) {
      x._value._value+=y; return x; }

    template<class T> inline
    ApproximateFloat<T>&
    operator*=(ApproximateFloat<T>& x, const int& y) {
      x._value._value*=y; return x; }


    template<class T> inline
    ApproximateFloat<T>&
    operator+=(ApproximateFloat<T>& x, const double& y) {
      x._value._value+=y; return x; }

    template<class T> inline
    ApproximateFloat<T>&
    operator-=(ApproximateFloat<T>& x, const double& y) {
      x._value._value-=y; return x; }

    template<class T> inline
    ApproximateFloat<T>&
    operator*=(ApproximateFloat<T>& x, const double& y) {
      x._value._value*=y; return x; }

    template<class T> inline
    ApproximateFloat<T>&
    operator/=(ApproximateFloat<T>& x, const double& y) {
      x._value._value/=y; return x; }


    template<class T> inline
    ApproximateFloat<T>&
    operator+=(ApproximateFloat<T>& x, const ApproximateFloat<T>& y) {
      x._value=add_approx(x._value,y._value); return x; }

    template<class T> inline
    ApproximateFloat<T>&
    operator-=(ApproximateFloat<T>& x, const ApproximateFloat<T>& y) {
      x._value=sub_approx(x._value,y._value); return x; }

    template<class T> inline
    ApproximateFloat<T>&
    operator*=(ApproximateFloat<T>& x, const ApproximateFloat<T>& y) {
      x._value=mul_approx(x._value,y._value); return x; }

    template<class T> inline
    ApproximateFloat<T>&
    operator/=(ApproximateFloat<T>& x, const ApproximateFloat<T>& y) {
      x._value=div_approx(x._value,y._value); return x; }



    template<class T> inline
    ApproximateFloat<T>
    operator*(int n, const ApproximateFloat<T>& x) {
      return ApproximateFloat<T>(n)*x; }

    template<class T> inline
    ApproximateFloat<T>
    operator*(const ApproximateFloat<T>& x, int n) {
      return x*ApproximateFloat<T>(n); }

    template<class T> inline
    ApproximateFloat<T>
    operator/(int n, const ApproximateFloat<T>& x) {
      return ApproximateFloat<T>(n)/x; }

    template<class T> inline
    ApproximateFloat<T>
    operator/(const ApproximateFloat<T>& x, int n) {
      return x/ApproximateFloat<T>(n); }

    template<class T> inline
    ApproximateFloat<T>
    operator/(double d, const ApproximateFloat<T>& x) {
      return ApproximateFloat<T>(d)/x; }


    template<class T> inline
    Interval< Float<T> > 
    operator+(const ApproximateFloat<T>& x, const Interval< Float<T> >& y) {
      return Float<T>(x._value)+y; }

    template<class T> inline
    Interval< Float<T> > 
    operator+(const Interval< Float<T> >& x, const ApproximateFloat<T>& y) {
      return x+Float<T>(y._value); }


    template<class T> inline
    Interval< Float<T> > 
    operator-(const ApproximateFloat<T>& x, const Interval< Float<T> >& y) {
      return Float<T>(x._value)-y; }

    template<class T> inline
    Interval< Float<T> > 
    operator-(const Interval< Float<T> >& x, const ApproximateFloat<T>& y) {
      return x-Float<T>(y._value); }


    template<class T> inline
    Interval< Float<T> > 
    operator*(const ApproximateFloat<T>& x, const Interval< Float<T> >& y) {
      return Float<T>(x._value)*y; }

    template<class T> inline
    Interval< Float<T> > 
    operator*(const Interval< Float<T> >& x, const ApproximateFloat<T>& y) {
      return x*Float<T>(y._value); }



    template<class T> inline
    Interval< Float<T> > 
    operator*(const ApproximateFloat<T>& x, const Float<T>& y) {
      return Float<T>(x._value)*y; }

    template<class T> inline
    Interval< Float<T> > 
    operator*(const Float<T>& x, const ApproximateFloat<T>& y) {
      return x*Float<T>(y._value); }


    template<class T> inline
    ApproximateFloat<T>
    sqrt(const ApproximateFloat<T>& y) {
      return std::sqrt(y._value._value); }

    template<class T> inline
    ApproximateFloat<T>
    exp(const ApproximateFloat<T>& y) {
      return std::exp(y._value._value); }

    template<class T> inline
    ApproximateFloat<T>
    log(const ApproximateFloat<T>& y) {
      return std::log(y._value._value); }

    template<class T> inline
    ApproximateFloat<T>
    sin(const ApproximateFloat<T>& y) {
      return std::sin(y._value._value); }

    template<class T> inline
    ApproximateFloat<T>
    cos(const ApproximateFloat<T>& y) {
      return std::cos(y._value._value); }

    template<class T> inline
    ApproximateFloat<T>
    tan(const ApproximateFloat<T>& y) {
      return std::tan(y._value._value); }

    template<class T> inline
    ApproximateFloat<T>
    asin(const ApproximateFloat<T>& y) {
      return std::asin(y._value._value); }

    template<class T> inline
    ApproximateFloat<T>
    acos(const ApproximateFloat<T>& y) {
      return std::acos(y._value._value); }

    template<class T> inline
    ApproximateFloat<T>
    atan(const ApproximateFloat<T>& y) {
      return std::atan(y._value._value); }

 
    template<class T> inline
    std::ostream& 
    operator<<(std::ostream& os, const ApproximateFloat<T>& x) {
      return os << x._value; }

    template<class T> inline
    std::istream& 
    operator>>(std::istream& is, ApproximateFloat<T>& x) {
      return is >> x._value; }


} // namespace Ariadne

#endif /* ARIADNE_NUMERIC_APPROXIMATE_FLOAT_H */
