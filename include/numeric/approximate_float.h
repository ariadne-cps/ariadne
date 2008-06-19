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

#include "float.h"
#include "interval.h"

namespace Ariadne {
  

    template<class T>
    class ApproximateFloat
    {
     public:
      Float<T> _value;
     public:
      ApproximateFloat() : _value(0) { }
      ApproximateFloat(const int& n) : _value(n) { }
      ApproximateFloat(const double& x) : _value(x) { }
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


    template<class T> inline
    bool
    operator==(const ApproximateFloat<T>& x, const ApproximateFloat<T>& y) {
      return x._value==y._value; }

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
    operator/(const int& x, const ApproximateFloat<T>& y) {
      return div_approx(x,y._value); }

    template<class T> inline
    ApproximateFloat<T> 
    operator/(const ApproximateFloat<T>& x, const int& y) {
      return div_approx(x._value,y); }

    template<class T> inline
    ApproximateFloat<T> 
    sqrt(const ApproximateFloat<T>& x) {
      return ApproximateFloat<T>(sqrt_approx(x._value)); }



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

  
} // namespace Ariadne

#endif /* ARIADNE_NUMERIC_APPROXIMATE_FLOAT_H */
