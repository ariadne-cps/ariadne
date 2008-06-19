/***************************************************************************
 *            numeric/error_float.h
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
 
/*! \file numeric/error_float.h
 *  \brief Floating-point numbers with error arithmetic. */

#ifndef ARIADNE_NUMERIC_ERROR_FLOAT_H
#define ARIADNE_NUMERIC_ERROR_FLOAT_H

#include "float.h"
#include "interval.h"

namespace Ariadne {
  

    template<class T>
    class ErrorFloat
    {
     public:
      Float<T> _value;
     public:
      ErrorFloat() : _value(0) { }
      ErrorFloat(const int& n) : _value(n) { }
      ErrorFloat(const double& x) : _value(x) { }
      explicit ErrorFloat(const Float<T>& x) : _value(x) { }
      explicit ErrorFloat(const Interval< Float<T> >& x) : _value(rad_up(x.lower(),x.upper())) { }
      ErrorFloat(const ErrorFloat<T>& x) : _value(x._value) { }
  
      ErrorFloat<T>& operator=(const int& n) { 
        _value=n; return *this; }
      ErrorFloat<T>& operator=(const double& x) { 
        _value=x; return *this; }
      ErrorFloat<T>& operator=(const ErrorFloat<T>& x) { 
        _value=x._value; return *this; }
    };


    template<class T> inline
    bool
    operator==(const ErrorFloat<T>& x, const ErrorFloat<T>& y) {
      return x._value==y._value; }

    template<class T> inline
    bool
    operator!=(const ErrorFloat<T>& x, const ErrorFloat<T>& y) {
      return x._value!=y._value; }

    template<class T> inline
    bool
    operator>=(const ErrorFloat<T>& x, const ErrorFloat<T>& y) {
      return x._value>=y._value; }

    template<class T> inline
    bool
    operator<=(const ErrorFloat<T>& x, const ErrorFloat<T>& y) {
      return x._value<=y._value; }

    template<class T> inline
    bool
    operator> (const ErrorFloat<T>& x, const ErrorFloat<T>& y) {
      return x._value> y._value; }

    template<class T> inline
    bool
    operator< (const ErrorFloat<T>& x, const ErrorFloat<T>& y) {
      return x._value< y._value; }


    template<class T> inline
    ErrorFloat<T> 
    rad(const Float<T>& x, const Float<T>& y) {
      return ErrorFloat<T>(div_up(sub_up(x._value,y._value),2)); }

    template<class T> inline
    ErrorFloat<T> 
    operator+(const ErrorFloat<T>& x) {
      return ErrorFloat<T>(x._value); }

    template<class T> inline
    ErrorFloat<T> 
    operator+(const ErrorFloat<T>& x, const ErrorFloat<T>& y) {
      return ErrorFloat<T>(add_up(x._value,y._value)); }
  
    template<class T> inline
    ErrorFloat<T> 
    operator*(const ErrorFloat<T>& x, const int& y) {
      return ErrorFloat<T>(mul_up(x._value,y)); }

    template<class T> inline
    ErrorFloat<T> 
    operator/(const ErrorFloat<T>& x, const int& y) {
      return div_up(x._value,y); }



    template<class T> inline
    ErrorFloat<T>&
    operator+=(ErrorFloat<T>& x, const ErrorFloat<T>& y) {
      x._value=add_up(x._value,y._value); return x; }

    template<class T> inline
    ErrorFloat<T>&
    operator*=(ErrorFloat<T>& x, const ErrorFloat<T>& y) {
      x._value=mul_up(x._value,y._value); return x; }

    template<class T> inline
    ErrorFloat<T>&
    operator*=(ErrorFloat<T>& x, const int& y) {
      x._value=mul_up(x._value,y); return x; }



} // namespace Ariadne


#endif /* ARIADNE_NUMERIC_ERROR_FLOAT_H */
