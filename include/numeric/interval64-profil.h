/***************************************************************************
 *            numeric/interval64-profil.h
 *
 *  Copyright  2008  Pieter Collins
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
 
/*! \file numeric/interval64-profil.h
 *  \brief Type definitions and conversion operators for 64-bit fixed precision floating point numbers using the PROFIL library.
 */

#ifndef ARIADNE_NUMERIC_INTERVAL64_PROFIL_H
#define ARIADNE_NUMERIC_INTERVAL64_PROFIL_H

#include <iosfwd>
#include <iostream>
#include <sstream>

#include <Interval.h>
#include <Functions.h>

namespace Ariadne {
namespace Numeric {

class Integer;
class Rational;
template<class T> class Float;
template<class R> class Interval;

typedef Float<double> Float64;
typedef Interval<Float64> Interval64;


/*!\ingroup Numeric
 * \brief A templated class representing an interval of real numbers.
 * 
 * An interval of real numbers with endpoints of type \a R.
 * All operations on an interval must be guarenteed to return an interval contining the exact result.
 * If \a val of real numbers with endpoints of type \a R.
 * All operations on an interval must be guarenteed to return an interval contining the exact result.
 * If \a T supports exact evaluation of a function, then the exact evaluation must be used.
 * If \a T is dense in the reals, e.g. dyadic or rational, then any approximate operations may be given a max_imum error of computation.
 *
 * Currently implemented as a wrapper around the boost::numeric::interval class template from the Boost C++ library.
 */
template<>
class Interval< Float64 >
{
  typedef Float64 R;
 public:
  INTERVAL _value;
  public:
  //@{
  //! \name Constructors and assignment operators
  /*! \brief Default constructer constructs empty interval. */
  Interval() : _value() { }
  /*! \brief Construct from a string literal. */
  Interval(const std::string& str) { std::stringstream ss(str); ss >> _value; }
  /*! \brief Construct a built-in integer. */
  Interval(const int& x) : _value(x) { }
  Interval(const uint& x) : _value(x) { }
  /*! \brief Construct a built-in double. */
  Interval(const double& x) : _value(x) { }
  /*! \brief Construct from lower and upper bounds. */
  Interval(const R& l, const R& u) : _value(l._value,u._value) { }
  Interval(const double& l, const double& u) : _value(l,u) { }
  /*! \brief Copy constructor. */
  Interval(const Interval<R>& ivl) : _value(ivl._value) { };
  Interval(const INTERVAL& ivl) : _value(ivl) { };
  
  /*! \brief Assign from a number. */
  Interval<R>& operator=(const R& x) { _value=x._value; return *this; } 
  /*! \brief Copy assignment operator. */
  Interval<R>& operator=(const Interval<R>& ivl) { _value=ivl._value; return *this; }
  
  /*! \brief Construct from an expression. */
  template<class E> Interval(const Expression<E>& e);
  /*! \brief Assign from an expression. */
  template<class E> Interval<R>& operator=(const Expression<E>& e);
  
  /*! \brief Construct a one-point interval. */
  template<class RX> Interval(const RX& x);
  /*! \brief Construct from lower and upper bounds. */
  template<class RL,class RU> Interval(const RL& l, const RU& u);
  /*! \brief Construct an interval with possibly different real type. */
  template<class RX> Interval(const Interval<RX>& ivl);
  /*! \brief Assign from a point of a possibly different type. */
  template<class RX> Interval<R>& operator=(const RX& x);
  /*! \brief Assign from an interval of possibly a different type. */
  template<class RX> Interval<R>& operator=(const Interval<RX>& ivl);
  //@}
  
  //@{
  //! \name Data access
  /*! \brief The lower bound. */
  const R lower() const { return ::Inf(_value); }
  /*! \brief The upper bound. */
  const R upper() const { return ::Sup(_value); }
  /*! \brief The midpoint of the interval, given by \f$(a+b)/2\f$. */
  R midpoint() const { return ::Mid(_value); }
  /*! \brief The radius of the interval, given by \f$(b-a)/2\f$. */
  R radius() const;
  /*! \brief The width of the interval, given by \f$b-a\f$. */
  R width() const { return Diam(_value); }
  //@}
  
  //@{
  //! \name Geometric operations
  
  /*! \brief Tests if the interval is empty. */
  bool empty() const { return this->lower() > this->upper(); }
  /*! \brief Tests if the interval consists of a single point. */
  bool singleton() const { return this->lower() == this->upper(); }
  /*! \brief Tests if the interval contains \a x. */
  template<class RX> bool encloses(const RX& x) const {
    return this->lower() <= x && x <= this->upper(); }
  /*! \brief Tests if the interval contains \a r. */
  template<class RX> bool refines(const Interval<RX>& ivl) const {
    return ivl.lower() <= this->lower() && this->upper() <= ivl->upper(); };
  //}
};
  
inline Interval64 pos(const Interval64& x) { 
  return x; }
inline Interval64 neg(const Interval64& x) { 
  return -x._value; }
 
inline Interval64 add(const Float64& x, const Float64& y) {
  return AddBounds(x._value,y._value); }
inline Interval64 sub(const Float64& x, const Float64& y) {
  return SubBounds(x._value,y._value); }
inline Interval64 mul(const Float64& x, const Float64& y) {
  return MulBounds(x._value,y._value); }
inline Interval64 div(const Float64& x, const Float64& y) {
  return DivBounds(x._value,y._value); }

inline Interval64 add(const Interval64& x, const Float64& y) {
  return operator+(x._value,y._value); }
inline Interval64 sub(const Interval64& x, const Float64& y) {
  return operator-(x._value,y._value); }
inline Interval64 mul(const Interval64& x, const Float64& y) {
  return operator*(x._value,y._value); }
inline Interval64 div(const Interval64& x, const Float64& y) {
  return operator/(x._value,y._value); }

inline Interval64 add(const Float64& x, const Interval64& y) {
  return operator+(x._value,y._value); }
inline Interval64 sub(const Float64& x, const Interval64& y) {
  return operator-(x._value,y._value); }
inline Interval64 mul(const Float64& x, const Interval64& y) {
  return operator*(x._value,y._value); }
inline Interval64 div(const Float64& x, const Interval64& y) {
  return operator/(x._value,y._value); }

inline Interval64 add(const Interval64& x, const Interval64& y) {
  return operator+(x._value,y._value); }
inline Interval64 sub(const Interval64& x, const Interval64& y) {
  return operator-(x._value,y._value); }
inline Interval64 mul(const Interval64& x, const Interval64& y) {
  return operator*(x._value,y._value); }
inline Interval64 div(const Interval64& x, const Interval64& y) {
  return operator/(x._value,y._value); }




inline Interval64 operator+(const Interval64& x) {
  return x; }
inline Interval64 operator-(const Interval64& x) {
  return -x._value; }

inline Interval64 operator+(const Float64& x, const Float64& y) {
  return AddBounds(x._value,y._value); }
inline Interval64 operator-(const Float64& x, const Float64& y) {
  return SubBounds(x._value,y._value); }
inline Interval64 operator*(const Float64& x, const Float64& y) {
  return MulBounds(x._value,y._value); }
inline Interval64 operator/(const Float64& x, const Float64& y) {
  return DivBounds(x._value,y._value); }

inline Interval64 operator+(const Interval64& x, const Float64& y) {
  return operator+(x._value,y._value); }
inline Interval64 operator-(const Interval64& x, const Float64& y) {
  return operator-(x._value,y._value); }
inline Interval64 operator*(const Interval64& x, const Float64& y) {
  return operator*(x._value,y._value); }
inline Interval64 operator/(const Interval64& x, const Float64& y) {
  return operator/(x._value,y._value); }

inline Interval64 operator+(const Float64& x, const Interval64& y) {
  return operator+(x._value,y._value); }
inline Interval64 operator-(const Float64& x, const Interval64& y) {
  return operator-(x._value,y._value); }
inline Interval64 operator*(const Float64& x, const Interval64& y) {
  return operator*(x._value,y._value); }
inline Interval64 operator/(const Float64& x, const Interval64& y) {
  return operator/(x._value,y._value); }

inline Interval64 operator+(const Interval64& x, const Interval64& y) {
  return operator+(x._value,y._value); }
inline Interval64 operator-(const Interval64& x, const Interval64& y) {
  return operator-(x._value,y._value); }
inline Interval64 operator*(const Interval64& x, const Interval64& y) {
  return operator*(x._value,y._value); }
inline Interval64 operator/(const Interval64& x, const Interval64& y) {
  return operator/(x._value,y._value); }

inline Interval64 operator+=(Interval64& x, const Float64& y) {
  x._value+=y._value; return x; }
inline Interval64 operator-=(Interval64& x, const Float64& y) {
  x._value-=y._value; return x; }
inline Interval64 operator*=(Interval64& x, const Float64& y) {
  x._value*=y._value; return x; }
inline Interval64 operator/=(Interval64& x, const Float64& y) {
  x._value/=y._value; return x; }

inline Interval64 operator+=(Interval64& x, const Interval64& y) {
  x._value+=y._value; return x; }
inline Interval64 operator-=(Interval64& x, const Interval64& y) {
  x._value-=y._value; return x; }
inline Interval64 operator*=(Interval64& x, const Interval64& y) {
  x._value*=y._value; return x; }
inline Interval64 operator/=(Interval64& x, const Interval64& y) {
  x._value/=y._value; return x; }

inline Interval64 abs(const Interval64& x) {
  return ::Abs(x._value); }
inline Interval64 pow(const Interval64& x, int n) {
  return ::Power(x._value,n); }
inline Interval64 sqrt(const Interval64& x) {
  return ::Sqrt(x._value); }
inline Interval64 exp(const Interval64& x) {
  return ::Exp(x._value); }
inline Interval64 log(const Interval64& x) {
  return ::Log(x._value); }
inline Interval64 sin(const Interval64& x) {
  return ::Sin(x._value); }
inline Interval64 cos(const Interval64& x) {
  return ::Cos(x._value); }
inline Interval64 tan(const Interval64& x) {
  return ::Tan(x._value); }
inline Interval64 asin(const Interval64& x) {
  return ::ArcSin(x._value); }
inline Interval64 acos(const Interval64& x) {
  return ::ArcCos(x._value); }
inline Interval64 atan(const Interval64& x) {
  return ::ArcTan(x._value); }
inline Interval64 sinh(const Interval64& x) {
  return ::Sinh(x._value); }
inline Interval64 cosh(const Interval64& x) {
  return ::Cosh(x._value); }
inline Interval64 tanh(const Interval64& x) {
  return ::Tanh(x._value); }
inline Interval64 asinh(const Interval64& x) {
  return ::ArSinh(x._value); }
inline Interval64 acosh(const Interval64& x) {
  return ::ArCosh(x._value); }
inline Interval64 atanh(const Interval64& x) {
  return ::ArTanh(x._value); }

inline Interval64 pow(const Float64& x, int n) {
  return pow(Interval64(x),n); }
inline Interval64 sqrt(const Float64& x) {
  return sqrt(Interval64(x)); }

inline tribool operator<(const Float64& x, const Interval64& y) {
  if(x<y.lower()) { return true; } 
  if(x>=y.upper()) { return false; }
  return indeterminate;
}

inline tribool operator<(const Interval64& x, const Float64& y) {
  if(x.upper()<y) { return true; } 
  if(x.lower()>=y) { return false; }
  return indeterminate;
}
inline tribool operator<(const Interval64& x, const Interval64& y) {
  if(x.upper()<y.lower()) { return true; } 
  if(x.lower()>=y.upper()) { return false; }
  return indeterminate;
}

inline std::ostream& operator<<(std::ostream& os, const Interval64& x) {
  return os << x._value; }
inline std::istream& operator>>(std::istream& is, Interval64& x) {
  return is >> x._value; }


template<class R> void instantiate_interval();

template<> inline
void instantiate_interval<Float64>() {
}



}
}

#endif /* ARIADNE_NUMERIC_INTERVAL64_PROFIL_H */
