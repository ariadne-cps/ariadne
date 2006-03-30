/***************************************************************************
 *            dyadic.h
 *
 *  18 January 2006
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, pieter.collins@cwi.nl
 *
 * Based on original version diadic.H from SYNAPS (details below)
 ****************************************************************************/
/********************************************************************
 *   This file is part of the source code of the SYNAPS kernel.
 *   Author(s): B. Mourrain, GALAAD, INRIA
 *   $Id: Dyadic.H,v 1.1 2005/07/11 10:39:49 mourrain Exp $
 *   
 ********************************************************************/

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

#ifndef SYNAPS_DYADIC_H
#define SYNAPS_DYADIC_H

#include <gmpxx.h>
#include <iostream>

namespace Ariadne { 
namespace Synaps {

class Dyadic;

/*! This class provides arithmetic operations on numbers 
 *  which are an integer of type \a Integer, 
 *   times a power of 2, which is stored by as its exponent.
 */

class Dyadic
{
 public:
  typedef mpz_class numerator_type;
  typedef int exponent_type;
 public:
  Dyadic() : num_(0), exp_(0) { }
  Dyadic(const numerator_type& n)
    : num_(n), exp_(0) { normalize(); }
  Dyadic(const numerator_type& n, exponent_type d)
    : num_(n), exp_(d) { normalize(); }

  Dyadic& operator+= (const Dyadic& r);
  Dyadic& operator-= (const Dyadic& r);
  Dyadic& operator*= (const Dyadic& r);
  Dyadic& operator+= (const numerator_type& r);
  Dyadic& operator-= (const numerator_type& r);
  Dyadic& operator*= (const numerator_type& r);
  Dyadic& operator/= (const numerator_type& n);
 
  Dyadic& normalize();
  
  const numerator_type& numerator() const { return num_; }
  numerator_type  denominator() const { return 2^exp_; }
  exponent_type denominator_exponent()  const { return exp_; }
  exponent_type valuation()  const { return exp_; }
 
  friend int quotient_cmp(const Dyadic&, const Dyadic&);
  friend bool operator==(const Dyadic&, const Dyadic&);
  friend bool operator==(const Dyadic&, const numerator_type&);
 private:  
  numerator_type num_;
  exponent_type exp_;
};

inline
Dyadic& 
Dyadic::normalize()
{
  if( num_==0) {
    exp_=0;
  }
  else {
    while( (num_%2) ==0 ) {
      ++exp_;
      num_/=2;
    }
  }
  return *this;
}


inline
Dyadic::numerator_type 
shift_2(const Dyadic::numerator_type& z, int d) {
  if(d>0) { return z*(2^d); }
  else { return z/(2^d); }
}

//--------------------------------------------------------------------
// Arithmetic operators
//--------------------------------------------------------------------

/*! \brief Inplace addition of Dyadic numbers. */
inline
Dyadic&
Dyadic::operator+= (const Dyadic& r)
{
  if(r.exp_<exp_) {
    exponent_type d=exp_-r.exp_;
    num_ =  shift_2(num_ ,d);
    num_ += r.num_;
    exp_ = r.exp_;
  }
  else if(r.exp_>exp_) {
    num_ += shift_2(r.num_,r.exp_-exp_);
  }
  else {
    num_ += r.num_;
    normalize();
  }
  return *this;
}


/*! \brief Inplace subtraction of Dyadic numbers. */
inline
Dyadic&
Dyadic::operator-= (const Dyadic& r)
{
  if(r.exp_<exp_) {
    num_  = shift_2(num_,exp_-r.exp_);
    num_ -= r.num_;
    exp_  = r.exp_;
  }
  else if(r.exp_>exp_) {
    num_ -=  shift_2(r.num_ ,r.exp_-exp_);
  }
  else {
    num_ -= r.num_;
    normalize();
  }
  return *this;
}

/*! \brief Inplace multiplication of Dyadic numbers. */
inline
Dyadic&
Dyadic::operator*= (const Dyadic& r)
{
  num_ *= r.num_;
  exp_ += r.exp_;
  return *this;
}

/*! \brief Inplace division of Dyadic numbers. 
 * 
 *  The numerator is divided inplace
 *  by the other denominators and the exponents are substracted.
 *  This corresponds to a euclidean division.
 */
/*
Dyadic&
Dyadic::operator/=(const Dyadic& r)
{
  assert( r.num_ != 0 );
  exp_ -= r.exp_;
  num_ /= r.num_;
  normalize();
  return *this;
}
*/ 

/*! \brief Inplace division of a Dyadic number by an machine integer @n@. 
 *  The number @n@ is supposed to be a power of two.
 */
inline
Dyadic&
Dyadic::operator/= (const numerator_type& n)
{
  assert( n != 0 );
  numerator_type d=n;
  while( (d%2) ==0 ) {
    --exp_;
    d/=2;
  }
  num_/=d;
  return *this;
}

/*! Inplace addition of a Dyadic number with an integer. */
inline
Dyadic&
Dyadic::operator+= (const numerator_type& r)
{
  if(exp_<0) {
    num_ += shift_2(r,exp_);
  }
  else if(exp_>0) {
    num_ = shift_2(num_,exp_);
    num_+=r;
  }
  else {
    num_+=r;
  }
  normalize();
  return *this;
}

/*! Inplace subtraction of a Dyadic number with an integer. */
inline
Dyadic&
Dyadic::operator-= (const numerator_type& r)
{
  if(exp_<0) {
    num_ -= shift_2(r,exp_);
  }
  else if(exp_>0) {
    num_ = shift_2(num_,exp_);
    num_-= r;
  }
  else {
    num_-= r;
  }
  normalize();
  return *this;
}

/*! Inplace multiplication of a Dyadic number with an integer. */
inline
Dyadic&
Dyadic::operator*= (const numerator_type& r)
{
  num_ *= r;
  normalize();
  return *this;
}

inline
int 
sign(const Dyadic::numerator_type& n) {
  if(n>0) { return +1; }
  else if(n<0) { return -1; }
  else { return 0; }
}

// Help function. Do not use it :-)
inline
int
quotient_cmp(const Dyadic& a, const Dyadic& b)
{
  int asign = sign(a.num_);
  int bsign = sign(b.num_);
  
  if ((asign == 0) && (bsign == 0)) { return 0; }
  if (asign < bsign) { return -1; }
  if (asign > bsign) { return +1; }

  // now a and b have the same sign
  Dyadic::numerator_type ta=a.num_;
  Dyadic::numerator_type tb=b.num_;
  Dyadic::exponent_type d=a.exp_-b.exp_;

  if(d>0) {
    ta = shift_2(ta,d);
  }
  else {
    tb = shift_2(tb,-d);
  }

  if(ta<tb) { return -1; }
  else if(tb<ta) { return +1; }
  else { return 0; }
}

/*! Returns the sign of y-x. */
inline
int
compare(const Dyadic& x, const Dyadic& y)
{ 
  return quotient_cmp(x, y); 
}


/*! Stream insertion operator. */
inline
std::ostream&
operator<<(std::ostream& os, const Dyadic& x)
{
  os << x.numerator();
  os<< "*2^(" << x.valuation() << ")";
  return os;
}

//--------------------------------------------------------------------
// Arithmetic with Dyadics
//--------------------------------------------------------------------

/*! \brief Addition of dyadic numbers. */
inline
Dyadic
operator+(const Dyadic& x, const Dyadic& y)
{
  return Dyadic(x)+=y;
}

/*! Substraction of dyadic numbers. */
inline
Dyadic
operator-(const Dyadic& x, const Dyadic& y)
{ 
  return Dyadic(x)-=y; 
}


/*! Multiplication of dyadic numbers. */
inline
Dyadic
operator*(const Dyadic& x, const Dyadic& y)
{
  return Dyadic(x)*=y;
}

/*! Division of dyadic numbers. */
/*
Dyadic
operator/(const Dyadic& x, const Dyadic& y)
{
  return Dyadic(x)/=y;
}
*/

/*! Division of a dyadic number by a power of 2. */
inline
Dyadic
operator/(const Dyadic& x, int d)
{
  return Dyadic(x)/=d;
}

/*! Negation of a dyadic number. */
inline
Dyadic
operator-(const Dyadic& x)
{ 
  return Dyadic(-x.numerator(),x.denominator_exponent()); 
}

/*! Addition of an integer and a dyadic number. */
inline
Dyadic
operator+(const Dyadic::numerator_type& x, const Dyadic& y)
{
  return Dyadic(y)+=x;
}

/*! Addition of a dyadic number and an integer. */
inline
Dyadic
operator+(const Dyadic& x, const Dyadic::numerator_type& y)
{
  return Dyadic(x)+=y;
}

/*! Subtraction of an integer and a dyadic number. */
inline
Dyadic
operator-(const Dyadic::numerator_type& x, const Dyadic& y)
{
  return Dyadic(y)-=x;
}

/*! Subtraction of a dyadic number and an integer. */
inline
Dyadic
operator-(const Dyadic& x, const Dyadic::numerator_type& y)
{
  return Dyadic(x)+=y;
}

/*! Multiplication of an integer and a dyadic number. */
inline
Dyadic
operator*(const Dyadic::numerator_type& x, const Dyadic& y)
{
  return Dyadic(y)*=x;
}
/*! Multiplication of a dyadic number and an integer. */
inline
Dyadic
operator*(const Dyadic& x, const Dyadic::numerator_type& y)
{
  return Dyadic(x)+=y;
}

//--------------------------------------------------------------------
// Boolean operations
//--------------------------------------------------------------------

/// 
inline
bool
operator==(const Dyadic& x, const Dyadic& y)
{ 
  return shift_2(x.num_ ,y.exp_) == shift_2(y.num_,x.exp_); 
}

///
inline
bool operator==(const Dyadic& x, const Dyadic::numerator_type& y)
{ 
  return shift_2(y, x.exp_) == x.num_; 
}

///
inline
bool
operator==(const Dyadic::numerator_type& x, const Dyadic& y)
{ 
  return y == x; 
}

///
inline
bool
operator!=(const Dyadic& x, const Dyadic& y)
{ 
  return ! (x == y); 
}

///
inline
bool
operator!=(const Dyadic& x, const Dyadic::numerator_type& y)
{ 
  return ! (x == y); 
}

///
inline
bool
operator!=(const Dyadic::numerator_type& x, const Dyadic& y)
{ 
  return ! (x == y); 
}

///
inline
bool
operator<(const Dyadic& x, const Dyadic& y)
{
  return quotient_cmp(x,y) == -1;
}

///
inline
bool
operator<(const Dyadic& x, const Dyadic::numerator_type& y)
{
  return quotient_cmp(x,Dyadic(y)) == -1;
}

///
inline
bool
operator<(const Dyadic::numerator_type& x, const Dyadic& y)
{
  return quotient_cmp(Dyadic(x),y) == -1;
}

///
inline
bool
operator>(const Dyadic& x, const Dyadic& y)
{ 
  return(y < x); 
}

///
inline
bool
operator>(const Dyadic& x, const Dyadic::numerator_type& y)
{ 
  return y < x; 
}

///
inline
bool
operator>(const Dyadic::numerator_type& x, const Dyadic& y)
{ 
  return y < x; 
}

///
inline
bool
operator<=(const Dyadic& x, const Dyadic& y)
{ 
  return ! (y < x); 
}

///
inline
bool
operator<=(const Dyadic& x, const Dyadic::numerator_type& y)
{ 
  return ! (y < x); 
}


inline
bool
operator<=(const Dyadic::numerator_type& x, const Dyadic& y)
{ 
  return ! (y < x); 
}

///
inline
bool
operator>=(const Dyadic& x, const Dyadic& y)
{
  return ! (x < y); 
}

///
inline
bool
operator>=(const Dyadic& x, const Dyadic::numerator_type& y)
{ 
  return ! (x < y); 
}

///
inline
bool
operator>=(const Dyadic::numerator_type& x, const Dyadic& y)
{ 
  return ! (x < y); 
}



/*! The denominator of a dyadic number. */
inline
Dyadic::numerator_type denominator(const Dyadic& q)
{ 
  return q.denominator(); 
}

/*! The numerator of a dyadic number. */
inline
const Dyadic::numerator_type&
numerator(const Dyadic& q)
{ 
  return q.numerator(); 
}

} // namespace Synaps
} // namespace Ariadne

#endif // SYNAPS_DYADIC_H
