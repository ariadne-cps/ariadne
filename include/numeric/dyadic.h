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

/*! \file dyadic.h
 *  \brief Dyadic numbers.
 */
 
#ifndef SYNAPS_DYADIC_H
#define SYNAPS_DYADIC_H

#include <gmpxx.h>
#include <iostream>


namespace Ariadne { 
namespace Synaps {

/*! \brief This class provides arithmetic operations on numbers 
 *  which are an integer of type \a Integer, 
 *   times a power of 2, which is stored by as its exponent.
 */
class dyadic
{
 public:
  typedef mpz_class numerator_type;
  typedef int exponent_type;
 public:
  dyadic() : num_(0), exp_(0) { }
  dyadic(const double& x) { mpq_class q=x; num_=q.get_num(); exp_=-log2_floor_(q.get_den()); normalize_(); }

  dyadic(const numerator_type& n)
    : num_(n), exp_(0) { normalize_(); }
  dyadic(const numerator_type& n, exponent_type d)
    : num_(n), exp_(d) { normalize_(); }

/* These definitions conflict with the template function denominator
#if __GNU_MP_VERSION_MINOR < 2 
  template<class T> dyadic(const __gmp_expr<__gmpz_value, T>& n) { 
    *this=dyadic(mpz_class(n)); }
#else 
  template<class T> dyadic(const __gmp_expr<T, T>& n) { 
    *this=dyadic(mpz_class(n)); }
#endif
*/
  
  dyadic& operator+= (const dyadic& r);
  dyadic& operator-= (const dyadic& r);
  dyadic& operator*= (const dyadic& r);
  dyadic& operator+= (const int& r);
  dyadic& operator-= (const int& r);
  dyadic& operator*= (const int& r);
  dyadic& operator/= (const int& n);
  dyadic& operator+= (const numerator_type& r);
  dyadic& operator-= (const numerator_type& r);
  dyadic& operator*= (const numerator_type& r);
  dyadic& operator/= (const numerator_type& n);
  dyadic& operator+= (const double& r);
  dyadic& operator-= (const double& r);
  dyadic& operator*= (const double& r);
  dyadic& operator/= (const double& r);
   
  exponent_type precision() const { return log2_floor_(num_); }

  const numerator_type& numerator() const { return num_; }
  numerator_type  denominator() const { return 2^exp_; }

  exponent_type exponent() const { return exp_; }
  dyadic mantissa() const { return dyadic(num_,exp_-log2_floor_(num_)); }

  operator mpq_class() const { return mpq_class(numerator(),denominator()); }
  double get_d() const { return mpq_class(numerator(),denominator()).get_d(); }

  friend int quotient_cmp(const dyadic&, const dyadic&);
  friend bool operator==(const dyadic&, const dyadic&);
  friend bool operator==(const dyadic&, const double&);
  friend bool operator==(const dyadic&, const mpz_class&);
  friend bool operator==(const dyadic&, const int&);
  friend std::ostream& operator<<(std::ostream&, const dyadic&);
  friend std::istream& operator>>(std::istream&, dyadic&);
 private:  
  void normalize_();
  static int log2_floor_(mpz_class n);
 private:  
  numerator_type num_;
  exponent_type exp_;
};

inline
int 
dyadic::log2_floor_(mpz_class n)
{
  if(n==0) { return -1; }
  if(n<0) { n=-n; }
  int result=0;
  while(n!=1) {
    ++result;
    n=n/2;
  }
  return result;
}


inline
void
dyadic::normalize_()
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
}


inline
dyadic::numerator_type 
shift_2(const dyadic::numerator_type& z, int d) {
  if(d>0) { return z*(2^d); }
  else { return z/(2^d); }
}

//--------------------------------------------------------------------
// Arithmetic operators
//--------------------------------------------------------------------

/*! \brief Inplace addition of dyadic numbers. */
inline
dyadic&
dyadic::operator+= (const dyadic& r)
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
    normalize_();
  }
  return *this;
}


/*! \brief Inplace subtraction of dyadic numbers. */
inline
dyadic&
dyadic::operator-= (const dyadic& r)
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
    normalize_();
  }
  return *this;
}

/*! \brief Inplace multiplication of dyadic numbers. */
inline
dyadic&
dyadic::operator*= (const dyadic& r)
{
  num_ *= r.num_;
  exp_ += r.exp_;
  return *this;
}


/*! \brief Inplace addition of a dyadic number with an integer. */
inline
dyadic&
dyadic::operator+= (const numerator_type& r)
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
  normalize_();
  return *this;
}

/*! \brief Inplace subtraction of a dyadic number with an integer. */
inline
dyadic&
dyadic::operator-= (const numerator_type& r)
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
  normalize_();
  return *this;
}

/*! \brief Inplace multiplication of a dyadic number with an integer. */
inline
dyadic&
dyadic::operator*= (const numerator_type& r)
{
  num_ *= r;
  normalize_();
  return *this;
}

/*! \brief Inplace division of a dyadic number by an integer \a n. 
 *  The number \a n must be a power of two.
 */
inline
dyadic&
dyadic::operator/= (const numerator_type& n)
{
  assert( n != 0 );
  numerator_type d=n;
  while( (d%2)==0 ) {
    --exp_;
    d/=2;
  }
  assert(d==1);
  return *this;
}

/*! \brief Inplace addition of a dyadic number with a machine integer \a n. */
inline
dyadic&
dyadic::operator+= (const int& n)
{
  return this->operator+=(numerator_type(n));
}

/*! \brief Inplace subtraction of a dyadic number with a machine integer \a n. */
inline
dyadic&
dyadic::operator-= (const int& n)
{
  return this->operator-=(numerator_type(n));
}

/*! \brief Inplace multiplication of a dyadic number with a machine integer \a n. */
inline
dyadic&
dyadic::operator*= (const int& n)
{
  return this->operator*=(numerator_type(n));
}

/*! \brief Inplace division of a dyadic number by a machine integer \a n
 *  The number \a n is supposed to be a power of two.
 */
inline
dyadic&
dyadic::operator/= (const int& n)
{
  return this->operator/=(numerator_type(n));
}


/*! \brief Inplace addition of a dyadic number with a double. */
inline
dyadic&
dyadic::operator+= (const double& r)
{
  return this->operator+=(dyadic(r));
}

/*! \brief Inplace subtraction of a dyadic number with a double. */
inline
dyadic&
dyadic::operator-= (const double& r)
{
  return this->operator-=(dyadic(r));
}

/*! \brief Inplace multiplication of a dyadic number with a double. */
inline
dyadic&
dyadic::operator*= (const double& r)
{
  return this->operator*=(dyadic(r));
}


inline
int 
sign(const dyadic::numerator_type& n) {
  if(n>0) { return +1; }
  else if(n<0) { return -1; }
  else { return 0; }
}


// Help function. Do not use it :-)
inline
int
quotient_cmp(const dyadic& a, const dyadic& b)
{
  int asign = sign(a.num_);
  int bsign = sign(b.num_);
  
  if ((asign == 0) && (bsign == 0)) { return 0; }
  if (asign < bsign) { return -1; }
  if (asign > bsign) { return +1; }

  // now a and b have the same sign
  dyadic::numerator_type ta=a.num_;
  dyadic::numerator_type tb=b.num_;
  dyadic::exponent_type d=a.exp_-b.exp_;

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

/*! \brief Returns the sign of y-x. */
inline
int
compare(const dyadic& x, const dyadic& y)
{ 
  return quotient_cmp(x, y); 
}


/*! \brief Stream insertion operator. */
inline
std::ostream&
operator<<(std::ostream& os, const dyadic& x)
{
  os << x.num_;
  os << "*2^" << x.exp_ << "";
  return os;
}

/*! \brief Stream extraction operator. */
inline
std::istream&
operator>>(std::istream& is, dyadic& x)
{
  double y;
  is >> y;
  x=dyadic(y);
  return is;
}

//--------------------------------------------------------------------
// Arithmetic with dyadics
//--------------------------------------------------------------------

/*! \brief Addition of dyadic numbers. */
inline
dyadic
operator+(const dyadic& x, const dyadic& y)
{
  return dyadic(x)+=y;
}

/*! \brief Substraction of dyadic numbers. */
inline
dyadic
operator-(const dyadic& x, const dyadic& y)
{ 
  return dyadic(x)-=y; 
}


/*! \brief Multiplication of dyadic numbers. */
inline
dyadic
operator*(const dyadic& x, const dyadic& y)
{
  return dyadic(x)*=y;
}

/*! \brief Division of dyadic numbers, which yields a rational. */
inline
mpq_class
operator/(const dyadic& x, const dyadic& y)
{
  return mpq_class(x.numerator()*y.denominator(),x.denominator()*y.numerator());
}

/*! \brief Negation of a dyadic number. */
inline
dyadic
operator-(const dyadic& x)
{ 
  return dyadic(-x.numerator(),x.exponent()); 
}

//--------------------------------------------------------------------
// Arithmetic with dyadics and mpz_class
//--------------------------------------------------------------------

/*! \brief Addition of an integer and a dyadic number. */
inline
dyadic
operator+(const mpz_class& n, const dyadic& x)
{
  return dyadic(n)+=x;
}

/*! \brief Addition of a dyadic number and an integer. */
inline
dyadic
operator+(const dyadic& x, const mpz_class& n)
{
  return dyadic(n)+=x;
}

/*! \brief Subtraction of an integer and a dyadic number. */
inline
dyadic
operator-(const mpz_class& n, const dyadic& x)
{
  return dyadic(n)-=x;
}

/*! \brief Subtraction of a dyadic number and an integer. */
inline
dyadic
operator-(const dyadic& x, const mpz_class& n)
{
  return dyadic(n)-=x;
}

/*! \brief Multiplication of an integer and a dyadic number. */
inline
dyadic
operator*(const mpz_class& n, const dyadic& x)
{
  return dyadic(n)*=x;
}

/*! \brief Multiplication of a dyadic number and an integer. */
inline
dyadic
operator*(const dyadic& x, const mpz_class& n)
{
  return dyadic(n)*=x;
}

/*! \brief Division of a dyadic number by an integer. */
inline
dyadic
operator/(const dyadic& x, const mpz_class& n)
{
  return dyadic(x)/=n;
}

//--------------------------------------------------------------------
// Arithmetic with dyadics and ints
//--------------------------------------------------------------------

inline
dyadic
operator+(const dyadic& x, const int& n)
{ 
  return dyadic(x)+=n;
}

inline
dyadic
operator+(const int& n, const dyadic& x)
{ 
  return dyadic(x)+=n;
}

inline
dyadic
operator-(const dyadic& x, const int& n)
{ 
  return dyadic(x)-=n;
}

inline
dyadic
operator-(const int& n, const dyadic& x)
{ 
  return dyadic(x)-=n;
}

inline
dyadic
operator*(const dyadic& x, const int& n)
{ 
  return dyadic(x)*=n;
}

inline
dyadic
operator*(const int& n, const dyadic& x)
{ 
  return dyadic(x)*=n;
}

inline
dyadic
operator/(const dyadic& x, const int& n)
{ 
  return dyadic(x)/=n;
}

//--------------------------------------------------------------------
// Arithmetic with dyadics and unsigned ints
//--------------------------------------------------------------------

inline
dyadic
operator+(const dyadic& x, const unsigned int& n)
{ 
  return dyadic(x)+=int(n);
}

inline
dyadic
operator+(const unsigned int& n, const dyadic& x)
{ 
  return dyadic(x)+=int(n);
}

inline
dyadic
operator-(const dyadic& x, const unsigned int& n)
{ 
  return dyadic(x)-=int(n);
}

inline
dyadic
operator-(const unsigned int& n, const dyadic& x)
{ 
  return dyadic(x)-=int(n);
}

inline
dyadic
operator*(const dyadic& x, const unsigned int& n)
{ 
  return dyadic(x)*=int(n);
}

inline
dyadic
operator*(const unsigned int& n, const dyadic& x)
{ 
  return dyadic(x)*=int(n);
}

inline
dyadic
operator/(const dyadic& x, const unsigned int& n)
{ 
  return dyadic(x)/=int(n);
}


//--------------------------------------------------------------------
// Arithmetic with dyadics and doubles
//--------------------------------------------------------------------

inline
dyadic
operator+(const dyadic& x, const double& y)
{ 
  return dyadic(y)+=x;
}

inline
dyadic
operator+(const double& x, const dyadic& y)
{ 
  return dyadic(x)+=y;
}

inline
dyadic
operator-(const dyadic& x, const double& y)
{ 
  return dyadic(x)-=y;
}

inline
dyadic
operator-(const double& x, const dyadic& y)
{ 
  return dyadic(x)-=y;
}

inline
dyadic
operator*(const dyadic& x, const double& y)
{ 
  return x*dyadic(y);
}

inline
dyadic
operator*(const double& x, const dyadic& y)
{ 
  return dyadic(y)*x;
}

//--------------------------------------------------------------------
// Arithmetic with dyadics and rationals
//--------------------------------------------------------------------

///
inline
mpq_class
operator+(const dyadic& x, const mpq_class& y)
{ 
  return mpq_class(x)+y;
}

///
inline
mpq_class
operator+(const mpq_class& x, const dyadic& y)
{ 
  return x+mpq_class(y);
}

///
inline
mpq_class
operator-(const dyadic& x, const mpq_class& y)
{ 
  return mpq_class(x)-y;
}

///
inline
mpq_class
operator-(const mpq_class& x, const dyadic& y)
{ 
  return x-mpq_class(y);
}

///
inline
mpq_class
operator*(const dyadic& x, const mpq_class& y)
{ 
  return mpq_class(x)*y;
}

///
inline
mpq_class
operator*(const mpq_class& x, const dyadic& y)
{ 
  return x*mpq_class(y);
}

inline
mpq_class
operator/(const dyadic& x, const mpq_class& y)
{ 
  return mpq_class(x)/y;
}

///
inline
mpq_class
operator/(const mpq_class& x, const dyadic& y)
{ 
  return x/mpq_class(y);
}

//--------------------------------------------------------------------
// Boolean operations
//--------------------------------------------------------------------

/// 
inline
bool
operator==(const dyadic& x, const dyadic& y)
{ 
  return shift_2(x.num_ ,y.exp_) == shift_2(y.num_,x.exp_); 
}

///
inline
bool operator==(const dyadic& x, const double& y)
{ 
  return mpq_class(x)==mpq_class(y);
}

///
inline
bool operator==(const double& x, const dyadic& y)
{ 
  return y==x;
}

///
inline
bool operator==(const dyadic& x, const int& n)
{ 
  return shift_2(n, x.exp_) == x.num_; 
}

///
inline
bool
operator==(const int& n, const dyadic& x)
{ 
  return x == n; 
}

///
inline
bool operator==(const dyadic& x, const mpz_class& n)
{ 
  return shift_2(n, x.exp_) == x.num_; 
}

///
inline
bool
operator==(const mpz_class& n, const dyadic& x)
{ 
  return x == n; 
}

///
inline
bool
operator!=(const dyadic& x, const dyadic& y)
{ 
  return ! (x == y); 
}

///
inline
bool
operator!=(const dyadic& x, const double& y)
{ 
  return ! (x == y); 
}

///
inline
bool
operator!=(const double& x, const dyadic& y)
{ 
  return ! (x == y); 
}

///
inline
bool
operator!=(const dyadic& x, const mpz_class& n)
{ 
  return ! (x == n); 
}

///
inline
bool
operator!=(const mpz_class& n, const dyadic& x)
{ 
  return ! (x == n); 
}

///
inline
bool
operator!=(const dyadic& x, const int& n)
{ 
  return ! (x == n); 
}

///
inline
bool
operator!=(const int& n, const dyadic& x)
{ 
  return ! (x == n); 
}

///
inline
bool
operator<(const dyadic& x, const dyadic& y)
{
  return quotient_cmp(x,y) == -1;
}

///
inline
bool
operator<(const dyadic& x, const mpz_class& y)
{
  return quotient_cmp(x,dyadic(y)) == -1;
}

///
inline
bool
operator<(const mpz_class& x, const dyadic& y)
{
  return quotient_cmp(dyadic(x),y) == -1;
}

inline
bool
operator<(const dyadic& x, const int& y)
{
  return quotient_cmp(x,dyadic(y)) == -1;
}

///
inline
bool
operator<(const int& x, const dyadic& y)
{
  return quotient_cmp(dyadic(x),y) == -1;
}

///
inline
bool
operator<(const dyadic& x, const double& y)
{
  return x<dyadic(y);
}

///
inline
bool
operator<(const double& x, const dyadic& y)
{
  return dyadic(x)<y;
}

///
inline
bool
operator>(const dyadic& x, const dyadic& y)
{ 
  return(y < x); 
}

///
inline
bool
operator>(const dyadic& x, const mpz_class& y)
{ 
  return y < x; 
}

///
inline
bool
operator>(const mpz_class& x, const dyadic& y)
{ 
  return y < x; 
}

///
inline
bool
operator>(const dyadic& x, const int& y)
{ 
  return y < x; 
}

///
inline
bool
operator>(const int& x, const dyadic& y)
{ 
  return y < x; 
}

///
inline
bool
operator>(const dyadic& x, const double& y)
{ 
  return y < x; 
}

///
inline
bool
operator>(const double& x, const dyadic& y)
{ 
  return y < x; 
}

///
inline
bool
operator<=(const dyadic& x, const dyadic& y)
{ 
  return ! (y < x); 
}

///
inline
bool
operator<=(const dyadic& x, const mpz_class& y)
{ 
  return ! (y < x); 
}


inline
bool
operator<=(const mpz_class& x, const dyadic& y)
{ 
  return ! (y < x); 
}

inline
bool
operator<=(const dyadic& x, const int& y)
{ 
  return ! (y < x); 
}


inline
bool
operator<=(const int& x, const dyadic& y)
{ 
  return ! (y < x); 
}

inline
bool
operator<=(const dyadic& x, const double& y)
{ 
  return ! (y < x); 
}


inline
bool
operator<=(const double& x, const dyadic& y)
{ 
  return ! (y < x); 
}

///
inline
bool
operator>=(const dyadic& x, const dyadic& y)
{
  return ! (x < y); 
}

///
inline
bool
operator>=(const dyadic& x, const mpz_class& y)
{ 
  return ! (x < y); 
}

///
inline
bool
operator>=(const mpz_class& x, const dyadic& y)
{ 
  return ! (x < y); 
}

///
inline
bool
operator>=(const dyadic& x, const int& y)
{ 
  return ! (x < y); 
}

///
inline
bool
operator>=(const int& x, const dyadic& y)
{ 
  return ! (x < y); 
}

///
inline
bool
operator>=(const dyadic& x, const double& y)
{ 
  return ! (x < y); 
}

///
inline
bool
operator>=(const double& x, const dyadic& y)
{ 
  return ! (x < y); 
}



} // namespace Synaps
} // namespace Ariadne

#endif // SYNAPS_DYADIC_H
