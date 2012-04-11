/***************************************************************************
 *            interval.h
 *
 *  Copyright 2008-10  Pieter Collins
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

/*! \file interval.h
 *  \brief Interval numeric class.
 */
#ifndef ARIADNE_INTERVAL_H
#define ARIADNE_INTERVAL_H

#ifdef HAVE_GMPXX_H
#include <gmpxx.h>
#endif // HAVE_GMPXX_H

#include <iostream>
#include <cassert>

#include "tribool.h"
#include "rounding.h"
#include "float.h"

// Simplifying typedefs for unsigned types
typedef unsigned int uint;

namespace Ariadne {

// Forward declarations
class Float;
class Interval;
class Real;


//! \ingroup NumericModule
//! \brief Intervals with floating-point endpoints supporting outwardly-rounded arithmetic.
//! \details
//! Note that <c>%Interval(3.3)</c> yields the singleton interval \f$[3.2999999999999998224,3.2999999999999998224]\f$ (the constant is first interpreted by the C++ compiler to give a C++ \c double, whereas <c>%Interval("3.3")</c> yields the interval \f$[3.2999999999999998224,3.3000000000000002665]\f$ enclosing \f$3.3\f$.
//!
//! Comparison tests on \c Interval use the idea that an interval represents a single number with an unknown value.
//! Hence the result is of type \c tribool, which can take values { \c True, \c False, \c Indeterminate }.
//! Hence a test \f$[l_1,u_1]\leq [l_2,u_2]\f$ returns \c True if \f$u_1\leq u_2\f$, since in this case \f$x_1\leq x_2\f$ whenever \f$x_1\in[l_1,u_2]\f$ and \f$x_2\in[l_2,u_2]\f$, \c False if \f$l_1>u_2\f$, since in this case we know \f$x_1>x_2\f$, and \c Indeterminate otherwise, since in this case we can find \f$x_1,x_2\f$ making the result either true or false.
//! In the case of equality, the comparison \f$[l_1,u_1]\f$==\f$[l_2,u_2]\f$ only returns \c True if both intervals are singletons, since otherwise we can find values making the result either true of false.
//!
//! To obtain the lower and upper bounds of an interval, use \c ivl.lower() and \c ivl.upper().
//! To obtain the midpoint and radius, use \c ivl.midpoint() and \c ivl.radius().
//! Alternatives \c midpoint(ivl) and \c radius(ivl) are also provided.
//! Note that \c midpoint and \c radius return approximations to the true midpoint and radius of the interval. If \f$m\f$ and \f$r\f$ are the returned midpoint and radius of the interval \f$[l,u]\f$, the using exact arithmetic, we guarentee \f$m-r\leq l\f$ and \f$m+r\geq u\f$
//!
//! To test if an interval contains a point or another interval, use \c encloses(Interval,Float) or \c encloses(Interval,Interval).
//! The test \c refines(Interval,Interval) can also be used.
//! \sa Float
//!
//! \par Python interface
//!
//! In the Python interface, %Ariadne intervals can be constructed from Python literals of the form \c {a:b} or (deprecated) \c [a,b] .
//! The former is preferred, as it cannot be confused with literals for other classes such as Vector and Array types.
//! Automatic conversion is used to convert Interval literals of the form \c {a,b} to an Interval in functions.
//!
//! Care must be taken when defining intervals using floating-point coefficients, since values are first converted to the nearest
//! representable value by the Python interpreter. <br><br>
//! \code
//!   Interval({1.1:2.3}) # Create the interval [1.1000000000000001, 2.2999999999999998]
//!   Interval({2.5:4.25}) # Create the interval [2.5, 4.25], which can be represented exactly
//!   Interval([2.5,4.25]) # Alternative syntax for creating the interval [2.5, 4.25]
//! \endcode
class Interval {
  public:
    typedef Interval NumericType;
  public:
    //! \brief Default constructor yields the singleton zero interval \a [0,0].
    Interval() : l(0.0), u(0.0) { }
    Interval(uint m) : l(m), u(m) { }
    Interval(int n) : l(n), u(n) { }
    //! \brief Convert from a builtin double-precision floating-point value. Yields the singleton interval \a [x,x].
    Interval(double x) : l(x), u(x) { }
    //! \brief Create from a floating-point value. Yields the singleton interval \a [x,x].
    //! Cannot be used in conversions since the \c %Interval class provides stronger accuracy guarantees than the \c %Float class.
    explicit Interval(const Float& x) : l(x), u(x) { }
    //! \brief Copy constructor.
    Interval(const Interval& i) : l(i.l), u(i.u) { }
    //! \brief Convert from a general real number. Yields an interval containing the exact value.
    Interval(const Real& x);
    //! \brief Convert from a floating-point number with an exact representation.
    Interval(const ExactFloat& x) : l(x.value()), u(x.value()) { }

    //! \brief Convert to a floating-point approximation.
    operator Float () const { return this->midpoint(); }

    //! \brief Create from explicitly given lower and upper bounds. Yields the interval \a [lower,upper].
    Interval(double lower, double upper) : l(lower), u(upper) { }
    //! \brief Create from explicitly given lower and upper bounds. Yields the interval \a [lower,upper].
    Interval(const Float& lower, const Float& upper) : l(lower), u(upper) { }
        // ARIADNE_ASSERT_MSG(lower<=upper, "lower = "<<lower<<", upper ="<<upper);
#ifdef HAVE_GMPXX_H
    typedef mpq_class Rational;
    Interval(const Rational& q);
    Interval& operator=(const Rational& q);
    Interval(const Rational& lower, const Rational& upper);
#endif // HAVE_GMPXX_H

    Interval& operator=(uint m) { l=m; u=m; return *this; }
    Interval& operator=(int n) { l=n; u=n; return *this; }
    Interval& operator=(double c) { l=c; u=c; return *this; }
    Interval& operator=(const Float& x) { l=x; u=x; return *this; }
    Interval& operator=(const Real& x);
    Interval& operator=(const ExactFloat& x) { l=x.value(); u=x.value(); return *this; };

    //! \brief The lower bound of the interval.
    const Float& lower() const { return l; }
    //! \brief The upper bound of the interval.
    const Float& upper() const { return u; }
    //! \brief An approximation to the midpoint of the interval.
    const Float midpoint() const { return add_approx(l,u)/2; }
    //! \brief An over-approximation to the radius of the interval.
    const Float radius() const { return sub_up(u,l)/2; }
    //! \brief An over-approximation to the width of the interval.
    const Float width() const { return sub_up(u,l); }

    //! \brief Tests if the interval is empty.
    bool empty() const { return l>u; }
    //! \brief Tests if the interval is a singleton.
    bool singleton() const { return l==u; }

    //! \brief Sets the interval to a "canonical" empty interval \a [1,0].
    void set_empty() { l=+std::numeric_limits< double >::infinity(); u=-std::numeric_limits< double >::infinity(); }
    void set_lower(const Float& lower) { l=lower; }
        // ARIADNE_ASSERT(lower<=this->u);
    void set_upper(const Float& upper) { u=upper; }
        // ARIADNE_ASSERT(this->l<=upper);
    void set(const Float& lower, const Float& upper) { l=lower; u=upper; }
        // ARIADNE_ASSERT(lower<=upper);
  public:
    //! \brief Extract a double-precision point approximation to the value represented by the interval.
    double get_d() const { return (this->l.get_d()+this->u.get_d())/2; }
  private:
    Float l, u;
};

std::ostream& operator<<(std::ostream& os, const Interval& ivl);

inline Float midpoint(Interval i) {
    return add_approx(i.lower(),i.upper())/2;
}

inline Float radius(Interval i) {
    return sub_up(i.upper(),i.lower())/2;
}

inline Float width(Interval i) {
    return sub_up(i.upper(),i.lower());
}

//! \related Interval \brief Test if the intervals are equal (as sets).
inline bool equal(Interval i1, Interval i2) {
    //std::cerr<<"equal(i1,i2) with i1="<<i1<<"; i2="<<i2<<std::endl;
    return i1.lower()==i2.lower() && i1.upper()==i2.upper();
}

//! \related Interval \brief Test if the interval is empty.
inline bool empty(Interval i) {
    return i.lower()>i.upper();
}

//! \related Interval \brief Test if the interval is bounded.
inline bool bounded(Interval i) {
    return i.lower()!=-inf && i.upper()!=+inf;
}

//! \related Interval \brief The intersection of two intervals.
inline Interval intersection(Interval i1, Interval i2) {
    return Interval(max(i1.lower(),i2.lower()),min(i1.upper(),i2.upper()));
}

//! \related Interval \brief The hull of two intervals, equal to the smallest interval containing both as subsets.
inline Interval hull(Interval i1, Interval i2) {
    assert(i1.lower()<=i1.upper() && i2.lower()<=i2.upper());
    return Interval(min(i1.lower(),i2.lower()),max(i1.upper(),i2.upper()));
}

//! \related Interval \brief The hull of an interval and a point, equal to the smallest interval containing both.
inline Interval hull(Interval i1, Float x2) {
    return Interval(min(i1.lower(),x2),max(i1.upper(),x2));
}

// An interval one ulp wider
//! \related Interval \brief An interval containing the given interval in its interior.
Interval widen(Interval i);
//! \related Interval \brief An interval contained in the interior of the given interval.
Interval narrow(Interval i);

// Over-approximate by an interval with float coefficients
//! \related Interval \brief Over-approximate the interval by one using builtin single-precision floating-point values as endpoints.
Interval trunc(Interval);
Interval trunc(Interval, uint eps);

//! \related Interval \brief The midpoint of the interval.
inline Float med(Interval i) { return (i.lower()+i.upper())/2; }
//! \related Interval \brief An over-approximation to the radius of the interval.
inline Float rad(Interval i) { return up((i.upper()-i.lower())/2); }
//! \related Interval \brief An over-approximation to the width of the interval.
inline Float diam(Interval i) { return up(i.upper()-i.lower()); }

//! \related Interval \brief The interval of possible maximum values. Yields the interval between \c i1.upper() and \c i2.upper().
inline Interval max(Interval i1,Interval i2);
//! \related Interval \brief The interval of possible minimum values. Yields the interval between \c i1.lower() and \c i2.lower().
inline Interval min(Interval,Interval);
//! \related Interval \brief The interval of possible absolute values. Yields \f$\{ |x| \mid x\in I\}\f$.
inline Interval abs(Interval);

//! \related Interval \brief Unary plus function. Yields the identity \f$I=\{+x | x\in I\}\f$.
inline Interval pos(Interval i);
//! \related Interval \brief Unary negation function. Yields the exact interval \f$\{-x | x\in I\}\f$.
inline Interval neg(Interval i);
//! \related Interval \brief Unary square function. Yields an over-approximation to \f$\{ x^2 \mid x\in I\}\f$.
//! Note that if \a I contains positive and negative values, \c sqr(I) is tighter than \c I*I .
Interval sqr(Interval i);
//! \related Interval \brief Unary reciprocal function. Yields an over-approximation to \f$\{ 1/x \mid x\in I\}\f$.
//! Yields \f$[-\infty,+\infty]\f$ if \a I contains \a 0 in its interior, and an interval containing \f$[1/u,+\infty]\f$ if \a I=[0,u] .
Interval rec(Interval i);

//! \related Interval \brief Binary addition function. Yields an over-approximation to \f$\{ x_1+x_2 \mid x_1\in I_1 \wedge x_2\in I_2\}\f$.
inline Interval add(Interval, Interval);
//! \related Interval \brief Subtraction function. Yields an over-approximation to \f$\{ x_1-x_2 \mid x_1\in I_1 \wedge x_2\in I_2\}\f$.
inline Interval sub(Interval, Interval);
//! \related Interval \brief Binary multiplication function. Yields an over-approximation to \f$\{ x_1\times x_2 \mid x_1\in I_1 \wedge x_2\in I_2\}\f$.
Interval mul(Interval, Interval);
//! \related Interval \brief Division function. Yields an over-approximation to \f$\{ x_1 \div x_2 \mid x_1\in I_1 \wedge x_2\in I_2\}\f$.
Interval div(Interval, Interval);

inline Interval add(Interval, Float);
inline Interval add(Float, Interval);
inline Interval sub(Interval, Float);
inline Interval sub(Float, Interval);
Interval mul(Interval, Float);
Interval mul(Float,Interval);
Interval div(Interval, Float);
Interval div(Float, Interval);

extern const Interval pi_ivl;

inline Interval neg_ivl(Float);
inline Interval rec_ivl(Float);
inline Interval add_ivl(Float, Float);
inline Interval sub_ivl(Float, Float);
inline Interval mul_ivl(Float, Float);
inline Interval div_ivl(Float, Float);

//! \related Interval \brief Positive integer power function. Yields an over-approximation to \f$\{ x^m \mid x\in I\}\f$.
Interval pow(Interval i, uint m);
//! \related Interval \brief %Integer power function. Yields an over-approximation to \f$\{ x^n \mid x\in I\}\f$.
Interval pow(Interval i, int n);

//! \related Interval \brief Square-root function. Yields an over-approximation to \f$\{ \sqrt{x} \mid x\in I\}\f$.
//! Requires \c I.lower()>=0 .
Interval sqrt(Interval);
//! \related Interval \brief Exponential function. Yields an over-approximation to \f$\{ \exp{x} \mid x\in I\}\f$.
Interval exp(Interval);
//! \related Interval \brief Natural logarithm function. Yields an over-approximation to \f$\{ \log{x} \mid x\in I\}\f$.
//! Requires \c I.lower()>0 .
Interval log(Interval);

//! \related Interval \brief Sine function. Yields an over-approximation to \f$\{ \sin{x} \mid x\in I\}\f$.
Interval sin(Interval);
//! \related Interval \brief Cosine function. Yields an over-approximation to \f$\{ \sin{x} \mid x\in I\}\f$.
Interval cos(Interval);
//! \related Interval \brief Tangent function. Yields an over-approximation to \f$\{ \sin{x} \mid x\in I\}\f$.
Interval tan(Interval);
Interval asin(Interval);
Interval acos(Interval);
Interval atan(Interval);


//! \related Interval \brief The magnitude of the interval \a I. Yields \f$ \max\{ |x|\,\mid\,x\in I \}\f$.
inline Float mag(Interval i) { return max(abs(i.lower()),abs(i.upper())); }
//! \related Interval \brief The mignitude of the interval \a I. Yields \f$ \min\{ |x|\,\mid\,x\in I \}\f$.
inline Float mig(Interval i) { return min(abs(i.lower()),abs(i.upper())); }

//! \related Interval \brief Test if the interval \a I contains the number \a x.
inline bool contains(Interval i, Float x) { return i.lower()<=x && x<=i.upper(); }

//! \related Interval \brief Test if the interval \a I1 is a subset of \a I2.
inline bool subset(Interval i1, Interval i2) { return i1.lower()>=i2.lower() && i1.upper()<=i2.upper(); }
//! \related Interval \brief Test if the interval \a I1 is a superset of \a I2.
inline bool superset(Interval i1, Interval i2) { return i1.lower()<=i2.lower() && i1.upper()>=i2.upper(); }
//! \related Interval \brief Test if the interval \a I1 is disjoint from \a I2. Returns \c false even if the two intervals only have an endpoint in common.
inline bool disjoint(Interval i1, Interval i2) { return i1.lower()>i2.upper() || i1.upper()<i2.lower(); }
//! \related Interval \brief Test if the interval \a I1 intersects \a I2. Returns \c true even if the two intervals only have an endpoint in common.
inline bool intersect(Interval i1, Interval i2) { return i1.lower()<=i2.upper() && i1.upper()>=i2.lower(); }

//! \related Interval \brief Test if the closed interval \a I1 is disjoint from the closed interval \a I2.
//! Returns \c false if the two intervals only have an endpoint in common.
inline bool separated(Interval i1, Interval i2) { return i1.lower()>i2.upper() || i1.upper()<i2.lower(); }
//! \related Interval \brief Test if the interval \a I1 overlaps \a I2.
//! Returns \c false if the two intervals only have an endpoint in common.
//! Returns \c true if one of the intervals is a singleton in the interior of the other.
inline bool overlap(Interval i1, Interval i2) { return i1.lower()<i2.upper() && i1.upper()>i2.lower(); }
//! \related Interval \brief Test if the (closed) interval \a I1 is a subset of the interior of \a I2.
inline bool inside(Interval i1, Interval i2) { return i1.lower()>i2.lower() && i1.upper()<i2.upper(); }
//! \related Interval \brief Test if the interior of the interval \a I1 is a superset of the (closed) interval \a I2.
inline bool covers(Interval i1, Interval i2) { return i1.lower()<i2.lower() && i1.upper()>i2.upper(); }

inline Interval max(Interval i1, Interval i2)
{
    return Interval(max(i1.lower(),i2.lower()),max(i1.upper(),i2.upper()));
}

inline Interval min(Interval i1, Interval i2)
{
    return Interval(min(i1.lower(),i2.lower()),min(i1.upper(),i2.upper()));
}


inline Interval abs(Interval i)
{
    if(i.lower()>=0) {
        return Interval(i.lower(),i.upper());
    } else if(i.upper()<=0) {
        return Interval(-i.upper(),-i.lower());
    } else {
        return Interval(static_cast<Float>(0.0),max(-i.lower(),i.upper()));
    }
}

inline Interval pos(Interval i)
{
    return Interval(+i.lower(),+i.upper());
}

inline Interval pos_ivl(Float x)
{
    return Interval(+x,+x);
}

inline Interval neg(Interval i)
{
    return Interval(-i.upper(),-i.lower());
}

inline Interval neg_ivl(Float x)
{
    return Interval(-x,-x);
}

inline Interval sqr_ivl(Float x)
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double xv=internal_cast<volatile double&>(x);
    set_rounding_mode(downward);
    volatile double rl=xv*xv;
    set_rounding_mode(upward);
    volatile double ru=xv*xv;
    set_rounding_mode(rnd);
    return Interval(rl,ru);
}

inline Interval rec_ivl(Float x)
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double xv=internal_cast<volatile double&>(x);
    set_rounding_mode(downward);
    volatile double rl=1.0/xv;
    set_rounding_mode(upward);
    volatile double ru=1.0/xv;
    set_rounding_mode(rnd);
    return Interval(rl,ru);
}



inline Interval add(Interval i1, Interval i2)
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double i1l=internal_cast<volatile double&>(i1.lower());
    volatile double i1u=internal_cast<volatile double&>(i1.upper());
    volatile double i2l=internal_cast<volatile double&>(i2.lower());
    volatile double i2u=internal_cast<volatile double&>(i2.upper());
    set_rounding_mode(downward);
    volatile double rl=i1l+i2l;
    set_rounding_mode(upward);
    volatile double ru=i1u+i2u;
    set_rounding_mode(rnd);
    return Interval(rl,ru);
}

inline Interval add(Interval i1, Float x2)
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double i1l=internal_cast<volatile double&>(i1.lower());
    volatile double i1u=internal_cast<volatile double&>(i1.upper());
    volatile double x2v=internal_cast<volatile double&>(x2);
    set_rounding_mode(downward);
    volatile double rl=i1l+x2v;
    set_rounding_mode(upward);
    volatile double ru=i1u+x2v;
    set_rounding_mode(rnd);
    return Interval(rl,ru);
}

inline Interval add(Float x1, Interval i2)
{
    return add(i2,x1);
}

inline Interval add_ivl(Float x1, Float x2)
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double x1v=internal_cast<volatile double&>(x1);
    volatile double x2v=internal_cast<volatile double&>(x2);
    set_rounding_mode(downward);
    volatile double rl=x1v+x2v;
    set_rounding_mode(upward);
    volatile double ru=x1v+x2v;
    set_rounding_mode(rnd);
    return Interval(rl,ru);
}

inline Interval sub(Interval i1, Interval i2)
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double i1l=internal_cast<volatile double&>(i1.lower());
    volatile double i1u=internal_cast<volatile double&>(i1.upper());
    volatile double i2l=internal_cast<volatile double&>(i2.lower());
    volatile double i2u=internal_cast<volatile double&>(i2.upper());
    set_rounding_mode(downward);
    volatile double rl=i1l-i2u;
    set_rounding_mode(upward);
    volatile double ru=i1u-i2l;
    set_rounding_mode(rnd);
    return Interval(rl,ru);
}

inline Interval sub(Interval i1, Float x2)
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double i1l=internal_cast<volatile double&>(i1.lower());
    volatile double i1u=internal_cast<volatile double&>(i1.upper());
    volatile double x2v=internal_cast<volatile double&>(x2);
    set_rounding_mode(downward);
    volatile double rl=i1l-x2v;
    set_rounding_mode(upward);
    volatile double ru=i1u-x2v;
    set_rounding_mode(rnd);
    return Interval(rl,ru);
}

inline Interval sub(Float x1, Interval i2)
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double x1v=internal_cast<volatile double&>(x1);
    volatile double i2l=internal_cast<volatile double&>(i2.lower());
    volatile double i2u=internal_cast<volatile double&>(i2.upper());
    set_rounding_mode(downward);
    volatile double rl=x1v-i2u;
    set_rounding_mode(upward);
    volatile double ru=x1v-i2l;
    set_rounding_mode(rnd);
    return Interval(rl,ru);
}

inline Interval sub_ivl(Float x1, Float x2)
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double x1v=internal_cast<volatile double&>(x1);
    volatile double x2v=internal_cast<volatile double&>(x2);
    set_rounding_mode(downward);
    volatile double rl=x1v-x2v;
    set_rounding_mode(upward);
    volatile double ru=x1v-x2v;
    set_rounding_mode(rnd);
    return Interval(rl,ru);
}

inline Interval mul_ivl(Float x1, Float x2)
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double x1v=internal_cast<volatile double&>(x1);
    volatile double x2v=internal_cast<volatile double&>(x2);
    set_rounding_mode(downward);
    volatile double rl=x1v*x2v;
    set_rounding_mode(upward);
    volatile double ru=x1v*x2v;
    set_rounding_mode(rnd);
    return Interval(rl,ru);
}

inline Interval div_ivl(Float x1, Float x2)
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double x1v=internal_cast<volatile double&>(x1);
    volatile double x2v=internal_cast<volatile double&>(x2);
    set_rounding_mode(downward);
    volatile double rl=x1v/x2v;
    set_rounding_mode(upward);
    volatile double ru=x1v/x2v;
    set_rounding_mode(rnd);
    return Interval(rl,ru);
}

inline Interval pow_ivl(Float x1, int n2)
{
    return pow(Interval(x1),n2);
}

inline Interval med_ivl(Float x1, Float x2)
{
    return add_ivl(x1/2,x2/2);
}

inline Interval rad_ivl(Float x1, Float x2)
{
    return sub_ivl(x2/2,x1/2);
}

inline Interval med_ivl(Interval i) {
    return add_ivl(i.lower()/2,i.upper()/2);
}

inline Interval rad_ivl(Interval i) {
    return sub_ivl(i.upper()/2,i.lower()/2);
}


//! \related Interval \brief Unary plus operator. Should be implemented exactly and yield \f$\{ +x \mid x\in I\}\f$.
inline Interval operator+(const Interval& i) { return Interval(i.lower(),i.upper()); }
//! \related Interval \brief Unary negation operator. Should be implemented exactly and yield \f$\{ -x \mid x\in I\}\f$.
inline Interval operator-(const Interval& i) { return Interval(-i.upper(),-i.lower()); }
//! \related Interval \brief Binary addition operator. Guaranteed to yield an over approximation to \f$\{ x_1+x_2 \mid x_1\in I_1 \wedge x_2\in I_2\}\f$.
inline Interval operator+(const Interval& i1, const Interval& i2) { return add(i1,i2); }
//! \related Interval \brief Binary addition operator. Guaranteed to yield an over approximation to \f$\{ x_1-x_2 \mid x_1\in I_1 \wedge x_2\in I_2\}\f$.
inline Interval operator-(const Interval& i1, const Interval& i2) { return sub(i1,i2); }
//! \related Interval \brief Binary addition operator. Guaranteed to yield an over approximation to \f$\{ x_1*x_2 \mid x_1\in I_1 \wedge x_2\in I_2\}\f$.
inline Interval operator*(const Interval& i1, const Interval& i2) { return mul(i1,i2); }
//! \related Interval \brief Binary addition operator. Guaranteed to yield an over approximation to \f$\{ x_1/x_2 \mid x_1\in I_1 \wedge x_2\in I_2\}\f$. Yields \f$[-\infty,+\infty]\f$ if \f$0\in I_2\f$.
inline Interval operator/(const Interval& i1, const Interval& i2) { return div(i1,i2); };

//! \related Interval \brief Inplace addition operator.
inline Interval& operator+=(Interval& i1, const Interval& i2) { i1=add(i1,i2); return i1; }
//! \related Interval \brief Inplace subtraction operator.
inline Interval& operator-=(Interval& i1, const Interval& i2) { i1=sub(i1,i2); return i1; }
//! \related Interval \brief Inplace multiplication operator.
inline Interval& operator*=(Interval& i1, const Interval& i2) { i1=mul(i1,i2); return i1; }
//! \related Interval \brief Inplace division operator.
inline Interval& operator/=(Interval& i1, const Interval& i2) { i1=div(i1,i2); return i1; }

inline Interval operator+(const Interval& i1, const Float& x2) { return add(i1,x2); }
inline Interval operator-(const Interval& i1, const Float& x2) { return sub(i1,x2); }
inline Interval operator*(const Interval& i1, const Float& x2) { return mul(i1,x2); }
inline Interval operator/(const Interval& i1, const Float& x2) { return div(i1,x2); }
inline Interval operator+(const Float& x1, const Interval& i2) { return add(i2,x1); }
inline Interval operator-(const Float& x1, const Interval& i2) { return sub(x1,i2); }
inline Interval operator*(const Float& x1, const Interval& i2) { return mul(i2,x1); }
inline Interval operator/(const Float& x1, const Interval& i2) { return div(x1,i2); }

inline Interval operator+(const ExactFloat& x1, const ExactFloat& x2) { return add_ivl(x1,x2); }
inline Interval operator-(const ExactFloat& x1, const ExactFloat& x2) { return sub_ivl(x1,x2); }
inline Interval operator*(const ExactFloat& x1, const ExactFloat& x2) { return mul_ivl(x1,x2); }
inline Interval operator/(const ExactFloat& x1, const ExactFloat& x2) { return div_ivl(x1,x2); }
inline Interval operator+(const Interval& i1, const ExactFloat& x2) { return add(i1,static_cast<Float>(x2)); }
inline Interval operator-(const Interval& i1, const ExactFloat& x2) { return sub(i1,static_cast<Float>(x2)); }
inline Interval operator*(const Interval& i1, const ExactFloat& x2) { return mul(i1,static_cast<Float>(x2)); }
inline Interval operator/(const Interval& i1, const ExactFloat& x2) { return div(i1,static_cast<Float>(x2)); }
inline Interval operator+(const ExactFloat& x1, const Interval& i2) { return add(static_cast<Float>(x1),i2); }
inline Interval operator-(const ExactFloat& x1, const Interval& i2) { return sub(static_cast<Float>(x1),i2); }
inline Interval operator*(const ExactFloat& x1, const Interval& i2) { return mul(static_cast<Float>(x1),i2); }
inline Interval operator/(const ExactFloat& x1, const Interval& i2) { return div(static_cast<Float>(x1),i2); }

inline Interval& operator+=(Interval& i1, const Float& x2) { i1=add(i1,x2); return i1; }
inline Interval& operator-=(Interval& i1, const Float& x2) { i1=sub(i1,x2); return i1; }
inline Interval& operator*=(Interval& i1, const Float& x2) { i1=mul(i1,x2); return i1; }
inline Interval& operator/=(Interval& i1, const Float& x2) { i1=div(i1,x2); return i1; }

inline Interval operator+(const Interval& i1, double x2) { return add(i1,static_cast<Float>(x2)); }
inline Interval operator+(double x1, const Interval& i2) { return add(i2,static_cast<Float>(x1)); }
inline Interval operator-(const Interval& i1, double x2) { return sub(i1,static_cast<Float>(x2)); }
inline Interval operator-(double x1, const Interval& i2) { return sub(static_cast<Float>(x1),i2); }
inline Interval operator*(const Interval& i1, double x2) { return mul(i1,static_cast<Float>(x2)); }
inline Interval operator*(double x1, const Interval& i2) { return mul(i2,static_cast<Float>(x1)); }
inline Interval operator/(const Interval& i1, double x2) { return div(i1,static_cast<Float>(x2)); }
inline Interval operator/(double x1, const Interval& i2) { return div(static_cast<Float>(x1),i2); }

inline Interval& operator+=(Interval& i1, double x2) { i1=add(i1,static_cast<Float>(x2)); return i1; }
inline Interval& operator-=(Interval& i1, double x2) { i1=sub(i1,static_cast<Float>(x2)); return i1; }
inline Interval& operator*=(Interval& i1, double x2) { i1=mul(i1,static_cast<Float>(x2)); return i1; }
inline Interval& operator/=(Interval& i1, double x2) { i1=div(i1,static_cast<Float>(x2)); return i1; }

//inline Interval operator/(const Interval& i1, int n2) { return div(i1,Float(n2)); }
//inline Interval operator/(const Interval& i1, double x2) { return div(i1,Float(x2)); }

// Standard equality operators
//! \related Interval \brief Equality operator. Tests equality of intervals as geometric objects, so \c [0,1]==[0,1] returns \c true.
inline bool operator==(const Interval& i1, const Interval& i2) { return i1.lower()==i2.lower() && i1.upper()==i2.upper(); }
//! \related Interval \brief Inequality operator. Tests equality of intervals as geometric objects, so \c [0,2]!=[1,3] returns \c true,
//! even though the intervals possibly represent the same exact real value.
inline bool operator!=(const Interval& i1, const Interval& i2) { return i1.lower()!=i2.lower() || i1.upper()!=i2.upper(); }

// Boost-style tribool (in)equality operators
//inline tribool operator==(const Interval& i1, const Interval& i2) {
//  if(i1.lower()>i2.upper() || i1.upper()<i2.lower()) { return false; } else if(i1.lower()==i2.upper() && i1.upper()==i2.lower()) { return true; } else { return indeterminate; } }
//inline tribool operator!=(const Interval& i1, const Interval& i2) { return !(i1==i2); }

//! \related Interval \brief Equality operator. Tests equality of represented real-point value.
//! Hence \c [0.0,2.0]==1.0 yields \c indeterminate since the interval may represent a real number other than \c 1.0 .
inline tribool operator==(const Interval& i1, const Float& x2) {
    if(i1.upper()<x2 || i1.lower()>x2) { return false; }
    else if(i1.lower()==x2 && i1.upper()==x2) { return true; }
    else { return indeterminate; }
}

//! \related Interval \brief Equality operator. Tests equality of represented real-point value.
//! Hence \c [0.0,2.0]!=1.0 yields \c indeterminate since the interval may represent a real number equal to \c 1.0 .
inline tribool operator!=(const Interval& i1, const Float& x2) {
    if(i1.upper()<x2 || i1.lower()>x2) { return true; }
    else if(i1.lower()==x2 && i1.upper()==x2) { return false; }
    else { return indeterminate; }
}

inline tribool operator> (const Interval& i1, const Float& x2) {
    if(i1.lower()> x2) { return true; }
    else if(i1.upper()<=x2) { return false; }
    else { return indeterminate; }
}

inline tribool operator< (const Interval& i1, const Float& x2) {
    if(i1.upper()< x2) { return true; }
    else if(i1.lower()>=x2) { return false; }
    else { return indeterminate; }
}

inline tribool operator>=(const Interval& i1, const Float& x2) {
    if(i1.lower()>=x2) { return true; }
    else if(i1.upper()< x2) { return false; }
    else { return indeterminate; }
}

inline tribool operator<=(const Interval& i1, const Float& x2) {
    if(i1.upper()<=x2) { return true; }
    else if(i1.lower()> x2) { return false; }
    else { return indeterminate; }
}

inline tribool operator==(const Interval& i1, double x2) { return i1==static_cast<Float>(x2); }
inline tribool operator!=(const Interval& i1, double x2) { return i1!=static_cast<Float>(x2); }
inline tribool operator<=(const Interval& i1, double x2) { return i1<=static_cast<Float>(x2); }
inline tribool operator>=(const Interval& i1, double x2) { return i1>=static_cast<Float>(x2); }
inline tribool operator< (const Interval& i1, double x2) { return i1< static_cast<Float>(x2); }
inline tribool operator> (const Interval& i1, double x2) { return i1> static_cast<Float>(x2); }



//! \related Interval \brief Strict greater-than comparison operator. Tests equality of represented real-point value.
//! Hence \c [1.0,3.0]>[0.0,2.0] yields \c indeterminate since the first interval could represent the number 1.25 and the second 1.75.
inline tribool operator> (Interval i1, Interval i2) {
    if(i1.lower()> i2.upper()) { return true; }
    else if(i1.upper()<=i2.lower()) { return false; }
    else { return indeterminate; }
}

//! \related Interval \brief Strict greater-than comparison operator. Tests equality of represented real-point value.
inline tribool operator< (Interval i1, Interval i2) {
    if(i1.upper()< i2.lower()) { return true; }
    else if(i1.lower()>=i2.upper()) { return false; }
    else { return indeterminate; }
}

//! \related Interval \brief Strict greater-than comparison operator. Tests equality of represented real-point value.
inline tribool operator>=(Interval i1, Interval i2) {
    if(i1.lower()>=i2.upper()) { return true; }
    else if(i1.upper()< i2.lower()) { return false; }
    else { return indeterminate; }
}

//! \related Interval \brief Strict greater-than comparison operator. Tests equality of represented real-point value.
inline tribool operator<=(Interval i1, Interval i2) {
    if(i1.upper()<=i2.lower()) { return true; }
    else if(i1.lower()> i2.upper()) { return false; }
    else { return indeterminate; }
}

#ifdef ARIADNE_ENABLE_SERIALIZATION
  template<class A> void serialize(A& a, Interval& ivl, const uint version) {
    a & ivl.lower() & ivl.upper(); }
#endif

std::ostream& operator<<(std::ostream&, const Interval&);
std::istream& operator>>(std::istream&, Interval&);


} // namespace Ariadne

#endif
