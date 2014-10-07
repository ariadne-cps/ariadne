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

#include "declarations.h"

#include "tribool.h"
#include "rounding.h"
#include "float.h"
#include "float-exact.h"
#include "float-validated.h"
#include "float-approximate.h"

// Simplifying typedefs for unsigned types
typedef unsigned int uint;

namespace Ariadne {

// Forward declarations
class Float;
class ApproximateFloat;
class ValidatedFloat;
class ExactFloat;

class Real;

class Integer;
class Rational;
class Dyadic;
class Decimal;

class Interval;
class UnitInterval;
class UpperInterval;

typedef UpperInterval NumericInterval;

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
//! To obtain the lower and upper bounds of an interval, use \c ivl.lower() and \c ivl.upper_raw().
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
    explicit Interval(const Real& x);
    //! \brief Convert from a floating-point number with an exact representation.
    explicit Interval(const ExactFloat& x) : l(x.raw()), u(x.raw()) { }
    //! \brief Convert from a floating-point number with an exact representation.
    explicit Interval(const ValidatedFloat& x) : l(x.lower()), u(x.upper_raw()) { }
    //! \brief Convert from a dyadic number.
    explicit Interval(const Dyadic& x);
    //! \brief Convert from a decimal number.
    explicit Interval(const Decimal& x);

    //! \brief Convert to a floating-point approximation.
    explicit operator Float () const { return half_exact(add_near(l.raw(),u.raw())); }
    //! \brief Convert from a floating-point number with an exact representation.
    explicit operator ValidatedFloat () const { return ValidatedFloat(this->lower(),this->upper()); }

    //! \brief Create from explicitly given lower and upper bounds. Yields the interval \a [lower,upper].
    Interval(double lower, double upper) : l(lower), u(upper) { }
    //! \brief Create from explicitly given lower and upper bounds. Yields the interval \a [lower,upper].
    Interval(const Float& lower, const Float& upper) : l(lower), u(upper) { }
    //! \brief Convert from a floating-point number with an exact representation.
    Interval(const ExactFloat& lower, const ExactFloat& upper) : l(lower.raw()), u(upper.raw()) { }
    //! \brief Construct an over-approximating interval.
    explicit Interval(const LowerFloat& lower, const UpperFloat& upper) : l(lower.raw()), u(upper.raw()) { }
    //! \brief Create from explicitly given lower and upper bounds. Yields the interval \a [lower,upper].
    explicit Interval(const Real& lower, const Real& upper);
#ifdef HAVE_GMPXX_H
    Interval(const Integer& z);
    Interval(const Rational& q);
    Interval& operator=(const Rational& q);
    Interval(const Rational& lower, const Rational& upper);
#endif // HAVE_GMPXX_H

    Interval& operator=(uint m) { l=m; u=m; return *this; }
    Interval& operator=(int n) { l=n; u=n; return *this; }
    Interval& operator=(double c) { l=c; u=c; return *this; }
    Interval& operator=(const Float& x) { l=x; u=x; return *this; }
    Interval& operator=(const Real& x);
    Interval& operator=(const ExactFloat& x) { l=x.raw(); u=x.raw(); return *this; };

    //! \brief The lower bound of the interval.
    const Float& lower_value() const { return l; }
    const Float& lower_raw() const { return l; }
    //! \brief The upper bound of the interval.
    const Float& upper_value() const { return u; }
    const Float& upper_raw() const { return u; }

    //! \brief The lower bound of the interval.
    const ExactFloatType& lower() const { return reinterpret_cast<ExactFloatType const&>(l); }
    //! \brief The upper bound of the interval.
    const ExactFloatType& upper() const { return reinterpret_cast<ExactFloatType const&>(u); }
    //! \brief The midpoint of the interval.
    const ValidatedFloatType midpoint() const { return half(this->lower()+this->upper()); }
    //! \brief The radius of the interval.
    const ValidatedFloatType radius() const { return half(this->upper()-this->lower()); }
    //! \brief An over-approximation to the width of the interval.
    const ErrorFloatType width() const { return ErrorFloatType(sub_up(u,l)); }

    //! \brief An approximation to the midpoint of the interval.
    const ExactFloatType centre() const { return ExactFloatType(half(add_near(l,u))); }
    //! \brief An over-approximation to the distance between the centre and the endpoints.
    const ErrorFloatType error() const { Float c=half(add_near(l,u)); return ErrorFloatType(max(sub_up(u,c),sub_up(c,l))); }

    //! \brief Tests if the interval is empty.
    bool empty() const { return l>u; }
    //! \brief Tests if the interval is a singleton.
    bool singleton() const { return l==u; }

    //! \brief Sets the interval to a "canonical" empty interval \a [1,0].
    void set_empty() { l=+std::numeric_limits< double >::infinity(); u=-std::numeric_limits< double >::infinity(); }
    void set_lower(const ExactFloatType& lower) { l=lower.raw(); }
        // ARIADNE_ASSERT(lower<=this->u);
    void set_upper(const ExactFloatType& upper) { u=upper.raw(); }
        // ARIADNE_ASSERT(this->l<=upper);
    void set(const ExactFloatType& lower, const ExactFloatType& upper) { l=lower.raw(); u=upper.raw(); }
        // ARIADNE_ASSERT(lower<=upper);
    void set_lower(const Float& lower) { l=lower; }
        // ARIADNE_ASSERT(lower<=this->u);
    void set_upper(const Float& upper) { u=upper; }
        // ARIADNE_ASSERT(this->l<=upper);
    void set(const Float& lower, const Float& upper) { l=lower; u=upper; }
        // ARIADNE_ASSERT(lower<=upper);
  public:
    //! \brief Extract a double-precision point approximation to the value represented by the interval.
    double get_d() const { return (this->l.get_d()+this->u.get_d())/2; }
    static uint output_precision;
    static void set_output_precision(uint p) { output_precision=p; }
  private:
    Float l, u;
};

std::ostream& operator<<(std::ostream& os, const Interval& ivl);

inline ValidatedFloatType midpoint(Interval i) {
    return i.midpoint();
}

inline ValidatedFloatType radius(Interval i) {
    return i.radius();
}

inline ErrorFloatType width(Interval i) {
    return i.width();
}

//! \related Interval \brief Test if the intervals are equal (as sets).
inline bool equal(Interval i1, Interval i2) {
    //std::cerr<<"equal(i1,i2) with i1="<<i1<<"; i2="<<i2<<std::endl;
    return i1.lower_raw()==i2.lower_raw() && i1.upper_raw()==i2.upper_raw();
}

//! \related Interval \brief Test if the interval is empty.
inline bool empty(Interval i) {
    return i.lower_raw()>i.upper_raw();
}

//! \related Interval \brief Test if the interval is bounded.
inline bool bounded(Interval i) {
    return i.lower_raw()!=-inf && i.upper_raw()!=+inf;
}

//! \related Interval \brief The intersection of two intervals.
inline Interval intersection(Interval i1, Interval i2) {
    return Interval(max(i1.lower_raw(),i2.lower_raw()),min(i1.upper_raw(),i2.upper_raw()));
}

//! \related Interval \brief The hull of two intervals, equal to the smallest interval containing both as subsets.
inline Interval hull(Interval i1, Interval i2) {
    assert(i1.lower_raw()<=i1.upper_raw() && i2.lower_raw()<=i2.upper_raw());
    return Interval(min(i1.lower_raw(),i2.lower_raw()),max(i1.upper_raw(),i2.upper_raw()));
}

//! \related Interval \brief The hull of an interval and a point, equal to the smallest interval containing both.
inline Interval hull(Interval i1, Float x2) {
    return Interval(min(i1.lower_raw(),x2),max(i1.upper_raw(),x2));
}

//! \related Interval \brief The hull of two points, equal to the smallest interval containing both.
inline Interval hull(Float x1, Float x2) {
    return Interval(min(x1,x2),max(x1,x2));
}


//! \related Interval \brief Test if the interval \a I contains the number \a x.
inline bool contains(Interval i, ExactFloatType x) { return i.lower_raw()<=x.raw() && x.raw()<=i.upper_raw(); }
inline bool contains(Interval i, ValidatedFloatType x) { return i.lower_raw()<=x.lower_raw() && x.upper_raw()<=i.upper_raw(); }
inline bool contains(Interval i, RawFloatType x) { return i.lower_raw()<=x && x<=i.upper_raw(); }

inline bool element(ExactFloatType x, Interval i) { return i.lower_raw()<=x.raw() && x.raw()<=i.upper_raw(); }
inline bool element(ValidatedFloatType x, Interval i) { return i.lower_raw()<=x.lower_raw() && x.upper_raw()<=i.upper_raw(); }
inline bool element(RawFloatType x, Interval i) { return i.lower_raw()<=x && x<=i.upper_raw(); }

//! \related Interval \brief Test if the interval \a I1 is a subset of \a I2.
inline bool subset(Interval i1, Interval i2) { return i1.lower_raw()>=i2.lower_raw() && i1.upper_raw()<=i2.upper_raw(); }
//! \related Interval \brief Test if the interval \a I1 is a superset of \a I2.
inline bool superset(Interval i1, Interval i2) { return i1.lower_raw()<=i2.lower_raw() && i1.upper_raw()>=i2.upper_raw(); }
//! \related Interval \brief Test if the interval \a I1 is disjoint from \a I2. Returns \c false even if the two intervals only have an endpoint in common.
inline bool disjoint(Interval i1, Interval i2) { return i1.lower_raw()>i2.upper_raw() || i1.upper_raw()<i2.lower_raw(); }
//! \related Interval \brief Test if the interval \a I1 intersects \a I2. Returns \c true even if the two intervals only have an endpoint in common.
inline bool intersect(Interval i1, Interval i2) { return i1.lower_raw()<=i2.upper_raw() && i1.upper_raw()>=i2.lower_raw(); }

//! \related Interval \brief Test if the closed interval \a I1 is disjoint from the closed interval \a I2.
//! Returns \c false if the two intervals only have an endpoint in common.
inline bool separated(Interval i1, Interval i2) { return i1.lower_raw()>i2.upper_raw() || i1.upper_raw()<i2.lower_raw(); }
//! \related Interval \brief Test if the interval \a I1 overlaps \a I2.
//! Returns \c false if the two intervals only have an endpoint in common.
//! Returns \c true if one of the intervals is a singleton in the interior of the other.
inline bool overlap(Interval i1, Interval i2) { return i1.lower_raw()<i2.upper_raw() && i1.upper_raw()>i2.lower_raw(); }
//! \related Interval \brief Test if the (closed) interval \a I1 is a subset of the interior of \a I2.
inline bool inside(Interval i1, Interval i2) { return i1.lower_raw()>i2.lower_raw() && i1.upper_raw()<i2.upper_raw(); }
//! \related Interval \brief Test if the interior of the interval \a I1 is a superset of the (closed) interval \a I2.
inline bool covers(Interval i1, Interval i2) { return i1.lower_raw()<i2.lower_raw() && i1.upper_raw()>i2.upper_raw(); }


//! \brief An over-approximation to an interval set.
class UpperInterval {
  public:
    typedef UpperInterval NumericType;

    //! \brief Default constructor yields the singleton zero interval \a [0,0].
    explicit UpperInterval() : l(0.0), u(0.0) { }
    //! \brief Construct a singleton interval.
    // FIXME: Should make explicit, but this interferes with role as a numeric type
    UpperInterval(Float point) : l(point), u(point) { }
    // FIXME: Should make explicit, but this interferes with role as a numeric type
    UpperInterval(int  point) : l(point), u(point) { }
    UpperInterval(double  point) : l(point), u(point) { }
    //! \brief Create from explicitly given lower and upper bounds. Yields the interval \a [lower,upper].
    explicit UpperInterval(Float lower, Float upper) : l(lower), u(upper) { }
    explicit UpperInterval(double lower, double upper) : l(lower), u(upper) { }

    //! \brief Construct an over-approximating interval.
    UpperInterval(LowerFloat lower, UpperFloat upper) : UpperInterval(lower.raw(),upper.raw()) { }
    //! \brief Convert from an exact interval.
    UpperInterval(Interval ivl) : UpperInterval(ivl.lower_raw(),ivl.upper_raw()) { }

    //! \brief Construct a singleton interval.
    UpperInterval(ExactFloat point) : l(point.raw()), u(point.raw()) { }
    UpperInterval(ValidatedFloat point) : l(point.lower_raw()), u(point.upper_raw()) { }

    //! \brief Set the lower bound of the interval.
    void set_lower(LowerFloat lower) { l=lower.raw(); }
    //! \brief Set the upper bound of the interval.
    void set_upper(UpperFloat upper) { u=upper.raw(); }

    //! \brief The lower bound of the interval.
    const Float& lower_raw() const { return l; }
    //! \brief The upper bound of the interval.
    const Float& upper_raw() const { return u; }

    //! \brief The lower bound of the interval.
    const LowerFloat& lower() const { return reinterpret_cast<LowerFloat const&>(l); }
    //! \brief The upper bound of the interval.
    const UpperFloat& upper() const { return reinterpret_cast<UpperFloat const&>(u); }
    //! \brief The midpoint of the interval.
    const ApproximateFloatType midpoint() const { return half(this->lower()+this->upper()); }
    const ExactFloatType centre() const { return make_exact(half(this->lower()+this->upper())); }
    //! \brief The radius of the interval.
    const PositiveUpperFloat radius() const { return half(this->upper()-this->lower()); }
    //! \brief An over-approximation to the width of the interval.
    const PositiveUpperFloat width() const { return this->upper()-this->lower(); }

    explicit operator ValidatedFloatType() const { return ValidatedFloatType(l,u); }

    tribool empty() const { return (l>u) || tribool(indeterminate); }

    friend const ApproximateFloatType midpoint(UpperInterval const& ivl) {
        return ivl.midpoint(); }
    friend bool empty(UpperInterval const& ivl) {
        return (ivl.l>ivl.u); }
    friend bool bounded(UpperInterval const& ivl) {
        return -inf<ivl.l && ivl.u<+inf; }
    friend bool contains(UpperInterval const& ivl, ExactNumber const& x) {
        return ivl.lower_raw() <= x.raw() && x.raw() <= ivl.upper_raw(); }
    friend tribool inside(UpperInterval const& ivl1, Interval const& ivl2) {
        return (ivl1.l>ivl2.lower_raw() && ivl1.u<ivl2.upper_raw()) || tribool(indeterminate); }
    friend tribool subset(UpperInterval const& ivl1, Interval const& ivl2) {
        return (ivl1.l>=ivl2.lower_raw() && ivl1.u<=ivl2.upper_raw()) || tribool(indeterminate); }
    friend bool equal(UpperInterval const& ivl1, UpperInterval const& ivl2) {
        return ivl1.l==ivl2.l && ivl1.u==ivl2.u; }
    friend bool refines(UpperInterval const& ivl1, UpperInterval const& ivl2) {
        return ivl1.l>=ivl2.l && ivl1.u<=ivl2.u; }
    friend bool models(UpperInterval const& ivl1, Interval const& ivl2) {
        return ivl1.l>=ivl2.lower_raw() && ivl1.u<=ivl2.upper_raw(); }
    friend tribool disjoint(UpperInterval const& ivl1, UpperInterval const& ivl2) {
        return (ivl1.u<ivl2.l || ivl2.u<ivl1.l) || tribool(indeterminate); }
    friend UpperInterval hull(UpperInterval const& ivl1, UpperInterval const& ivl2) {
        return UpperInterval(min(ivl1.l,ivl2.l),max(ivl1.u,ivl2.u)); }
    friend UpperInterval intersection(UpperInterval const& ivl1, UpperInterval const& ivl2) {
        return UpperInterval(max(ivl1.l,ivl2.l),min(ivl1.u,ivl2.u)); }
    friend UpperInterval widen(UpperInterval x) {
        return UpperInterval(widen(ValidatedFloat(x.lower_raw(),x.upper_raw()))); }
    friend std::ostream& operator<<(std::ostream& os, UpperInterval const& ivl) {
        return os << Interval(ivl.lower_raw(),ivl.upper_raw()); }
  private:
    Float l, u;
};

const ApproximateFloatType midpoint(UpperInterval const& ivl);
bool bounded(UpperInterval const& ivl);
bool contains(UpperInterval const& ivl, ExactNumber const& x);
bool refines(UpperInterval const& ivl1, UpperInterval const& ivl2);
bool models(UpperInterval const& ivl1, Interval const& ivl2);
tribool inside(UpperInterval const& ivl1, Interval const& ivl2);
tribool subset(UpperInterval const& ivl1, Interval const& ivl2);
tribool disjoint(UpperInterval const& ivl1, UpperInterval const& ivl2);
UpperInterval widen(UpperInterval i);

inline Interval make_exact_interval(UpperInterval ivl) { return Interval(ivl.lower_raw(),ivl.upper_raw()); }

// An interval one ulp wider
//! \related Interval \brief An interval containing the given interval in its interior.
Interval widen(Interval i);
//! \related Interval \brief An interval contained in the interior of the given interval.
Interval narrow(Interval i);

// Over-approximate by an interval with float coefficients
//! \related Interval \brief Over-approximate the interval by one using builtin single-precision floating-point values as endpoints.
Interval trunc(Interval);
Interval trunc(Interval, uint eps);

//! \related Interval \brief The nearest representable number to the midpoint of the interval.
inline Float med(Interval i) { return half_exact(add_near(i.lower_raw(),i.upper_raw())); }
//! \related Interval \brief An over-approximation to the radius of the interval.
inline Float rad(Interval i) { return half_exact(sub_up(i.upper_raw(),i.lower_raw())); }
//! \related Interval \brief An over-approximation to the width of the interval.
inline Float diam(Interval i) { return sub_up(i.upper_raw(),i.lower_raw()); }

//! \related UpperInterval \brief The interval of possible maximum values. Yields the interval between \c i1.upper_raw() and \c i2.upper_raw().
inline UpperInterval max(UpperInterval i1,UpperInterval i2);
//! \related UpperInterval \brief The interval of possible minimum values. Yields the interval between \c i1.lower() and \c i2.lower().
inline UpperInterval min(UpperInterval,UpperInterval);
//! \related UpperInterval \brief The interval of possible absolute values. Yields \f$\{ |x| \mid x\in I\}\f$.
inline UpperInterval abs(UpperInterval);

//! \related UpperInterval \brief Unary plus function. Yields the identity \f$I=\{+x | x\in I\}\f$.
inline UpperInterval pos(UpperInterval i);
//! \related UpperInterval \brief Unary negation function. Yields the exact interval \f$\{-x | x\in I\}\f$.
inline UpperInterval neg(UpperInterval i);
//! \related UpperInterval \brief Unary square function. Yields an over-approximation to \f$\{ x^2 \mid x\in I\}\f$.
//! Note that if \a I contains positive and negative values, \c sqr(I) is tighter than \c I*I .
UpperInterval sqr(UpperInterval i);
//! \related UpperInterval \brief Unary reciprocal function. Yields an over-approximation to \f$\{ 1/x \mid x\in I\}\f$.
//! Yields \f$[-\infty,+\infty]\f$ if \a I contains \a 0 in its interior, and an interval containing \f$[1/u,+\infty]\f$ if \a I=[0,u] .
UpperInterval rec(UpperInterval i);

//! \related UpperInterval \brief Binary addition function. Yields an over-approximation to \f$\{ x_1+x_2 \mid x_1\in I_1 \wedge x_2\in I_2\}\f$.
inline UpperInterval add(UpperInterval, UpperInterval);
//! \related UpperInterval \brief Subtraction function. Yields an over-approximation to \f$\{ x_1-x_2 \mid x_1\in I_1 \wedge x_2\in I_2\}\f$.
inline UpperInterval sub(UpperInterval, UpperInterval);
//! \related UpperInterval \brief Binary multiplication function. Yields an over-approximation to \f$\{ x_1\times x_2 \mid x_1\in I_1 \wedge x_2\in I_2\}\f$.
UpperInterval mul(UpperInterval, UpperInterval);
//! \related UpperInterval \brief Division function. Yields an over-approximation to \f$\{ x_1 \div x_2 \mid x_1\in I_1 \wedge x_2\in I_2\}\f$.
UpperInterval div(UpperInterval, UpperInterval);

inline UpperInterval add(UpperInterval, Float);
inline UpperInterval add(Float, UpperInterval);
inline UpperInterval sub(UpperInterval, Float);
inline UpperInterval sub(Float, UpperInterval);
UpperInterval mul(UpperInterval, Float);
UpperInterval mul(Float,UpperInterval);
UpperInterval div(UpperInterval, Float);
UpperInterval div(Float, UpperInterval);

extern const UpperInterval pi_ivl;

inline UpperInterval neg_ivl(Float);
inline UpperInterval rec_ivl(Float);
inline UpperInterval add_ivl(Float, Float);
inline UpperInterval sub_ivl(Float, Float);
inline UpperInterval mul_ivl(Float, Float);
inline UpperInterval div_ivl(Float, Float);

//! \related UpperInterval \brief Positive integer power function. Yields an over-approximation to \f$\{ x^m \mid x\in I\}\f$.
UpperInterval pow(UpperInterval i, uint m);
//! \related UpperInterval \brief %Integer power function. Yields an over-approximation to \f$\{ x^n \mid x\in I\}\f$.
UpperInterval pow(UpperInterval i, int n);

//! \related UpperInterval \brief Square-root function. Yields an over-approximation to \f$\{ \sqrt{x} \mid x\in I\}\f$.
//! Requires \c I.lower()>=0 .
UpperInterval sqrt(UpperInterval);
//! \related UpperInterval \brief Exponential function. Yields an over-approximation to \f$\{ \exp{x} \mid x\in I\}\f$.
UpperInterval exp(UpperInterval);
//! \related UpperInterval \brief Natural logarithm function. Yields an over-approximation to \f$\{ \log{x} \mid x\in I\}\f$.
//! Requires \c I.lower()>0 .
UpperInterval log(UpperInterval);

//! \related UpperInterval \brief Sine function. Yields an over-approximation to \f$\{ \sin{x} \mid x\in I\}\f$.
UpperInterval sin(UpperInterval);
//! \related UpperInterval \brief Cosine function. Yields an over-approximation to \f$\{ \sin{x} \mid x\in I\}\f$.
UpperInterval cos(UpperInterval);
//! \related UpperInterval \brief Tangent function. Yields an over-approximation to \f$\{ \sin{x} \mid x\in I\}\f$.
UpperInterval tan(UpperInterval);
UpperInterval asin(UpperInterval);
UpperInterval acos(UpperInterval);
UpperInterval atan(UpperInterval);


//! \related UpperInterval \brief The magnitude of the interval \a I. Yields \f$ \max\{ |x|\,\mid\,x\in I \}\f$.
inline Float mag(UpperInterval i) { return max(abs(i.lower_raw()),abs(i.upper_raw())); }
//! \related UpperInterval \brief The mignitude of the interval \a I. Yields \f$ \min\{ |x|\,\mid\,x\in I \}\f$.
inline Float mig(UpperInterval i) { return min(abs(i.lower_raw()),abs(i.upper_raw())); }

inline UpperInterval max(UpperInterval i1, UpperInterval i2)
{
    return UpperInterval(max(i1.lower_raw(),i2.lower_raw()),max(i1.upper_raw(),i2.upper_raw()));
}

inline UpperInterval min(UpperInterval i1, UpperInterval i2)
{
    return UpperInterval(min(i1.lower_raw(),i2.lower_raw()),min(i1.upper_raw(),i2.upper_raw()));
}


inline UpperInterval abs(UpperInterval i)
{
    if(i.lower_raw()>=0) {
        return UpperInterval(i.lower_raw(),i.upper_raw());
    } else if(i.upper_raw()<=0) {
        return UpperInterval(-i.upper_raw(),-i.lower_raw());
    } else {
        return UpperInterval(static_cast<Float>(0.0),max(-i.lower_raw(),i.upper_raw()));
    }
}

inline UpperInterval pos(UpperInterval i)
{
    return UpperInterval(+i.lower_raw(),+i.upper_raw());
}

inline UpperInterval pos_ivl(Float x)
{
    return UpperInterval(+x,+x);
}

inline UpperInterval neg(UpperInterval i)
{
    return UpperInterval(-i.upper_raw(),-i.lower_raw());
}

inline UpperInterval neg_ivl(Float x)
{
    return UpperInterval(-x,-x);
}

inline UpperInterval sqr_ivl(Float x)
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double xv=internal_cast<volatile double&>(x);
    set_rounding_mode(downward);
    volatile double rl=xv*xv;
    set_rounding_mode(upward);
    volatile double ru=xv*xv;
    set_rounding_mode(rnd);
    return UpperInterval(rl,ru);
}

inline UpperInterval rec_ivl(Float x)
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double xv=internal_cast<volatile double&>(x);
    set_rounding_mode(downward);
    volatile double rl=1.0/xv;
    set_rounding_mode(upward);
    volatile double ru=1.0/xv;
    set_rounding_mode(rnd);
    return UpperInterval(rl,ru);
}



inline UpperInterval add(UpperInterval i1, UpperInterval i2)
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double i1l=internal_cast<volatile double&>(i1.lower_raw());
    volatile double i1u=internal_cast<volatile double&>(i1.upper_raw());
    volatile double i2l=internal_cast<volatile double&>(i2.lower_raw());
    volatile double i2u=internal_cast<volatile double&>(i2.upper_raw());
    set_rounding_mode(downward);
    volatile double rl=i1l+i2l;
    set_rounding_mode(upward);
    volatile double ru=i1u+i2u;
    set_rounding_mode(rnd);
    return UpperInterval(rl,ru);
}

inline UpperInterval add(UpperInterval i1, Float x2)
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double i1l=internal_cast<volatile double&>(i1.lower_raw());
    volatile double i1u=internal_cast<volatile double&>(i1.upper_raw());
    volatile double x2v=internal_cast<volatile double&>(x2);
    set_rounding_mode(downward);
    volatile double rl=i1l+x2v;
    set_rounding_mode(upward);
    volatile double ru=i1u+x2v;
    set_rounding_mode(rnd);
    return UpperInterval(rl,ru);
}

inline UpperInterval add(Float x1, UpperInterval i2)
{
    return add(i2,x1);
}

inline UpperInterval add_ivl(Float x1, Float x2)
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double x1v=internal_cast<volatile double&>(x1);
    volatile double x2v=internal_cast<volatile double&>(x2);
    set_rounding_mode(downward);
    volatile double rl=x1v+x2v;
    set_rounding_mode(upward);
    volatile double ru=x1v+x2v;
    set_rounding_mode(rnd);
    return UpperInterval(rl,ru);
}

inline UpperInterval sub(UpperInterval i1, UpperInterval i2)
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double i1l=internal_cast<volatile double&>(i1.lower_raw());
    volatile double i1u=internal_cast<volatile double&>(i1.upper_raw());
    volatile double i2l=internal_cast<volatile double&>(i2.lower_raw());
    volatile double i2u=internal_cast<volatile double&>(i2.upper_raw());
    set_rounding_mode(downward);
    volatile double rl=i1l-i2u;
    set_rounding_mode(upward);
    volatile double ru=i1u-i2l;
    set_rounding_mode(rnd);
    return UpperInterval(rl,ru);
}

inline UpperInterval sub(UpperInterval i1, Float x2)
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double i1l=internal_cast<volatile double&>(i1.lower_raw());
    volatile double i1u=internal_cast<volatile double&>(i1.upper_raw());
    volatile double x2v=internal_cast<volatile double&>(x2);
    set_rounding_mode(downward);
    volatile double rl=i1l-x2v;
    set_rounding_mode(upward);
    volatile double ru=i1u-x2v;
    set_rounding_mode(rnd);
    return UpperInterval(rl,ru);
}

inline UpperInterval sub(Float x1, UpperInterval i2)
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double x1v=internal_cast<volatile double&>(x1);
    volatile double i2l=internal_cast<volatile double&>(i2.lower_raw());
    volatile double i2u=internal_cast<volatile double&>(i2.upper_raw());
    set_rounding_mode(downward);
    volatile double rl=x1v-i2u;
    set_rounding_mode(upward);
    volatile double ru=x1v-i2l;
    set_rounding_mode(rnd);
    return UpperInterval(rl,ru);
}

inline UpperInterval sub_ivl(Float x1, Float x2)
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double x1v=internal_cast<volatile double&>(x1);
    volatile double x2v=internal_cast<volatile double&>(x2);
    set_rounding_mode(downward);
    volatile double rl=x1v-x2v;
    set_rounding_mode(upward);
    volatile double ru=x1v-x2v;
    set_rounding_mode(rnd);
    return UpperInterval(rl,ru);
}

inline UpperInterval mul_ivl(Float x1, Float x2)
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double x1v=internal_cast<volatile double&>(x1);
    volatile double x2v=internal_cast<volatile double&>(x2);
    set_rounding_mode(downward);
    volatile double rl=x1v*x2v;
    set_rounding_mode(upward);
    volatile double ru=x1v*x2v;
    set_rounding_mode(rnd);
    return UpperInterval(rl,ru);
}

inline UpperInterval div_ivl(Float x1, Float x2)
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double x1v=internal_cast<volatile double&>(x1);
    volatile double x2v=internal_cast<volatile double&>(x2);
    set_rounding_mode(downward);
    volatile double rl=x1v/x2v;
    set_rounding_mode(upward);
    volatile double ru=x1v/x2v;
    set_rounding_mode(rnd);
    return UpperInterval(rl,ru);
}

inline UpperInterval pow_ivl(Float x1, int n2)
{
    return pow(UpperInterval(x1),n2);
}

inline UpperInterval med_ivl(Float x1, Float x2)
{
    return add_ivl(half(x1),half(x2));
}

inline UpperInterval rad_ivl(Float x1, Float x2)
{
    return sub_ivl(half(x2),half(x1));
}

inline UpperInterval med_ivl(UpperInterval i) {
    return add_ivl(half(i.lower_raw()),half(i.upper_raw()));
}

inline UpperInterval rad_ivl(UpperInterval i) {
    return sub_ivl(half(i.upper_raw()),half(i.lower_raw()));
}


//! \related UpperInterval \brief Unary plus operator. Should be implemented exactly and yield \f$\{ +x \mid x\in I\}\f$.
inline UpperInterval operator+(const UpperInterval& i) { return UpperInterval(i.lower_raw(),i.upper_raw()); }
//! \related UpperInterval \brief Unary negation operator. Should be implemented exactly and yield \f$\{ -x \mid x\in I\}\f$.
inline UpperInterval operator-(const UpperInterval& i) { return UpperInterval(-i.upper_raw(),-i.lower_raw()); }
//! \related UpperInterval \brief Binary addition operator. Guaranteed to yield an over approximation to \f$\{ x_1+x_2 \mid x_1\in I_1 \wedge x_2\in I_2\}\f$.
inline UpperInterval operator+(const UpperInterval& i1, const UpperInterval& i2) { return add(i1,i2); }
//! \related UpperInterval \brief Binary addition operator. Guaranteed to yield an over approximation to \f$\{ x_1-x_2 \mid x_1\in I_1 \wedge x_2\in I_2\}\f$.
inline UpperInterval operator-(const UpperInterval& i1, const UpperInterval& i2) { return sub(i1,i2); }
//! \related UpperInterval \brief Binary addition operator. Guaranteed to yield an over approximation to \f$\{ x_1*x_2 \mid x_1\in I_1 \wedge x_2\in I_2\}\f$.
inline UpperInterval operator*(const UpperInterval& i1, const UpperInterval& i2) { return mul(i1,i2); }
//! \related UpperInterval \brief Binary addition operator. Guaranteed to yield an over approximation to \f$\{ x_1/x_2 \mid x_1\in I_1 \wedge x_2\in I_2\}\f$. Yields \f$[-\infty,+\infty]\f$ if \f$0\in I_2\f$.
inline UpperInterval operator/(const UpperInterval& i1, const UpperInterval& i2) { return div(i1,i2); };

//! \related UpperInterval \brief Inplace addition operator.
inline UpperInterval& operator+=(UpperInterval& i1, const UpperInterval& i2) { i1=add(i1,i2); return i1; }
//! \related UpperInterval \brief Inplace subtraction operator.
inline UpperInterval& operator-=(UpperInterval& i1, const UpperInterval& i2) { i1=sub(i1,i2); return i1; }
//! \related UpperInterval \brief Inplace multiplication operator.
inline UpperInterval& operator*=(UpperInterval& i1, const UpperInterval& i2) { i1=mul(i1,i2); return i1; }
//! \related UpperInterval \brief Inplace division operator.
inline UpperInterval& operator/=(UpperInterval& i1, const UpperInterval& i2) { i1=div(i1,i2); return i1; }

inline UpperInterval operator+(const UpperInterval& i1, const Float& x2) { return add(i1,x2); }
inline UpperInterval operator-(const UpperInterval& i1, const Float& x2) { return sub(i1,x2); }
inline UpperInterval operator*(const UpperInterval& i1, const Float& x2) { return mul(i1,x2); }
inline UpperInterval operator/(const UpperInterval& i1, const Float& x2) { return div(i1,x2); }
inline UpperInterval operator+(const Float& x1, const UpperInterval& i2) { return add(i2,x1); }
inline UpperInterval operator-(const Float& x1, const UpperInterval& i2) { return sub(x1,i2); }
inline UpperInterval operator*(const Float& x1, const UpperInterval& i2) { return mul(i2,x1); }
inline UpperInterval operator/(const Float& x1, const UpperInterval& i2) { return div(x1,i2); }

inline UpperInterval operator+(const UpperInterval& i1, const ExactFloat& x2) { return add(i1,static_cast<UpperInterval>(x2)); }
inline UpperInterval operator-(const UpperInterval& i1, const ExactFloat& x2) { return sub(i1,static_cast<UpperInterval>(x2)); }
inline UpperInterval operator*(const UpperInterval& i1, const ExactFloat& x2) { return mul(i1,static_cast<UpperInterval>(x2)); }
inline UpperInterval operator/(const UpperInterval& i1, const ExactFloat& x2) { return div(i1,static_cast<UpperInterval>(x2)); }
inline UpperInterval operator+(const ExactFloat& x1, const UpperInterval& i2) { return add(static_cast<UpperInterval>(x1),i2); }
inline UpperInterval operator-(const ExactFloat& x1, const UpperInterval& i2) { return sub(static_cast<UpperInterval>(x1),i2); }
inline UpperInterval operator*(const ExactFloat& x1, const UpperInterval& i2) { return mul(static_cast<UpperInterval>(x1),i2); }
inline UpperInterval operator/(const ExactFloat& x1, const UpperInterval& i2) { return div(static_cast<UpperInterval>(x1),i2); }

inline UpperInterval& operator+=(UpperInterval& i1, const ValidatedFloat& x2) { i1=add(i1,UpperInterval(x2)); return i1; }
inline UpperInterval& operator-=(UpperInterval& i1, const ValidatedFloat& x2) { i1=sub(i1,UpperInterval(x2)); return i1; }
inline UpperInterval& operator*=(UpperInterval& i1, const ValidatedFloat& x2) { i1=mul(i1,UpperInterval(x2)); return i1; }
inline UpperInterval& operator/=(UpperInterval& i1, const ValidatedFloat& x2) { i1=div(i1,UpperInterval(x2)); return i1; }

inline UpperInterval& operator+=(UpperInterval& i1, const Float& x2) { i1=add(i1,x2); return i1; }
inline UpperInterval& operator-=(UpperInterval& i1, const Float& x2) { i1=sub(i1,x2); return i1; }
inline UpperInterval& operator*=(UpperInterval& i1, const Float& x2) { i1=mul(i1,x2); return i1; }
inline UpperInterval& operator/=(UpperInterval& i1, const Float& x2) { i1=div(i1,x2); return i1; }

inline UpperInterval operator+(const UpperInterval& i1, double x2) { return add(i1,static_cast<Float>(x2)); }
inline UpperInterval operator+(double x1, const UpperInterval& i2) { return add(i2,static_cast<Float>(x1)); }
inline UpperInterval operator-(const UpperInterval& i1, double x2) { return sub(i1,static_cast<Float>(x2)); }
inline UpperInterval operator-(double x1, const UpperInterval& i2) { return sub(static_cast<Float>(x1),i2); }
inline UpperInterval operator*(const UpperInterval& i1, double x2) { return mul(i1,static_cast<Float>(x2)); }
inline UpperInterval operator*(double x1, const UpperInterval& i2) { return mul(i2,static_cast<Float>(x1)); }
inline UpperInterval operator/(const UpperInterval& i1, double x2) { return div(i1,static_cast<Float>(x2)); }
inline UpperInterval operator/(double x1, const UpperInterval& i2) { return div(static_cast<Float>(x1),i2); }

inline UpperInterval& operator+=(UpperInterval& i1, double x2) { i1=add(i1,static_cast<Float>(x2)); return i1; }
inline UpperInterval& operator-=(UpperInterval& i1, double x2) { i1=sub(i1,static_cast<Float>(x2)); return i1; }
inline UpperInterval& operator*=(UpperInterval& i1, double x2) { i1=mul(i1,static_cast<Float>(x2)); return i1; }
inline UpperInterval& operator/=(UpperInterval& i1, double x2) { i1=div(i1,static_cast<Float>(x2)); return i1; }

//inline UpperInterval operator/(const UpperInterval& i1, int n2) { return div(i1,Float(n2)); }
//inline UpperInterval operator/(const UpperInterval& i1, double x2) { return div(i1,Float(x2)); }

// Standard equality operators
//! \related UpperInterval \brief Equality operator. Tests equality of intervals as geometric objects, so \c [0,1]==[0,1] returns \c true.
inline bool operator==(const UpperInterval& i1, const UpperInterval& i2) { return i1.lower_raw()==i2.lower_raw() && i1.upper_raw()==i2.upper_raw(); }
//! \related UpperInterval \brief Inequality operator. Tests equality of intervals as geometric objects, so \c [0,2]!=[1,3] returns \c true,
//! even though the intervals possibly represent the same exact real value.
inline bool operator!=(const UpperInterval& i1, const UpperInterval& i2) { return i1.lower_raw()!=i2.lower_raw() || i1.upper_raw()!=i2.upper_raw(); }

// Boost-style tribool (in)equality operators
//inline tribool operator==(const UpperInterval& i1, const UpperInterval& i2) {
//  if(i1.lower_raw()>i2.upper_raw() || i1.upper_raw()<i2.lower_raw()) { return false; } else if(i1.lower_raw()==i2.upper_raw() && i1.upper_raw()==i2.lower_raw()) { return true; } else { return indeterminate; } }
//inline tribool operator!=(const UpperInterval& i1, const UpperInterval& i2) { return !(i1==i2); }

//! \related UpperInterval \brief Equality operator. Tests equality of represented real-point value.
//! Hence \c [0.0,2.0]==1.0 yields \c indeterminate since the interval may represent a real number other than \c 1.0 .
inline tribool operator==(const UpperInterval& i1, const Float& x2) {
    if(i1.upper_raw()<x2 || i1.lower_raw()>x2) { return false; }
    else if(i1.lower_raw()==x2 && i1.upper_raw()==x2) { return true; }
    else { return indeterminate; }
}

//! \related UpperInterval \brief Equality operator. Tests equality of represented real-point value.
//! Hence \c [0.0,2.0]!=1.0 yields \c indeterminate since the interval may represent a real number equal to \c 1.0 .
inline tribool operator!=(const UpperInterval& i1, const Float& x2) {
    if(i1.upper_raw()<x2 || i1.lower_raw()>x2) { return true; }
    else if(i1.lower_raw()==x2 && i1.upper_raw()==x2) { return false; }
    else { return indeterminate; }
}

inline tribool operator> (const UpperInterval& i1, const Float& x2) {
    if(i1.lower_raw()> x2) { return true; }
    else if(i1.upper_raw()<=x2) { return false; }
    else { return indeterminate; }
}

inline tribool operator< (const UpperInterval& i1, const Float& x2) {
    if(i1.upper_raw()< x2) { return true; }
    else if(i1.lower_raw()>=x2) { return false; }
    else { return indeterminate; }
}

inline tribool operator>=(const UpperInterval& i1, const Float& x2) {
    if(i1.lower_raw()>=x2) { return true; }
    else if(i1.upper_raw()< x2) { return false; }
    else { return indeterminate; }
}

inline tribool operator<=(const UpperInterval& i1, const Float& x2) {
    if(i1.upper_raw()<=x2) { return true; }
    else if(i1.lower_raw()> x2) { return false; }
    else { return indeterminate; }
}

inline tribool operator==(const UpperInterval& i1, double x2) { return i1==static_cast<Float>(x2); }
inline tribool operator!=(const UpperInterval& i1, double x2) { return i1!=static_cast<Float>(x2); }
inline tribool operator<=(const UpperInterval& i1, double x2) { return i1<=static_cast<Float>(x2); }
inline tribool operator>=(const UpperInterval& i1, double x2) { return i1>=static_cast<Float>(x2); }
inline tribool operator< (const UpperInterval& i1, double x2) { return i1< static_cast<Float>(x2); }
inline tribool operator> (const UpperInterval& i1, double x2) { return i1> static_cast<Float>(x2); }



//! \related UpperInterval \brief Strict greater-than comparison operator. Tests equality of represented real-point value.
//! Hence \c [1.0,3.0]>[0.0,2.0] yields \c indeterminate since the first interval could represent the number 1.25 and the second 1.75.
inline tribool operator> (UpperInterval i1, UpperInterval i2) {
    if(i1.lower_raw()> i2.upper_raw()) { return true; }
    else if(i1.upper_raw()<=i2.lower_raw()) { return false; }
    else { return indeterminate; }
}

//! \related UpperInterval \brief Strict greater-than comparison operator. Tests equality of represented real-point value.
inline tribool operator< (UpperInterval i1, UpperInterval i2) {
    if(i1.upper_raw()< i2.lower_raw()) { return true; }
    else if(i1.lower_raw()>=i2.upper_raw()) { return false; }
    else { return indeterminate; }
}

//! \related UpperInterval \brief Strict greater-than comparison operator. Tests equality of represented real-point value.
inline tribool operator>=(UpperInterval i1, UpperInterval i2) {
    if(i1.lower_raw()>=i2.upper_raw()) { return true; }
    else if(i1.upper_raw()< i2.lower_raw()) { return false; }
    else { return indeterminate; }
}

//! \related UpperInterval \brief Strict greater-than comparison operator. Tests equality of represented real-point value.
inline tribool operator<=(UpperInterval i1, UpperInterval i2) {
    if(i1.upper_raw()<=i2.lower_raw()) { return true; }
    else if(i1.lower_raw()> i2.upper_raw()) { return false; }
    else { return indeterminate; }
}

#ifdef ARIADNE_ENABLE_SERIALIZATION
  template<class A> void serialize(A& a, Interval& ivl, const uint version) {
    a & ivl.lower_raw() & ivl.upper_raw(); }
#endif

std::ostream& operator<<(std::ostream&, const Interval&);
std::istream& operator>>(std::istream&, Interval&);

inline ValidatedNumberType make_singleton(Interval const& ivl) {
    return ValidatedNumberType(ivl.lower_raw(),ivl.upper_raw());
}

inline ValidatedNumberType make_singleton(UpperInterval const& ivl) {
    return ValidatedNumberType(ivl.lower_raw(),ivl.upper_raw());
}

//! \brief An over-approximation to an interval set.
class ApproximateInterval {
  public:
    explicit ApproximateInterval() : l(0.0), u(0.0) { }
    explicit ApproximateInterval(Float point) : l(point), u(point) { }
    explicit ApproximateInterval(Float lower, Float upper) : l(lower), u(upper) { }
    explicit ApproximateInterval(ApproximateFloat point) : l(point.raw()), u(point.raw()) { }
    explicit ApproximateInterval(ApproximateFloat lower, ApproximateFloat upper) : l(lower.raw()), u(upper.raw()) { }
    ApproximateInterval(Interval ivl) : l(ivl.lower_raw()), u(ivl.upper_raw()) { }
    ApproximateInterval(UpperInterval ivl) : l(ivl.lower_raw()), u(ivl.upper_raw()) { }
    Float const& lower_raw() const { return l; }
    Float const& upper_raw() const { return l; }
    ApproximateFloat lower() const { return ApproximateFloat(l); }
    ApproximateFloat upper() const { return ApproximateFloat(u); }
    ApproximateFloat midpoint() const { return ApproximateFloat((l+u)/2); }
    ApproximateFloat radius() const { return ApproximateFloat((u-l)/2); }
    ApproximateFloat width() const { return ApproximateFloat(u-l); }
    friend bool contains(ApproximateInterval const& ivl, ApproximateFloat const& x) {
        return ivl.lower_raw()<=x.raw() && x.raw()<=ivl.upper_raw(); }
    friend std::ostream& operator<<(std::ostream& os, const ApproximateInterval& ivl) {
        return os << Interval(ivl.lower_raw(),ivl.upper_raw()); }
  private:
    Float l, u;
};

class UnitInterval
    : public Interval
{
  public:
    UnitInterval() : Interval(-1,+1) { }
};

} // namespace Ariadne

#endif
