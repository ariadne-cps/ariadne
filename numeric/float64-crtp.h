/***************************************************************************
 *            numeric/float64-crtp.h
 *
 *  Copyright 2013-17  Pieter Collins
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

/*! \file numeric/float64-crtp.h
 *  \brief 
 */



#ifndef ARIADNE_FLOAT64_H
#define ARIADNE_FLOAT64_H

#include "utility/module.h"
#include "utility/metaprogramming.h"
#include "numeric/paradigm.h"

namespace Ariadne {

/************ Float64 ********************************************************/

struct ExactTag { };

class Boolean;
class ValidatedKleenean;
class ValidatedSierpinskian;
class Fuzzy;

class Float64Approximation;
class Float64LowerBound;
class Float64UpperBound;
class ValidFloat64;
class Float64Bounds;
class Float64Error;
class Float64Value;

class Rational;

Float64Value operator"" _x(long double);
Float64Error operator"" _e(long double);
Float64LowerBound operator"" _l(long double);
Float64UpperBound operator"" _u(long double);
Float64Approximation operator"" _a(long double);

Void Float64::set_rounding_to_nearest();
Void Float64::set_rounding_downward();
Void Float64::set_rounding_upward();
Void Float64::set_rounding_toward_zero();

template<class N> inline N& operator+=(NumberObject<N>& n1, const NumberObject<N>& n2) {
    n1.upcast()=n1.upcast()+n2.upcast(); return n1.upcast(); }
template<class N> inline N& operator*=(NumberObject<N>& n1, const NumberObject<N>& n2) {
    return n1.upcast()=n1.upcast()*n2.upcast(); return n1.upcast(); }

enum class Sign { NEGATIVE=-1, ZERO=0, POSITIVE=+1 };
enum class Comparison { LESS=-1, EQUAL=0, GREATER=+1 };

class Float64 : public ScalarObject<Float64> {
    volatile double d;
  public:
    typedef Void Paradigm;
    Float64() : d(0.0) { }
    Float64(double dbl) : d(dbl) { }
    explicit operator double() const { return d; }
    Float64 get_flt() const { return *this; }
    double get_d() const { return d; }
    friend Float64 operator+(Float64 x);
    friend Float64 operator-(Float64 x);
    friend Float64 operator+(Float64 x1, Float64 x2);
    friend Float64 operator-(Float64 x1, Float64 x2);
    friend Float64 operator*(Float64 x1, Float64 x2);
    friend Float64 operator/(Float64 x1, Float64 x2);
    friend Float64 operator/(Float64 x1, Int n2);
    friend Float64& operator+=(Float64& x1, Float64 x2);
    friend Float64& operator-=(Float64& x1, Float64 x2);
    friend Float64& operator*=(Float64& x1, Float64 x2);
    friend Float64& operator/=(Float64& x1, Float64 x2);
    friend Float64 max(Float64 x1, Float64 x2);
    friend Float64 neg(Float64 x);
    friend Float64 abs(Float64 x);
    friend Int32 floor(Float64 x);
    friend Int32 ceil(Float64 x);
    friend Sign sng(Float64 x);
    friend Comparison cmp(Float64 x1, Float64 x2);
    friend Bool operator==(Float64 x1, Float64 x2);
    friend Bool operator<=(Float64 x1, Float64 x2);
    friend OutputStream& operator<<(OutputStream& os, Float64 x);
};
inline Float64 operator+(Float64 x) { return Float64{+x.d}; }
inline Float64 operator-(Float64 x) { return Float64{-x.d}; }
inline Float64 operator+(Float64 x1, Float64 x2) { return Float64{x1.d+x2.d}; }
inline Float64 operator-(Float64 x1, Float64 x2) { return Float64{x1.d-x2.d}; }
inline Float64 operator*(Float64 x1, Float64 x2) { return Float64{x1.d*x2.d}; }
inline Float64 operator/(Float64 x1, Float64 x2) { return Float64{x1.d/x2.d}; }
inline Float64 operator/(Float64 x1, Int n2) { return Float64{x1.d/n2}; }
inline Float64& operator+=(Float64& x1, Float64 x2) { x1.d+=x2.d; return x1; }
inline Float64& operator-=(Float64& x1, Float64 x2) { x1.d-=x2.d; return x1; }
inline Float64& operator*=(Float64& x1, Float64 x2) { x1.d*=x2.d; return x1; }
inline Float64& operator/=(Float64& x1, Float64 x2) { x1.d/=x2.d; return x1; }
inline Comparison cmp(Float64 x1, Float64 x2) { return (x1.d==x2.d)?Comparison::EQUAL:(x1.d>x2.d)?Comparison::GREATER:Comparison::LESS; }
Bool isnan(Float64 x);

Float64 abs(Float64 x);
inline Float64 mig(Float64 x) { return abs(x); }
inline Float64 mag(Float64 x) { return abs(x); }

inline Float64 neg(Float64 x) { return Float64{-x.d}; }

inline Bool operator==(Float64 x1, Float64 x2) { return x1.d==x2.d; }
inline Bool operator!=(Float64 x1, Float64 x2) { return !(x1==x2); }
inline Bool operator<=(Float64 x1, Float64 x2) { return x1.d<=x2.d; }
inline Bool operator>=(Float64 x1, Float64 x2) { return (x2<=x1); }
inline Bool operator< (Float64 x1, Float64 x2) { return !(x2<=x1); }
inline Bool operator> (Float64 x1, Float64 x2) { return !(x1<=x2); }

/************ Float64Value ***************************************************/

extern const Float64Value inf;

class Float64Value
    : public NumberObject<Float64Value>
{
    friend class Float64Bounds;
    friend class Float64Approximation;
    volatile double v;
    static Int output_places;
  public:
    static Void set_output_places(Int);
  public:
    typedef ExactTag Paradigm;
    Float64Value();
    template<class N, typename std::enable_if<std::is_integral<N>::value>::type> Float64Value(N);
    Float64Value(Nat);
    Float64Value(Int);
    explicit Float64Value(double);
    explicit Float64Value(Float64);
    operator Float64Bounds () const;
    operator ValidFloat64 () const;
    operator Float64LowerBound () const;
    operator Float64UpperBound () const;
    operator Float64Approximation () const;
    Float64 get_flt() const;
    double get_d() const;
    friend Float64Value operator+(Float64Value);
    friend Float64Value operator-(Float64Value);
    friend Float64Bounds operator+(Float64Value x1, Float64Value x2);
    friend Float64Bounds operator-(Float64Value x1, Float64Value x2);
    friend Float64Bounds operator*(Float64Value x1, Float64Value x2);
    friend Float64Bounds operator/(Float64Value x1, Float64Value x2);
    friend OutputStream& operator<<(OutputStream& os, Float64Value const&);
  private:
    Float64Value(long long int n, std::nullptr_t);
};
Boolean operator==(Float64Value x1, Float64Value x2);
Boolean operator!=(Float64Value x1, Float64Value x2);
Boolean operator<=(Float64Value x1, Float64Value x2);
Boolean operator>=(Float64Value x1, Float64Value x2);
Boolean operator< (Float64Value x1, Float64Value x2);
Boolean operator> (Float64Value x1, Float64Value x2);
Float64Value min(Float64Value x1, Float64Value x2);
Float64Value max(Float64Value x1, Float64Value x2);

class Float64Error : public NumberObject<Float64UpperBound> {
    volatile double e;
    static Int output_places;
  public:
    static Void set_output_places(Int);
  public:
    Float64Error();
    template<class M, EnableIf<IsSame<M,Nat>> = dummy> Float64Error(M m);
    template<class X, EnableIf<IsSame<X,double>> = dummy> explicit Float64Error(X x);
    explicit Float64Error(double);
    explicit Float64Error(Float64);
    explicit operator Float64UpperBound () const;
    Float64 get_flt() const;
    Float64& get_flt();
    double get_d() const;
    explicit operator double() const;
    friend Float64Error operator+(Float64Error);
    friend Float64Error operator+(Float64Error,Float64Error);
    friend Float64Error operator*(Float64Error,Float64Error);
    friend Float64Error pow(Float64Error,Nat);
    friend Float64UpperBound operator/(Float64Error x1, Float64Value x2);
    friend OutputStream& operator<<(OutputStream& os, Float64Error const&);
    friend Fuzzy operator==(Float64Error, Float64);
};

Float64Error max(Float64Error x1, Float64Error x2);

template<class M, typename std::enable_if<std::is_unsigned<M>::value,Int>::type=0> inline Float64Error operator/(Float64Error x1, M m2) {
    return Float64Error(x1.get_d()/m2); }

class ValidFloat64 : public NumberObject<ValidFloat64> {
    volatile double v; volatile double e;
  private:
    ValidFloat64(double d, ExactTag);
  public:
    typedef ValidatedTag Paradigm;
    ValidFloat64();
    template<class N, typename std::enable_if<std::is_integral<N>::value,Int>::type = 0>
        ValidFloat64(N n);
    explicit ValidFloat64(double);
    explicit ValidFloat64(double,double);
    ValidFloat64(Float64Value, Float64Error);
    operator Float64Approximation () const;
    Float64Value value() const;
    Float64Error error() const;
    Float64UpperBound upper() const;
    Float64LowerBound lower() const;
    friend ValidFloat64 operator+(ValidFloat64);
    friend ValidFloat64 operator-(ValidFloat64);
    friend ValidFloat64 operator+(ValidFloat64,ValidFloat64);
    friend ValidFloat64 operator-(ValidFloat64,ValidFloat64);
    friend ValidFloat64 operator*(ValidFloat64,ValidFloat64);
    friend ValidFloat64 operator/(ValidFloat64,ValidFloat64);
    friend ValidFloat64 sqr(ValidFloat64);
    friend ValidFloat64 rec(ValidFloat64 x);
    friend Bool operator==(ValidFloat64,Int);
    friend OutputStream& operator<<(OutputStream& os, ValidFloat64 const&);
};
ValidFloat64 min(ValidFloat64 x1, ValidFloat64 x2);
ValidFloat64 max(ValidFloat64 x1, ValidFloat64 x2);


Bool same(ValidFloat64 x1, ValidFloat64 x2);

template<class N, typename std::enable_if<std::is_integral<N>::value,Int>::type> inline
ValidFloat64::ValidFloat64(N n) : ValidFloat64(double(n),ExactTag())
{ }

class Float64Bounds : public NumberObject<Float64Bounds> {
    friend class Float64Approximation; friend class ValidFloat64;
    volatile double l; volatile double u;
    static Int output_places;
  public:
    static Void set_output_places(Int);
  public:
    typedef ValidatedTag Paradigm;
    Float64Bounds();
    Float64Bounds(Int);
    Float64Bounds(const Rational&);
    explicit Float64Bounds(double);
    explicit Float64Bounds(double,double);
    Float64Bounds(Float64LowerBound, Float64UpperBound);
    operator Float64Approximation () const;
    operator Float64UpperBound () const;
    operator Float64LowerBound () const;
    operator ValidFloat64 () const;
    friend Float64Bounds operator+(Float64Bounds);
    friend Float64Bounds operator-(Float64Bounds);
    friend Float64Bounds operator+(Float64Bounds,Float64Bounds);
    friend Float64Bounds operator-(Float64Bounds,Float64Bounds);
    friend Float64Bounds operator*(Float64Bounds,Float64Bounds);
    friend Float64Bounds operator/(Float64Bounds,Float64Bounds);
    friend Float64Bounds sqr(Float64Bounds);
    friend Float64Bounds max(Float64Bounds,Float64Bounds);
    friend Float64Bounds min(Float64Bounds,Float64Bounds);
    friend Bool operator==(Float64Bounds,Int);
    friend ValidatedKleenean operator> (Float64Bounds,Float64Bounds);
    friend OutputStream& operator<<(OutputStream& os, Float64Bounds const&);

    Float64LowerBound lower() const;
    Float64UpperBound upper() const;
    Float64Error error() const;
    Float64Error width() const;
    Float64Error radius() const;

};
Float64Bounds min(Float64Bounds x1, Float64Bounds x2);
Float64Bounds max(Float64Bounds x1, Float64Bounds x2);
Float64Error mag(Float64Bounds x);
Float64Bounds sqr(Float64Bounds);

Bool same(Float64Bounds x1, Float64Bounds x2);

ValidatedKleenean operator==(Float64Bounds x1, Float64Bounds x2);
ValidatedKleenean operator!=(Float64Bounds x1, Float64Bounds x2);
ValidatedKleenean operator<=(Float64Bounds x1, Float64Bounds x2);
ValidatedKleenean operator>=(Float64Bounds x1, Float64Bounds x2);
ValidatedKleenean operator< (Float64Bounds x1, Float64Bounds x2);
ValidatedKleenean operator> (Float64Bounds x1, Float64Bounds x2);

class Float64LowerBound : public NumberObject<Float64LowerBound> {
    volatile double l;
  public:
    Float64LowerBound();
    Float64LowerBound(Int);
    explicit Float64LowerBound(double);
    operator Float64Approximation () const;
    Float64 get_flt() const;
    friend Float64LowerBound operator+(Float64LowerBound);
    friend Float64UpperBound operator-(Float64LowerBound);
    friend Float64LowerBound operator-(Float64UpperBound);
    friend Float64LowerBound operator+(Float64LowerBound,Float64LowerBound);
    friend Float64LowerBound operator-(Float64LowerBound,Float64UpperBound);
    friend Float64UpperBound operator-(Float64UpperBound,Float64LowerBound);
    friend Float64UpperBound min(Float64UpperBound,Float64UpperBound);
    friend Float64UpperBound max(Float64UpperBound,Float64UpperBound);
    friend OutputStream& operator<<(OutputStream& os, Float64LowerBound const&);
};
Float64LowerBound min(Float64LowerBound,Float64LowerBound);
Float64LowerBound max(Float64LowerBound,Float64LowerBound);

Float64UpperBound rec(Float64LowerBound);

class Float64UpperBound : public NumberObject<Float64UpperBound> {
    volatile double u;
  public:
    Float64UpperBound();
    Float64UpperBound(Int);
    explicit Float64UpperBound(double);
    operator Float64Approximation () const;
    Float64 get_flt() const;
    friend Float64UpperBound operator+(Float64UpperBound);
    friend Float64UpperBound operator-(Float64LowerBound);
    friend Float64LowerBound operator-(Float64UpperBound);
    friend Float64UpperBound operator+(Float64UpperBound,Float64UpperBound);
    friend Float64UpperBound operator-(Float64UpperBound,Float64LowerBound);
    friend Float64LowerBound operator-(Float64LowerBound,Float64UpperBound);
    friend Float64UpperBound operator*(Float64UpperBound,Float64UpperBound);
    friend Float64LowerBound min(Float64LowerBound,Float64LowerBound);
    friend Float64LowerBound max(Float64LowerBound,Float64LowerBound);
    friend OutputStream& operator<<(OutputStream& os, Float64UpperBound const&);
};
Float64UpperBound min(Float64UpperBound,Float64UpperBound);
Float64UpperBound max(Float64UpperBound,Float64UpperBound);
Float64LowerBound rec(Float64UpperBound);

Float64UpperBound operator+(Float64UpperBound,Float64Error);
Float64LowerBound operator-(Float64LowerBound,Float64Error);

class Float64Approximation : public NumberObject<Float64Approximation> {
    volatile double a;
  public:
    typedef ApproximateTag Paradigm;
    Float64Approximation();
    Float64Approximation(double);
    Float64Approximation(const Rational&);
    Float64 get_flt() const;
    double get_d() const;
    friend Float64Approximation operator+(Float64Approximation);
    friend Float64Approximation operator-(Float64Approximation);
    friend Float64Approximation operator+(Float64Approximation,Float64Approximation);
    friend Float64Approximation operator-(Float64Approximation,Float64Approximation);
    friend Float64Approximation operator*(Float64Approximation,Float64Approximation);
    friend Float64Approximation operator/(Float64Approximation,Float64Approximation);
    friend Float64Approximation& operator+=(Float64Approximation&,Float64Approximation);
    friend Float64Approximation neg(Float64Approximation);
    friend Float64Approximation sqr(Float64Approximation);
    friend Float64Approximation min(Float64Approximation,Float64Approximation);
    friend Float64Approximation max(Float64Approximation,Float64Approximation);
    friend Fuzzy operator==(Float64Approximation,Float64Approximation);
    friend OutputStream& operator<<(OutputStream& os, Float64Approximation const&);
};

Fuzzy operator==(Float64Approximation,Float64Approximation);
Fuzzy operator!=(Float64Approximation,Float64Approximation);
Fuzzy operator<=(Float64Approximation,Float64Approximation);
Fuzzy operator>=(Float64Approximation,Float64Approximation);
Fuzzy operator< (Float64Approximation,Float64Approximation);
Fuzzy operator> (Float64Approximation,Float64Approximation);
Float64Approximation neg(Float64Approximation);
Float64Approximation abs(Float64Approximation);
Float64Approximation sqr(Float64Approximation);


class TwoExp {
    Int _n;
  public:
    TwoExp(Int n) : _n(n) { }
    double get_d() const { if(_n>=0) { return 1<<_n; } else { return 1.0/(1<<(-_n)); } }
    operator Float64Value () const { return Float64Value(this->get_d()); }
    operator Float64Error () const { return Float64Error(this->get_d()); }
    operator ValidFloat64 () const { return ValidFloat64(this->get_d()); }
    operator Float64Bounds () const { return Float64Bounds(this->get_d()); }
};
inline TwoExp two_exp(Int n) { return TwoExp(n); }


Float64Value cast_exact(const Float64&);
Float64Error make_error(const Float64&);
ValidFloat64 make_valid(const Float64&);
Float64 cast_raw(const Float64Approximation&);

} // namespace Ariadne

#endif
