/***************************************************************************
 *            numeric/float64-crtp.h
 *
 *  Copyright 2013-14  Pieter Collins
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
class Tribool;
class Sierpinski;
class Fuzzy;

class ApprxFloat64;
class LowerFloat64;
class UpperFloat64;
class ValidFloat64;
class BoundFloat64;
class ErrorFloat64;
class ExactFloat64;

class Rational;

ExactFloat64 operator"" _x(long double);
ErrorFloat64 operator"" _e(long double);
LowerFloat64 operator"" _l(long double);
UpperFloat64 operator"" _u(long double);
ApprxFloat64 operator"" _a(long double);

Void set_rounding_to_nearest();
Void set_rounding_downward();
Void set_rounding_upward();
Void set_rounding_toward_zero();

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

/************ ExactFloat64 ***************************************************/

extern const ExactFloat64 inf;

class ExactFloat64
    : public NumberObject<ExactFloat64>
{
    friend class BoundFloat64;
    friend class ApprxFloat64;
    volatile double v;
    static Int output_precision;
  public:
    static Void set_output_precision(Int);
  public:
    typedef Exact Paradigm;
    ExactFloat64();
    template<class N, typename std::enable_if<std::is_integral<N>::value>::type> ExactFloat64(N);
    ExactFloat64(Nat);
    ExactFloat64(Int);
    explicit ExactFloat64(double);
    explicit ExactFloat64(Float64);
    operator BoundFloat64 () const;
    operator ValidFloat64 () const;
    operator LowerFloat64 () const;
    operator UpperFloat64 () const;
    operator ApprxFloat64 () const;
    Float64 get_flt() const;
    double get_d() const;
    friend ExactFloat64 operator+(ExactFloat64);
    friend ExactFloat64 operator-(ExactFloat64);
    friend BoundFloat64 operator+(ExactFloat64 x1, ExactFloat64 x2);
    friend BoundFloat64 operator-(ExactFloat64 x1, ExactFloat64 x2);
    friend BoundFloat64 operator*(ExactFloat64 x1, ExactFloat64 x2);
    friend BoundFloat64 operator/(ExactFloat64 x1, ExactFloat64 x2);
    friend OutputStream& operator<<(OutputStream& os, ExactFloat64 const&);
  private:
    ExactFloat64(long long int n, std::nullptr_t);
};
Boolean operator==(ExactFloat64 x1, ExactFloat64 x2);
Boolean operator!=(ExactFloat64 x1, ExactFloat64 x2);
Boolean operator<=(ExactFloat64 x1, ExactFloat64 x2);
Boolean operator>=(ExactFloat64 x1, ExactFloat64 x2);
Boolean operator< (ExactFloat64 x1, ExactFloat64 x2);
Boolean operator> (ExactFloat64 x1, ExactFloat64 x2);
ExactFloat64 min(ExactFloat64 x1, ExactFloat64 x2);
ExactFloat64 max(ExactFloat64 x1, ExactFloat64 x2);

class ErrorFloat64 : public NumberObject<UpperFloat64> {
    volatile double e;
    static Int output_precision;
  public:
    static Void set_output_precision(Int);
  public:
    ErrorFloat64();
    template<class M, EnableIf<IsSame<M,Nat>> = dummy> ErrorFloat64(M m);
    template<class X, EnableIf<IsSame<X,double>> = dummy> explicit ErrorFloat64(X x);
    explicit ErrorFloat64(double);
    explicit ErrorFloat64(Float64);
    explicit operator UpperFloat64 () const;
    Float64 get_flt() const;
    Float64& get_flt();
    double get_d() const;
    explicit operator double() const;
    friend ErrorFloat64 operator+(ErrorFloat64);
    friend ErrorFloat64 operator+(ErrorFloat64,ErrorFloat64);
    friend ErrorFloat64 operator*(ErrorFloat64,ErrorFloat64);
    friend ErrorFloat64 pow(ErrorFloat64,Nat);
    friend UpperFloat64 operator/(ErrorFloat64 x1, ExactFloat64 x2);
    friend OutputStream& operator<<(OutputStream& os, ErrorFloat64 const&);
    friend Fuzzy operator==(ErrorFloat64, Float64);
};

ErrorFloat64 max(ErrorFloat64 x1, ErrorFloat64 x2);

template<class M, typename std::enable_if<std::is_unsigned<M>::value,Int>::type=0> inline ErrorFloat64 operator/(ErrorFloat64 x1, M m2) {
    return ErrorFloat64(x1.get_d()/m2); }

class ValidFloat64 : public NumberObject<ValidFloat64> {
    volatile double v; volatile double e;
  private:
    ValidFloat64(double d, ExactTag);
  public:
    typedef Validated Paradigm;
    ValidFloat64();
    template<class N, typename std::enable_if<std::is_integral<N>::value,Int>::type = 0>
        ValidFloat64(N n);
    explicit ValidFloat64(double);
    explicit ValidFloat64(double,double);
    ValidFloat64(ExactFloat64, ErrorFloat64);
    operator ApprxFloat64 () const;
    ExactFloat64 value() const;
    ErrorFloat64 error() const;
    UpperFloat64 upper() const;
    LowerFloat64 lower() const;
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

class BoundFloat64 : public NumberObject<BoundFloat64> {
    friend class ApprxFloat64; friend class ValidFloat64;
    volatile double l; volatile double u;
    static Int output_precision;
  public:
    static Void set_output_precision(Int);
  public:
    typedef Validated Paradigm;
    BoundFloat64();
    BoundFloat64(Int);
    BoundFloat64(const Rational&);
    explicit BoundFloat64(double);
    explicit BoundFloat64(double,double);
    BoundFloat64(LowerFloat64, UpperFloat64);
    operator ApprxFloat64 () const;
    operator UpperFloat64 () const;
    operator LowerFloat64 () const;
    operator ValidFloat64 () const;
    friend BoundFloat64 operator+(BoundFloat64);
    friend BoundFloat64 operator-(BoundFloat64);
    friend BoundFloat64 operator+(BoundFloat64,BoundFloat64);
    friend BoundFloat64 operator-(BoundFloat64,BoundFloat64);
    friend BoundFloat64 operator*(BoundFloat64,BoundFloat64);
    friend BoundFloat64 operator/(BoundFloat64,BoundFloat64);
    friend BoundFloat64 sqr(BoundFloat64);
    friend BoundFloat64 max(BoundFloat64,BoundFloat64);
    friend BoundFloat64 min(BoundFloat64,BoundFloat64);
    friend Bool operator==(BoundFloat64,Int);
    friend Tribool operator> (BoundFloat64,BoundFloat64);
    friend OutputStream& operator<<(OutputStream& os, BoundFloat64 const&);

    LowerFloat64 lower() const;
    UpperFloat64 upper() const;
    ErrorFloat64 error() const;
    ErrorFloat64 width() const;
    ErrorFloat64 radius() const;

};
BoundFloat64 min(BoundFloat64 x1, BoundFloat64 x2);
BoundFloat64 max(BoundFloat64 x1, BoundFloat64 x2);
ErrorFloat64 mag(BoundFloat64 x);
BoundFloat64 sqr(BoundFloat64);

Bool same(BoundFloat64 x1, BoundFloat64 x2);

Tribool operator==(BoundFloat64 x1, BoundFloat64 x2);
Tribool operator!=(BoundFloat64 x1, BoundFloat64 x2);
Tribool operator<=(BoundFloat64 x1, BoundFloat64 x2);
Tribool operator>=(BoundFloat64 x1, BoundFloat64 x2);
Tribool operator< (BoundFloat64 x1, BoundFloat64 x2);
Tribool operator> (BoundFloat64 x1, BoundFloat64 x2);

class LowerFloat64 : public NumberObject<LowerFloat64> {
    volatile double l;
  public:
    LowerFloat64();
    LowerFloat64(Int);
    explicit LowerFloat64(double);
    operator ApprxFloat64 () const;
    Float64 get_flt() const;
    friend LowerFloat64 operator+(LowerFloat64);
    friend UpperFloat64 operator-(LowerFloat64);
    friend LowerFloat64 operator-(UpperFloat64);
    friend LowerFloat64 operator+(LowerFloat64,LowerFloat64);
    friend LowerFloat64 operator-(LowerFloat64,UpperFloat64);
    friend UpperFloat64 operator-(UpperFloat64,LowerFloat64);
    friend UpperFloat64 min(UpperFloat64,UpperFloat64);
    friend UpperFloat64 max(UpperFloat64,UpperFloat64);
    friend OutputStream& operator<<(OutputStream& os, LowerFloat64 const&);
};
LowerFloat64 min(LowerFloat64,LowerFloat64);
LowerFloat64 max(LowerFloat64,LowerFloat64);

UpperFloat64 rec(LowerFloat64);

class UpperFloat64 : public NumberObject<UpperFloat64> {
    volatile double u;
  public:
    UpperFloat64();
    UpperFloat64(Int);
    explicit UpperFloat64(double);
    operator ApprxFloat64 () const;
    Float64 get_flt() const;
    friend UpperFloat64 operator+(UpperFloat64);
    friend UpperFloat64 operator-(LowerFloat64);
    friend LowerFloat64 operator-(UpperFloat64);
    friend UpperFloat64 operator+(UpperFloat64,UpperFloat64);
    friend UpperFloat64 operator-(UpperFloat64,LowerFloat64);
    friend LowerFloat64 operator-(LowerFloat64,UpperFloat64);
    friend UpperFloat64 operator*(UpperFloat64,UpperFloat64);
    friend LowerFloat64 min(LowerFloat64,LowerFloat64);
    friend LowerFloat64 max(LowerFloat64,LowerFloat64);
    friend OutputStream& operator<<(OutputStream& os, UpperFloat64 const&);
};
UpperFloat64 min(UpperFloat64,UpperFloat64);
UpperFloat64 max(UpperFloat64,UpperFloat64);
LowerFloat64 rec(UpperFloat64);

UpperFloat64 operator+(UpperFloat64,ErrorFloat64);
LowerFloat64 operator-(LowerFloat64,ErrorFloat64);

class ApprxFloat64 : public NumberObject<ApprxFloat64> {
    volatile double a;
  public:
    typedef Approximate Paradigm;
    ApprxFloat64();
    ApprxFloat64(double);
    ApprxFloat64(const Rational&);
    Float64 get_flt() const;
    double get_d() const;
    friend ApprxFloat64 operator+(ApprxFloat64);
    friend ApprxFloat64 operator-(ApprxFloat64);
    friend ApprxFloat64 operator+(ApprxFloat64,ApprxFloat64);
    friend ApprxFloat64 operator-(ApprxFloat64,ApprxFloat64);
    friend ApprxFloat64 operator*(ApprxFloat64,ApprxFloat64);
    friend ApprxFloat64 operator/(ApprxFloat64,ApprxFloat64);
    friend ApprxFloat64& operator+=(ApprxFloat64&,ApprxFloat64);
    friend ApprxFloat64 neg(ApprxFloat64);
    friend ApprxFloat64 sqr(ApprxFloat64);
    friend ApprxFloat64 min(ApprxFloat64,ApprxFloat64);
    friend ApprxFloat64 max(ApprxFloat64,ApprxFloat64);
    friend Fuzzy operator==(ApprxFloat64,ApprxFloat64);
    friend OutputStream& operator<<(OutputStream& os, ApprxFloat64 const&);
};

Fuzzy operator==(ApprxFloat64,ApprxFloat64);
Fuzzy operator!=(ApprxFloat64,ApprxFloat64);
Fuzzy operator<=(ApprxFloat64,ApprxFloat64);
Fuzzy operator>=(ApprxFloat64,ApprxFloat64);
Fuzzy operator< (ApprxFloat64,ApprxFloat64);
Fuzzy operator> (ApprxFloat64,ApprxFloat64);
ApprxFloat64 neg(ApprxFloat64);
ApprxFloat64 abs(ApprxFloat64);
ApprxFloat64 sqr(ApprxFloat64);


class TwoExp {
    Int _n;
  public:
    TwoExp(Int n) : _n(n) { }
    double get_d() const { if(_n>=0) { return 1<<_n; } else { return 1.0/(1<<(-_n)); } }
    operator ExactFloat64 () const { return ExactFloat64(this->get_d()); }
    operator ErrorFloat64 () const { return ErrorFloat64(this->get_d()); }
    operator ValidFloat64 () const { return ValidFloat64(this->get_d()); }
    operator BoundFloat64 () const { return BoundFloat64(this->get_d()); }
};
inline TwoExp two_exp(Int n) { return TwoExp(n); }


ExactFloat64 make_exact(const Float64&);
ErrorFloat64 make_error(const Float64&);
ValidFloat64 make_valid(const Float64&);
Float64 make_raw(const ApprxFloat64&);

} // namespace Ariadne

#endif
