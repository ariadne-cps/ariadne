/***************************************************************************
 *            numeric/floatmp.h
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
 *  You should have received _a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

/*! \file numeric/floatmp.h
 *  \brief
 */



#ifndef ARIADNE_FLOATMP_H
#define ARIADNE_FLOATMP_H

#include "paradigm.h"
#include "number.h"
#include "mixins.h"
#include <mpfr.h>

namespace Ariadne {

template<class F> struct NumberObject;
template<class F> struct FltMPObject : NumberObject<F> { };

/************ FloatMP ********************************************************/

//enum class RoundingModeMP : mpfr_rnd_t { NEAREST=MPFR_RNDN, UPWARD=MPFR_RNDU, DOWNWARD=MPFR_RNDD };
enum class RoundingModeMP : Nat { NEAREST=MPFR_RNDN, UPWARD=MPFR_RNDU, DOWNWARD=MPFR_RNDD };

template<class F> struct FltMPExpression { };

template<class N> class NumberObject { };

struct PrecisionMP {
    mpfr_prec_t prec;
    PrecisionMP(mpfr_prec_t pr) : prec(pr) { }
    operator mpfr_prec_t () const { return prec; }
};

struct NoInit { };
struct RawPtr { };

//! \ingroup FltMPSubModule
//! \brief Multiple-precision floating-point numbers.
//! Currently defined as _a wrapper around \c mpfr_t from the MPFE library.
//! Default arithmetic operations are approximate, and comparisons are exact, so this class is \em unsafe.
class FloatMP {
  private: public:
    mpfr_t _mpfr;
    typedef decltype(_mpfr[0]) MpfrReference;
  public:
    typedef Raw Paradigm;
    typedef FloatMP NumericType;
    typedef PrecisionMP PrecisionType;
    typedef mpfr_rnd_t RoundingModeType;
    static Void set_rounding_mode(RoundingModeMP rnd);
    static Void set_default_precision(PrecisionMP prec);
    static PrecisionMP get_default_precision();
    static mpfr_rnd_t current_rounding_mode;
  public:
    ~FloatMP();
    explicit FloatMP();
    explicit FloatMP(NoInit);
    explicit FloatMP(PrecisionMP);
    explicit FloatMP(PrecisionMP, NoInit);
    explicit FloatMP(Int32);
    explicit FloatMP(Int32, PrecisionMP);
    explicit FloatMP(Integer const&, PrecisionMP, RoundingModeType);
    explicit FloatMP(Rational const&, PrecisionMP, RoundingModeType);
    FloatMP(double);
    FloatMP(double, PrecisionMP);
    FloatMP(Float64);
//    FloatMP(const mpfr_t);
    FloatMP(const FloatMP&);
    FloatMP(FloatMP&&);
    template<class N, EnableIf<IsIntegral<N>> =dummy> FloatMP& operator=(N n);
    FloatMP& operator=(const FloatMP&);
    FloatMP& operator=(FloatMP&&);
    explicit operator Rational() const;

    PrecisionMP precision() const;
    Void set_precision(PrecisionMP);

    friend FloatMP next_up(FloatMP x);
    friend FloatMP next_down(FloatMP x);

    friend FloatMP operator+(FloatMP const& x);
    friend FloatMP operator-(FloatMP const& x);
    friend FloatMP operator+(FloatMP const& x1, FloatMP const& x2);
    friend FloatMP operator-(FloatMP const& x1, FloatMP const& x2);
    friend FloatMP operator*(FloatMP const& x1, FloatMP const& x2);
    friend FloatMP operator/(FloatMP const& x1, FloatMP const& x2);
    friend FloatMP& operator+=(FloatMP& x1, FloatMP const& x2);
    friend FloatMP& operator-=(FloatMP& x1, FloatMP const& x2);
    friend FloatMP& operator*=(FloatMP& x1, FloatMP const& x2);
    friend FloatMP& operator/=(FloatMP& x1, FloatMP const& x2);

    //friend Integer floor(FloatMP const& x);
    //friend Integer ceil(FloatMP const& x);
    friend FloatMP floor(FloatMP const& x);
    friend FloatMP ceil(FloatMP const& x);

    friend FloatMP nul(FloatMP const& x);
    friend FloatMP pos(FloatMP const& x);
    friend FloatMP neg(FloatMP const& x);
    friend FloatMP half(FloatMP const& x);
    friend FloatMP mag(FloatMP const& x);
    friend FloatMP abs(FloatMP const& x);
    friend FloatMP min(FloatMP const& x1, FloatMP const& x2);
    friend FloatMP max(FloatMP const& x1, FloatMP const& x2);

    friend FloatMP sqr(FloatMP const& x);
    friend FloatMP rec(FloatMP const& x);
    friend FloatMP add(FloatMP const& x1, FloatMP const& x2);
    friend FloatMP sub(FloatMP const& x1, FloatMP const& x2);
    friend FloatMP mul(FloatMP const& x1, FloatMP const& x2);
    friend FloatMP div(FloatMP const& x1, FloatMP const& x2);
    friend FloatMP pow(FloatMP const& x, Int n);
    friend FloatMP sqrt(FloatMP const& x);
    friend FloatMP exp(FloatMP const& x);
    friend FloatMP log(FloatMP const& x);

    friend FloatMP pos(FloatMP const& x, RoundingModeType);
    friend FloatMP neg(FloatMP const& x, RoundingModeType);
    friend FloatMP rec(FloatMP const& x, RoundingModeType);
    friend FloatMP add(FloatMP const& x1, FloatMP const& x2, RoundingModeType);
    friend FloatMP sub(FloatMP const& x1, FloatMP const& x2, RoundingModeType);
    friend FloatMP mul(FloatMP const& x1, FloatMP const& x2, RoundingModeType);
    friend FloatMP div(FloatMP const& x1, FloatMP const& x2, RoundingModeType);
    friend FloatMP exp(FloatMP const& x, RoundingModeType);

    friend FloatMP pos_exact(FloatMP const& x);
    friend FloatMP neg_exact(FloatMP const& x);
    friend FloatMP abs_exact(FloatMP const& x);
    friend FloatMP half_exact(FloatMP const& x);
    friend FloatMP add_approx(FloatMP const& x1, FloatMP const& x2);
    friend FloatMP add_near(FloatMP const& x1, FloatMP const& x2);
    friend FloatMP add_down(FloatMP const& x1, FloatMP const& x2);
    friend FloatMP add_up(FloatMP const& x1, FloatMP const& x2);
    friend FloatMP sub_approx(FloatMP const& x1, FloatMP const& x2);
    friend FloatMP sub_near(FloatMP const& x1, FloatMP const& x2);
    friend FloatMP sub_down(FloatMP const& x1, FloatMP const& x2);
    friend FloatMP sub_up(FloatMP const& x1, FloatMP const& x2);
    friend FloatMP mul_approx(FloatMP const& x1, FloatMP const& x2);
    friend FloatMP mul_near(FloatMP const& x1, FloatMP const& x2);
    friend FloatMP mul_down(FloatMP const& x1, FloatMP const& x2);
    friend FloatMP mul_up(FloatMP const& x1, FloatMP const& x2);
    friend FloatMP div_approx(FloatMP const& x1, FloatMP const& x2);
    friend FloatMP div_near(FloatMP const& x1, FloatMP const& x2);
    friend FloatMP div_down(FloatMP const& x1, FloatMP const& x2);
    friend FloatMP div_up(FloatMP const& x1, FloatMP const& x2);
    friend FloatMP pow_approx(FloatMP const& x, Int n);
    friend FloatMP pow_down(FloatMP const& x, Int n);
    friend FloatMP pow_up(FloatMP const& x, Int n);

    friend FloatMP sqrt_approx(FloatMP const& x);
    friend FloatMP exp_approx(FloatMP const& x);
    friend FloatMP log_approx(FloatMP const& x);
    friend FloatMP sin_approx(FloatMP const& x);
    friend FloatMP cos_approx(FloatMP const& x);
    friend FloatMP tan_approx(FloatMP const& x);
    friend FloatMP asin_approx(FloatMP const& x);
    friend FloatMP acos_approx(FloatMP const& x);
    friend FloatMP atan_approx(FloatMP const& x);


    friend FloatMP add_rnd(FloatMP const& x1, FloatMP const& x2);
    friend FloatMP sub_rnd(FloatMP const& x1, FloatMP const& x2);
    friend FloatMP mul_rnd(FloatMP const& x1, FloatMP const& x2);
    friend FloatMP div_rnd(FloatMP const& x1, FloatMP const& x2);
    friend FloatMP pow_rnd(FloatMP const& x, Int n);
    friend FloatMP sqrt_rnd(FloatMP const& x);
    friend FloatMP exp_rnd(FloatMP const& x);
    friend FloatMP log_rnd(FloatMP const& x);
    friend FloatMP sin_rnd(FloatMP const& x);
    friend FloatMP cos_rnd(FloatMP const& x);
    friend FloatMP abs(FloatMP const& x, RoundingModeType);

    friend Comparison cmp(FloatMP const& x1, FloatMP const& x2);
    friend Bool operator==(FloatMP const& x1, FloatMP const& x2);
    friend Bool operator!=(FloatMP const& x1, FloatMP const& x2);
    friend Bool operator<=(FloatMP const& x1, FloatMP const& x2);
    friend Bool operator>=(FloatMP const& x1, FloatMP const& x2);
    friend Bool operator< (FloatMP const& x1, FloatMP const& x2);
    friend Bool operator> (FloatMP const& x1, FloatMP const& x2);

    friend Comparison cmp(FloatMP const& x1, double x2);
    friend Bool operator==(FloatMP const& x1, double x2);
    friend Bool operator!=(FloatMP const& x1, double x2);
    friend Bool operator<=(FloatMP const& x1, double x2);
    friend Bool operator>=(FloatMP const& x1, double x2);
    friend Bool operator< (FloatMP const& x1, double x2);
    friend Bool operator> (FloatMP const& x1, double x2);

    template<class FE> FloatMP(const FltMPExpression<FE>&);
    FloatMP& raw() { return *this; }
    FloatMP const& raw() const { return *this; }
    friend OutputStream& operator<<(OutputStream& os, FloatMP const&);
    friend InputStream& operator>>(InputStream& is, FloatMP&);
  public:
    FloatMP& operator=(double d);
    MpfrReference get_mpfr();
    MpfrReference get_mpfr() const;
    double get_d() const;
};

template<class R, class A> inline R integer_cast(const A& a);
template<> inline Int integer_cast(const FloatMP& x) { return static_cast<Int>(x.get_d()); }
template<> inline Nat integer_cast(const FloatMP& x) { return static_cast<Nat>(x.get_d()); }

template<class N, EnableIf<IsIntegral<N>>> inline FloatMP& FloatMP::operator=(N n) { double x=n; return *this=x; }

template<class P> class FloatMPTemplate;
using ApprxFloatMP = FloatMPTemplate<Apprx>;
using LowerFloatMP = FloatMPTemplate<Lower>;
using UpperFloatMP = FloatMPTemplate<Upper>;
using BoundFloatMP = FloatMPTemplate<Bound>;
using MetrcFloatMP = FloatMPTemplate<Metrc>;
using ExactFloatMP = FloatMPTemplate<Exact>;
using ErrorFloatMP = FloatMPTemplate<Error>;

template<> class FloatMPTemplate<Error>
    : public NumberObject<ErrorFloatMP>
{
    FloatMP _e;
  public:
    typedef Error Paradigm;
    FloatMPTemplate<Error>(Nat);
    explicit FloatMPTemplate<Error>(double);
    explicit FloatMPTemplate<Error>(double,PrecisionMP);
    explicit FloatMPTemplate<Error>(FloatMP);
    FloatMP const& get_flt() const;
    friend ErrorFloatMP operator+(ErrorFloatMP);
    friend ErrorFloatMP operator+(ErrorFloatMP x1, ErrorFloatMP x2);
    friend ErrorFloatMP operator*(ErrorFloatMP x1, ErrorFloatMP x2);
    friend OutputStream& operator<<(OutputStream& os, ErrorFloatMP const&);
};

template<> class FloatMPTemplate<Exact>
    : public NumberObject<ExactFloatMP>
{
    friend class FloatMPTemplate<Metrc>;
    friend class FloatMPTemplate<Bound>;
    friend class FloatMPTemplate<Apprx>;
    FloatMP _v;
  public:
    typedef Exact Paradigm;
    template<class N, EnableIf<IsIntegral<N>> = dummy>
        FloatMPTemplate<Exact>(N n) : _v(n) { }
    template<class N, EnableIf<IsIntegral<N>> = dummy>
        FloatMPTemplate<Exact>(N n, PrecisionMP pr) : _v(n,pr) { }
    template<class N, EnableIf<IsIntegral<N>> = dummy>
        FloatMPTemplate<Exact>& operator=(N n) { _v=n; return *this; }
    explicit FloatMPTemplate<Exact>(double);
    explicit FloatMPTemplate<Exact>(FloatMP);
    explicit FloatMPTemplate<Exact>(Integer const&, PrecisionMP);
    operator FloatMPTemplate<Metrc>() const;
    operator FloatMPTemplate<Bound>() const;
    operator FloatMPTemplate<Apprx>() const;
    Void set_precision(PrecisionMP pr);
    PrecisionMP precision() const;
    MetrcFloatMP pm(FloatMPTemplate<Error>) const;
    friend ExactFloatMP operator+(ExactFloatMP);
    friend ExactFloatMP operator-(ExactFloatMP);
    friend BoundFloatMP operator+(ExactFloatMP x1, ExactFloatMP x2);
    friend BoundFloatMP operator-(ExactFloatMP x1, ExactFloatMP x2);
    friend BoundFloatMP operator*(ExactFloatMP x1, ExactFloatMP x2);
    friend BoundFloatMP operator/(ExactFloatMP x1, ExactFloatMP x2);
    friend ExactFloatMP pos(ExactFloatMP x);
    friend ExactFloatMP neg(ExactFloatMP x);
    friend BoundFloatMP rec(ExactFloatMP x);
    friend BoundFloatMP abs(ExactFloatMP x);
    friend BoundFloatMP add(ExactFloatMP x1, ExactFloatMP x2);
    friend BoundFloatMP sub(ExactFloatMP x1, ExactFloatMP x2);
    friend BoundFloatMP mul(ExactFloatMP x1, ExactFloatMP x2);
    friend BoundFloatMP div(ExactFloatMP x1, ExactFloatMP x2);
    friend BoundFloatMP const& min(ExactFloatMP const& x1, ExactFloatMP const& x2);
    friend BoundFloatMP const& max(ExactFloatMP const& x1, ExactFloatMP const& x2);
    FloatMP const& get_flt() const;
    friend OutputStream& operator<<(OutputStream& os, ExactFloatMP const&);
};

template<> class FloatMPTemplate<Metrc> : public NumberObject<MetrcFloatMP>
    , public DeclareAnalyticOperations<MetrcFloatMP>
    , public DeclareOrderedOperations<MetrcFloatMP>
    , public ProvideFieldOperators<MetrcFloatMP>
{
    friend class FloatMPTemplate<Apprx>;
    FloatMP _v; FloatMP _e;
  public:
    typedef Metrc Paradigm;
    FloatMPTemplate<Metrc>();
    template<class N, EnableIf<IsIntegral<N>> =dummy> FloatMPTemplate<Metrc>(N n) : FloatMPTemplate<Metrc>(Int32(n)) { }
    explicit FloatMPTemplate<Metrc>(Integer const&, PrecisionMP);
    explicit FloatMPTemplate<Metrc>(Rational const&, PrecisionMP);
    explicit FloatMPTemplate<Metrc>(Real const&, PrecisionMP);
//    explicit FloatMPTemplate<Metrc>(double);
//    explicit FloatMPTemplate<Metrc>(double,double);
    explicit FloatMPTemplate<Metrc>(FloatMP v);
    explicit FloatMPTemplate<Metrc>(FloatMP v, FloatMP e);
    explicit operator FloatMPTemplate<Bound> () const;
    operator FloatMPTemplate<Upper> () const;
    operator FloatMPTemplate<Lower> () const;
    operator FloatMPTemplate<Apprx> () const;
    ExactFloatMP const& value() const;
    ErrorFloatMP const& error() const;
    Void set_precision(PrecisionMP pr);
    PrecisionMP precision() const;
    double get_d() const;
    friend MetrcFloatMP nul(MetrcFloatMP x);
    friend MetrcFloatMP pos(MetrcFloatMP x);
    friend MetrcFloatMP sqr(MetrcFloatMP x);
    friend MetrcFloatMP neg(MetrcFloatMP x);
    friend MetrcFloatMP rec(MetrcFloatMP x);
    friend MetrcFloatMP add(MetrcFloatMP x1, MetrcFloatMP x2);
    friend MetrcFloatMP sub(MetrcFloatMP x1, MetrcFloatMP x2);
    friend MetrcFloatMP mul(MetrcFloatMP x1, MetrcFloatMP x2);
    friend MetrcFloatMP div(MetrcFloatMP x1, MetrcFloatMP x2);
    friend MetrcFloatMP abs(MetrcFloatMP x);
    friend MetrcFloatMP max(MetrcFloatMP x1, MetrcFloatMP x2);
    friend MetrcFloatMP min(MetrcFloatMP x1, MetrcFloatMP x2);
    friend Bool operator==(MetrcFloatMP,Int);
    friend OutputStream& operator<<(OutputStream& os, MetrcFloatMP const&);
 private:
    friend ApprxFloatMP operator+(ApprxFloatMP,ApprxFloatMP);
    FloatMPTemplate<Metrc>(Int32);
};

template<> class FloatMPTemplate<Bound>
    : public NumberObject<BoundFloatMP>
    , public DeclareAnalyticOperations<BoundFloatMP>
    , public DeclareOrderedOperations<BoundFloatMP>
    , public ProvideFieldOperators<BoundFloatMP>
{
    friend class FloatMPTemplate<Apprx>;
    FloatMP _l; FloatMP _u;
  public:
    typedef Bound Paradigm;
    FloatMPTemplate<Bound>();
    template<class N, EnableIf<IsIntegral<N>> =dummy> FloatMPTemplate<Bound>(N n) : FloatMPTemplate<Bound>(Int32(n)) { }
//    explicit FloatMPTemplate<Metrc>(double);
//    explicit FloatMPTemplate<Metrc>(double,double);
    explicit FloatMPTemplate<Bound>(Integer const&, PrecisionMP);
    explicit FloatMPTemplate<Bound>(Rational const&, PrecisionMP);
    explicit FloatMPTemplate<Bound>(Real const&, PrecisionMP);
    explicit FloatMPTemplate<Bound>(FloatMP);
    explicit FloatMPTemplate<Bound>(FloatMP,FloatMP);
    operator FloatMPTemplate<Metrc> () const;
    operator FloatMPTemplate<Upper> () const;
    operator FloatMPTemplate<Lower> () const;
    operator FloatMPTemplate<Apprx> () const;
    LowerFloatMP const& lower() const;
    UpperFloatMP const& upper() const;
    ExactFloatMP value() const;
    ErrorFloatMP error() const;
    Void set_precision(PrecisionMP pr);
    PrecisionMP precision() const;
    friend BoundFloatMP nul(BoundFloatMP);
    friend BoundFloatMP pos(BoundFloatMP);
    friend BoundFloatMP sqr(BoundFloatMP);
    friend BoundFloatMP neg(BoundFloatMP);
    friend BoundFloatMP rec(BoundFloatMP);
    friend BoundFloatMP add(BoundFloatMP,BoundFloatMP);
    friend BoundFloatMP sub(BoundFloatMP,BoundFloatMP);
    friend BoundFloatMP mul(BoundFloatMP,BoundFloatMP);
    friend BoundFloatMP div(BoundFloatMP,BoundFloatMP);
    friend BoundFloatMP pow(BoundFloatMP,Int);
    friend BoundFloatMP sqrt(BoundFloatMP);
    friend BoundFloatMP exp(BoundFloatMP);
    friend BoundFloatMP log(BoundFloatMP);
    friend BoundFloatMP sin(BoundFloatMP);
    friend BoundFloatMP cos(BoundFloatMP);
    friend BoundFloatMP tan(BoundFloatMP);
    friend BoundFloatMP atan(BoundFloatMP);
    friend BoundFloatMP abs(BoundFloatMP);
    friend BoundFloatMP min(BoundFloatMP,BoundFloatMP);
    friend BoundFloatMP max(BoundFloatMP,BoundFloatMP);
    double get_d() const;
    friend Bool operator==(BoundFloatMP,Int);
    friend OutputStream& operator<<(OutputStream& os, BoundFloatMP const&);
 private:
    friend ApprxFloatMP operator+(ApprxFloatMP,ApprxFloatMP);
    FloatMPTemplate<Bound>(Int32);
};

template<> class FloatMPTemplate<Upper>
    : public NumberObject<UpperFloatMP>
{
    FloatMP _u;
  public:
    typedef Upper Paradigm;
    FloatMPTemplate<Upper>();
    FloatMPTemplate<Upper>(FloatMP);
    FloatMP const& get_flt() const;
    explicit FloatMPTemplate<Upper>(Real const&, PrecisionMP);
    operator FloatMPTemplate<Apprx> () const;
    friend UpperFloatMP operator+(UpperFloatMP);
    friend LowerFloatMP operator-(UpperFloatMP);
    friend UpperFloatMP operator+(UpperFloatMP,UpperFloatMP);
    friend UpperFloatMP operator-(UpperFloatMP,LowerFloatMP);
    friend LowerFloatMP operator-(LowerFloatMP,UpperFloatMP);
    friend UpperFloatMP operator*(UpperFloatMP,UpperFloatMP);
    friend UpperFloatMP operator/(UpperFloatMP,LowerFloatMP);
    friend LowerFloatMP operator/(LowerFloatMP,UpperFloatMP);
    friend UpperFloatMP pos(UpperFloatMP);
    friend LowerFloatMP neg(UpperFloatMP);
    friend UpperFloatMP add(UpperFloatMP,UpperFloatMP);
    friend UpperFloatMP sub(UpperFloatMP,LowerFloatMP);
    friend LowerFloatMP sub(LowerFloatMP,UpperFloatMP);
    friend UpperFloatMP mul(UpperFloatMP,UpperFloatMP);
    friend UpperFloatMP div(UpperFloatMP,LowerFloatMP);
    friend LowerFloatMP div(LowerFloatMP,UpperFloatMP);
    friend OutputStream& operator<<(OutputStream& os, UpperFloatMP const&);
};

template<> class FloatMPTemplate<Lower>
    : public NumberObject<LowerFloatMP>
{
    FloatMP _l;
  public:
    typedef Lower Paradigm;
    FloatMPTemplate<Lower>();
    FloatMPTemplate<Lower>(FloatMP);
    explicit FloatMPTemplate<Lower>(Real const&, PrecisionMP);
    FloatMP const& get_flt() const;
    operator FloatMPTemplate<Apprx> () const;
    friend LowerFloatMP operator+(LowerFloatMP);
    friend UpperFloatMP operator-(LowerFloatMP);
    friend LowerFloatMP operator+(LowerFloatMP,LowerFloatMP);
    friend LowerFloatMP operator-(LowerFloatMP,UpperFloatMP);
    friend UpperFloatMP operator-(UpperFloatMP,LowerFloatMP);
    friend LowerFloatMP operator*(LowerFloatMP,LowerFloatMP);
    friend LowerFloatMP operator/(LowerFloatMP,UpperFloatMP);
    friend UpperFloatMP operator/(UpperFloatMP,LowerFloatMP);
    friend LowerFloatMP pos(LowerFloatMP);
    friend UpperFloatMP neg(LowerFloatMP);
    friend LowerFloatMP add(LowerFloatMP,LowerFloatMP);
    friend LowerFloatMP sub(LowerFloatMP,UpperFloatMP);
    friend UpperFloatMP sub(UpperFloatMP,LowerFloatMP);
    friend LowerFloatMP mul(LowerFloatMP,LowerFloatMP);
    friend LowerFloatMP div(LowerFloatMP,UpperFloatMP);
    friend UpperFloatMP div(UpperFloatMP,LowerFloatMP);
    friend OutputStream& operator<<(OutputStream& os, LowerFloatMP const&);
};

template<> class FloatMPTemplate<Apprx>
    : public NumberObject<ApprxFloatMP>
    , public DeclareAnalyticOperations<ApprxFloatMP>
    , public DeclareOrderedOperations<ApprxFloatMP>
    , public ProvideFieldOperators<ApprxFloatMP>
{
    FloatMP _a;
  public:
    typedef Apprx Paradigm;
    FloatMPTemplate<Apprx>();
    FloatMPTemplate<Apprx>(double);
    FloatMPTemplate<Apprx>(FloatMP);
    explicit FloatMPTemplate<Apprx>(Integer const&, PrecisionMP);
    explicit FloatMPTemplate<Apprx>(Rational const&, PrecisionMP);
    explicit FloatMPTemplate<Apprx>(Real const&, PrecisionMP);
    FloatMP const& get_flt() const;
    Void set_precision(PrecisionMP pr);
    PrecisionMP precision() const;
    friend ApprxFloatMP nul(ApprxFloatMP);
    friend ApprxFloatMP sqr(ApprxFloatMP);
    friend ApprxFloatMP pos(ApprxFloatMP);
    friend ApprxFloatMP neg(ApprxFloatMP);
    friend ApprxFloatMP rec(ApprxFloatMP);
    friend ApprxFloatMP add(ApprxFloatMP,ApprxFloatMP);
    friend ApprxFloatMP sub(ApprxFloatMP,ApprxFloatMP);
    friend ApprxFloatMP mul(ApprxFloatMP,ApprxFloatMP);
    friend ApprxFloatMP div(ApprxFloatMP,ApprxFloatMP);
    friend ApprxFloatMP abs(ApprxFloatMP);
    friend ApprxFloatMP max(ApprxFloatMP,ApprxFloatMP);
    friend ApprxFloatMP min(ApprxFloatMP,ApprxFloatMP);
    friend OutputStream& operator<<(OutputStream& os, ApprxFloatMP const&);
};


// Concrete Real type
template<class P> inline auto
    operator+(FloatMPTemplate<P> x, Real const& y) -> decltype(x+declval<ExactFloatMP>()) { return x+BoundFloatMP(y,x.precision()); }
template<class P> inline auto
    operator-(FloatMPTemplate<P> x, Real const& y) -> decltype(x-declval<ExactFloatMP>()) { return x-BoundFloatMP(y,x.precision()); }
template<class P> inline auto
    operator*(FloatMPTemplate<P> x, Real const& y) -> decltype(x*declval<ExactFloatMP>()) { return x*BoundFloatMP(y,x.precision()); }
template<class P> inline auto
    operator/(FloatMPTemplate<P> x, Real const& y) -> decltype(x/declval<ExactFloatMP>()) { return x/BoundFloatMP(y,x.precision()); }
template<class P> inline auto
    operator+(Real const& y, FloatMPTemplate<P> x) -> decltype(declval<ExactFloatMP>()+x) { return BoundFloatMP(y,x.precision())+x; }
template<class P> inline auto
    operator-(Real const& y, FloatMPTemplate<P> x) -> decltype(declval<ExactFloatMP>()-x) { return BoundFloatMP(y,x.precision())-x; }
template<class P> inline auto
    operator*(Real const& y, FloatMPTemplate<P> x) -> decltype(declval<ExactFloatMP>()*x) { return BoundFloatMP(y,x.precision())*x; }
template<class P> inline auto
    operator/(Real const& y, FloatMPTemplate<P> x) -> decltype(declval<ExactFloatMP>()/x) { return BoundFloatMP(y,x.precision())/x; }

} // namespace Ariadne

#endif
