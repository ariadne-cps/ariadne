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

/************ FltMP ********************************************************/

//enum class RoundingModeMP : mpfr_rnd_t { NEAREST=MPFR_RNDN, UPWARD=MPFR_RNDU, DOWNWARD=MPFR_RNDD };
enum class RoundingModeMP : uint { NEAREST=MPFR_RNDN, UPWARD=MPFR_RNDU, DOWNWARD=MPFR_RNDD };

template<class F> struct FltMPExpression { };

template<class N> class NumberObject { };

struct PrecisionMP {
    mpfr_prec_t prec;
    PrecisionMP(mpfr_prec_t pr) : prec(pr) { }
    operator mpfr_prec_t () const { return prec; }
};

struct NoInit { };

//! \ingroup FltMPSubModule
//! \brief Multiple-precision floating-point numbers.
//! Currently defined as _a wrapper around \c mpfr_t from the MPFE library.
//! Default arithmetic operations are approximate, and comparisons are exact, so this class is \em unsafe.
class FltMP {
  private: public:
    mpfr_t _mpfr;
    typedef decltype(_mpfr[0]) MpfrReference;
  public:
    typedef PrecisionMP PrecisionType;
    typedef mpfr_rnd_t RoundingModeType;
    static Void set_rounding_mode(RoundingModeMP rnd);
    static Void set_default_precision(PrecisionMP prec);
    static PrecisionMP get_default_precision();
    static mpfr_rnd_t current_rounding_mode;
  public:
    ~FltMP();
    explicit FltMP();
    explicit FltMP(NoInit);
    explicit FltMP(PrecisionMP);
    explicit FltMP(PrecisionMP, NoInit);
    explicit FltMP(Int32);
    explicit FltMP(Int32, PrecisionMP);
    explicit FltMP(Integer const&, PrecisionMP, RoundingModeType);
    explicit FltMP(Rational const&, PrecisionMP, RoundingModeType);
    FltMP(double);
    FltMP(double, PrecisionMP);
    FltMP(const mpfr_t);
    FltMP(const FltMP&);
    FltMP(FltMP&&);
    template<class N, EnableIf<IsIntegral<N>> =dummy> FltMP& operator=(N n);
    FltMP& operator=(const FltMP&);
    FltMP& operator=(FltMP&&);
    explicit operator Rational() const;

    PrecisionMP precision() const;
    Void set_precision(PrecisionMP);

    friend FltMP next_up(FltMP x);
    friend FltMP next_down(FltMP x);

    friend FltMP operator+(FltMP const& x);
    friend FltMP operator-(FltMP const& x);
    friend FltMP operator+(FltMP const& x1, FltMP const& x2);
    friend FltMP operator-(FltMP const& x1, FltMP const& x2);
    friend FltMP operator*(FltMP const& x1, FltMP const& x2);
    friend FltMP operator/(FltMP const& x1, FltMP const& x2);

    //friend Integer floor(FltMP const& x);
    //friend Integer ceil(FltMP const& x);
    friend FltMP floor(FltMP const& x);
    friend FltMP ceil(FltMP const& x);

    friend FltMP nul(FltMP const& x);
    friend FltMP pos(FltMP const& x);
    friend FltMP neg(FltMP const& x);
    friend FltMP abs(FltMP const& x);
    friend FltMP half(FltMP&& x);

    friend FltMP pos(FltMP const& x, RoundingModeType);
    friend FltMP neg(FltMP const& x, RoundingModeType);
    friend FltMP rec(FltMP const& x, RoundingModeType);
    friend FltMP add(FltMP const& x1, FltMP const& x2, RoundingModeType);
    friend FltMP sub(FltMP const& x1, FltMP const& x2, RoundingModeType);
    friend FltMP mul(FltMP const& x1, FltMP const& x2, RoundingModeType);
    friend FltMP div(FltMP const& x1, FltMP const& x2, RoundingModeType);
    friend FltMP exp(FltMP const& x, RoundingModeType);

    friend FltMP abs(FltMP const& x, RoundingModeType);

    friend Comparison cmp(FltMP const& x1, FltMP const& x2);
    friend Bool operator==(FltMP const& x1, FltMP const& x2);
    friend Bool operator!=(FltMP const& x1, FltMP const& x2);
    friend Bool operator<=(FltMP const& x1, FltMP const& x2);
    friend Bool operator>=(FltMP const& x1, FltMP const& x2);
    friend Bool operator< (FltMP const& x1, FltMP const& x2);
    friend Bool operator> (FltMP const& x1, FltMP const& x2);

    friend Comparison cmp(FltMP const& x1, double x2);
    friend Bool operator==(FltMP const& x1, double x2);
    friend Bool operator!=(FltMP const& x1, double x2);
    friend Bool operator<=(FltMP const& x1, double x2);
    friend Bool operator>=(FltMP const& x1, double x2);
    friend Bool operator< (FltMP const& x1, double x2);
    friend Bool operator> (FltMP const& x1, double x2);

    template<class FE> FltMP(const FltMPExpression<FE>&);
    friend OutputStream& operator<<(OutputStream& os, FltMP const&);
  public:
    FltMP& operator=(double d);
    MpfrReference get_mpfr();
    MpfrReference get_mpfr() const;
    double get_d() const;
};

template<class N, EnableIf<IsIntegral<N>>> inline FltMP& FltMP::operator=(N n) { double x=n; return *this=x; }

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
    FltMP _e;
  public:
    typedef Error Paradigm;
    FloatMPTemplate<Error>(uint);
    explicit FloatMPTemplate<Error>(double);
    explicit FloatMPTemplate<Error>(double,PrecisionMP);
    explicit FloatMPTemplate<Error>(FltMP);
    FltMP const& get_flt() const;
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
    FltMP _v;
  public:
    typedef Exact Paradigm;
    template<class N, EnableIf<IsIntegral<N>> = dummy>
        FloatMPTemplate<Exact>(N n) : _v(n) { }
    template<class N, EnableIf<IsIntegral<N>> = dummy>
        FloatMPTemplate<Exact>(N n, PrecisionMP pr) : _v(n,pr) { }
    template<class N, EnableIf<IsIntegral<N>> = dummy>
        FloatMPTemplate<Exact>& operator=(N n) { _v=n; return *this; }
    explicit FloatMPTemplate<Exact>(double);
    explicit FloatMPTemplate<Exact>(FltMP);
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
    FltMP const& get_flt() const;
    friend OutputStream& operator<<(OutputStream& os, ExactFloatMP const&);
};

template<> class FloatMPTemplate<Metrc> : public NumberObject<MetrcFloatMP>
    , public DeclareAnalyticOperations<MetrcFloatMP>
    , public DeclareOrderedOperations<MetrcFloatMP>
    , public ProvideFieldOperators<MetrcFloatMP>
{
    friend class FloatMPTemplate<Apprx>;
    FltMP _v; FltMP _e;
  public:
    typedef Metrc Paradigm;
    FloatMPTemplate<Metrc>();
    template<class N, EnableIf<IsIntegral<N>> =dummy> FloatMPTemplate<Metrc>(N n) : FloatMPTemplate<Metrc>(Int32(n)) { }
    explicit FloatMPTemplate<Metrc>(Integer const&, PrecisionMP);
    explicit FloatMPTemplate<Metrc>(Rational const&, PrecisionMP);
    explicit FloatMPTemplate<Metrc>(Real const&, PrecisionMP);
//    explicit FloatMPTemplate<Metrc>(double);
//    explicit FloatMPTemplate<Metrc>(double,double);
    explicit FloatMPTemplate<Metrc>(FltMP v);
    explicit FloatMPTemplate<Metrc>(FltMP v, FltMP e);
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
    friend bool operator==(MetrcFloatMP,int);
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
    FltMP _l; FltMP _u;
  public:
    typedef Bound Paradigm;
    FloatMPTemplate<Bound>();
    template<class N, EnableIf<IsIntegral<N>> =dummy> FloatMPTemplate<Bound>(N n) : FloatMPTemplate<Bound>(Int32(n)) { }
//    explicit FloatMPTemplate<Metrc>(double);
//    explicit FloatMPTemplate<Metrc>(double,double);
    explicit FloatMPTemplate<Bound>(Integer const&, PrecisionMP);
    explicit FloatMPTemplate<Bound>(Rational const&, PrecisionMP);
    explicit FloatMPTemplate<Bound>(Real const&, PrecisionMP);
    explicit FloatMPTemplate<Bound>(FltMP);
    explicit FloatMPTemplate<Bound>(FltMP,FltMP);
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
    friend bool operator==(BoundFloatMP,int);
    friend OutputStream& operator<<(OutputStream& os, BoundFloatMP const&);
 private:
    friend ApprxFloatMP operator+(ApprxFloatMP,ApprxFloatMP);
    FloatMPTemplate<Bound>(Int32);
};

template<> class FloatMPTemplate<Upper>
    : public NumberObject<UpperFloatMP>
{
    FltMP _u;
  public:
    typedef Upper Paradigm;
    FloatMPTemplate<Upper>();
    FloatMPTemplate<Upper>(FltMP);
    FltMP const& get_flt() const;
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
    FltMP _l;
  public:
    typedef Lower Paradigm;
    FloatMPTemplate<Lower>();
    FloatMPTemplate<Lower>(FltMP);
    explicit FloatMPTemplate<Lower>(Real const&, PrecisionMP);
    FltMP const& get_flt() const;
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
    FltMP _a;
  public:
    typedef Apprx Paradigm;
    FloatMPTemplate<Apprx>();
    FloatMPTemplate<Apprx>(double);
    FloatMPTemplate<Apprx>(FltMP);
    explicit FloatMPTemplate<Apprx>(Integer const&, PrecisionMP);
    explicit FloatMPTemplate<Apprx>(Rational const&, PrecisionMP);
    explicit FloatMPTemplate<Apprx>(Real const&, PrecisionMP);
    FltMP const& get_flt() const;
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
