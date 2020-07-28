/***************************************************************************
 *            numeric/float_ball.hpp
 *
 *  Copyright  2008-20  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

/*! \file numeric/float_ball.hpp
 *  \brief Floating-point approximations to real numbers.
 */

#ifndef ARIADNE_FLOAT_BALL_HPP
#define ARIADNE_FLOAT_BALL_HPP

#include "../utility/macros.hpp"

#include "number.decl.hpp"
#include "float.decl.hpp"

#include "float_operations.hpp"
#include "float_traits.hpp"

namespace Ariadne {

template<class PRE, class PR, EnableIf<IsDefaultConstructible<PRE>> = dummy> inline
    PRE _error_precision(PR const&) { return PRE(); }
template<class PRE, class PR, DisableIf<IsDefaultConstructible<PRE>> = dummy, EnableIf<IsConstructible<PRE,PR>> = dummy> inline
    PRE _error_precision(PR const& pr) { return PRE(pr); }

//! \ingroup NumericModule
//! \brief Floating point approximations to a real number with guaranteed error bounds.
template<class F, class FE> class Ball
    : public DefineConcreteGenericOperators<Ball<F,FE>>
    , public DefineFieldOperators<Ball<F,FE>>
    , public DefineComparisonOperators<Ball<F,FE>,LessTrait<Ball<F,FE>>,EqualsTrait<Ball<F,FE>>>
    , public ProvideConvertedFieldOperations<Bounds<F>,Ball<F,FE>>
    , public ProvideConvertedFieldOperations<Ball<F,FE>,Value<F>>
{
    typedef ValidatedTag P; typedef typename F::PrecisionType PR; typedef typename FE::PrecisionType PRE;
    static_assert(IsConstructible<PR,PRE>::value or IsDefaultConstructible<PRE>::value,"");
  public:
    typedef P Paradigm;
    typedef Ball<F,FE> NumericType;
    typedef Number<P> GenericType;
    typedef F RawType;
    typedef FE RawErrorType;
    typedef PR PrecisionType;
    typedef PRE ErrorPrecisionType;
    typedef PR PropertiesType;
  public:
    explicit Ball<F,FE>(PrecisionType pr) : _v(0.0_x,pr), _e(0.0_x,_error_precision<PRE>(pr)) { }
    explicit Ball<F,FE>(PrecisionType pr, ErrorPrecisionType pre) : _v(0.0_x,pr), _e(0.0_x,pre) { }
    explicit Ball<F,FE>(F const& v) : _v(v), _e(0.0_x,_error_precision<PRE>(v.precision())) { }
    explicit Ball<F,FE>(F const& v, PRE pre) : _v(v), _e(0.0_x,pre) { }
    explicit Ball<F,FE>(F const& v, FE const& e) : _v(v), _e(e) { }
    Ball<F,FE>(Value<F> const& value, Error<FE> const& error);
    Ball<F,FE>(Value<F> const& x, PRE pre);
    Ball<F,FE>(Bounds<F> const& x, PRE pre);
    Ball<F,FE>(LowerBound<F> const& lower, UpperBound<F> const& upper) = delete;

    Ball<F,FE>(const ExactDouble& d, PR pr);
        Ball<F,FE>(const TwoExp& t, PR pr);
        Ball<F,FE>(const Integer& z, PR pr);
        Ball<F,FE>(const Dyadic& w, PR pr);
        Ball<F,FE>(const Decimal& d, PR pr);
        Ball<F,FE>(const Rational& q, PR pr);
        Ball<F,FE>(const Real& r, PR pr);
        Ball<F,FE>(const Ball<F,FE>& x, PR pr);
    Ball<F,FE>(const ValidatedNumber& y, PR pr);

    // FIXME: Constructors for other types
        Ball<F,FE>(const Dyadic& w, PR pr, PRE pre);
        Ball<F,FE>(const Rational& q, PR pr, PRE pre);
        Ball<F,FE>(const Real& q, PR pr, PRE pre);
    Ball<F,FE>(const ValidatedNumber& y, PR pr, PRE pre);

    explicit Ball<F,FE>(Bounds<F> const& x);
    Ball<F,FE>(Value<F> const& x);

    Ball<F,FE>& operator=(const ValidatedNumber& y) { return *this=Ball<F,FE>(y,this->precision(),this->error_precision()); }

    operator ValidatedNumber () const;

    Ball<F,FE> create(const ValidatedNumber& y) const { return Ball<F,FE>(y,this->precision(),this->error_precision()); }

    LowerBound<F> const lower() const;
    UpperBound<F> const upper() const;
    Value<F> const value() const;
    Error<FE> const error() const;

    friend Value<F> value(Ball<F,FE> const& x) { return x.value(); }
    friend Error<FE> error(Ball<F,FE> const& x) { return x.error(); }

    RawType const lower_raw() const { return sub(down,_v,_e); }
    RawType const upper_raw() const { return add(up,_v,_e); }
    RawType const& value_raw() const { return _v; }
    RawErrorType const& error_raw() const { return _e; }
    double get_d() const { return _v.get_d(); }

    PrecisionType precision() const { return _v.precision(); }
    ErrorPrecisionType error_precision() const { return _e.precision(); }
    PropertiesType properties() const { return _v.precision(); }
    GenericType generic() const { return this->operator GenericType(); }
    Ball<F,FE> pm(Error<FE> const& e) const;
    friend Approximation<F> round(Approximation<F> const& x);
  public:
    friend Bool is_nan(Ball<F,FE> const& x) {
        return is_nan(x._v) || is_nan(x._e); }

    friend Ball<F,FE> max(Ball<F,FE> const& x1, Ball<F,FE> const& x2) {
        return Operations<Ball<F,FE>>::_max(x1,x2); }
    friend Ball<F,FE> min(Ball<F,FE> const& x1, Ball<F,FE> const& x2) {
        return Operations<Ball<F,FE>>::_min(x1,x2); }
    friend Ball<F,FE> abs(Ball<F,FE> const& x) {
        return Operations<Ball<F,FE>>::_abs(x); }
    friend PositiveLowerBound<F> mig(Ball<F,FE> const& x) {
        return Operations<Ball<F,FE>>::_mig(x); }
    friend PositiveUpperBound<F> mag(Ball<F,FE> const& x) {
        return Operations<Ball<F,FE>>::_mag(x); }

    friend Ball<F,FE> nul(Ball<F,FE> const& x) {
        return Ball<F,FE>(nul(x._v),nul(x._e)); }
    friend Ball<F,FE> pos(Ball<F,FE> const& x) {
        return Ball<F,FE>(pos(x._v),x._e); }
    friend Ball<F,FE> neg(Ball<F,FE> const& x) {
        return Ball<F,FE>(neg(x._v),x._e); }
    friend Ball<F,FE> hlf(Ball<F,FE> const& x) {
        return Ball<F,FE>(hlf(x._v),hlf(x._e)); }
    friend Ball<F,FE> sqr(Ball<F,FE> const& x) {
        return Operations<Ball<F,FE>>::_sqr(x); }
    friend Ball<F,FE> rec(Ball<F,FE> const& x) {
        return Operations<Ball<F,FE>>::_rec(x); }

    friend Ball<F,FE> operator+(Ball<F,FE> const& x1, Ball<F,FE> const& x2);

    friend Ball<F,FE> add(Ball<F,FE> const& x1, Ball<F,FE> const& x2) {
        return Operations<Ball<F,FE>>::_add(x1,x2); }
    friend Ball<F,FE> sub(Ball<F,FE> const& x1, Ball<F,FE> const& x2) {
        return Operations<Ball<F,FE>>::_sub(x1,x2); }
    friend Ball<F,FE> mul(Ball<F,FE> const& x1, Ball<F,FE> const& x2) {
        return Operations<Ball<F,FE>>::_mul(x1,x2); }
    friend Ball<F,FE> div(Ball<F,FE> const& x1, Ball<F,FE> const& x2) {
        return Operations<Ball<F,FE>>::_div(x1,x2); }
    friend Ball<F,FE> fma(Ball<F,FE> const& x1, Ball<F,FE> const& x2, Ball<F,FE> const& x3) {
        return Operations<Ball<F,FE>>::_fma(x1,x2,x3); }
    friend Ball<F,FE> pow(Ball<F,FE> const& x, Int n) {
        return Operations<Ball<F,FE>>::_pow(x,n); }
    friend Ball<F,FE> pow(Ball<F,FE> const& x, Nat m) {
        return Operations<Ball<F,FE>>::_pow(x,m); }

    friend Ball<F,FE> sqrt(Ball<F,FE> const& x) {
        return Ball<F,FE>(sqrt(Bounds<F>(x)),x.error_precision()); }
    friend Ball<F,FE> exp(Ball<F,FE> const& x) {
        return Ball<F,FE>(exp(Bounds<F>(x)),x.error_precision()); }
    friend Ball<F,FE> log(Ball<F,FE> const& x) {
        return Ball<F,FE>(log(Bounds<F>(x)),x.error_precision()); }
    friend Ball<F,FE> sin(Ball<F,FE> const& x) {
        return Ball<F,FE>(sin(Bounds<F>(x)),x.error_precision()); }
    friend Ball<F,FE> cos(Ball<F,FE> const& x) {
        return Ball<F,FE>(cos(Bounds<F>(x)),x.error_precision()); }
    friend Ball<F,FE> tan(Ball<F,FE> const& x) {
        return Ball<F,FE>(tan(Bounds<F>(x)),x.error_precision()); }
    friend Ball<F,FE> asin(Ball<F,FE> const& x) {
        return Ball<F,FE>(asin(Bounds<F>(x)),x.error_precision()); }
    friend Ball<F,FE> acos(Ball<F,FE> const& x) {
        return Ball<F,FE>(acos(Bounds<F>(x)),x.error_precision()); }
    friend Ball<F,FE> atan(Ball<F,FE> const& x) {
        return Ball<F,FE>(atan(Bounds<F>(x)),x.error_precision()); }

    //! \related Ball<F,FE> \brief Strict greater-than comparison operator. Tests equality of represented real-point value.
    friend LogicalType<ValidatedTag> eq(Ball<F,FE> const& x1, Ball<F,FE> const& x2) {
        return Operations<Ball<F,FE>>::_eq(x1,x2); }
    //! \related Ball<F,FE> \brief Strict greater-than comparison operator. Tests equality of represented real-point value.
    friend LogicalType<ValidatedTag> lt(Ball<F,FE> const& x1, Ball<F,FE> const& x2) {
        return Operations<Ball<F,FE>>::_lt(x1,x2); }


    friend Ball<F,FE> round(Ball<F,FE> const& x) {
        return Ball<F,FE>(round(x.lower_raw()),round(x.upper_raw())); }
    friend Ball<F,FE> widen(Ball<F,FE> const& x) {
        const F m=std::numeric_limits<float>::min(); return Ball<F,FE>(sub(down,x._l,m),add(up,x._u,m)); }
    friend Ball<F,FE> narrow(Ball<F,FE> const& x) {
        const F m=std::numeric_limits<float>::min(); return Ball<F,FE>(add(up,x._l,m),add(down,x._u,m)); }
    friend Ball<F,FE> trunc(Ball<F,FE> const& x) {
        return Operations<Ball<F,FE>>::_trunc(x); }
    friend Ball<F,FE> trunc(Ball<F,FE> const& x, Nat n) {
        return Operations<Ball<F,FE>>::_trunc(x,n); }

    friend auto is_zero(Ball<F,FE> const& x) -> LogicalType<ValidatedTag> {
        if(x.lower_raw()>0.0 || x.upper_raw()<0.0) { return false; }
        else if(x.lower_raw()==0.0 && x.upper_raw()==0.0) { return true; }
        else { return indeterminate; } }
    friend auto is_positive(Ball<F,FE> const& x) -> LogicalType<ValidatedTag> {
        if(x.lower_raw()>=0.0) { return true; } else if(x.upper_raw()<0.0) { return false; } else { return indeterminate; } }

    friend Bool same(Ball<F,FE> const& x1, Ball<F,FE> const& x2) {
        return x1._v==x2._v && x1._e==x2._e; }
    friend Bool models(Ball<F,FE> const& x1, Value<F> const& x2) {
        return x1._l<=x2._v && x1._u >= x2._v; }
    friend Bool consistent(Ball<F,FE> const& x1, Ball<F,FE> const& x2) {
        return x1._l<=x2._u && x1._u >= x2._l; }
    friend Bool inconsistent(Ball<F,FE> const& x1, Ball<F,FE> const& x2) {
        return x1._l>x2._u || x1._u < x2._l; }
    friend Bool refines(Ball<F,FE> const& x1, Ball<F,FE> const& x2) {
        return x1._l>=x2._l && x1._u <= x2._u; }
    friend Ball<F,FE> refinement(Ball<F,FE> const& x1, Ball<F,FE> const& x2) {
        return Ball<F,FE>(refinement(Bounds<F>(x1),Bounds<F>(x2)),x1.error_precision()); }

    friend Integer cast_integer(Ball<F,FE> const& x) {
        return Operations<Ball<F,FE>>::_cast_integer(x); }

    friend OutputStream& operator<<(OutputStream& os, Ball<F,FE> const& x) { return Operations<Ball<F,FE>>::_write(os,x); }
    friend InputStream& operator>>(InputStream& is, Ball<F,FE>& x) { return Operations<Ball<F,FE>>::_read(is,x); }
  private: public:
    F _v; FE _e;
};

template<class PR> Ball(ValidatedNumber, PR) -> Ball<RawFloatType<PR>>;
template<class PR, class PRE> Ball(ValidatedNumber, PR, PRE) -> Ball<RawFloatType<PR>,RawFloatType<PRE>>;
template<class F, class FE> Ball(F,FE) -> Ball<F,FE>;

template<class F, class FE> inline FloatFactory<PrecisionType<F>> factory(Ball<F,FE> const& flt) {
    return FloatFactory<PrecisionType<F>>(flt.precision());
}

template<class F, class FE> class Positive<Ball<F,FE>> : public Ball<F,FE> {
    using PR = typename Ball<F,FE>::PrecisionType;
    using PRE = typename Ball<F,FE>::ErrorPrecisionType;
  public:
    Positive<Ball<F,FE>>() : Bounds<F>() { }
    template<class M, EnableIf<IsBuiltinUnsignedIntegral<M>> =dummy>
        Positive<Ball<F,FE>>(M m, PRE pre) : Ball<F,FE>(m,pre) { }
    explicit Positive<Ball<F,FE>>(PR pr, PRE pre) : Ball<F,FE>(pr,pre) { }
    explicit Positive<Ball<F,FE>>(Ball<F,FE> const& x) : Ball<F,FE>(x) { }
};

template<class F, class FE> inline PositiveBall<F,FE> cast_positive(Ball<F,FE> const& x) {
    return PositiveBall<F,FE>(x); }



template<class F, class FE> struct Operations<Ball<F,FE>> {
    typedef typename FE::PrecisionType PRE;

    static FE _make_error(F const& e) {
        static_assert(IsSame<F,FE>::value or IsDefaultConstructible<PRE>::value);
        if constexpr (IsSame<F,FE>()) { return e; }
        else if constexpr (IsDefaultConstructible<PRE>()) { return FE(e,up,PRE()); }
    }

    static Ball<F,FE> _nul(Ball<F,FE> const& x) {
        return Ball<F,FE>(nul(x._v),nul(x._e));
    }

    static Ball<F,FE> _pos(Ball<F,FE> const& x) {
        return Ball<F,FE>(pos(x._v),x._e);
    }

    static Ball<F,FE> _neg(Ball<F,FE> const& x) {
        return Ball<F,FE>(neg(x._v),x._e);
    }

    static Ball<F,FE> _hlf(Ball<F,FE> const& x) {
        return Ball<F,FE>(hlf(x._v),hlf(x._e));
    }

    static Ball<F,FE> _sqr(Ball<F,FE> const& x) {
        Ball<F,FE> r=x*x;
        if(r._e>r._v) {
            r._e=hlf(add(up,r._e,_make_error(r._v)));
            r._v=F(Dyadic(r._e),upward,x.precision());
        }
        return r;
    }

    static Ball<F,FE> _rec(Ball<F,FE> const& x) {
        // Use this code to find value same as reciprocal value
        auto rv=rec(approx,x._v);
        auto ru=rec(up,sub(down,x._v,x._e));
        auto rl=rec(down,add(up,x._v,x._e));
        auto re=max(sub(up,ru,rv),sub(up,rv,rl));
        return Ball<F,FE>(rv,_make_error(re));
    #ifdef ARIADNE_UNDEFINED
        // Use this code to get same result as interval computation
        auto ru=rec(up,sub(down,x._v,x._e));
        auto rl=rec(down,add(up,x._v,x._e));
        auto re=hlf(sub(up,ru,rl));
        auto rv=hlf(add(near,rl,ru));
        return Ball<F,FE>(rv,re);
    #endif
    }

    static Ball<F,FE> _add(Ball<F,FE> const& x, Ball<F,FE> const& y) {
        auto rv=add(near,x._v,y._v);
        auto ru=add(up,x._v,y._v);
        auto rl=add(down,x._v,y._v);
        auto se=add(up,x._e,y._e);
        auto ae=hlf(FE(sub(up,ru,rl),up,se.precision()));
        auto re=add(up,ae,se);
        return Ball<F,FE>(rv,re);
    }

    static Ball<F,FE> _sub(Ball<F,FE> const& x, Ball<F,FE> const& y) {
        auto rv=sub(near,x._v,y._v);
        auto ru=sub(up,x._v,y._v);
        auto rl=sub(down,x._v,y._v);
        auto ae=_make_error(hlf(sub(up,ru,rl)));
        auto re=add(up,ae,add(up,x._e,y._e));
        return Ball<F,FE>(rv,re);
    }

    static Ball<F,FE> _mul(Ball<F,FE> const& x, Ball<F,FE> const& y) {
        auto rv=mul(near,x._v,y._v);
        auto ru=mul(up,x._v,y._v);
        auto rl=mul(down,x._v,y._v);
        auto re0=_make_error(hlf(sub(up,ru,rl)));
        auto re1=add(up,re0,mul(up,x._e,y._e));
        auto re2=add(up,mul(up,_make_error(abs(x._v)),y._e),mul(up,x._e,_make_error(abs(y._v))));
        auto re=add(up,re1,re2);
        return Ball<F,FE>(rv,re);
    }

    static Ball<F,FE> _div(Ball<F,FE> const& x, Ball<F,FE> const& y) {
        return x*rec(y);
    }

    static Ball<F,FE> _add(Value<F> const& x, Value<F> const& y, PRE pre) {
        auto rv=add(near,x._v,y._v);
        auto ru=add(up,x._v,y._v);
        auto rl=add(down,x._v,y._v);
        auto re=FE(hlf(sub(up,ru,rl)),up,pre);
        return Ball<F,FE>(rv,re);
    }

    static Ball<F,FE> _sub(Value<F> const& x, Value<F> const& y, PRE pre) {
        auto rv=sub(near,x._v,y._v);
        auto ru=sub(up,x._v,y._v);
        auto rl=sub(down,x._v,y._v);
        auto re=FE(hlf(sub(up,ru,rl)),up,pre);
        return Ball<F,FE>(rv,re);
    }

    static Ball<F,FE> _mul(Value<F> const& x, Value<F> const& y, PRE pre) {
        auto rv=mul(near,x._v,y._v);
        auto ru=mul(up,x._v,y._v);
        auto rl=mul(down,x._v,y._v);
        auto re=FE(hlf(sub(up,ru,rl)),up,pre);
        return Ball<F,FE>(rv,re);
    }

    static Ball<F,FE> _div(Value<F> const& x, Value<F> const& y, PRE pre) {
        auto rv=div(near,x._v,y._v);
        auto ru=div(up,x._v,y._v);
        auto rl=div(down,x._v,y._v);
        auto re=FE(hlf(sub(up,ru,rl)),up,pre);
        return Ball<F,FE>(rv,re);
    }

    static Ball<F,FE> _pow(Ball<F,FE> const& x, Nat m) {
        return Ball<F,FE>(pow(Bounds<F>(x),m));
    }

    static Ball<F,FE> _pow(Ball<F,FE> const& x, Int n) {
        return Ball<F,FE>(pow(Bounds<F>(x),n));
    }

    static Ball<F,FE> _sqrt(Ball<F,FE> const& x) {
        return Ball<F,FE>(sqrt(Bounds<F>(x)));
    }

    static Ball<F,FE> _exp(Ball<F,FE> const& x) {
        return Ball<F,FE>(exp(Bounds<F>(x)));
    }

    static Ball<F,FE> _log(Ball<F,FE> const& x) {
        return Ball<F,FE>(log(Bounds<F>(x)));
    }

    static Ball<F,FE> _sin(Ball<F,FE> const& x) {
        return Ball<F,FE>(sin(Bounds<F>(x)));
    }

    static Ball<F,FE> _cos(Ball<F,FE> const& x) {
        return Ball<F,FE>(cos(Bounds<F>(x)));
    }

    static Ball<F,FE> _tan(Ball<F,FE> const& x) {
        return Ball<F,FE>(tan(Bounds<F>(x)));
    }

    static Ball<F,FE> _asin(Ball<F,FE> const& x) {
        return Ball<F,FE>(asin(Bounds<F>(x)));
    }

    static Ball<F,FE> _acos(Ball<F,FE> const& x) {
        return Ball<F,FE>(acos(Bounds<F>(x)));
    }

    static Ball<F,FE> _atan(Ball<F,FE> const& x) {
        return Ball<F,FE>(atan(Bounds<F>(x)));
    }


    static Ball<F,FE> _abs(Ball<F,FE> const& x) {
        if(x._e<abs(x._v)) { return x; }
        else { auto rv=hlf(add(up,abs(x._v),x._e)); return Ball<F,FE>(rv,_make_error(rv)); }
    }

    static Ball<F,FE> _max(Ball<F,FE> const& x1, Ball<F,FE> const& x2) {
        return hlf((x1+x2)+abs(x1-x2));
    }

    static Ball<F,FE> _min(Ball<F,FE> const& x1, Ball<F,FE> const& x2) {
        return hlf((x1+x2)-abs(x1-x2));
    }

    static Error<F> _mag(Ball<F,FE> const& x) {
        return PositiveUpperBound<F>(add(up,abs(x._v),x._e));
    }

    static PositiveLowerBound<F> _mig(Ball<F,FE> const& x) {
        return PositiveLowerBound<F>(max(F(0,x._v.precision()),sub(down,abs(x._v),x._e)));
    }

    //! \related Bounds<F> \brief Strict greater-than comparison operator. Tests equality of represented real-point value.
    static ValidatedKleenean _eq(Ball<F,FE> const& x1, Ball<F,FE> const& x2) {
        return Bounds<F>(x1) == Bounds<F>(x2);
    }

    //! \related Bounds<F> \brief Strict greater-than comparison operator. Tests equality of represented real-point value.
    static ValidatedKleenean _lt(Ball<F,FE> const& x1, Ball<F,FE> const& x2) {
        return Bounds<F>(x1) <  Bounds<F>(x2);
    }

    static Bool _same(Ball<F,FE> const& x1, Ball<F,FE> const& x2) {
        return x1._v==x2._v && x1._e==x2._e;
    }

    static Bool _models(Ball<F,FE> const& x1, Value<F> const& x2) {
        return (x1._v>=x2._v ? sub(up,x1._v,x2._v) : sub(up,x2._v,x1._v)) <= x1._e;
    }

    static Bool _consistent(Ball<F,FE> const& x1, Ball<F,FE> const& x2) {
        return consistent(Bounds<F>(x1),Bounds<F>(x2));
    }

    static Bool _inconsistent(Ball<F,FE> const& x1, Ball<F,FE> const& x2) {
        return inconsistent(Bounds<F>(x1),Bounds<F>(x2));
    }

    static Bool _refines(Ball<F,FE> const& x1, Ball<F,FE> const& x2) {
        return (x1._v>=x2._v ? sub(up,x1._v,x2._v) : sub(up,x2._v,x1._v)) <= sub(down,x2._e, x1._e);
    }

    static Ball<F,FE> _refinement(Ball<F,FE> const& x1, Ball<F,FE> const& x2) {
        return Ball<F,FE>(refinement(Bounds<F>(x1),Bounds<F>(x2)));
    }

    static Integer _cast_integer(Ball<F,FE> const& x);

    static OutputStream& _write(OutputStream& os, Ball<F,FE> const& x);

    static InputStream& _read(InputStream& is, Ball<F,FE>& x);
};



}

#endif
