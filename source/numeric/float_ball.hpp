/***************************************************************************
 *            float_ball.hpp
 *
 *  Copyright 2008-17  Pieter Collins
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

/*! \file float_ball.hpp
 *  \brief Floating-point approximations to real numbers.
 */

#ifndef ARIADNE_FLOAT_BALL_HPP
#define ARIADNE_FLOAT_BALL_HPP

#include "../utility/macros.hpp"

#include "number.decl.hpp"
#include "float.decl.hpp"

#include "float_operations.hpp"

namespace Ariadne {

template<class FE, class FLT, DisableIf<IsSame<FE,FLT>> =dummy> inline FE _make_error(FLT const& x) {
    typename FE::PrecisionType pre; return FE(Dyadic(x),upward,pre); }
template<class FE, class FLT, EnableIf<IsSame<FE,FLT>> =dummy> inline FE _make_error(FLT const& x) {
    return x; }
template<class FE, class FLT, class PRE> inline FE _make_error(FLT const& x, PRE pre) {
    return FE(Dyadic(x),upward,pre); }

template<class F, class FE> struct NumericTraits<Ball<F,FE>> {
    typedef ValidatedNumber GenericType;
    typedef Ball<F,FE> OppositeType;
    typedef PositiveBall<F,FE> PositiveType;
    typedef ValidatedKleenean LessType;
    typedef ValidatedKleenean EqualsType;
};

template<class PRE, class PR, EnableIf<IsDefaultConstructible<PRE>> = dummy> PRE _error_precision(PR const&) { return PRE(); }
template<class PRE, class PR, DisableIf<IsDefaultConstructible<PRE>> = dummy, EnableIf<IsConstructible<PRE,PR>> = dummy> PRE _error_precision(PR const& pr) { return PRE(pr); }

//! \ingroup NumericModule
//! \brief Floating point approximations to a real number with guaranteed error bounds.
template<class F, class FE> class Ball
    : public DefineConcreteGenericOperators<Ball<F,FE>>
    , public DefineFieldOperators<Ball<F,FE>>
    , public DefineComparisonOperators<Ball<F,FE>,LessTrait<Ball<F,FE>>,EqualsTrait<Ball<F,FE>>>
//    : public DispatchFloatOperations<Ball<F,FE>>
//    , public ProvideConvertedFieldOperations<Bounds<F>,Ball<F,FE>>
//    , public ProvideConvertedFieldOperations<Ball<F,FE>,Value<F>>
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
    Ball<F,FE>() : _v(0.0), _e(0.0) { }
    explicit Ball<F,FE>(PrecisionType pr) : _v(0.0,pr), _e(0.0,_error_precision<PRE>(pr)) { }
    explicit Ball<F,FE>(PrecisionType pr, ErrorPrecisionType pre) : _v(0.0,pr), _e(0.0,pre) { }
    explicit Ball<F,FE>(F const& v) : _v(v), _e(0.0,_error_precision<PRE>(v.precision())) { }
    explicit Ball<F,FE>(F const& v, PRE pre) : _v(v), _e(0.0,pre) { }
    explicit Ball<F,FE>(F const& v, FE const& e) : _v(v), _e(e) { }
    Ball<F,FE>(Value<F> const& value, Error<FE> const& error);
    Ball<F,FE>(Bounds<F> const& x, PRE pre);
    Ball<F,FE>(LowerBound<F> const& lower, UpperBound<F> const& upper) = delete;

    Ball<F,FE>(ExactDouble d, PR pr);
        Ball<F,FE>(TwoExp t, PR pr);
        Ball<F,FE>(const Integer& z, PR pr);
        Ball<F,FE>(const Dyadic& w, PR pr);
        Ball<F,FE>(const Decimal& d, PR pr);
        Ball<F,FE>(const Rational& q, PR pr);
        Ball<F,FE>(const Real& r, PR pr);
        Ball<F,FE>(const Ball<F,FE>& x, PR pr);
    Ball<F,FE>(const ValidatedNumber& y, PR pr);

    // FIXME: Constructors for other types
        Ball<F,FE>(const Rational& q, PR pr, PRE pre);
        Ball<F,FE>(const Real& q, PR pr, PRE pre);
    Ball<F,FE>(const ValidatedNumber& y, PR pr, PRE pre);

    explicit Ball<F,FE>(Bounds<F> const& x);
    Ball<F,FE>(Value<F> const& x);

    Ball<F,FE>& operator=(const ValidatedNumber& y);

    operator ValidatedNumber () const { return ValidatedNumber(new NumberWrapper<Ball<F,FE>>(*this)); }

    Ball<F,FE> create(const ValidatedNumber& y) const;

    LowerBound<F> const lower() const;
    UpperBound<F> const upper() const;
    Value<F> const value() const;
    Error<FE> const error() const;

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
    friend Ball<F,FE> nul(Ball<F,FE> const& x) {
        return Ball<F,FE>(nul(x._v),nul(x._e)); }
    friend Ball<F,FE> pos(Ball<F,FE> const& x) {
        return Ball<F,FE>(pos(x._v),x._e); }
    friend Ball<F,FE> neg(Ball<F,FE> const& x) {
        return Ball<F,FE>(neg(x._v),x._e); }
    friend Ball<F,FE> hlf(Ball<F,FE> const& x) {
        return Ball<F,FE>(hlf(x._v),hlf(x._e)); }
    friend Ball<F,FE> sqr(Ball<F,FE> const& x) {
        Ball<F,FE> r=x*x; if(r._e>r._v) { r._e=hlf(add(up,r._e,_make_error<FE>(r._v))); r._v=F(Dyadic(r._e),upward,x.precision()); } return r; }
    friend Ball<F,FE> rec(Ball<F,FE> const& x) {
        return Operations<Ball<F,FE>>::_rec(x); }

    friend Ball<F,FE> add(Ball<F,FE> const& x1, Ball<F,FE> const& x2) {
        return Operations<Ball<F,FE>>::_add(x1,x2); }
    friend Ball<F,FE> sub(Ball<F,FE> const& x1, Ball<F,FE> const& x2) {
        return Operations<Ball<F,FE>>::_sub(x1,x2); }
    friend Ball<F,FE> mul(Ball<F,FE> const& x1, Ball<F,FE> const& x2) {
        return Operations<Ball<F,FE>>::_mul(x1,x2); }
    friend Ball<F,FE> div(Ball<F,FE> const& x1, Ball<F,FE> const& x2) {
        return Operations<Ball<F,FE>>::_div(x1,x2); }

    friend Ball<F,FE> pow(Ball<F,FE> const& x, Nat m) {
        return Ball<F,FE>(pow(Bounds<F>(x),m),x.error_precision()); }
    friend Ball<F,FE> pow(Ball<F,FE> const& x, Int n) {
        return Ball<F,FE>(pow(Bounds<F>(x),n),x.error_precision()); }
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

    friend Ball<F,FE> max(Ball<F,FE> const& x1, Ball<F,FE> const& x2) {
        return hlf((x1+x2)+abs(x1-x2)); }
    friend Ball<F,FE> min(Ball<F,FE> const& x1, Ball<F,FE> const& x2) {
        return hlf((x1+x2)-abs(x1-x2)); }
    friend Ball<F,FE> abs(Ball<F,FE> const& x) {
        if(x._e<abs(x._v)) { return x; } else { auto rv=hlf(add(up,abs(x._v),x._e)); return Ball<F,FE>(rv,_make_error<FE>(rv)); } }
    friend PositiveUpperBound<F> mag(Ball<F,FE> const& x) {
        return PositiveUpperBound<F>(add(up,abs(x._v),x._e)); }
    friend PositiveLowerBound<F> mig(Ball<F,FE> const& x) {
        return PositiveLowerBound<F>(max(0,sub(down,abs(x._v),x._e))); }

    //! \related Bounds<F> \brief Strict greater-than comparison operator. Tests equality of represented real-point value.
    friend ValidatedKleenean eq(Ball<F,FE> const& x1, Ball<F,FE> const& x2) {
        return Bounds<F>(x1) == Bounds<F>(x2); }
    //! \related Bounds<F> \brief Strict greater-than comparison operator. Tests equality of represented real-point value.
    friend ValidatedKleenean lt(Ball<F,FE> const& x1, Ball<F,FE> const& x2) {
        return Bounds<F>(x1) <  Bounds<F>(x2); }

    friend Bool same(Ball<F,FE> const& x1, Ball<F,FE> const& x2) {
        return x1._v==x2._v && x1._e==x2._e; }
    friend Bool models(Ball<F,FE> const& x1, Value<F> const& x2) {
        return (x1._v>=x2._v ? sub(up,x1._v,x2._v) : sub(up,x2._v,x1._v)) <= x1._e; }
    friend Bool consistent(Ball<F,FE> const& x1, Ball<F,FE> const& x2) {
        return consistent(Bounds<F>(x1),Bounds<F>(x2)); }
    friend Bool inconsistent(Ball<F,FE> const& x1, Ball<F,FE> const& x2) {
        return inconsistent(Bounds<F>(x1),Bounds<F>(x2)); }
    friend Bool refines(Ball<F,FE> const& x1, Ball<F,FE> const& x2) {
        return (x1._v>=x2._v ? sub(up,x1._v,x2._v) : sub(up,x2._v,x1._v)) <= sub(down,x2._e, x1._e); }
    friend Ball<F,FE> refinement(Ball<F,FE> const& x1, Ball<F,FE> const& x2) {
        return Ball<F,FE>(refinement(Bounds<F>(x1),Bounds<F>(x2))); }

    friend OutputStream& operator<<(OutputStream& os, Ball<F,FE> const& x) {
        return Operations<Ball<F,FE>>::_write(os,x); }
    friend InputStream& operator>>(InputStream& is, Ball<F,FE>& x) {
        return Operations<Ball<F,FE>>::_read(is,x); }

  private: public:
    F _v; FE _e;
};

template<class F, class FE> inline FloatFactory<PrecisionType<F>> factory(Ball<F,FE> const& flt) {
    return FloatFactory<PrecisionType<F>>(flt.precision());
}

template<class F, class FE> class Positive<Ball<F,FE>> : public Ball<F,FE>
    , public DefineArithmeticOperators<Positive<Ball<F,FE>>>
    , public DefineConcreteGenericOperators<Positive<Ball<F,FE>>>
    , public ProvidePositiveFieldOperations<Positive<Ball<F,FE>>>
{
    using typename Ball<F,FE>::PRE;
  public:
    Positive<Ball<F,FE>>() : Bounds<F>() { }
    template<class M, EnableIf<IsBuiltinUnsignedIntegral<M>> =dummy>
        Positive<Ball<F,FE>>(M m, PRE pre) : Ball<F,FE>(m,pre) { }
    explicit Positive<Ball<F,FE>>(Ball<F,FE> const& x) : Ball<F,FE>(x) { }
};

template<class F, class FE> inline PositiveBall<F,FE> cast_positive(Ball<F,FE> const& x) {
    return PositiveBall<F,FE>(x); }



template<class F, class FE> struct Operations<Ball<F,FE>> {
    static Ball<F,FE> _rec(Ball<F,FE> const& x);
    static Ball<F,FE> _add(Ball<F,FE> const& x, Ball<F,FE> const& y);
    static Ball<F,FE> _sub(Ball<F,FE> const& x, Ball<F,FE> const& y);
    static Ball<F,FE> _mul(Ball<F,FE> const& x, Ball<F,FE> const& y);
    static Ball<F,FE> _div(Ball<F,FE> const& x, Ball<F,FE> const& y);
    static OutputStream& _write(OutputStream& os, Ball<F,FE> const& x);
    static InputStream& _read(InputStream& is, Ball<F,FE>& x);

};

} // namespace Ariadne

#include "float_ball.inl.hpp"

#endif
