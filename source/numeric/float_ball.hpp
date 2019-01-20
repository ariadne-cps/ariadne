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
    : public DispatchFloatOperations<Ball<F,FE>>
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
    friend PositiveUpperBound<F> mag(Ball<F,FE> const&);
    friend PositiveLowerBound<F> mig(Ball<F,FE> const&);
    friend Bool same(Ball<F,FE> const&, Ball<F,FE> const&);
    friend Bool models(Ball<F,FE> const&, Value<F> const&);
    friend Bool consistent(Ball<F,FE> const&, Ball<F,FE> const&);
    friend Bool inconsistent(Ball<F,FE> const&, Ball<F,FE> const&);
    friend Bool refines(Ball<F,FE> const&, Ball<F,FE> const&);
    friend Ball<F,FE> refinement(Ball<F,FE> const&, Ball<F,FE> const&);
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
    using typename Ball<F,FE>::PRE;
  public:
    Positive<Ball<F,FE>>() : Bounds<F>() { }
    template<class M, EnableIf<IsBuiltinUnsignedIntegral<M>> =dummy>
        Positive<Ball<F,FE>>(M m, PRE pre) : Ball<F,FE>(m,pre) { }
    explicit Positive<Ball<F,FE>>(Ball<F,FE> const& x) : Ball<F,FE>(x) { }
};

template<class F, class FE> inline PositiveBall<F,FE> cast_positive(Ball<F,FE> const& x) {
    return PositiveBall<F,FE>(x); }


}

#endif
