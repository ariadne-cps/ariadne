/***************************************************************************
 *            numeric/float_error.hpp
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

/*! \file numeric/float_error.hpp
 *  \brief Floating-point error bounds for metric spaces.
 */

#ifndef ARIADNE_FLOAT_ERROR_HPP
#define ARIADNE_FLOAT_ERROR_HPP

#include "../utility/macros.hpp"
#include "../utility/metaprogramming.hpp"

#include "logical.decl.hpp"
#include "number.decl.hpp"
#include "float.decl.hpp"

#include "float_traits.hpp"
#include "float_operations.hpp"

#include "float_upper_bound.hpp"

namespace Ariadne {

//! \ingroup NumericModule
//! \brief Floating-point upper bounds for positive real numbers, suitable for use as an upper bound for an error in a metric space.
template<class F> class Error
    : public DefineDirectedGroupOperators<UpperBound<F>,LowerBound<F>>
    , public DefineDirectedGroupOperators<LowerBound<F>,UpperBound<F>>
    , public DefineDirectedComparisonOperators<UpperBound<F>,LowerBound<F>,LessTrait<UpperBound<F>>,EqualsTrait<UpperBound<F>>>
    , public DefineDirectedComparisonOperators<LowerBound<F>,UpperBound<F>,LessTrait<LowerBound<F>>,EqualsTrait<LowerBound<F>>>
    , public DefineConcreteGenericOperators<UpperBound<F>>

    , public DeclarePositiveDirectedNumericOperations<PositiveUpperBound<F>,PositiveLowerBound<F>>
    , public ProvideConcreteGenericDirectedSemiFieldOperations<PositiveUpperBound<F>,PositiveLowerBound<F>,Nat,Nat>
{
    using PR=typename F::PrecisionType;
  private: public:
    F _e;
  public:
    typedef PositiveValidatedUpperNumber GenericType;
    typedef PR PrecisionType;
    typedef PR PropertiesType;
  public:
    Error<F>(PositiveValue<F> const& x);
    Error<F>(PositiveBounds<F> const& x);
    Error<F>(PositiveUpperBound<F> const& x) : _e(x._u) { }
    operator PositiveUpperBound<F> const& () const { return reinterpret_cast<PositiveUpperBound<F>const&>(*this); }
    operator PositiveUpperBound<F>& () { return reinterpret_cast<PositiveUpperBound<F>&>(*this); }
  public:
    Error<F>() : _e() { }
    explicit Error<F>(PR const& pr) : _e(pr) { }
    // Check error not being negative to allow for NaN as a valid input.
    explicit Error<F>(F const& x) : _e(x) { ARIADNE_PRECONDITION_MSG(!(this->_e<0),"e="<<*this); }
    template<class M, EnableIf<IsBuiltinUnsignedIntegral<M>> =dummy> Error<F>(M m, PR pr) : _e(m,pr) { }
    explicit Error<F>(UpperBound<F> const& x) : Error<F>(x._u) { }
    explicit Error<F>(ValidatedUpperNumber const& y, PR pr) : Error<F>(UpperBound<F>(y,pr)) { }
    explicit Error<F>(const ExactDouble& d, PR pr) : Error<F>(UpperBound<F>(d,pr)) { }
    explicit Error<F>(const TwoExp& t, PR pr) : Error<F>(UpperBound<F>(t,pr)) { }
    Error<F>& operator=(Nat m) { reinterpret_cast<UpperBound<F>&>(*this)=m; return *this; }
    Error<F>& operator=(ValidatedErrorNumber y) { return *this=cast_positive(y.get(this->precision())); }
    operator ValidatedErrorNumber() const;
  public:
    ValidatedErrorNumber generic() const { return this->operator ValidatedErrorNumber(); }
    PrecisionType precision() const { return _e.precision(); }
    PropertiesType properties() const { return _e.precision(); }
    F const& raw() const { return _e; }
    F& raw() { return _e; }
  public:
    friend Error<F> operator+(Error<F> const& x1, Error<F> const& x2) { return Error<F>(add(up,x1._e,x2._e)); }
    friend Error<F> operator*(Error<F> const& x1, Error<F> const& x2) { return Error<F>(mul(up,x1._e,x2._e)); }

    friend Error<F> nul(Error<F> const& x) { return Error<F>(0u,x.precision()); }
    friend UpperBound<F> pos(Error<F> const& x) { return UpperBound<F>(pos(x._e)); }
    friend LowerBound<F> neg(Error<F> const& x) { return LowerBound<F>(neg(x._e)); }
    friend Error<F> add(Error<F> const& x1, Error<F> const& x2) { return Error<F>(add(up,x1._e,x2._e)); }
    friend Error<F> mul(Error<F> const& x1, Error<F> const& x2) { return Error<F>(mul(up,x1._e,x2._e)); }
    friend Error<F> sqr(Error<F> const& x) { return Error<F>(sqr(up,x._e)); }
    friend Error<F> pow(Error<F> const& x, Nat m) { return Error<F>(pow(up,x._e,static_cast<Int>(m))); }
    friend Error<F> exp(Error<F> const& x) { return Error<F>(exp(up,x._e)); }
    friend UpperBound<F> log(Error<F> const& x) { return UpperBound<F>(log(up,x._e)); }
    friend Error<F> max(Error<F> const& x1, Error<F> const& x2) { return Error<F>(max(x1._e,x2._e)); }
    friend Error<F> min(Error<F> const& x1, Error<F> const& x2) { return Error<F>(min(x1._e,x2._e)); }
    friend Error<F> abs(Error<F> const& x) { return x; }
    friend Error<F> mag(Error<F> const& x) { return x; }

        friend Error<F> operator+(Error<F> const& x1, PositiveValue<F> const& x2) { return Error<F>(add(up,x1._e,x2._v)); }
        friend Error<F> operator+(PositiveValue<F> const& x1, Error<F> const& x2) { return Error<F>(add(up,x1._v,x2._e)); }
        friend Error<F> operator*(Error<F> const& x1, PositiveValue<F> const& x2) { return Error<F>(mul(up,x1._e,x2._v)); }
        friend Error<F> operator*(PositiveValue<F> const& x1, Error<F> const& x2) { return Error<F>(mul(up,x1._v,x2._e)); }

        friend Error<F> operator+(Error<F> const& x1, PositiveBounds<F> const& x2) { return Error<F>(add(up,x1._e,x2._u)); }
        friend Error<F> operator+(PositiveBounds<F> const& x1, Error<F> const& x2) { return Error<F>(add(up,x1._u,x2._e)); }
        friend Error<F> operator*(Error<F> const& x1, PositiveBounds<F> const& x2) { return Error<F>(mul(up,x1._e,x2._u)); }
        friend Error<F> operator*(PositiveBounds<F> const& x1, Error<F> const& x2) { return Error<F>(mul(up,x1._u,x2._e)); }

        friend Error<F> operator+(Error<F> const& x1, PositiveUpperBound<F> const& x2) { return Error<F>(add(up,x1._e,x2._u)); }
        friend Error<F> operator+(PositiveUpperBound<F> const& x1, Error<F> const& x2) { return Error<F>(add(up,x1._u,x2._e)); }
        friend Error<F> operator*(Error<F> const& x1, PositiveUpperBound<F> const& x2) { return Error<F>(mul(up,x1._e,x2._u)); }
        friend Error<F> operator*(PositiveUpperBound<F> const& x1, Error<F> const& x2) { return Error<F>(mul(up,x1._u,x2._e)); }

        friend Error<F> operator/(Error<F> const& x1, PositiveLowerBound<F> const& x2) { return Error<F>(div(up,x1._e,x2._l)); }

        friend Error<F> operator+(Error<F> const& x1, Nat m) { return Error<F>(add(up,x1._e,F(m,up,x1.precision()))); }
        friend Error<F> operator*(Error<F> const& x1, Nat m) { return Error<F>(mul(up,x1._e,F(m,up,x1.precision()))); }
        friend Error<F> operator/(Error<F> const& x1, Nat m) { return Error<F>(div(up,x1._e,F(m,down,x1.precision()))); }

        friend Approximation<F> operator-(Error<F> const& x1, Error<F> const& x2) { return Approximation<F>(sub(up,x1._e,x2._e)); }
        friend UpperBound<F> operator-(UpperBound<F> const& x1, LowerBound<F> const& x2);
        friend LowerBound<F> operator-(LowerBound<F> const& x1, UpperBound<F> const& x2);

    /*
    friend PositiveUpperBound<F> operator+(PositiveUpperBound<F> const& x1, PositiveUpperBound<F> const& x2);
    friend PositiveUpperBound<F> operator*(PositiveUpperBound<F> const& x1, PositiveUpperBound<F> const& x2);
    friend PositiveUpperBound<F> add(PositiveUpperBound<F> const& x1, PositiveUpperBound<F> const& x2);
    friend PositiveUpperBound<F> mul(PositiveUpperBound<F> const& x1, PositiveUpperBound<F> const& x2);
    friend PositiveUpperBound<F> max(PositiveUpperBound<F> const& x1, PositiveUpperBound<F> const& x2);
    friend PositiveUpperBound<F> min(PositiveUpperBound<F> const& x1, PositiveUpperBound<F> const& x2);
*/
    friend Error<F>& operator+=(Error<F>& x1, Error<F> const& x2) { return x1=x1+x2; }
        friend Error<F>& operator+=(Error<F>& x1, PositiveUpperBound<F> const& x2) { return x1=x1+x2; }
        friend Error<F>& operator+=(Error<F>& x1, PositiveBounds<F> const& x2) { return x1=x1+x2; }
        friend Error<F>& operator+=(Error<F>& x1, PositiveValue<F> const& x2) { return x1=x1+x2; }
    friend Error<F>& operator*=(Error<F>& x1, Error<F> const& x2) { return x1=x1*x2; }
        friend Error<F>& operator*=(Error<F>& x1, PositiveUpperBound<F> const& x2) { return x1=x1*x2; }
        friend Error<F>& operator*=(Error<F>& x1, PositiveBounds<F> const& x2) { return x1=x1*x2; }
        friend Error<F>& operator*=(Error<F>& x1, PositiveValue<F> const& x2) { return x1=x1*x2; }

  public:
    friend UpperBound<F> operator+(Error<F> const& x) { return UpperBound<F>(+x._e); }
    friend LowerBound<F> operator-(Error<F> const& x) { return LowerBound<F>(-x._e); }
    friend UpperBound<F> operator+(Value<F> const& x1, Error<F> const& x2) { return UpperBound<F>(add(up,x1._v,x2._e)); }
    friend LowerBound<F> operator-(Value<F> const& x1, Error<F> const& x2) { return LowerBound<F>(sub(down,x1._v,x2._e)); }
    friend UpperBound<F> log2(Error<F> const& x) {
        return log(x)/cast_positive(log(Bounds<F>(2u,x.precision()))); }

    friend Bounds<F> pm(Error<F> const& x) { return Bounds<F>(-x._e,+x._e); }

    friend Bool same(Error<F> const& x1, Error<F> const& x2) { return x1._e==x2._e; }
    friend Bool same(Error<F> const& x1, Dyadic const& x2) { return x1._e==x2; }
    friend Bool refines(Error<F> const& x1, Error<F> const& x2) { return x1._e<=x2._e; }
    friend Error<F> refinement(Error<F> const& x1, Error<F> const& x2) { return Error<F>(min(x1._e,x2._e)); }
  public:
    friend Error<F>const& cast_positive(Error<F> const& x) { return x; }
  public:
    friend OutputStream& operator<<(OutputStream& os, Error<F> const& x) {
        return write(os,x.raw(),DecimalPrecision{Error<F>::output_places},upward); }
    friend InputStream& operator>>(InputStream& is, Error<F>& x) {
        UpperBound<F> xu; is >> xu; x=Error<F>(xu); return is; }
  public:
    static Nat output_places;
    static Void set_output_places(Nat p) { output_places=p; }
};

template<class PR> Error(ValidatedUpperNumber, PR) -> Error<RawFloatType<PR>>;
template<class F> Error(F) -> Error<F>;

extern template Ariadne::Nat Ariadne::Error<Ariadne::FloatDP>::output_places;
extern template Ariadne::Nat Ariadne::Error<Ariadne::FloatMP>::output_places;


}

#endif
