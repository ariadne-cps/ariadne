/***************************************************************************
 *            float_error.hpp
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

/*! \file float_error.hpp
 *  \brief Floating-point error bounds for metric spaces.
 */

#ifndef ARIADNE_FLOAT_ERROR_HPP
#define ARIADNE_FLOAT_ERROR_HPP

#include "../utility/macros.hpp"

#include "number.decl.hpp"
#include "float.decl.hpp"

namespace Ariadne {

template<class F> struct NumericTraits<Error<F>> {
    typedef PositiveValidatedUpperNumber GenericType;
    typedef Error<F> PositiveType;
    typedef PositiveLowerBound<F> OppositeType;
    typedef ValidatedLowerKleenean LessType;
    typedef ValidatedNegatedSierpinskian EqualsType;
};

//! \ingroup NumericModule
//! \brief Floating-point upper bounds for positive real numbers, suitable for use as an upper bound for an error in a metric space.
template<class F> class Error
    : public DispatchDirectedFloatOperations<UpperBound<F>>
    , public DispatchPositiveDirectedNumericOperations<PositiveUpperBound<F>,PositiveLowerBound<F>>
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
    Error<F>(PositiveBounds<F> const& x) : _e(x._u) { }
    Error<F>(PositiveUpperBound<F> const& x) : _e(x._u) { }
    operator PositiveUpperBound<F> const& () const { return reinterpret_cast<PositiveUpperBound<F>const&>(*this); }
    operator PositiveUpperBound<F>& () { return reinterpret_cast<PositiveUpperBound<F>&>(*this); }
  public:
    Error<F>() : _e() { }
    explicit Error<F>(PR const& pr) : _e(pr) { }
    explicit Error<F>(F const& x) : _e(x) { ARIADNE_PRECONDITION_MSG((this->_e>=0),"e="<<*this); }
    template<class M, EnableIf<IsBuiltinUnsignedIntegral<M>> =dummy> Error<F>(M m, PR pr) : _e(m,pr) { }
    explicit Error<F>(UpperBound<F> const& x) : Error<F>(x._u) { }
    explicit Error<F>(ValidatedUpperNumber const& y, PR pr) : Error<F>(UpperBound<F>(y,pr)) { }
    explicit Error<F>(const TwoExp& t, PR pr) : Error<F>(UpperBound<F>(t,pr)) { }
    Error<F>(PositiveValue<F> const& x) : _e(x._v) { }
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
    friend Error<F> mag(Error<F> const& x) { return x; }
    friend UpperBound<F> operator+(Error<F> const& x) { return UpperBound<F>(+x._e); }
    friend LowerBound<F> operator-(Error<F> const& x) { return LowerBound<F>(-x._e); }
    friend UpperBound<F> operator+(Value<F> const& x1, Error<F> const& x2) { return UpperBound<F>(add(up,x1._v,x2._e)); }
    friend LowerBound<F> operator-(Value<F> const& x1, Error<F> const& x2) { return LowerBound<F>(sub(down,x1._v,x2._e)); }
    friend UpperBound<F> log2(Error<F> const& x) {
        return log(x)/cast_positive(log(Bounds<F>(2u,x.precision()))); }

    friend Bounds<F> pm(Error<F> const& x) { return Bounds<F>(-x._e,+x._e); }
    
    friend Bool same(Error<F> const& x1, Error<F> const& x2) { return x1._e==x2._e; }
    friend Bool refines(Error<F> const& x1, Error<F> const& x2) { return x1._e<=x2._e; }
    friend Error<F> refinement(Error<F> const& x1, Error<F> const& x2) { return Error<F>(min(x1._e,x2._e)); }
    friend OutputStream& operator<<(OutputStream& os, Error<F> const& x) { return Operations<Error<F>>::_write(os,x); }
  public:
    static Nat output_places;
    static Void set_output_places(Nat p) { output_places=p; }
};

}

#endif
