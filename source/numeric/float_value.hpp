/***************************************************************************
 *            numeric/float_value.hpp
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

/*! \file numeric/float_value.hpp
 *  \brief Exact Floating-point representations of real numbers.
 */

#ifndef ARIADNE_FLOAT_VALUE_HPP
#define ARIADNE_FLOAT_VALUE_HPP

#include "logical.decl.hpp"
#include "number.decl.hpp"
#include "float.decl.hpp"

#include "float_traits.hpp"
#include "float_operations.hpp"
#include "float_factory.hpp"

#include "integer.hpp"
#include "dyadic.hpp"

#include "floatdp.hpp"
#include "floatmp.hpp"
#include "float_bounds.hpp"

namespace Ariadne {



template<ARawFloat F> class Positive<F> : public F {
    using PR = typename F::PrecisionType;
  public:
    Positive() : Value<F>() { }
    explicit Positive(PR const& pr) : Value<F>(pr) { }
    template<BuiltinUnsignedIntegral M> Positive(M m, PR pr) : Value<F>(m,pr) { }
    Positive(TwoExp const& ex, PR pr) : Value<F>(ex,pr) { }
    explicit Positive(Dyadic const& w, PR pr) : F(w,pr) { }
    explicit Positive(F const& x) : F(x) { }
  public:
    friend Positive<Bounds<F>> operator+(Positive<F> const& x1, Positive<F> const& x2) {
        return cast_positive(cast_unsigned(x1)+cast_unsigned(x2)); }
    friend Positive<Bounds<F>> operator*(Positive<F> const& x1, Positive<F> const& x2) {
        return cast_positive(cast_unsigned(x1)*cast_unsigned(x2)); }
    friend Positive<Bounds<F>> operator/(Positive<F> const& x1, Positive<F> const& x2) {
        return cast_positive(cast_unsigned(x1)/cast_unsigned(x2)); }
    friend Positive<F> max(Positive<F> const& x1, Positive<F> const& x2) {
        return cast_positive(max(cast_unsigned(x1),cast_unsigned(x2))); }
    friend Positive<F> max(Positive<F> const& x1, F const& x2) {
        return cast_positive(max(cast_unsigned(x1),x2)); }
    friend Positive<F> max(F const& x1, Positive<F> const& x2) {
        return cast_positive(max(x1,cast_unsigned(x2))); }
    friend Positive<F> min(Positive<F> const& x1, Positive<F> const& x2) {
        return cast_positive(min(cast_unsigned(x1),cast_unsigned(x2))); }

    friend Positive<F> nul(Positive<F> const& x) { return Positive<F>(nul(static_cast<F const&>(x))); }
    friend Positive<F> pos(Positive<F> const& x) { return Positive<F>(pos(static_cast<F const&>(x))); }
    friend Positive<F> hlf(Positive<F> const& x) { return Positive<F>(hlf(static_cast<F const&>(x))); }
    friend Positive<Bounds<F>> pow(Positive<F> const& x, Nat m) {
        return pow(Positive<Bounds<F>>(x),m); }
    friend Positive<Bounds<F>> pow(Positive<F> const& x, Int n) {
        return pow(Positive<Bounds<F>>(x),n); }

    friend Error<F> operator+(Positive<F> const& v1, Error<F> const& e2);
};

template<class F> inline Positive<F> cast_positive(F const& x) {
    return Positive<F>(x);
}
template<class F> inline Positive<F> cast_positive(Positive<F> const& x) {
    return x;
}
template<class F> inline F const& cast_unsigned(Positive<F> const& x) {
    return static_cast<F const&>(x);
}

}


#endif
