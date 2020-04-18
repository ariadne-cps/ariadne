/***************************************************************************
 *            numeric/builtin.hpp
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

/*! \file numeric/builtin.hpp
 *  \brief Inclusion header for wrapper for builtin numbers.
 */

#ifndef ARIADNE_BUILTIN_HPP
#define ARIADNE_BUILTIN_HPP

#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>

#include "../utility/metaprogramming.hpp"
#include "../numeric/number.decl.hpp"


namespace Ariadne {

struct ExactTag;
enum class Comparison : char;

class ApproximateDouble {
    double _d;
  public:
    typedef ApproximateTag Paradigm;
    ApproximateDouble() : _d() { }
    ApproximateDouble(int n) : _d(n) { }
    template<class X, EnableIf<IsBuiltinArithmetic<X>> =dummy> ApproximateDouble(X const& x) : _d(x) { }
    template<class X, DisableIf<IsBuiltinArithmetic<X>> =dummy> ApproximateDouble(X const& x) : _d(x.get_d()) { }
    friend ApproximateDouble operator+(ApproximateDouble x) { return ApproximateDouble(+x._d); }
    friend ApproximateDouble operator-(ApproximateDouble x) { return ApproximateDouble(-x._d); }
    operator double() const { return _d; }
};

//! \ingroup NumericModule
//! \brief A wrapper around a builtin double-precision floating-point number,
//! indicating that the stored value is the \em exact value of a real quantity.
class ExactDouble {
    double _d;
  public:
    typedef ExactTag Paradigm;
    double get_d() const { return this->_d; }
    ExactDouble() : _d() { }
    template<class N, EnableIf<IsBuiltinIntegral<N>> =dummy> ExactDouble(N n) : _d(n) { assert(_d==n); }
    template<class X, EnableIf<IsBuiltinFloatingPoint<X>> =dummy> explicit ExactDouble(X const& x) : _d(x) { assert(std::isnan(_d) || (_d==x)); }
    static ExactDouble infinity() { return ExactDouble(std::numeric_limits<double>::infinity()); }
    operator ExactNumber() const;
    friend ExactDouble operator+(ExactDouble x) { return ExactDouble(+x._d); }
    friend ExactDouble operator-(ExactDouble x) { return ExactDouble(-x._d); }
    friend Comparison cmp(ExactDouble const& x1, ExactDouble const& x2) {
        return Comparison( (x1.get_d()==x2.get_d()) ? 0 : (x1.get_d()<x2.get_d()) ? -1 : +1 ); }
    friend ExactDouble operator"" _x (long double lx) { double x=lx; assert(x==lx); return ExactDouble(x); }
    friend OutputStream& operator<<(OutputStream& os, ExactDouble x) { return os << std::setprecision(18) << x.get_d(); }
};
inline ExactDouble operator"" _x (long double lx);

} // namespace Ariadne

#endif
