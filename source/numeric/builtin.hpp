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

#include "twoexp.hpp"

namespace Ariadne {

struct ExactTag;
enum class Comparison : char;

class ApproximateDouble {
    double _d;
  public:
    typedef ApproximateTag Paradigm;
    ApproximateDouble() : _d() { }
    ApproximateDouble(double d) : _d(d) { }
    ApproximateDouble(ExactDouble const& x);
    ApproximateDouble(Real const& r0);
    friend ApproximateDouble operator+(ApproximateDouble x) { return ApproximateDouble(+x._d); }
    friend ApproximateDouble operator-(ApproximateDouble x) { return ApproximateDouble(-x._d); }
    friend ApproximateDouble operator*(ApproximateDouble x, TwoExp p) { return ApproximateDouble(x.get_d()*p.get_d()); }
    friend ApproximateDouble operator/(ApproximateDouble x, TwoExp p) { return ApproximateDouble(x.get_d()/p.get_d()); }
    friend ApproximateDouble& operator+=(ApproximateDouble& x1, ApproximateDouble x2) { x1._d+=x2._d; return x1; }
    operator double() const { return _d; }
    double get_d() const { return this->_d; }
    friend ApproximateDouble operator"" _a (long double lx) { double x=lx; return ApproximateDouble(x); }
};
inline ApproximateDouble operator"" _a (long double lx);

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
    friend ExactDouble nul(ExactDouble x) { return ExactDouble(0.0); }
    friend ExactDouble abs(ExactDouble x) { return ExactDouble(std::abs(x._d)); }
    friend ExactDouble operator+(ExactDouble x) { return ExactDouble(+x._d); }
    friend ExactDouble operator-(ExactDouble x) { return ExactDouble(-x._d); }
//    friend ApproximateDouble operator-(ApproximateDouble x1, ApproximateDouble x2);
    friend ExactDouble operator*(ExactDouble x, TwoExp p) { return ExactDouble(x.get_d()*p.get_d()); }
    friend ExactDouble operator/(ExactDouble x, TwoExp p) { return ExactDouble(x.get_d()/p.get_d()); }
    friend Comparison cmp(ExactDouble const& x1, ExactDouble const& x2) {
        return Comparison( (x1._d==x2._d) ? 0 : (x1._d<x2._d) ? -1 : +1 ); }
    friend Bool operator==(ExactDouble const& x1, ExactDouble const& x2) { return x1._d==x2._d; }
    friend Bool operator!=(ExactDouble const& x1, ExactDouble const& x2) { return x1._d!=x2._d; }
    friend Bool operator>=(ExactDouble const& x1, ExactDouble const& x2) { return x1._d>=x2._d; }
    friend Bool operator<=(ExactDouble const& x1, ExactDouble const& x2) { return x1._d<=x2._d; }
    friend Bool operator> (ExactDouble const& x1, ExactDouble const& x2) { return x1._d> x2._d; }
    friend Bool operator< (ExactDouble const& x1, ExactDouble const& x2) { return x1._d< x2._d; }
    friend ExactDouble operator"" _x (long double lx) { double x=lx; if (x!=lx) { std::cerr<<"lx="<<lx; assert(x==lx); } return ExactDouble(x); }
    friend ExactDouble operator"" _pr (long double lx) { double x=lx; return ExactDouble(x); }
    friend OutputStream& operator<<(OutputStream& os, ExactDouble x) { return os << std::setprecision(18) << x.get_d(); }
};
inline ExactDouble operator"" _x (long double lx);
inline ExactDouble operator"" _pr (long double lx);
inline ExactDouble cast_exact(double d) { return ExactDouble(d); }
inline ExactDouble cast_exact(ApproximateDouble ax) { return ExactDouble(ax.get_d()); }

static const ExactDouble inf = ExactDouble(std::numeric_limits<double>::infinity());


template<class X> class Positive;
template<> class Positive<ExactDouble> : public ExactDouble {
  public:
    explicit Positive<ExactDouble>(ExactDouble x) : ExactDouble(x) { }
};
inline Positive<ExactDouble> cast_positive(ExactDouble x) { return Positive<ExactDouble>(x); }

inline ApproximateDouble::ApproximateDouble(ExactDouble const& x) : _d(x.get_d()) { }

} // namespace Ariadne

#endif
