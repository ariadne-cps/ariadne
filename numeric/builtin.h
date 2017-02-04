/***************************************************************************
 *            builtin.h
 *
 *  Copyright 2008-15  Pieter Collins
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
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

/*! \file builtin.h
 *  \brief Inclusion header for wrapper for builtin numbers.
 */

#ifndef ARIADNE_BUILTIN_H
#define ARIADNE_BUILTIN_H

#include <cassert>
#include <iostream>
#include <iomanip>
#include <limits>

#include "utility/metaprogramming.h"
#include "numeric/number.decl.h"


namespace Ariadne {

struct ExactTag;
enum class Comparison : char;

class ApproximateDouble {
    double _d;
  public:
    typedef ApproximateTag Paradigm;
    ApproximateDouble(int n) : _d(n) { }
    template<class X, EnableIf<IsBuiltinArithmetic<X>> =dummy> ApproximateDouble(X const& x) : _d(x) { }
    template<class X, DisableIf<IsBuiltinArithmetic<X>> =dummy> ApproximateDouble(X const& x) : _d(x.get_d()) { }
    friend ApproximateDouble operator+(ApproximateDouble x) { return ApproximateDouble(+x._d); }
    friend ApproximateDouble operator-(ApproximateDouble x) { return ApproximateDouble(-x._d); }
    operator double() const { return _d; }
};

class ExactDouble {
    double _d;
  public:
    typedef ExactTag Paradigm;
    double get_d() const { return this->_d; }
    template<class N, EnableIf<IsBuiltinIntegral<N>> =dummy> ExactDouble(N n) : _d(n) { assert(_d==n); }
    template<class X, EnableIf<IsBuiltinFloatingPoint<X>> =dummy> explicit ExactDouble(X const& x) : _d(x) { assert(_d==x); }
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
