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

namespace Ariadne {

struct ExactTag;

using ApproximateDouble = double;

class ExactDouble {
    double _d;
  public:
    typedef ExactTag Paradigm;
    double get_d() const { return this->_d; }
    template<class N, EnableIf<IsIntegral<N>> =dummy> explicit ExactDouble(N n) : _d(n) { assert(_d==n); }
    explicit ExactDouble(double d) : _d(d) { }
    operator ExactNumber() const;
    friend ExactDouble operator+(ExactDouble x) { return ExactDouble(+x._d); }
    friend ExactDouble operator-(ExactDouble x) { return ExactDouble(-x._d); }
    friend ExactDouble operator"" _x (long double lx) { double x=lx; assert(x==lx); return ExactDouble(x); }
};
inline ExactDouble operator"" _x (long double lx);

} // namespace Ariadne

#endif
