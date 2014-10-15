/***************************************************************************
 *            utility/builtin.h
 *
 *  Copyright 2013-14  Pieter Collins
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

/*! \file utility/builtin.h
 *  \brief 
 */



#ifndef ARIADNE_BUILTIN_H
#define ARIADNE_BUILTIN_H

#include "metaprogramming.h"

namespace Ariadne {

class builtin_uint {
    long unsigned int _m;
  public:
    template<class M> builtin_uint(M m, EnableIf<IsUnsigned<M>> =dummy) : _m(m) { }
    explicit operator long unsigned int() const { return _m; }
    long unsigned int get_ui() const { return _m; }
};

class builtin_int {
    long int _n;
  public:
    template<class N> builtin_int(N n, EnableIf<IsIntegral<N>> =dummy) : _n(n) { }
    explicit operator long int () const { return _n; }
    long int get_si() const { return _n; }
};

class builtin_double {
    double _x;
  public:
    template<class X> builtin_double(X x, EnableIf<IsFloatingPoint<X>> =dummy) : _x(x) { }
    explicit operator double () const { return _x; }
    double get_d() const { return _x; }
};

class exact_double {
    double _x;
  private:
    exact_double(double x) : _x(x) { }
  public:
    friend exact_double make_exact(double);
    exact_double() : _x(0.0) { }
    explicit operator double () const { return _x; }
    double get_d() const { return _x; }
};

exact_double make_exact(double x) { return exact_double(x); }

} // namespace Ariadne

#endif

