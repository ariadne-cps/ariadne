/***************************************************************************
 *            float_error.h
 *
 *  Copyright 2008-16  Pieter Collins
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

/*! \file float_error.h
 *  \brief Floating-point error bounds for metric spaces.
 */

#ifndef ARIADNE_FLOAT_ERROR_H
#define ARIADNE_FLOAT_ERROR_H

#include "utility/macros.h"

#include "number.decl.h"
#include "float.decl.h"

namespace Ariadne {

//! \ingroup NumericModule
//! \brief Floating-point upper bounds for positive real numbers, suitable for use as an upper bound for an error in a metric space.
template<class PR> class FloatError
    : public DispatchDirectedFloatOperations<FloatUpperBound<PR>>
    , public DispatchPositiveDirectedNumericOperations<PositiveFloatUpperBound<PR>,PositiveFloatLowerBound<PR>>
    , public ProvideConcreteGenericDirectedSemiFieldOperations<PositiveFloatUpperBound<PR>,PositiveFloatLowerBound<PR>,Nat,Nat>
{
  private: public:
    RawFloat<PR> _e;
  public:
    typedef PR PrecisionType;
  public:
    FloatError<PR>(PositiveFloatUpperBound<PR> const& x) : _e(x._u) { }
    operator PositiveFloatUpperBound<PR> const& () const { return reinterpret_cast<PositiveFloatUpperBound<PR>const&>(*this); }
    operator PositiveFloatUpperBound<PR>& () { return reinterpret_cast<PositiveFloatUpperBound<PR>&>(*this); }
  public:
    FloatError<PR>() : _e() { }
    explicit FloatError<PR>(PR const& pr) : _e(pr) { }
    explicit FloatError<PR>(RawFloat<PR> const& x) : _e(x) { ARIADNE_PRECONDITION_MSG((this->_e>=0),"e="<<*this); }
    template<class M, EnableIf<IsUnsignedIntegral<M>> =dummy> FloatError<PR>(M m, PR pr) : _e(m,pr) { }
    explicit FloatError<PR>(FloatUpperBound<PR> const& x) : FloatError<PR>(x._u) { }
    explicit FloatError<PR>(ValidatedUpperNumber const& y, PR pr) : FloatError(FloatUpperBound<PR>(y,pr)) { }
    FloatError<PR>(PositiveFloatValue<PR> const& x) : _e(x._v) { }
    FloatError<PR>& operator=(Nat m) { reinterpret_cast<FloatUpperBound<PR>&>(*this)=m; return *this; }
  public:
    PrecisionType precision() const { return _e.precision(); }
    RawFloat<PR> const& raw() const { return _e; }
    RawFloat<PR>& raw() { return _e; }
  public:
    friend FloatError<PR> mag(FloatError<PR> const& x) { return x; }
    friend FloatUpperBound<PR> operator+(FloatError<PR> const& x) { return FloatUpperBound<PR>(+x._e); }
    friend FloatLowerBound<PR> operator-(FloatError<PR> const& x) { return FloatLowerBound<PR>(-x._e); }
    friend FloatUpperBound<PR> operator+(FloatValue<PR> const& x1, FloatError<PR> const& x2) { return FloatUpperBound<PR>(add_up(x1._v,x2._e)); }
    friend FloatLowerBound<PR> operator-(FloatValue<PR> const& x1, FloatError<PR> const& x2) { return FloatLowerBound<PR>(sub_down(x1._v,x2._e)); }

    friend Bool same(FloatError<PR> const& x1, FloatError<PR> const& x2) { return x1._e==x2._e; }
    friend Bool refines(FloatError<PR> const& x1, FloatError<PR> const& x2) { return x1._e<=x2._e; }
    friend FloatError<PR> refinement(FloatError<PR> const& x1, FloatError<PR> const& x2) { return FloatError<PR>(min(x1._e,x2._e)); }
    friend OutputStream& operator<<(OutputStream& os, FloatError<PR> const& x) { return Operations<FloatError<PR>>::_write(os,x); }
  public:
    static Nat output_places;
    static Void set_output_places(Nat p) { output_places=p; }
};

}

#endif
