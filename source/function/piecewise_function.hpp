/***************************************************************************
 *            piecewise_function.hpp
 *
 *  Copyright 2018  Pieter Collins
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

/*! \file piecewise_function.hpp
 *  \brief Piecewise-defined functions
 */

#ifndef ARIADNE_PIECEWISE_FUNCTION_HPP
#define ARIADNE_PIECEWISE_FUNCTION_HPP

#include "geometry/interval.hpp"

namespace Ariadne {


template<class B, class F> class UnivariatePiecewiseFunction {
    typedef typename F::Paradigm P;
  private:
    List<B> _bounds;
    List<F> _pieces;
  public:
    UnivariatePiecewiseFunction(List<B> b, List<F> p) : _bounds(b), _pieces(p) { ARIADNE_PRECONDITION(b.size()==p.size()+1u); }
    Void append(F const& f, B const& ub) { ARIADNE_PRECONDITION(_bounds.back()<ub); _bounds.append(ub); _pieces.append(f); }

    Interval<B> const& bounds(SizeType i) const { return Interval<B>(_bounds[i],_bounds[i+1]); }
    B const& bound(SizeType i) const { return _bounds[i]; }
    F const& piece(SizeType i) const { return _pieces[i]; }

    template<class X> ResultOf<F(X)> operator() (X const& x) const { return this->_evaluate(x); }

    friend OutputStream& operator<<(OutputStream& os, UnivariatePiecewiseFunction<B,F> const& f) { return f._write(os); }
  private:
    OutputStream& _write(OutputStream& os) const;
    template<class X> ResultOf<F(X)> _evaluate(X const& x) const;
    template<class X> SizeType _piece(X const& x) const;
};


template<class B, class X> SizeType find_piece(List<B> const& bounds, X const& x, ExactTag);
template<class B, class X> SizeType find_piece(List<B> const& bounds, X const& x, ValidatedTag);
template<class B, class X> SizeType find_piece(List<B> const& bounds, X const& x, ValidatedTag);

template<class B, class F> template<class X> inline ResultOf<F(X)> UnivariatePiecewiseFunction<B,F>::_evaluate(X const& x) const {
    SizeType i=_piece(x); return _pieces[i](x); }

template<class B, class F> template<class X> inline SizeType UnivariatePiecewiseFunction<B,F>::_piece(X const& x) const {
    using P=typename X::Paradigm;  find_piece(_bounds,x,P()); }

template<class B, class X> SizeType find_piece(List<B> const& bounds, X const& x, ExactTag) {
    ARIADNE_PRECONDITION(bounds.front()<=x && x<=bounds.back());
    // TODO: Use binary search
    SizeType i=0;
    while (x>=bounds[i+1u]) { ++i; }
    return i;
}

template<class B, class X> SizeType find_piece(List<B> const& bounds, X const& x, ValidatedTag) {
    ARIADNE_PRECONDITION(definitely(bounds.front()<=x) && definitely(x<=bounds.back()));
    // TODO: Use binary search
    // FIXME: Allow multiple intervals
    SizeType i=0;
    while (definitely(x>=bounds[i+1u])) { ++i; }
    return i;
}

} // namespace Ariadne

#endif
