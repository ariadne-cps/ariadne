/***************************************************************************
 *            numeric/sequence.hpp
 *
 *  Copyright 2013-17  Pieter Collins
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

/*! \file numeric/sequence.hpp
 *  \brief
 */

#ifndef ARIADNE_SEQUENCE_HPP
#define ARIADNE_SEQUENCE_HPP

#include <functional>

namespace Ariadne {

class Natural;

template<class X> class Sequence {
    std::function<X(Natural)> _fn;
  public:
    Sequence<X>(std::function<X(Natural)> fn) : _fn(fn) { }
    X operator[](Natural const& n) const { return _fn(n); }
};

template<class Y> struct CompletionTypedef;
template<> struct CompletionTypedef<Dyadic> { typedef Real Type; };
template<> struct CompletionTypedef<Rational> { typedef Real Type; };
template<> struct CompletionTypedef<Real> { typedef Real Type; };
template<class Y> using CompletionType = typename CompletionTypedef<Y>::Type;

template<class X> class ConvergentSequence : public Sequence<X> {
  public:
    ConvergentSequence(std::function<X(Natural)> fn) : Sequence<X>(fn) { }
    ConvergentSequence(Sequence<X> const& seq) : Sequence<X>(seq) { }
};
template<class X> class StrongCauchySequence : public Sequence<X> {
  public:
    StrongCauchySequence(std::function<X(Natural)> fn) : Sequence<X>(fn) { }
    StrongCauchySequence(Sequence<X> const& seq) : Sequence<X>(seq) { }
    friend CompletionType<X> limit(StrongCauchySequence<X> const&);
};

} // namespace Ariadne

#endif
