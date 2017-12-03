/***************************************************************************
 *            evaluate.hpp
 *
 *  Copyright 2008-17  Pieter Collins
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

/*! \file evaluate.hpp
 *  \brief Functions for evaluating polynomial expansions.
 */

#ifndef ARIADNE_EVALUATE_HPP
#define ARIADNE_EVALUATE_HPP

namespace Ariadne {

template<class I, class X> class Expansion;
template<class A> class Vector;
class ReverseGradedLexicographicLess;
class ReverseLexicographicKeyLess;

//! \ingroup FunctionModule
//! \brief Evaluate a power series on an algebra using Horner's rule.
template<class X, class A> A horner_evaluate(const Expansion<MultiIndex,X>& e, const Vector<A>& y);

//! \ingroup FunctionModule
//! \brief Evaluate a power series on an algebra using termwise expansion (slow).
template<class X, class A> A power_evaluate(const Expansion<MultiIndex,X>& e, const Vector<A>& y);

template<class X, class A> A evaluate(const Expansion<MultiIndex,X>& e, const Vector<A>& y);

template<class X, class A> A simple_evaluate(const Expansion<MultiIndex,X>& e, const Vector<A>& y);


template<class X, class A> Vector<A> evaluate(const Vector< Expansion<MultiIndex,X> >& e, const Vector<A>& y);

} // namespace Ariadne

#include "evaluate.tpl.hpp"

#endif /* ARIADNE_EVALUATE_HPP */
