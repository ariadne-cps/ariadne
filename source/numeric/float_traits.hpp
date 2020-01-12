/***************************************************************************
 *            numeric/float_traits.hpp
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

/*! \file numeric/float_traits.hpp
 *  \brief Type traits for floating-point numbers.
 */

#ifndef ARIADNE_FLOAT_TRAITS_HPP
#define ARIADNE_FLOAT_TRAITS_HPP

#include "float.decl.hpp"

namespace Ariadne {

template<class F> struct NumericTraits<Approximation<F>> {
    typedef ApproximateNumber GenericType;
    typedef PositiveApproximation<F> PositiveType;
    typedef Approximation<F> OppositeType;
    typedef Fuzzy LessType;
    typedef Fuzzy EqualsType;
};

template<class F> struct NumericTraits<LowerBound<F>> {
    typedef ValidatedLowerNumber GenericType;
    typedef UpperBound<F> OppositeType;
    typedef PositiveLowerBound<F> PositiveType;
    typedef ValidatedUpperKleenean LessType;
    typedef ValidatedNegatedSierpinskian EqualsType;
};

template<class F> struct NumericTraits<UpperBound<F>> {
    typedef ValidatedUpperNumber GenericType;
    typedef LowerBound<F> OppositeType;
    typedef PositiveUpperBound<F> PositiveType;
    typedef ValidatedLowerKleenean LessType;
    typedef ValidatedNegatedSierpinskian EqualsType;
};

template<class F> struct NumericTraits<Bounds<F>> {
    typedef ValidatedNumber GenericType;
    typedef Bounds<F> OppositeType;
    typedef PositiveBounds<F> PositiveType;
    typedef ValidatedKleenean LessType;
    typedef ValidatedKleenean EqualsType;
};

template<class F, class FE> struct NumericTraits<Ball<F,FE>> {
    typedef ValidatedNumber GenericType;
    typedef Ball<F,FE> OppositeType;
    typedef PositiveBall<F,FE> PositiveType;
    typedef ValidatedKleenean LessType;
    typedef ValidatedKleenean EqualsType;
};

template<class F> struct NumericTraits<Value<F>> {
    typedef ExactNumber GenericType;
    typedef Value<F> OppositeType;
    typedef PositiveValue<F> PositiveType;
    typedef Boolean LessType;
    typedef Boolean EqualsType;
};

template<class F> struct NumericTraits<Error<F>> {
    typedef PositiveValidatedUpperNumber GenericType;
    typedef Error<F> PositiveType;
    typedef PositiveLowerBound<F> OppositeType;
    typedef ValidatedLowerKleenean LessType;
    typedef ValidatedNegatedSierpinskian EqualsType;
};

} // namespace Ariadne

#endif
