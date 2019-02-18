/***************************************************************************
 *            polynomial.cpp
 *
 *  Copyright 2008--17  Pieter Collins
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

#include "../numeric/numeric.hpp"
#include "../config.hpp"

#include "../geometry/interval.hpp"
#include "../function/polynomial.hpp"
#include "../function/polynomial.tpl.hpp"

namespace Ariadne {

template class MultivariatePolynomial<FloatDP>;
template class MultivariatePolynomial<FloatDPApproximation>;
template class MultivariatePolynomial<FloatDPBounds>;
template class MultivariatePolynomial<UpperIntervalType>;

template struct AlgebraOperations<MultivariatePolynomial<FloatDP>>;
template struct AlgebraOperations<MultivariatePolynomial<FloatDPApproximation>>;
template struct AlgebraOperations<MultivariatePolynomial<FloatDPBounds>>;
// template struct AlgebraOperations<MultivariatePolynomial<UpperIntervalType>>;

template<> Void MultivariatePolynomial<FloatDPValue>::cleanup() { }

template MultivariatePolynomial<FloatDPValue>::MultivariatePolynomial(SizeType);
template Expansion<MultiIndex,FloatDPValue>& MultivariatePolynomial<FloatDPValue>::expansion();
template OutputStream& MultivariatePolynomial<FloatDPValue>::_write(OutputStream&) const;
template OutputStream& MultivariatePolynomial<FloatDPValue>::_write(OutputStream&, List<String> const&) const;


template class MultivariatePolynomial<FloatMPApproximation>;
template class MultivariatePolynomial<FloatMPBounds>;
template struct AlgebraOperations<MultivariatePolynomial<FloatMPApproximation>>;
template struct AlgebraOperations<MultivariatePolynomial<FloatMPBounds>>;

} //namespace Ariadne


