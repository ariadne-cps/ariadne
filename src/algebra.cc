/***************************************************************************
 *            algebra.cc
 *
 *  Copyright 2011  Pieter Collins
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

#include "algebra_mixin.tcc"

namespace Ariadne {

template GradedAlgebra<Float> compose(const Series<Float>&, const GradedAlgebra<Float>&);
template GradedAlgebra<Interval> compose(const Series<Interval>&, const GradedAlgebra<Interval>&);
template GradedAlgebra<Real> compose(const Series<Real>&, const GradedAlgebra<Real>&);

template NormedAlgebra<Float> rec(const NormedAlgebra<Float>&);
template NormedAlgebra<Interval> rec(const NormedAlgebra<Interval>&);
template NormedAlgebra<Real> rec(const NormedAlgebra<Real>&);

template NormedAlgebra<Float> sqrt(const NormedAlgebra<Float>&);
template NormedAlgebra<Interval> sqrt(const NormedAlgebra<Interval>&);
template NormedAlgebra<Real> sqrt(const NormedAlgebra<Real>&);

template NormedAlgebra<Float> exp(const NormedAlgebra<Float>&);
template NormedAlgebra<Interval> exp(const NormedAlgebra<Interval>&);
template NormedAlgebra<Real> exp(const NormedAlgebra<Real>&);

template NormedAlgebra<Float> log(const NormedAlgebra<Float>&);
template NormedAlgebra<Interval> log(const NormedAlgebra<Interval>&);
template NormedAlgebra<Real> log(const NormedAlgebra<Real>&);

template NormedAlgebra<Float> sin(const NormedAlgebra<Float>&);
template NormedAlgebra<Interval> sin(const NormedAlgebra<Interval>&);
template NormedAlgebra<Real> sin(const NormedAlgebra<Real>&);

template NormedAlgebra<Float> cos(const NormedAlgebra<Float>&);
template NormedAlgebra<Interval> cos(const NormedAlgebra<Interval>&);
template NormedAlgebra<Real> cos(const NormedAlgebra<Real>&);

template NormedAlgebra<Float> tan(const NormedAlgebra<Float>&);
template NormedAlgebra<Interval> tan(const NormedAlgebra<Interval>&);
template NormedAlgebra<Real> tan(const NormedAlgebra<Real>&);

/*
template NormedAlgebra<Float> asin(const NormedAlgebra<Float>&);
template NormedAlgebra<Interval> asin(const NormedAlgebra<Interval>&);
template NormedAlgebra<Real> asin(const NormedAlgebra<Real>&);

template NormedAlgebra<Float> acos(const NormedAlgebra<Float>&);
template NormedAlgebra<Interval> acos(const NormedAlgebra<Interval>&);
template NormedAlgebra<Real> acos(const NormedAlgebra<Real>&);

template NormedAlgebra<Float> atan(const NormedAlgebra<Float>&);
template NormedAlgebra<Interval> atan(const NormedAlgebra<Interval>&);
template NormedAlgebra<Real> atan(const NormedAlgebra<Real>&);
*/

} // namespace Ariadne
