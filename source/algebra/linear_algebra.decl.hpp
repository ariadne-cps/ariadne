/***************************************************************************
 *            algebra/linear_algebra.decl.hpp
 *
 *  Copyright  2013-20  Pieter Collins
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

/*! \file algebra/linear_algebra.decl.hpp
 *  \brief
 */

#include "../numeric/number.decl.hpp"
#include "../numeric/float.decl.hpp"

#ifndef ARIADNE_LINEAR_ALGEBRA_DECL_HPP
#define ARIADNE_LINEAR_ALGEBRA_DECL_HPP

namespace Ariadne {

template<class X> using Scalar = X;
template<class X> class Vector;
template<class X> class Covector;
template<class X> class Matrix;
template<class X> class DiagonalMatrix;
template<class X> class SymmetricMatrix;

//@{
//! \relates Vector
//! \name Type synonyms
using DyadicVector = Vector<Dyadic>; //!< .
using RationalVector = Vector<Rational>; //!< .
using RealVector = Vector<Real>; //!< .
using FloatDPVector = Vector<FloatDP>; //!< .
using RawFloatDPVector = Vector<RawFloatDP>; //!< .
using FloatDPApproximationVector = Vector<FloatDPApproximation>; //!< .
using FloatDPBoundsVector = Vector<FloatDPBounds>; //!< .
using FloatDPBallVector = Vector<FloatDPBall>; //!< .
using FloatDPValueVector = Vector<FloatDPValue>; //!< .
using FloatMPVector = Vector<FloatMP>; //!< .
using RawFloatMPVector = Vector<RawFloatMP>; //!< .
using FloatMPApproximationVector = Vector<FloatMPApproximation>; //!< .
using FloatMPBoundsVector = Vector<FloatMPBounds>; //!< .
using FloatMPBallVector = Vector<FloatMPBall>; //!< .
using FloatMPValueVector = Vector<FloatMPValue>; //!< .
using FloatMPDPBallVector = Vector<FloatMPDPBall>; //!< .
//@}

//@{
//! \relates Covector
//! \name Type synonyms
using DyadicCovector = Covector<Dyadic>; //!< .
using RationalCovector = Covector<Rational>; //!< .
using RealCovector = Covector<Real>; //!< .
using FloatDPCovector = Covector<FloatDP>; //!< .
using RawFloatDPCovector = Covector<RawFloatDP>; //!< .
using FloatDPApproximationCovector = Covector<FloatDPApproximation>; //!< .
using FloatDPBoundsCovector = Covector<FloatDPBounds>; //!< .
using FloatDPBallCovector = Covector<FloatDPBall>; //!< .
using FloatDPValueCovector = Covector<FloatDPValue>; //!< .
using FloatMPCovector = Covector<FloatMP>; //!< .
using RawFloatMPCovector = Covector<RawFloatMP>; //!< .
using FloatMPApproximationCovector = Covector<FloatMPApproximation>; //!< .
using FloatMPBoundsCovector = Covector<FloatMPBounds>; //!< .
using FloatMPBallCovector = Covector<FloatMPBall>; //!< .
using FloatMPValueCovector = Covector<FloatMPValue>; //!< .
using FloatMPDPBallCovector = Covector<FloatMPDPBall>; //!< .
//@}


//@{
//! \relates Matrix
//! \name Type synonyms
using DyadicMatrix = Matrix<Dyadic>; //!< .
using RationalMatrix = Matrix<Rational>; //!< .
using RealMatrix = Matrix<Real>; //!< .
using FloatDPMatrix = Matrix<FloatDP>; //!< .
using RawFloatDPMatrix = Matrix<RawFloatDP>; //!< .
using FloatDPApproximationMatrix = Matrix<FloatDPApproximation>; //!< .
using FloatDPBoundsMatrix = Matrix<FloatDPBounds>; //!< .
using FloatDPBallMatrix = Matrix<FloatDPBall>; //!< .
using FloatDPValueMatrix = Matrix<FloatDPValue>; //!< .
using FloatMPApproximationMatrix = Matrix<FloatMPApproximation>; //!< .
using FloatMPBoundsMatrix = Matrix<FloatMPBounds>; //!< .
using FloatMPBallMatrix = Matrix<FloatMPBall>; //!< .
using FloatMPValueMatrix = Matrix<FloatMPValue>; //!< .
using FloatMPDPBallMatrix = Matrix<FloatMPDPBall>; //!< .
//@}


//@{
//! \relates DiagonalMatrix
//! \name Type synonyms
using FloatDPApproximationDiagonalMatrix = DiagonalMatrix<FloatDPApproximation>; //!< .
using FloatDPBoundsDiagonalMatrix = DiagonalMatrix<FloatDPBounds>; //!< .
using FloatDPBallDiagonalMatrix = DiagonalMatrix<FloatDPBall>; //!< .
using FloatDPValueDiagonalMatrix = DiagonalMatrix<FloatDPValue>; //!< .
using FloatMPApproximationDiagonalMatrix = DiagonalMatrix<FloatMPApproximation>; //!< .
using FloatMPBoundsDiagonalMatrix = DiagonalMatrix<FloatMPBounds>; //!< .
using FloatMPBallDiagonalMatrix = DiagonalMatrix<FloatMPBall>; //!< .
using FloatMPValueDiagonalMatrix = DiagonalMatrix<FloatMPValue>; //!< .
using FloatMPDPBallDiagonalMatrix = DiagonalMatrix<FloatMPDPBall>; //!< .
//@}

} // namespace Ariadne

#endif
