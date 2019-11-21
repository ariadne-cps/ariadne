/***************************************************************************
 *            pivot_matrix.hpp
 *
 *  Copyright 2005-19  Alberto Casagrande, Pieter Collins
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

#ifndef ARIADNE_PIVOT_MATRIX_HPP
#define ARIADNE_PIVOT_MATRIX_HPP

namespace Ariadne {

//! \ingroup LinearAlgebraModule
//! \brief A matrix describing a permutation in terms of successive pivots.
//! \details The matrix is formed by an array \f$p_i\f$ such that
//! the product \f$P A\f$ is formed by successively swapping
//! the \f$i\f$-th row of \f$A\f$ with the \f$p_i\f$-th row.
//! \relates Matrix<X>
class PivotMatrix {
    Array<SizeType> _ary;
  public:
    PivotMatrix(SizeType n=0u) : _ary(n) {
        for(SizeType i=0; i!=n; ++i) { _ary[i]=i; } }
    PivotMatrix(Array<SizeType> pv) : _ary(pv) { }
    SizeType size() const { return _ary.size(); }
    SizeType const& operator[](SizeType i) const { return _ary[i]; }
    SizeType& operator[](SizeType i) { return _ary[i]; }
    template<class X> operator Matrix<X> () const;
    friend OutputStream& operator<<(OutputStream& os, const PivotMatrix& pv);

    template<class X> friend Vector<X> operator*(PivotMatrix, Vector<X>);
    template<class X> friend Matrix<X> operator*(PivotMatrix, Matrix<X>);
    template<class X> friend Matrix<X> operator*(Matrix<X>, PivotMatrix);
};

} // namespace Ariadne

#endif
