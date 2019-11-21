/***************************************************************************
 *            pivot_matrix.tpl.hpp
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

#include "pivot_matrix.hpp"

namespace Ariadne {

template<class X> Matrix<X> pivot_matrix(const Array<SizeType>& pv)
{
    const SizeType n=pv.size();
    Array<SizeType> perm(n); for(Nat i=0; i!=n; ++i) { perm[i]=i; }
    for(SizeType i=0; i!=n; ++i) {
        std::swap(perm[i],perm[pv[i]]);
    }
    Matrix<X> P(n,n);
    for(SizeType i=0; i!=n; ++i) {
        P[i][perm[i]]=1;
    }
    return P;
}

template<class X> PivotMatrix::operator Matrix<X> () const {
    return pivot_matrix<X>(this->_ary);
}

template<class X> Matrix<X> operator*(PivotMatrix P, Matrix<X> A) {
    ARIADNE_PRECONDITION(P.size()==A.row_size());
    for(SizeType i=0; i!=A.row_size(); ++i) {
        for(SizeType j=0; j!=A.column_size(); ++j) {
            std::swap(A[i][j],A[P[i]][j]);
        }
    }
    return A;
}

template<class X> Matrix<X> operator*(Matrix<X> A, PivotMatrix P) {
    ARIADNE_PRECONDITION(A.column_size()==P.size());
    for(SizeType j=0; j!=A.column_size(); ++j) {
        for(SizeType i=0; i!=A.row_size(); ++i) {
            std::swap(A[i][j],A[i][P[j]]);
        }
    }
    return A;
}


} // namespace Ariadne
