/***************************************************************************
 *            matrix.cc
 *
 *  Copyright 2008  Alberto Casagrande, Pieter Collins
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
 
#include <numeric.h>
#include "vector.h"
#include "matrix.h"

template class boost::numeric::ublas::matrix<Ariadne::Float>;
template class boost::numeric::ublas::matrix<Ariadne::Interval>;

namespace Ariadne {

Matrix<Float> 
midpoint(const Matrix<Interval>& A) {
    Matrix<Float> R(A.row_size(),A.column_size());
    for(size_t i=0; i!=A.row_size(); ++i) {
        for(size_t j=0; j!=A.row_size(); ++j) {
            R[i][j]=A[i][j].midpoint();
        }
    }
    return R;
}


template class Matrix<Float>;
template class Matrix<Interval>;

} // namespace Ariadne

