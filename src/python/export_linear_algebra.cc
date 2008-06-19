/***************************************************************************
 *            python/export_linear_algebra.cc
 *
 *  Copyright  2007 Pieter Collins
 *   Pieter.Collins@cwi.nl
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

#include <utility>  //for std::pair

#include "numeric/rational.h"
#include "numeric/interval.h"
#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "linear_algebra/matrix_function.h"
#include "linear_algebra/qr_matrix.h"

#include "python/utilities.h"
#include "python/float.h"
#include "python/read_scalar.h"

using namespace Ariadne;

using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Python;

#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>
using namespace boost::python;

template<class R>
boost::python::tuple
qr(const Matrix<R>& A)
{
  const Matrix<double>& dA=reinterpret_cast<const Matrix<double>&>(A);
  const QRMatrix<double> dqr(dA);
  Matrix<double> dq=dqr.Q();
  Matrix<double> dr=dqr.R();
  const Matrix<R>& q=reinterpret_cast<const Matrix<R>&>(dq);
  const Matrix<R>& r=reinterpret_cast<const Matrix<R>&>(dr);

  return boost::python::make_tuple(q,r);
}

template<class R>
void export_linear_algebra()
{
  def("qr",qr<R>);
}

template void export_linear_algebra<FloatPy>();
