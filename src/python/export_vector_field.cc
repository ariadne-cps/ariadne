/***************************************************************************
 *            python/export_vector_field.cc
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is diself_ns::stributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include "python/float.h"

#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "differentiation/taylor_derivative.h"
#include "geometry/point.h"
#include "system/vector_field.h"
#include "system/affine_vector_field.h"

using namespace Ariadne;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;


template<class R>
void export_vector_field() 
{
  class_<VectorField<R> >("VectorField", init<const FunctionInterface<R>&>())
    .def(init<const VectorField<R>&>())
    .def("dimension", &VectorField<R>::dimension)
    .def("smoothness", &VectorField<R>::smoothness)
    .def("evaluate", &VectorField<R>::evaluate)
    .def("jacobian", &VectorField<R>::jacobian)
  ;
 
  class_< AffineVectorField<R>, bases< VectorField<R> > >("AffineVectorField", init<const Matrix<R>&, const Vector<R>&>())
    .def(self_ns::str(self));
}

template void export_vector_field<FloatPy>();
