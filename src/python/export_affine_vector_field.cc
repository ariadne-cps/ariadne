/***************************************************************************
 *            python/export_affine_vector_field.cc
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

#include "real_typedef.h"

#include "numeric/interval.h"
#include "linear_algebra/matrix.h"
#include "geometry/point.h"
#include "geometry/rectangle.h"
#include "geometry/parallelotope.h"
#include "geometry/simplex.h"
#include "geometry/zonotope.h"
#include "geometry/polytope.h"
#include "system/affine_vector_field.h"


#include "python/python_utilities.h"
using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Geometry;
using namespace Ariadne::System;

#include <boost/python.hpp>
using namespace boost::python;

template<typename R>
void export_affine_vector_field() 
{
  typedef typename numerical_traits<R>::arithmetic_type F;
  typedef Interval<R> I;
  
  class_< AffineVectorField<R>, bases< VectorField<R> > >("AffineVectorField",init< Matrix<R>, Vector<R> >())
    .def("dimension", &AffineVectorField<R>::dimension)
    .def("__call__", (Vector<F>(AffineVectorField<R>::*)(const Point<R>&)const)(&AffineVectorField<R>::operator()))
    .def("__call__", (Vector<I>(AffineVectorField<R>::*)(const Rectangle<R>&)const)(&AffineVectorField<R>::operator()))
    .def("jacobian", (Matrix<F>(AffineVectorField<R>::*)(const Point<R>&)const)(&AffineVectorField<R>::jacobian))
    .def("jacobian", (Matrix<I>(AffineVectorField<R>::*)(const Point<R>&)const)(&AffineVectorField<R>::jacobian))
    .def(self_ns::str(self))
  ;
}
