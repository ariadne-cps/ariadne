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


#include "numeric/interval.h"
#include "geometry/point.h"
#include "geometry/rectangle.h"
#include "geometry/parallelotope.h"
#include "geometry/simplex.h"
#include "geometry/zonotope.h"
#include "geometry/polyhedron.h"
#include "evaluation/affine_vector_field.h"


#include "python/typedefs.h"
#include "python/python_utilities.h"
using namespace Ariadne;

#include <boost/python.hpp>
using namespace boost::python;

typedef RVector (RAffineVectorField::*AffVectorFieldApplyPoint) (const RPoint&) const;
typedef RIntervalVector (RAffineVectorField::*AffVectorFieldApplyRectangle) (const RRectangle&) const;

//typedef RParallelotope (RAffineVectorField::*AffVectorFieldApplyParallelotope) (const RParallelotope&) const;
//typedef RZonotope (RAffineVectorField::*AffVectorFieldApplyZonotope) (const RZonotope&) const;
//typedef RSimplex (RAffineVectorField::*AffVectorFieldApplySimplex) (const RSimplex&) const;
//typedef RPolyhedron (RAffineVectorField::*AffVectorFieldApplyPolyhedron) (const RPolyhedron&) const;

typedef RMatrix (RAffineVectorField::*AffVectorFieldDerivativePoint) (const RPoint&) const;
typedef RIntervalMatrix (RAffineVectorField::*AffVectorFieldDerivativeRectangle) (const RRectangle&) const;

AffVectorFieldApplyPoint affine_vf_apply_point=&RAffineVectorField::apply;
AffVectorFieldApplyRectangle affine_vf_apply_rectangle=&RAffineVectorField::apply;
//AffVectorFieldApplyParallelotope affine_vf_apply_parallelotope=&RAffineVectorField::operator();
//AffVectorFieldApplyZonotope affine_vf_apply_zonotope=&RAffineVectorField::operator();
//AffVectorFieldApplySimplex affine_vf_apply_simplex=&RAffineVectorField::operator();
//AffVectorFieldApplyPolyhedron affine_vf_apply_polyhedron=&RAffineVectorField::operator();
AffVectorFieldDerivativePoint affine_vf_derivative_point=&RAffineVectorField::derivative;
AffVectorFieldDerivativeRectangle affine_vf_derivative_rectangle=&RAffineVectorField::derivative;

void export_affine_vector_field() {

  class_< RAffineVectorField, bases<RVectorFieldBase> >("AffineVectorField",init<RMatrix,RVector>())
    .def("dimension", &RAffineVectorField::dimension)
    .def("__call__", affine_vf_apply_point)
    .def("__call__", affine_vf_apply_rectangle)
//    .def("__call__", affine_vf_apply_parallelotope)
//    .def("__call__", affine_vf_apply_zonotope)
//    .def("__call__", affine_vf_apply_polyhedron)
//    .def("__call__", affine_vf_apply_simplex)
    .def("derivative", affine_vf_derivative_point)
    .def("derivative", affine_vf_derivative_rectangle)
//    .def(self_ns::str(self))    // __self_ns::str__
  ;
}
