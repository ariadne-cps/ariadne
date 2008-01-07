/***************************************************************************
 *            standard_integrator.h
 *
 *  Copyright  2006-7  Alberto Casagrande, Pieter Collins
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
 

#ifndef ARIADNE_STANDARD_INTEGRATOR_H
#define ARIADNE_STANDARD_INTEGRATOR_H

#include "numeric/interval.h"
#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "geometry/point.h"
#include "geometry/box.h"
#include "system/vector_field.h"

namespace {

using namespace Ariadne;

}



namespace Ariadne { 
  namespace Evaluation { 

    template<class R> 
    std::pair< Numeric::Rational, Geometry::Box<R> >
    standard_flow_bounds(const System::VectorField<R>& vf, 
                         const Geometry::Box<R>& bx,
                         const Numeric::Rational& t);

    template<class R>
    Geometry::Point< Numeric::Interval<R> >
    standard_flow_step(const System::VectorField<R>& vector_field, 
                       const Geometry::Point< Numeric::Interval<R> >& initial_point, 
                       const Numeric::Interval<R>& step_size, 
                       const Geometry::Box<R>& bounding_box);

    template<class R>
    LinearAlgebra::Matrix< Numeric::Interval<R> >
    standard_flow_step_jacobian(const System::VectorField<R>& vector_field, 
                                const Geometry::Point< Numeric::Interval<R> >& initial_point, 
                                const Numeric::Interval<R>& step_size, 
                                const Geometry::Box<R>& bounding_box);


  }
}

#endif /* STANDARD_INTEGRATOR_H */
