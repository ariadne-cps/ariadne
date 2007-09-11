/***************************************************************************
 *            python/export_vector_field_evolver.cc
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

#include "python/python_float.h"

#include "geometry/rectangle.h"
#include "geometry/parallelotope.h"
#include "geometry/zonotope.h"
#include "geometry/grid_set.h"
#include "geometry/list_set.h"

#include "system/vector_field.h"
#include "system/affine_vector_field.h"

#include "evaluation/evolution_parameters.h"
#include "evaluation/vector_field_evolver.h"
#include "evaluation/integrator_interface.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Geometry;
using namespace Ariadne::System;
using namespace Ariadne::Evaluation;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;

template<class R>
void export_vector_field_evolver() 
{

  class_< VectorFieldEvolver<R> >("VectorFieldEvolver",init<const EvolutionParameters<R>&,const IntegratorInterface<R>&>())
    .def("integrate",(ListSet< Rectangle<R> >(VectorFieldEvolver<R>::*)(const VectorFieldInterface<R>&,const ListSet< Rectangle<R> >&,const time_type&)const)
                                    (&VectorFieldEvolver<R>::integrate))
    .def("reach",(ListSet< Rectangle<R> >(VectorFieldEvolver<R>::*)(const VectorFieldInterface<R>&,const ListSet< Rectangle<R> >&,const time_type&)const)
                                    (&VectorFieldEvolver<R>::reach))
    .def("integrate",(GridMaskSet<R>(VectorFieldEvolver<R>::*)(const VectorFieldInterface<R>&,const GridMaskSet<R>&,const GridMaskSet<R>&,const time_type&)const)
                                    (&VectorFieldEvolver<R>::integrate))
    .def("reach",(GridMaskSet<R>(VectorFieldEvolver<R>::*)(const VectorFieldInterface<R>&,const GridMaskSet<R>&,const GridMaskSet<R>&,const time_type&)const)
                              (&VectorFieldEvolver<R>::reach))
    .def("chainreach",(GridMaskSet<R>(VectorFieldEvolver<R>::*)(const VectorFieldInterface<R>&,const GridMaskSet<R>&,const GridMaskSet<R>&)const)
         (&VectorFieldEvolver<R>::chainreach))
    .def("viable",(GridMaskSet<R>(VectorFieldEvolver<R>::*)(const VectorFieldInterface<R>&,const GridMaskSet<R>&)const)
         (&VectorFieldEvolver<R>::viable))
    .def("verify",(tribool(VectorFieldEvolver<R>::*)(const VectorFieldInterface<R>&,const GridMaskSet<R>&,const GridMaskSet<R>&)const)
         (&VectorFieldEvolver<R>::verify))
  ;

}

template void export_vector_field_evolver<Float>();
