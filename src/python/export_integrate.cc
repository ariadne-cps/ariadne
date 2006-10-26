/***************************************************************************
 *            python/export_integrator.cc
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

#include "geometry/rectangle.h"
#include "geometry/parallelotope.h"
#include "geometry/zonotope.h"
#include "geometry/grid_set.h"
#include "geometry/list_set.h"

#include "system/vector_field.h"

#include "evaluation/integrator.h"
#include "evaluation/lohner_integrator.h"


using namespace Ariadne;
using namespace Ariadne::Geometry;
using namespace Ariadne::System;
using namespace Ariadne::Evaluation;

#include <boost/python.hpp>
using namespace boost::python;

template<class R> inline 
Parallelotope<R> integrate_parallelotope(
  const C1LohnerIntegrator<R> li, const VectorField<R>& vf, const Parallelotope<R>& p, const time_type& t)
{
  return li.integrate(vf,p,t);
}

template<class R> inline 
Zonotope<R> integrate_zonotope(
  const C1LohnerIntegrator<R> li, const VectorField<R>& vf, const Zonotope<R>& z, const time_type& t)
{
  return li.integrate(vf,z,t);
}

template<class R>
void export_integrate() 
{
  class_< C1LohnerIntegrator<R> >("C1LohnerIntegrator",init<R,R,R>())
    .def(init<double,double,double>()) 
    .def_readwrite("maximum_step_size", &C1LohnerIntegrator<R>::maximum_step_size)
    .def_readwrite("maximum_basic_set_radius", &C1LohnerIntegrator<R>::maximum_basic_set_radius)
    .def("integration_step",(Rectangle<R>(C0LohnerIntegrator<R>::*)(const VectorField<R>&,const Rectangle<R>&,time_type&)const)
                              (&C0LohnerIntegrator<R>::integration_step))
    .def("integration_step",(Parallelotope<R>(C1LohnerIntegrator<R>::*)(const VectorField<R>&,const Parallelotope<R>&,time_type&)const)
                              (&C1LohnerIntegrator<R>::integration_step))
    .def("reach_step",(Rectangle<R>(C0LohnerIntegrator<R>::*)(const VectorField<R>&,const Rectangle<R>&,time_type&)const)
                              (&C0LohnerIntegrator<R>::reachability_step))
    .def("reach_step",(Zonotope<R>(C1LohnerIntegrator<R>::*)(const VectorField<R>&,const Zonotope<R>&,time_type&)const)
                              (&C1LohnerIntegrator<R>::reachability_step))

    .def("integrate",(Rectangle<R>(C0LohnerIntegrator<R>::*)(const VectorField<R>&,const Rectangle<R>&,const time_type&)const)
                              (&C0LohnerIntegrator<R>::integrate))
    .def("integrate",(Parallelotope<R>(C1LohnerIntegrator<R>::*)(const VectorField<R>&,const Parallelotope<R>&,const time_type&)const)
                              (&C1LohnerIntegrator<R>::integrate))
    .def("integrate",(Zonotope<R>(C1LohnerIntegrator<R>::*)(const VectorField<R>&,const Zonotope<R>&,const time_type&)const)
                              (&C1LohnerIntegrator<R>::integrate))
    .def("integrate",(ListSet<R,Parallelotope>(C1LohnerIntegrator<R>::*)(const VectorField<R>&,const ListSet<R,Parallelotope>&,const time_type&)const)
                              (&C1LohnerIntegrator<R>::integrate))
    .def("integrate",(ListSet<R,Zonotope>(C1LohnerIntegrator<R>::*)(const VectorField<R>&,const ListSet<R,Zonotope>&,const time_type&)const)
                              (&C1LohnerIntegrator<R>::integrate))
    .def("integrate",(GridMaskSet<R>(C1LohnerIntegrator<R>::*)(const VectorField<R>&,const GridMaskSet<R>&,const GridMaskSet<R>&,const time_type&)const)
                              (&C1LohnerIntegrator<R>::integrate))

    .def("reach",(ListSet<R,Parallelotope>(C1LohnerIntegrator<R>::*)(const VectorField<R>&,const Parallelotope<R>&,const time_type&)const)
                              (&C1LohnerIntegrator<R>::reach))

    .def("reach",(ListSet<R,Zonotope>(C1LohnerIntegrator<R>::*)(const VectorField<R>&,const Zonotope<R>&,const time_type&)const)
                              (&C1LohnerIntegrator<R>::reach))

    .def("reach",(ListSet<R,Parallelotope>(C1Integrator<R>::*)(const VectorField<R>&,const ListSet<R,Parallelotope>&,const time_type&)const)
                              (&C1Integrator<R>::reach))
    .def("reach",(ListSet<R,Zonotope>(C1Integrator<R>::*)(const VectorField<R>&,const ListSet<R,Zonotope>&,const time_type&)const)
                              (&C1Integrator<R>::reach))
    .def("reach",(GridMaskSet<R>(C1LohnerIntegrator<R>::*)(const VectorField<R>&,const GridMaskSet<R>&,const GridMaskSet<R>&,const time_type&)const)
                              (&C1LohnerIntegrator<R>::reach))

    .def("chainreach",(GridMaskSet<R>(Integrator<R>::*)(const VectorField<R>&,const GridMaskSet<R>&,const GridMaskSet<R>&)const)
                              (&Integrator<R>::reach),"chain reach of a set")
    .def("verify",(bool(Integrator<R>::*)(const VectorField<R>&,const GridMaskSet<R>&,const GridMaskSet<R>&)const)
                              (&Integrator<R>::reach),"verify a safety property")

;
}

template void export_integrate<Real>();
