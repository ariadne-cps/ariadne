/***************************************************************************
 *            python/export_apply.cc
 *
 *  Copyright  2005-6  Alberto Casagrande, Pieter Collins
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
#include "geometry/list_set.h"
#include "geometry/grid_set.h"

#include "system/map.h"

#include "evaluation/applicator.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::Geometry;
using namespace Ariadne::System;
using namespace Ariadne::Evaluation;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;

template<class R>
void export_apply() 
{
  typedef Interval<R> I;

  class_< Applicator<R> >("Applicator",init<const R&, const R&>())
    .def(init<>())
    .def("image",(SetInterface<R>*(Applicator<R>::*)(const MapInterface<R>&,const SetInterface<R>&)const)
                   (&Applicator<R>::image),return_value_policy<manage_new_object>(),"Compute the image of a set under a map" )
    .def("preimage",(SetInterface<R>*(Applicator<R>::*)(const MapInterface<R>&,const SetInterface<R>&,const SetInterface<R>&)const)
                   (&Applicator<R>::preimage),return_value_policy<manage_new_object>(),"Compute the preimage of a set under a map" )
    .def("reach",(SetInterface<R>*(Applicator<R>::*)(const MapInterface<R>&,const SetInterface<R>&)const)
                   (&Applicator<R>::reach),return_value_policy<manage_new_object>(),"Compute the reachable set under a map" )
    .def("chainreach",(SetInterface<R>*(Applicator<R>::*)(const MapInterface<R>&,const SetInterface<R>&,const SetInterface<R>&)const)
                        (&Applicator<R>::chainreach),return_value_policy<manage_new_object>(), "Compute the chain reachable set")
    .def("viable",(SetInterface<R>*(Applicator<R>::*)(const MapInterface<R>&,const SetInterface<R>&)const)
                        (&Applicator<R>::viable),return_value_policy<manage_new_object>(), "Compute the viability kernel")
    .def("verify",(tribool(Applicator<R>::*)(const MapInterface<R>&,const SetInterface<R>&,const SetInterface<R>&)const)
                    (&Applicator<R>::verify), "Verify that the reachable set lies within a safe set")

    .def("evaluate",(Rectangle<R>(Applicator<R>::*)(const MapInterface<R>&,const Rectangle<R>&)const)
                   (&Applicator<R>::evaluate),"Compute the image of a rectangle under a map")
    .def("evaluate",(Zonotope<R>(Applicator<R>::*)(const MapInterface<R>&,const Zonotope<R>&)const)
                   (&Applicator<R>::evaluate),"Compute the image of a zonotope under a map" )
    .def("evaluate",(Zonotope<Interval<R>,R>(Applicator<R>::*)(const MapInterface<R>&,const Zonotope<Interval<R>,R>&)const)
                   (&Applicator<R>::evaluate),"Compute the image of a zonotope under a map" )
    .def("evaluate",(Zonotope<Interval<R> >(Applicator<R>::*)(const MapInterface<R>&,const Zonotope< Interval<R> >&)const)
                   (&Applicator<R>::evaluate),"Compute the image of an interval zonotope under a map" )
    .def("image",(ListSet< Zonotope<R> >(Applicator<R>::*)(const MapInterface<R>&,const ListSet< Zonotope<R> >&)const)
                   (&Applicator<R>::image),"Compute the image of a list of zonotopes under a map" )
    .def("image",(ListSet< Zonotope<Interval<R>,R> >(Applicator<R>::*)(const MapInterface<R>&,const ListSet< Zonotope<I,R> >&)const)
                        (&Applicator<R>::image), "Compute the image of a list of zonotopes under a map")
    .def("image",(GridMaskSet<R>(Applicator<R>::*)(const MapInterface<R>&,const GridMaskSet<R>&,const GridMaskSet<R>&)const)
                   (&Applicator<R>::image),"Compute the image of a set under a map" )
    .def("preimage",(GridMaskSet<R>(Applicator<R>::*)(const MapInterface<R>&,const GridMaskSet<R>&,const GridMaskSet<R>&)const)
                   (&Applicator<R>::image),"Compute the preimage of a set under a map" )
    .def("reach",(ListSet< Zonotope<I,R> >(Applicator<R>::*)(const MapInterface<R>&,const ListSet< Zonotope<I,R> >&)const)
                        (&Applicator<R>::reach), "Compute the reachable set")
    .def("chainreach",(GridMaskSet<R>(Applicator<R>::*)(const MapInterface<R>&,const GridMaskSet<R>&,const GridMaskSet<R>&)const)
                        (&Applicator<R>::chainreach), "Compute the chain reachable set")
    .def("viable",(GridMaskSet<R>(Applicator<R>::*)(const MapInterface<R>&,const GridMaskSet<R>&)const)
                        (&Applicator<R>::viable), "Compute the viability kernel")
    .def("verify",(tribool(Applicator<R>::*)(const MapInterface<R>&,const GridMaskSet<R>&,const GridMaskSet<R>&)const)
                    (&Applicator<R>::verify), "Verify that the reachable set lies within a safe set")
  ;
 
}

template void export_apply<Float>();
