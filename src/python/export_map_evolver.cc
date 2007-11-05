/***************************************************************************
 *            python/export_map_evolver.cc
 *
 *  Copyright  2005-7  Alberto Casagrande, Pieter Collins
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

#include "evaluation/map_evolver.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::Geometry;
using namespace Ariadne::System;
using namespace Ariadne::Evaluation;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;

template<class R>
void export_map_evolver() 
{
  typedef Interval<R> I;
  typedef Zonotope<I,R> BS;

  class_< MapEvolver<R> > evolver_class("MapEvolver",init< EvolutionParameters<R> >());
  evolver_class.def(init<const EvolutionParameters<R>&,const ApplicatorInterface<BS>&>());
  evolver_class.def("image",(SetInterface<R>*(MapEvolver<R>::*)(const MapInterface<R>&,const SetInterface<R>&)const)
                    (&MapEvolver<R>::image),return_value_policy<manage_new_object>(),"Compute the image of a set under a map" );
  evolver_class.def("preimage",(SetInterface<R>*(MapEvolver<R>::*)(const MapInterface<R>&,const SetInterface<R>&,const SetInterface<R>&)const)
                    (&MapEvolver<R>::preimage),return_value_policy<manage_new_object>(),"Compute the preimage of a set under a map" );
  evolver_class.def("iterate",(SetInterface<R>*(MapEvolver<R>::*)(const MapInterface<R>&,const SetInterface<R>&,const Integer&)const)
                    (&MapEvolver<R>::iterate),return_value_policy<manage_new_object>(),"Compute an approximation to the iterate of a set under a map" );
  evolver_class.def("reach",(SetInterface<R>*(MapEvolver<R>::*)(const MapInterface<R>&,const SetInterface<R>&,const Integer&)const)
                    (&MapEvolver<R>::reach),return_value_policy<manage_new_object>(),"Compute an approximation to the timed reachable set under a map" );
  evolver_class.def("lower_reach",(SetInterface<R>*(MapEvolver<R>::*)(const MapInterface<R>&,const SetInterface<R>&)const)
                    (&MapEvolver<R>::lower_reach),return_value_policy<manage_new_object>(),"Compute a lower-approximation to the reachable set under a map" );
  evolver_class.def("chainreach",(SetInterface<R>*(MapEvolver<R>::*)(const MapInterface<R>&,const SetInterface<R>&,const SetInterface<R>&)const)
                    (&MapEvolver<R>::chainreach),return_value_policy<manage_new_object>(), "Compute an outer-approximation to the chain reachable set");
  evolver_class.def("viable",(SetInterface<R>*(MapEvolver<R>::*)(const MapInterface<R>&,const SetInterface<R>&)const)
                    (&MapEvolver<R>::viable),return_value_policy<manage_new_object>(), "Compute the viability kernel");
  evolver_class.def("verify",(tribool(MapEvolver<R>::*)(const MapInterface<R>&,const SetInterface<R>&,const SetInterface<R>&)const)
                    (&MapEvolver<R>::verify), "Verify that the reachable set lies within a safe set");


}

template void export_map_evolver<Float>();
