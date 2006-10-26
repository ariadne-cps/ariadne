/***************************************************************************
 *            python/export_apply.cc
 *
 *  6 February 2006
 *  Copyright  2005  Alberto Casagrande, Pieter Collins
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
#include "geometry/list_set.h"
#include "geometry/grid_set.h"

#include "system/map.h"

#include "evaluation/apply.h"

using namespace Ariadne;
using namespace Ariadne::Geometry;
using namespace Ariadne::System;
using namespace Ariadne::Evaluation;

#include <boost/python.hpp>
using namespace boost::python;

template<class R>
void export_apply() 
{
  class_< C1Applicator<R> >("Applicator",init<>())
    .def("apply",(Rectangle<R>(C0Applicator<R>::*)(const Map<R>&,const Rectangle<R>&)const)
                   (&C0Applicator<R>::apply),"apply the image of a map to a set")
    .def("apply",(Parallelotope<R>(C1Applicator<R>::*)(const Map<R>&,const Parallelotope<R>&)const)
                   (&C1Applicator<R>::apply),"apply the image of a map to a set" )
    .def("apply",(ListSet<R,Parallelotope>(C1Applicator<R>::*)(const Map<R>&,const ListSet<R,Parallelotope>&)const)
                   (&C1Applicator<R>::apply),"apply the image of a map to a set" )
    .def("apply",(GridMaskSet<R>(C1Applicator<R>::*)(const Map<R>&,const GridMaskSet<R>&,const GridMaskSet<R>&)const)
                   (&C1Applicator<R>::apply),"apply the image of a map to a set" )
    .def("chainreach",(GridMaskSet<R>(C1Applicator<R>::*)(const Map<R>&,const GridMaskSet<R>&,const GridMaskSet<R>&)const)
                        (&C1Applicator<R>::chainreach), "Compute the chain reachable set")
    .def("verify",(bool(C1Applicator<R>::*)(const Map<R>&,const GridMaskSet<R>&,const GridMaskSet<R>&)const)
                    (&C1Applicator<R>::verify), "Verify that the reachable set lies within a safe set")
  ;
 
  def("apply", (Rectangle<R>(*)(const Map<R>&,const Rectangle<R>&))(&apply), "apply the image of a map to a set" );
  def("apply", (Parallelotope<R>(*)(const Map<R>&,const Parallelotope<R>&))(&apply), "apply the image of a map to a set" );
  def("apply", (ListSet<R,Parallelotope>(*)(const Map<R>&,const ListSet<R,Parallelotope>&))(&apply), "apply the image of a map to a set" );
  def("apply", (GridMaskSet<R>(*)(const Map<R>&,const GridMaskSet<R>&,const GridMaskSet<R>&))(&apply), "apply the image of a map to a set" );
  def("chainreach", (GridMaskSet<R>(*)(const Map<R>&,const GridMaskSet<R>&,const GridMaskSet<R>&))(&chainreach), "Compute the chain reachable set" );
  def("verify", (bool(*)(const Map<R>&,const GridMaskSet<R>&,const GridMaskSet<R>&))(&verify), "Compute the chain reachable set" );
  
}

template void export_apply<Real>();
