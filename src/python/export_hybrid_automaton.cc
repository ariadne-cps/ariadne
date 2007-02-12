/***************************************************************************
 *            python/export_hybrid_automaton.cc
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

#include "system/discrete_mode.h"
#include "system/discrete_transition.h"
#include "system/hybrid_automaton.h"

using namespace Ariadne;
using namespace Ariadne::Geometry;
using namespace Ariadne::System;

#include <boost/python.hpp>
using namespace boost::python;
typedef return_value_policy<copy_const_reference> return_copy_const_reference;
typedef return_value_policy<reference_existing_object> return_reference_existing_object;


template<class R>
void export_hybrid_automaton() 
{
  init<const std::string&> hybrid_automaton_init;
  
  class_< DiscreteMode<R> >("DiscreteMode",no_init)
    .def("id",&DiscreteMode<R>::id)
    .def("name",&DiscreteMode<R>::name,return_copy_const_reference())
    .def("dimension",&DiscreteMode<R>::dimension)
    .def("dynamic",&DiscreteMode<R>::dynamic,return_reference_existing_object())
    .def("invariant",&DiscreteMode<R>::invariant,return_reference_existing_object())
    .def(self_ns::str(self))
  ;
  

  class_< DiscreteTransition<R> >("DiscreteTransition",no_init)
    .def("source",&DiscreteTransition<R>::source,return_reference_existing_object())
    .def("destination",&DiscreteTransition<R>::destination,return_reference_existing_object())
    .def("activation",&DiscreteTransition<R>::activation,return_reference_existing_object())
    .def("reset",&DiscreteTransition<R>::reset,return_reference_existing_object())
    .def(self_ns::str(self))
  ;


  class_< HybridAutomaton<R> >("HybridAutomaton",hybrid_automaton_init)
    .def("new_mode",(DiscreteMode<R>&(HybridAutomaton<R>::*)(const VectorField<R>&,const Geometry::Set<R>&))
           (&HybridAutomaton<R>::new_mode),
         return_reference_existing_object())
    .def("new_transition",
         (DiscreteTransition<R>&(HybridAutomaton<R>::*)
             (const Map<R>&,const Geometry::Set<R>&,const DiscreteMode<R>&,const DiscreteMode<R>&))
           (&HybridAutomaton<R>::new_transition),
         return_reference_existing_object())
    .def("new_transition",
         (DiscreteTransition<R>&(HybridAutomaton<R>::*)(const Map<R>&,const Geometry::Set<R>&,const id_type&,const id_type&))
           (&HybridAutomaton<R>::new_transition),
         return_reference_existing_object())
    .def("name",&HybridAutomaton<R>::name,return_copy_const_reference())
    .def("mode",&HybridAutomaton<R>::mode,return_copy_const_reference())
    .def("transition",&HybridAutomaton<R>::transition,return_copy_const_reference())
    .def(self_ns::str(self))
  ;
}

template void export_hybrid_automaton<Real>();
