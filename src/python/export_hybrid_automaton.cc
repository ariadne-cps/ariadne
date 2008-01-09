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

#include "python/float.h"

#include "geometry/set_reference.h"
#include "geometry/hybrid_space.h"
#include "geometry/hybrid_set.h"
#include "system/hybrid_automaton.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::Geometry;
using namespace Ariadne::System;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;
typedef return_value_policy<copy_const_reference> return_copy_const_reference;
typedef return_value_policy<reference_existing_object> return_reference_existing_object;


template<class R>
void export_hybrid_automaton() 
{
  init<const std::string&> hybrid_automaton_init;
  
  class_< DiscreteMode<R> >("DiscreteMode",no_init)
    .def("discrete_state",&DiscreteMode<R>::discrete_state)
    .def("dimension",&DiscreteMode<R>::dimension)
    .def("dynamic",&DiscreteMode<R>::dynamic,return_reference_existing_object())
    .def("invariant",&DiscreteMode<R>::invariant,return_reference_existing_object())
    .def(self_ns::str(self))
  ;
  

  class_< DiscreteTransition<R> >("DiscreteTransition",no_init)
    .def("discrete_event",&DiscreteTransition<R>::discrete_event)
    .def("source",&DiscreteTransition<R>::source,return_reference_existing_object())
    .def("destination",&DiscreteTransition<R>::destination,return_reference_existing_object())
    .def("activation",&DiscreteTransition<R>::activation,return_reference_existing_object())
    .def("reset",&DiscreteTransition<R>::reset,return_reference_existing_object())
    .def(self_ns::str(self))
  ;


  class_< HybridAutomaton<R> >("HybridAutomaton",hybrid_automaton_init)
    .def("new_mode",(const DiscreteMode<R>&(HybridAutomaton<R>::*)(DiscreteState, const VectorField<R>&,const Geometry::ConstraintSet<R>&))
           (&HybridAutomaton<R>::new_mode),
         return_reference_existing_object())
    .def("new_transition",
         (const DiscreteTransition<R>&(HybridAutomaton<R>::*)
             (DiscreteEvent,const DiscreteMode<R>&,const DiscreteMode<R>&,const Map<R>&,const Geometry::ConstraintSet<R>&))
           (&HybridAutomaton<R>::new_transition),
         return_reference_existing_object())
    .def("new_transition",
         (const DiscreteTransition<R>&(HybridAutomaton<R>::*)
             (DiscreteEvent,DiscreteState,DiscreteState,const Map<R>&,const Geometry::ConstraintSet<R>&))
           (&HybridAutomaton<R>::new_transition),
         return_reference_existing_object())
    .def("name",&HybridAutomaton<R>::name,return_copy_const_reference())
    .def("has_mode",&HybridAutomaton<R>::has_mode)
    .def("has_transition",&HybridAutomaton<R>::has_transition)
    .def("mode",&HybridAutomaton<R>::mode,return_copy_const_reference())
    .def("transition",&HybridAutomaton<R>::transition,return_copy_const_reference())
    .def("locations",&HybridAutomaton<R>::locations)
    .def("invariant",&HybridAutomaton<R>::invariant)
    .def(self_ns::str(self))
  ;
}

template void export_hybrid_automaton<FloatPy>();
