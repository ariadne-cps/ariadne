/***************************************************************************
 *            python/export_set_based_hybrid_automaton.cc
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
#include "system/set_based_hybrid_automaton.h"

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
void export_set_based_hybrid_automaton() 
{
  init<const std::string&> hybrid_automaton_init;
  
  class_< SetBasedDiscreteMode<R> >("SetBasedDiscreteMode",no_init)
    .def("discrete_state",&SetBasedDiscreteMode<R>::discrete_state)
    .def("dimension",&SetBasedDiscreteMode<R>::dimension)
    .def("dynamic",&SetBasedDiscreteMode<R>::dynamic,return_reference_existing_object())
    .def("invariant",&SetBasedDiscreteMode<R>::invariant,return_reference_existing_object())
    .def(self_ns::str(self))
  ;
  

  class_< SetBasedDiscreteTransition<R> >("SetBasedDiscreteTransition",no_init)
    .def("discrete_event",&SetBasedDiscreteTransition<R>::discrete_event)
    .def("source",&SetBasedDiscreteTransition<R>::source,return_reference_existing_object())
    .def("destination",&SetBasedDiscreteTransition<R>::destination,return_reference_existing_object())
    .def("activation",&SetBasedDiscreteTransition<R>::activation,return_reference_existing_object())
    .def("reset",&SetBasedDiscreteTransition<R>::reset,return_reference_existing_object())
    .def(self_ns::str(self))
  ;


  class_< SetBasedHybridAutomaton<R> >("SetBasedHybridAutomaton",hybrid_automaton_init)
    .def("new_mode",(const SetBasedDiscreteMode<R>&(SetBasedHybridAutomaton<R>::*)(DiscreteState, const VectorField<R>&,const Geometry::ConstraintSet<R>&))
           (&SetBasedHybridAutomaton<R>::new_mode),
         return_reference_existing_object())
    .def("new_transition",
         (const SetBasedDiscreteTransition<R>&(SetBasedHybridAutomaton<R>::*)
             (DiscreteEvent,const SetBasedDiscreteMode<R>&,const SetBasedDiscreteMode<R>&,const Map<R>&,const Geometry::ConstraintSet<R>&))
           (&SetBasedHybridAutomaton<R>::new_transition),
         return_reference_existing_object())
    .def("new_transition",
         (const SetBasedDiscreteTransition<R>&(SetBasedHybridAutomaton<R>::*)
             (DiscreteEvent,DiscreteState,DiscreteState,const Map<R>&,const Geometry::ConstraintSet<R>&))
           (&SetBasedHybridAutomaton<R>::new_transition),
         return_reference_existing_object())
    .def("name",&SetBasedHybridAutomaton<R>::name,return_copy_const_reference())
    .def("has_mode",&SetBasedHybridAutomaton<R>::has_mode)
    .def("has_transition",&SetBasedHybridAutomaton<R>::has_transition)
    .def("mode",&SetBasedHybridAutomaton<R>::mode,return_copy_const_reference())
    .def("transition",&SetBasedHybridAutomaton<R>::transition,return_copy_const_reference())
    .def("locations",&SetBasedHybridAutomaton<R>::locations)
    .def("invariant",&SetBasedHybridAutomaton<R>::invariant)
    .def(self_ns::str(self))
  ;
}

template void export_set_based_hybrid_automaton<FloatPy>();
