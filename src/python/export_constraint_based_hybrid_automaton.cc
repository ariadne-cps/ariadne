/***************************************************************************
 *            python/export_constraint_based_hybrid_automaton.cc
 *
 *  Copyright  2007 Pieter Collins
 *  Pieter.Collins@cwi.nl
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

#include "geometry/constraint_interface.h"
#include "system/map_interface.h"
#include "system/vector_field_interface.h"

#include "geometry/hybrid_set.h"
#include "system/constraint_based_hybrid_automaton.h"

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
void export_constraint_based_hybrid_automaton() 
{
  init<const std::string&> constraint_based_hybrid_automaton_init;
  
  class_< ConstraintBasedDiscreteMode<R> >("HybridAutomatonDiscreteMode",no_init)
    .def("id",&ConstraintBasedDiscreteMode<R>::id)
    .def(self_ns::str(self))
  ;
  

  class_< ConstraintBasedDiscreteTransition<R> >("HybridAutomatonDiscreteTransition",no_init)
    .def("id",&ConstraintBasedDiscreteTransition<R>::id)
    .def(self_ns::str(self))
  ;


  class_< ConstraintBasedHybridAutomaton<R> >("ConstraintBasedHybridAutomaton",init<const std::string&>())
    .def("new_dynamic",
         (const VectorFieldInterface<R>&(ConstraintBasedHybridAutomaton<R>::*)(id_type, const VectorFieldInterface<R>&)) &ConstraintBasedHybridAutomaton<R>::new_dynamic,
         return_reference_existing_object())
    .def("new_reset",
         (const MapInterface<R>&(ConstraintBasedHybridAutomaton<R>::*)(id_type, const MapInterface<R>&)) &ConstraintBasedHybridAutomaton<R>::new_dynamic,
         return_reference_existing_object())
    .def("new_constraint",
         (const ConstraintInterface<R>&(ConstraintBasedHybridAutomaton<R>::*)(id_type, const ConstraintInterface<R>&))&ConstraintBasedHybridAutomaton<R>::new_dynamic,
         return_reference_existing_object())

    .def("new_mode",
         (const ConstraintBasedDiscreteMode<R>&(ConstraintBasedHybridAutomaton<R>::*)(id_type,id_type)) &ConstraintBasedHybridAutomaton<R>::new_mode,
         return_reference_existing_object())
    .def("new_invariant",
         (const ConstraintBasedDiscreteTransition<R>&(ConstraintBasedHybridAutomaton<R>::*)(id_type,id_type,id_type)) &ConstraintBasedHybridAutomaton<R>::new_invariant,
          return_reference_existing_object())
     .def("new_transition",
          (const ConstraintBasedDiscreteTransition<R>&(ConstraintBasedHybridAutomaton<R>::*)(id_type,id_type,id_type,id_type,id_type,bool)) &ConstraintBasedHybridAutomaton<R>::new_transition,
          return_reference_existing_object())
    .def("new_forced_transition",
         (const ConstraintBasedDiscreteTransition<R>&(ConstraintBasedHybridAutomaton<R>::*)(id_type,id_type,id_type,id_type,id_type)) &ConstraintBasedHybridAutomaton<R>::new_forced_transition,
         return_reference_existing_object())
    .def("new_unforced_transition",
         (const ConstraintBasedDiscreteTransition<R>&(ConstraintBasedHybridAutomaton<R>::*)(id_type,id_type,id_type,id_type,id_type)) &ConstraintBasedHybridAutomaton<R>::new_unforced_transition,
          return_reference_existing_object())

    .def("new_mode",(const ConstraintBasedDiscreteMode<R>&(ConstraintBasedHybridAutomaton<R>::*)(id_type,const VectorFieldInterface<R>&))
           (&ConstraintBasedHybridAutomaton<R>::new_mode),
         return_reference_existing_object())
    .def("new_invariant",(const ConstraintBasedDiscreteTransition<R>&(ConstraintBasedHybridAutomaton<R>::*)(id_type,id_type,const ConstraintInterface<R>&))
           (&ConstraintBasedHybridAutomaton<R>::new_invariant),
         return_reference_existing_object())
     .def("new_transition",
         (const ConstraintBasedDiscreteTransition<R>&(ConstraintBasedHybridAutomaton<R>::*)(id_type,id_type,id_type,const MapInterface<R>&,const ConstraintInterface<R>&,bool))
           (&ConstraintBasedHybridAutomaton<R>::new_transition),
         return_reference_existing_object())
    .def("new_forced_transition",
         (const ConstraintBasedDiscreteTransition<R>&(ConstraintBasedHybridAutomaton<R>::*)(id_type,id_type,id_type,const MapInterface<R>&,const ConstraintInterface<R>&))
           (&ConstraintBasedHybridAutomaton<R>::new_forced_transition),
         return_reference_existing_object())
    .def("new_unforced_transition",
         (const ConstraintBasedDiscreteTransition<R>&(ConstraintBasedHybridAutomaton<R>::*)(id_type,id_type,id_type,const MapInterface<R>&,const ConstraintInterface<R>&))
           (&ConstraintBasedHybridAutomaton<R>::new_unforced_transition),
         return_reference_existing_object())

    .def("name",&ConstraintBasedHybridAutomaton<R>::name,return_copy_const_reference())
    .def("has_mode",&ConstraintBasedHybridAutomaton<R>::has_mode)
    .def("has_transition",&ConstraintBasedHybridAutomaton<R>::has_transition)
    .def("mode",&ConstraintBasedHybridAutomaton<R>::mode,return_copy_const_reference())
    .def("transition",&ConstraintBasedHybridAutomaton<R>::transition,return_copy_const_reference())
    .def(self_ns::str(self))
  ;
}

template void export_constraint_based_hybrid_automaton<Float>();
