/***************************************************************************
 *            python/export_constraint_hybrid_automaton.cc
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

#include "geometry/set_reference.h"
#include "geometry/hybrid_set.h"
#include "system/constraint_hybrid_automaton.h"
#include "system/map.h"
#include "system/vector_field.h"

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
void export_constraint_hybrid_automaton() 
{
  init<const std::string&> hybrid_automaton_init;
  
  typedef ConstraintHybridAutomaton<R> HybridAutomaton;
  typedef typename ConstraintHybridAutomaton<R>::mode_type DiscreteMode;
  typedef typename ConstraintHybridAutomaton<R>::transition_type DiscreteTransition;

  class_< DiscreteMode >("ConstraintHybridAutomatonDiscreteMode",no_init)
    .def("id",&DiscreteMode::id)
    .def(self_ns::str(self))
  ;
  

  class_< DiscreteTransition >("ConstraintHybridAutomatonDiscreteTransition",no_init)
    .def("id",&DiscreteTransition::id)
    .def(self_ns::str(self))
  ;


  class_< HybridAutomaton >("ConstraintHybridAutomaton",init<const std::string&>())
    .def("new_mode",(const DiscreteMode&(HybridAutomaton::*)(id_type, const VectorFieldInterface<R>&,const ConstraintInterface<R>&))
           (&HybridAutomaton::new_mode),
         return_reference_existing_object())
    .def("new_transition",
         (const DiscreteTransition&(HybridAutomaton::*)(id_type,id_type,id_type,const MapInterface<R>&,const ConstraintInterface<R>&))
           (&HybridAutomaton::new_transition),
         return_reference_existing_object())
    .def("name",&HybridAutomaton::name,return_copy_const_reference())
    .def("has_mode",&HybridAutomaton::has_mode)
    .def("has_transition",&HybridAutomaton::has_transition)
    .def("mode",&HybridAutomaton::mode,return_copy_const_reference())
    .def("transition",&HybridAutomaton::transition,return_copy_const_reference())
    .def(self_ns::str(self))
  ;
}

template void export_constraint_hybrid_automaton<Float>();
