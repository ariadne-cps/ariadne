/***************************************************************************
 *            system_submodule.cpp
 *
 *  Copyright 2009--17  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include "boost_python.hpp"
#include "utilities.hpp"

#include <iostream>
#include <iomanip>
#include <functional>

#include "utility/tribool.hpp"
#include "numeric/numeric.hpp"
#include "function/function.hpp"
#include "symbolic/expression.hpp"
#include "symbolic/space.hpp"
#include "hybrid/discrete_event.hpp"
#include "hybrid/hybrid_automata.hpp"
#include "hybrid/hybrid_time.hpp"
#include "hybrid/hybrid_set.hpp"

using namespace boost::python;
using namespace Ariadne;


namespace Ariadne {

template<class T> Nat __hash__(const T&);
template<> Nat __hash__<DiscreteEvent>(const DiscreteEvent& e) {
    return reinterpret_cast<const ushort&>(e.name().c_str()[0]); }
template<> Nat __hash__<DiscreteLocation>(const DiscreteLocation& q) {
    return reinterpret_cast<const ushort&>(to_string(q).c_str()[0]); }

template<>
struct from_python< Set<DiscreteEvent> > {
    from_python() { converter::registry::push_back(&convertible,&construct,type_id< Set<DiscreteEvent> >()); }
    static Void* convertible(PyObject* obj_ptr) { if (!PyList_Check(obj_ptr)) { return 0; } return obj_ptr; }
    static Void construct(PyObject* obj_ptr,converter::rvalue_from_python_stage1_data* data) {
        Void* storage = ((converter::rvalue_from_python_storage< Set<DiscreteEvent> >*)data)->storage.bytes;
        boost::python::list elements=boost::python::extract<boost::python::list>(obj_ptr);
        Set<DiscreteEvent>* evnts_ptr = new (storage) Set<DiscreteEvent>();
        for(Int i=0; i!=len(elements); ++i) {
            boost::python::extract<String> xs(elements[i]);
            if(xs.check()) { evnts_ptr->insert(DiscreteEvent(xs())); }
            else { DiscreteEvent e=boost::python::extract< DiscreteEvent >(elements[i]); evnts_ptr->insert(e); }
        }
        data->convertible = storage;
    }
};

}

Void export_hybrid_automaton()
{
    // Don't use return_value_policy<copy_const_reference> since reference lifetime should not exceed automaton lifetime

    to_python< Set<DiscreteEvent> >();
    to_python< Set<DiscreteLocation> >();

    from_python< List<DiscreteEvent> >();

    class_<DiscreteLocation> discrete_state_class("DiscreteLocation",init<DiscreteLocation>());
    discrete_state_class.def("__eq__", &__eq__<Bool,DiscreteLocation,DiscreteLocation>);
    discrete_state_class.def("__ne__", &__ne__<Bool,DiscreteLocation,DiscreteLocation>);
    discrete_state_class.def("__hash__", &__hash__<DiscreteLocation>);
    discrete_state_class.def(self_ns::str(self));

    class_<DiscreteEvent> discrete_event_class("DiscreteEvent",init<DiscreteEvent>());
    discrete_event_class.def("__eq__", &__eq__<Bool,DiscreteEvent,DiscreteEvent>);
    discrete_event_class.def("__ne__", &__ne__<Bool,DiscreteEvent,DiscreteEvent>);
    discrete_event_class.def("__hash__", &__hash__<DiscreteEvent>);
    discrete_event_class.def(self_ns::str(self));
    implicitly_convertible<Int,DiscreteEvent>();
    implicitly_convertible<StringType,DiscreteEvent>();


    class_<HybridTime> hybrid_time_class("HybridTime",init<double,Int>());
    hybrid_time_class.def("continuous_time",&HybridTime::continuous_time,return_value_policy<copy_const_reference>());
    hybrid_time_class.def("discrete_time",&HybridTime::discrete_time,return_value_policy<copy_const_reference>());

    class_<HybridAutomaton> hybrid_automaton_class("HybridAutomaton",init<>());
    hybrid_automaton_class.def("locations", &HybridAutomaton::locations);
    hybrid_automaton_class.def("events", &HybridAutomaton::events);
    hybrid_automaton_class.def("event_kind", &HybridAutomaton::event_kind);
    hybrid_automaton_class.def("dynamic_function", &HybridAutomaton::dynamic_function);
    hybrid_automaton_class.def("guard_function", &HybridAutomaton::guard_function);
    hybrid_automaton_class.def("reset_function", &HybridAutomaton::reset_function);
    hybrid_automaton_class.def("new_mode",(Void(HybridAutomaton::*)(DiscreteLocation,List<DottedRealAssignment> const&)) &HybridAutomaton::new_mode, return_value_policy<reference_existing_object>());
    hybrid_automaton_class.def("new_invariant",(Void(HybridAutomaton::*)(DiscreteLocation,ContinuousPredicate const&,DiscreteEvent)) &HybridAutomaton::new_invariant);
    hybrid_automaton_class.def("new_transition",(Void(HybridAutomaton::*)(DiscreteLocation,DiscreteEvent,DiscreteLocation,List<PrimedRealAssignment> const&,ContinuousPredicate const&,EventKind)) &HybridAutomaton::new_transition);
    hybrid_automaton_class.def(self_ns::str(self));

}



Void system_submodule() {
    export_hybrid_automaton();
}

