/***************************************************************************
 *            system_submodule.cpp
 *
 *  Copyright  2009-20  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "pybind11.hpp"
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

using namespace Ariadne;


namespace Ariadne {

template<class T> Nat __hash__(const T&);
template<> Nat __hash__<DiscreteEvent>(const DiscreteEvent& e) {
    return reinterpret_cast<const unsigned short&>(e.name().c_str()[0]); }
template<> Nat __hash__<DiscreteLocation>(const DiscreteLocation& q) {
    return reinterpret_cast<const unsigned short&>(to_string(q).c_str()[0]); }

}

Void export_hybrid_automaton(pybind11::module& module)
{
    // Don't use return_value_policy<copy_const_reference> since reference lifetime should not exceed automaton lifetime

    pybind11::class_<DiscreteLocation> discrete_state_class(module,"DiscreteLocation");
    discrete_state_class.def(pybind11::init<DiscreteLocation>());
    discrete_state_class.def("__eq__", &__eq__<DiscreteLocation,DiscreteLocation , Return<Bool> >);
    discrete_state_class.def("__ne__", &__ne__<DiscreteLocation,DiscreteLocation , Return<Bool> >);
    discrete_state_class.def("__hash__", &__hash__<DiscreteLocation>);
    discrete_state_class.def("__str__",&__cstr__<DiscreteLocation>);

    pybind11::class_<DiscreteEvent> discrete_event_class(module,"DiscreteEvent");
    discrete_event_class.def(pybind11::init<DiscreteEvent>());
    discrete_event_class.def("__eq__", &__eq__<DiscreteEvent,DiscreteEvent , Return<Bool> >);
    discrete_event_class.def("__ne__", &__ne__<DiscreteEvent,DiscreteEvent , Return<Bool> >);
    discrete_event_class.def("__hash__", &__hash__<DiscreteEvent>);
    discrete_event_class.def("__str__", &__cstr__<DiscreteEvent>);
    pybind11::implicitly_convertible<Int,DiscreteEvent>();
    pybind11::implicitly_convertible<StringType,DiscreteEvent>();


    pybind11::class_<HybridTime> hybrid_time_class(module,"HybridTime");
    hybrid_time_class.def(pybind11::init<double,Int>());
    hybrid_time_class.def("continuous_time",&HybridTime::continuous_time);
    hybrid_time_class.def("discrete_time",&HybridTime::discrete_time);

    pybind11::class_<HybridAutomaton> hybrid_automaton_class(module,"HybridAutomaton");
    hybrid_automaton_class.def(pybind11::init<>());
    hybrid_automaton_class.def("locations", &HybridAutomaton::locations);
    hybrid_automaton_class.def("events", &HybridAutomaton::events);
    hybrid_automaton_class.def("event_kind", &HybridAutomaton::event_kind);
    hybrid_automaton_class.def("dynamic_function", &HybridAutomaton::dynamic_function);
    hybrid_automaton_class.def("guard_function", &HybridAutomaton::guard_function);
    hybrid_automaton_class.def("reset_function", &HybridAutomaton::reset_function);
    hybrid_automaton_class.def("new_mode",(Void(HybridAutomaton::*)(DiscreteLocation,List<DottedRealAssignment> const&)) &HybridAutomaton::new_mode, pybind11::return_value_policy::reference_internal);
    hybrid_automaton_class.def("new_invariant",(Void(HybridAutomaton::*)(DiscreteLocation,ContinuousPredicate const&,DiscreteEvent)) &HybridAutomaton::new_invariant);
    hybrid_automaton_class.def("new_transition",(Void(HybridAutomaton::*)(DiscreteLocation,DiscreteEvent,DiscreteLocation,List<PrimedRealAssignment> const&,ContinuousPredicate const&,EventKind)) &HybridAutomaton::new_transition);
    hybrid_automaton_class.def("__str__", &__cstr__<HybridAutomaton>);

}



Void system_submodule(pybind11::module& module) {
    export_hybrid_automaton(module);
}

