/***************************************************************************
 *            evolution_submodule.cpp
 *
 *  Copyright  2008-20  Pieter Collins
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
#include "function-utilities.hpp"

#include "evolution/orbit.hpp"
#include "evolution/vector_field_simulator.hpp"
#include "evolution/iterated_map_evolver.hpp"
#include "evolution/vector_field_evolver.hpp"
#include "evolution/hybrid_time.hpp"
#include "evolution/hybrid_enclosure.hpp"
#include "evolution/hybrid_simulator.hpp"
#include "system/hybrid_automaton_interface.hpp"


using namespace Ariadne;

Void export_labelled_drawable_2d_interface(pybind11::module& module) {
    pybind11::class_<LabelledDrawable2dInterface> labelled_drawable_interface_class(module,"LabelledDrawable2dInterface");
}

Void export_semantics(pybind11::module& module) {
    pybind11::enum_<Semantics> semantics_enum(module,"Semantics");
    semantics_enum.value("UPPER", Semantics::UPPER);
    semantics_enum.value("LOWER", Semantics::LOWER);
    //    semantics_enum.def("__repr__" , &__cstr__<Semantics>);
}

Void export_enclosure(pybind11::module& module) {
    pybind11::class_<Enclosure,pybind11::bases<Drawable2dInterface>> enclosure_class(module,"Enclosure");
    enclosure_class.def("bounding_box",&Enclosure::bounding_box);
    enclosure_class.def("state_set",&Enclosure::state_set);
    enclosure_class.def("state_auxiliary_set",&Enclosure::state_auxiliary_set);
    enclosure_class.def("state_time_auxiliary_set",&Enclosure::state_time_auxiliary_set);
    enclosure_class.def("__str__", &__cstr__<Enclosure>);
    pybind11::class_<LabelledEnclosure,pybind11::bases<LabelledDrawable2dInterface,Enclosure>> labelled_enclosure_class(module,"LabelledEnclosure");
    labelled_enclosure_class.def("bounding_box",&LabelledEnclosure::bounding_box);
    labelled_enclosure_class.def("__str__", &__cstr__<LabelledEnclosure>);
}

template<class T> Void export_list_set(pybind11::module& module, const char* name) {
    pybind11::class_<ListSet<T>> list_set_class(module,name);
    list_set_class.def("size",&ListSet<T>::size);
    list_set_class.def("__len__",&ListSet<T>::size);
    list_set_class.def("__getitem__", &__getitem__<ListSet<T>,SizeType>);
    list_set_class.def("__iter__", [](ListSet<T> const& l){return pybind11::make_iterator(l);});
    list_set_class.def("__str__", &__cstr__<ListSet<T>>);
}


template<class ORB>
Void export_orbit(pybind11::module& module, const char* name)
{
    pybind11::class_<ORB,pybind11::bases<LabelledDrawable2dInterface>> orbit_class(module,name);
    orbit_class.def("reach", &ORB::reach);
    orbit_class.def("evolve", &ORB::final);
    orbit_class.def("final", &ORB::final);
    orbit_class.def("__str__", &__cstr__<ORB>);
}

template<class SIM> Void export_simulator(pybind11::module& module, const char* name);

template<> Void export_simulator<VectorFieldSimulator>(pybind11::module& module, const char* name)
{
    typedef VectorFieldSimulator::TerminationType TerminationType;
    typedef VectorFieldSimulator::ApproximateListPointType ApproximateListPointType;
    typedef VectorFieldSimulator::OrbitType OrbitType;
    typedef VectorFieldSimulator::RealBoxType RealBoxType;

    auto const& reference_internal = pybind11::return_value_policy::reference_internal;

    pybind11::class_<LabelledInterpolatedCurve,pybind11::bases<LabelledDrawable2dInterface>> labelled_interpolated_curve_class(module,"LabelledInterpolatedCurve");

    pybind11::class_<OrbitType,pybind11::bases<LabelledDrawable2dInterface>> simulator_orbit_class(module,"ApproximatePointOrbit");
    simulator_orbit_class.def("curves", &OrbitType::curves);

    pybind11::class_<VectorFieldSimulator> simulator_class(module,name);
    simulator_class.def(pybind11::init<VectorFieldSimulator::SystemType const&>());
    simulator_class.def("configuration",pybind11::overload_cast<>(&VectorFieldSimulator::configuration),reference_internal);
    simulator_class.def("orbit", (OrbitType(VectorFieldSimulator::*)(const ApproximateListPointType&, const TerminationType&)const) &VectorFieldSimulator::orbit);
    simulator_class.def("orbit", pybind11::overload_cast<ApproximateListPointType const&,TerminationType const&>(&VectorFieldSimulator::orbit,pybind11::const_));
    simulator_class.def("orbit", pybind11::overload_cast<RealExpressionBoundedConstraintSet const&,TerminationType const&>(&VectorFieldSimulator::orbit,pybind11::const_));
    simulator_class.def("orbit", pybind11::overload_cast<RealBoxType const&,TerminationType const&>(&VectorFieldSimulator::orbit,pybind11::const_));

    typedef typename VectorFieldSimulator::ConfigurationType ConfigurationType;
    pybind11::class_<ConfigurationType> simulator_configuration_class(module,"VectorFieldSimulatorConfiguration");
    simulator_configuration_class.def("set_step_size", &ConfigurationType::set_step_size);
    simulator_configuration_class.def("set_discretisation_type", &ConfigurationType::set_discretisation_type);
    simulator_configuration_class.def("set_num_subdivisions", &ConfigurationType::set_num_subdivisions);
    simulator_configuration_class.def("__repr__",&__cstr__<ConfigurationType>);
}

template<class EV>
Void export_evolver_interface(pybind11::module& module, const char* name)
{
    pybind11::class_<EV> evolver_interface_class(module,name);
}

template<class EV, class... Params>
Void export_evolver(pybind11::module& module, const char* name)
{
    typedef EV Evolver;
    typedef typename EV::Interface Interface;
    typedef typename EV::EnclosureType EnclosureType;
    typedef typename EV::TerminationType TerminationType;
    typedef typename EV::ConfigurationType ConfigurationType;
    typedef typename EV::OrbitType OrbitType;

    auto const& reference_internal = pybind11::return_value_policy::reference_internal;

    pybind11::class_<Evolver,pybind11::bases<Interface>> evolver_class(module,name);
    evolver_class.def(pybind11::init<Params...>());
    evolver_class.def("orbit",(OrbitType(Evolver::*)(const EnclosureType&,const TerminationType&,Semantics)const) &Evolver::orbit);
    evolver_class.def("orbit",(OrbitType(Evolver::*)(const RealExpressionBoundedConstraintSet&,const TerminationType&,Semantics)const) &Evolver::orbit);
    evolver_class.def("configuration",(ConfigurationType&(Evolver::*)())&Evolver::configuration,reference_internal);
}

Void export_vector_field_evolver_configuration(pybind11::module& module) {

    typedef typename VectorFieldEvolver::ConfigurationType ConfigurationType;

    pybind11::class_<ConfigurationType> vector_field_evolver_configuration_class(module,"VectorFieldEvolverConfiguration");
    vector_field_evolver_configuration_class.def("set_maximum_step_size", &ConfigurationType::set_maximum_step_size);
    vector_field_evolver_configuration_class.def("set_maximum_enclosure_radius", &ConfigurationType::set_maximum_enclosure_radius);
    vector_field_evolver_configuration_class.def("set_maximum_spacial_error", &ConfigurationType::set_maximum_spacial_error);
    vector_field_evolver_configuration_class.def("set_enable_reconditioning", &ConfigurationType::set_enable_reconditioning);
    vector_field_evolver_configuration_class.def("set_enable_subdivisions", &ConfigurationType::set_enable_subdivisions);
    vector_field_evolver_configuration_class.def("__repr__",&__cstr__<ConfigurationType>);
}

Void evolution_submodule(pybind11::module& module)
{
    export_labelled_drawable_2d_interface(module);

    export_semantics(module);

    export_enclosure(module);

    export_list_set<LabelledEnclosure>(module,"LabelledEnclosureListSet");

    export_simulator<VectorFieldSimulator>(module,"VectorFieldSimulator");

    export_evolver_interface<IteratedMapEvolver::Interface>(module, "IteratedMapEvolverInterface");
    export_evolver_interface<VectorFieldEvolver::Interface>(module,"VectorFieldEvolverInterface");

    export_orbit<VectorFieldEvolver::OrbitType>(module,"LabelledEnclosureOrbit");

    export_evolver<IteratedMapEvolver, IteratedMap>(module, "IteratedMapEvolver");
    export_evolver<VectorFieldEvolver, VectorField, IntegratorInterface const&>(module,"VectorFieldEvolver");

    export_vector_field_evolver_configuration(module);

    // Hybrid evolution belongs here, not to a separate hybrid binding layer.
    pybind11::class_<HybridTime>(module,"HybridTime")
        .def(pybind11::init<Real,Integer>())
        .def("continuous_time", &HybridTime::continuous_time)
        .def("discrete_time", &HybridTime::discrete_time);
    pybind11::class_<HybridEnclosure>(module,"HybridEnclosure")
        .def("location", &HybridEnclosure::location)
        .def("state_space", &HybridEnclosure::state_space);
    pybind11::class_<HybridSimulator>(module,"HybridSimulator")
        .def(pybind11::init<HybridSimulator::SystemType const&>());

}
