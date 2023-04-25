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
#include "utilities.hpp"

#include "dynamics/orbit.hpp"
#include "dynamics/vector_field_simulator.hpp"
#include "dynamics/iterated_map_evolver.hpp"
#include "dynamics/vector_field_evolver.hpp"
#include "dynamics/reachability_analyser.hpp"
#include "pronest/configuration_interface.hpp"
#include "pronest/configurable.hpp"
#include "pronest/configurable.tpl.hpp"
#include "pronest/configuration_property.tpl.hpp"
#include "pexplore/task_runner_interface.hpp"

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

Void export_storage(pybind11::module& module) {
    pybind11::class_<Storage,pybind11::bases<Drawable2dInterface>> storage_class(module,"Storage");
    pybind11::class_<LabelledStorage,pybind11::bases<LabelledDrawable2dInterface>> labelled_storage_class(module,"LabelledStorage");
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
    list_set_class.def("__iter__", [](ListSet<T> const& l){return pybind11::make_iterator(l.begin(),l.end());});
    list_set_class.def("size",&ListSet<T>::size);
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

template<class EV, class... Params>
Void export_runnable_evolver(pybind11::module& module, const char* name)
{
    typedef EV Evolver;
    typedef typename EV::Interface Interface;
    typedef typename EV::EnclosureType EnclosureType;
    typedef typename EV::TerminationType TerminationType;
    typedef typename EV::OrbitType OrbitType;

    char* runnable_name = new char[100];
    std::strcpy(runnable_name,name);
    strcat(runnable_name,"TaskRunnable");

    pybind11::class_<TaskRunnable<EV>> runnable_class(module,runnable_name);

    pybind11::class_<Evolver,pybind11::bases<Interface,TaskRunnable<EV>>> evolver_class(module,name);
    evolver_class.def(pybind11::init<Params...>());
    evolver_class.def("orbit",(OrbitType(Evolver::*)(const EnclosureType&,const TerminationType&,Semantics)const) &Evolver::orbit);
    evolver_class.def("orbit",(OrbitType(Evolver::*)(const RealExpressionBoundedConstraintSet&,const TerminationType&,Semantics)const) &Evolver::orbit);
}

Void export_vector_field_evolver_configuration(pybind11::module& module) {

    typedef typename VectorFieldEvolver::ConfigurationType ConfigurationType;
    auto const& reference_internal = pybind11::return_value_policy::reference_internal;
    using pybind11::overload_cast;

    pybind11::class_<ConfigurationType> vector_field_evolver_configuration_class(module,"VectorFieldEvolverConfiguration");
    vector_field_evolver_configuration_class.def("set_maximum_step_size", overload_cast<double const&>(&ConfigurationType::set_maximum_step_size),reference_internal);
    vector_field_evolver_configuration_class.def("set_maximum_enclosure_radius", overload_cast<double const&>(&ConfigurationType::set_maximum_enclosure_radius),reference_internal);
    vector_field_evolver_configuration_class.def("set_maximum_spacial_error", overload_cast<double const&>(&ConfigurationType::set_maximum_spacial_error),reference_internal);
    vector_field_evolver_configuration_class.def("set_enable_reconditioning", &ConfigurationType::set_enable_reconditioning);
    vector_field_evolver_configuration_class.def("set_enable_subdivisions", &ConfigurationType::set_enable_subdivisions);
    vector_field_evolver_configuration_class.def("__repr__",&__cstr__<ConfigurationType>);
}

template<class RA> Void export_safety_certificate(pybind11::module& module, const char* name) {
    typedef typename RA::StorageType StorageType;
    typedef typename RA::StateSpaceType StateSpaceType;
    typedef SafetyCertificate<StateSpaceType> SafetyCertificateType;

    pybind11::class_<SafetyCertificateType> safety_certificate_class(module,name);
    safety_certificate_class.def(pybind11::init<ValidatedSierpinskian,StorageType,StorageType >());
    safety_certificate_class.def_readonly("is_safe", &SafetyCertificateType::is_safe);
    safety_certificate_class.def_readonly("chain_reach_set", &SafetyCertificateType::chain_reach_set);
    safety_certificate_class.def_readonly("safe_set", &SafetyCertificateType::safe_set);
}


template<class RA, class... Params>
Void export_reachability_analyser(pybind11::module& module, const char* name)
{
    typedef typename RA::ConfigurationType Configuration;
    typedef typename RA::StorageType StorageType;
    typedef typename RA::OvertSetInterfaceType OvertSetType;
    typedef typename RA::OpenSetInterfaceType OpenSetType;
    typedef typename RA::CompactSetInterfaceType CompactSetType;
    typedef typename RA::SafetyCertificateType SafetyCertificateType;
    typedef typename RA::TimeType TimeType;

    auto const& reference_internal = pybind11::return_value_policy::reference_internal;

    pybind11::class_<RA> reachability_analyser_class(module,name);
    reachability_analyser_class.def(pybind11::init<Params...>());
    reachability_analyser_class.def("configuration",(Configuration&(RA::*)())&RA::configuration,reference_internal);
    reachability_analyser_class.def("evolver",&RA::evolver,reference_internal);
    reachability_analyser_class.def("lower_reach",(StorageType(RA::*)(OvertSetType const&,TimeType const&)const) &RA::lower_reach);
    reachability_analyser_class.def("upper_reach", (StorageType(RA::*)(CompactSetType const&,TimeType const&)const) &RA::upper_reach);
    reachability_analyser_class.def("outer_chain_reach", (StorageType(RA::*)(CompactSetType const&)const)&RA::outer_chain_reach);
    reachability_analyser_class.def("verify_safety",(SafetyCertificateType(RA::*)(CompactSetType const&, OpenSetType const&)const) &RA::verify_safety);

    pybind11::class_<Configuration> reachability_analyser_configuration_class(module,"ReachabilityAnalyserConfiguration");
    reachability_analyser_configuration_class.def("set_transient_time", &Configuration::set_transient_time);
    reachability_analyser_configuration_class.def("set_maximum_grid_fineness", &Configuration::set_maximum_grid_fineness);
    reachability_analyser_configuration_class.def("set_lock_to_grid_time", &Configuration::set_lock_to_grid_time);
    reachability_analyser_configuration_class.def("set_maximum_grid_extent", &Configuration::set_maximum_grid_extent);
    reachability_analyser_configuration_class.def("set_bounding_domain", &Configuration::set_bounding_domain);
    reachability_analyser_configuration_class.def("__repr__",&__cstr__<Configuration>);
}


Void evolution_submodule(pybind11::module& module)
{
    export_labelled_drawable_2d_interface(module);

    export_semantics(module);

    export_storage(module);
    export_enclosure(module);

    export_list_set<LabelledEnclosure>(module,"LabelledEnclosureListSet");

    export_simulator<VectorFieldSimulator>(module,"VectorFieldSimulator");

    export_evolver_interface<IteratedMapEvolver::Interface>(module, "IteratedMapEvolverInterface");
    export_evolver_interface<VectorFieldEvolver::Interface>(module,"VectorFieldEvolverInterface");

    export_orbit<VectorFieldEvolver::OrbitType>(module,"LabelledEnclosureOrbit");

    export_vector_field_evolver_configuration(module);

    export_evolver<IteratedMapEvolver, IteratedMap>(module, "IteratedMapEvolver");
    export_runnable_evolver<VectorFieldEvolver, VectorField, Configuration<VectorFieldEvolver>, IntegratorInterface const&>(module,"VectorFieldEvolver");

    export_safety_certificate<ContinuousReachabilityAnalyser>(module,"SafetyCertificate");
    export_reachability_analyser<ContinuousReachabilityAnalyser,VectorFieldEvolver>(module,"ContinuousReachabilityAnalyser");
}
