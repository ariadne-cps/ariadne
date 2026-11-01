#include "pybind11.hpp"
#include "function-utilities.hpp"
#include "reasoning/reachability_analyser.hpp"
#include "reasoning/storage.hpp"
#include "reasoning/hybrid_storage.hpp"
#include "reasoning/hybrid_reachability_analyser.hpp"
#include "evolution/vector_field_evolver.hpp"

using namespace Ariadne;

Void export_storage(pybind11::module& module) {
    pybind11::class_<Storage,pybind11::bases<Drawable2dInterface>> storage_class(module,"Storage");
    pybind11::class_<LabelledStorage,pybind11::bases<LabelledDrawable2dInterface>> labelled_storage_class(module,"LabelledStorage");
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



void reasoning_submodule(pybind11::module& module) {
    export_storage(module);
    pybind11::class_<HybridStorage>(module,"HybridStorage");
    pybind11::class_<HybridReachabilityAnalyser>(module,"HybridReachabilityAnalyser");
    export_safety_certificate<ContinuousReachabilityAnalyser>(module,"SafetyCertificate");
    export_reachability_analyser<ContinuousReachabilityAnalyser,VectorFieldEvolver>(module,"ContinuousReachabilityAnalyser");
}
