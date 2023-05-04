/***************************************************************************
 *            dynamics/vector_field_evolver.hpp
 *
 *  Copyright  2007-20  Alberto Casagrande, Pieter Collins
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

/*! \file dynamics/vector_field_evolver.hpp
 *  \brief Evolver for vector_field systems.
 */

#ifndef ARIADNE_VECTOR_FIELD_EVOLVER_HPP
#define ARIADNE_VECTOR_FIELD_EVOLVER_HPP

#include <string>
#include <vector>
#include <list>
#include <iostream>


#include "utility/tuple.hpp"

#include "dynamics/vector_field.hpp"
#include "function/function_interface.hpp"
#include "solvers/configuration_interface.hpp"
#include "solvers/integrator.hpp"
#include "dynamics/evolver_interface.hpp"

#include "betterthreads/workload.hpp"
#include "pronest/configuration_interface.hpp"
#include "pronest/searchable_configuration.hpp"
#include "pronest/configurable.hpp"
#include "pronest/configuration_property.hpp"
#include "pexplore/task_runner_interface.hpp"
#include "pexplore/task.tpl.hpp"
#include "conclog/logging.hpp"

using namespace ConcLog;

namespace Ariadne {
class VectorFieldEvolver;
}

namespace ProNest {
template<> struct Configuration<Ariadne::VectorFieldEvolver>;
}

namespace Ariadne {

class VectorField;
template<class ES> class Orbit;

using Mutex = std::mutex;
template<class T> using LockGuard = std::lock_guard<T>;
using BetterThreads::DynamicWorkload;
using ProNest::Configuration;
using pExplore::TaskRunnable;
class VectorFieldEvolverConfiguration;

class RealExpressionBoundedConstraintSet;

//! \brief A class for computing the evolution of a vector_field system.
//!
//! The actual evolution steps are performed by the Integrator class.
class VectorFieldEvolver
    : public EvolverInterface<VectorField,LabelledEnclosure,typename VectorField::TimeType>,
      public TaskRunnable<VectorFieldEvolver>
{
  public:
    typedef EvolverInterface<VectorField,LabelledEnclosure,typename VectorField::TimeType> Interface;
    typedef Configuration<VectorFieldEvolver> ConfigurationType;
    typedef VectorField SystemType;
    typedef IntegratorInterface IntegratorType;
    typedef typename VectorField::TimeType TimeType;
    typedef Dyadic TimeStepType;
    typedef TimeType TerminationType;
    typedef LabelledEnclosure EnclosureType;
    typedef ListSet<EnclosureType> EnclosureListType;
    typedef Pair<TimeStepType,EnclosureType> TimedEnclosureType;
    typedef Orbit<EnclosureType> OrbitType;
    typedef ValidatedFunctionPatchFactoryInterface FunctionFactoryType;
  private:
    //! \brief Synchronised wrapping of orbit to allow concurrent adjoining
    struct SynchronisedOrbit : public OrbitType {
        SynchronisedOrbit(const EnclosureType& initial) : OrbitType(initial) { }
        Void adjoin_reach(const EnclosureType& set) override { LockGuard<Mutex> lock(_mux); OrbitType::adjoin_reach(set); }
        Void adjoin_intermediate(const EnclosureType& set) override { LockGuard<Mutex> lock(_mux); OrbitType::adjoin_intermediate(set); }
        Void adjoin_final(const EnclosureType& set) override { LockGuard<Mutex> lock(_mux); OrbitType::adjoin_final(set); }
        SizeType reach_size() { LockGuard<Mutex> lock(_mux); return OrbitType::reach().size(); }
      private:
        Mutex _mux;
    };
    typedef DynamicWorkload<TimedEnclosureType,TimeType const&,Semantics,SharedPointer<SynchronisedOrbit>> WorkloadType;
  public:

    //! \brief Construct from \a system, \a configuration and \a integrator.
    VectorFieldEvolver(SystemType const& system, ConfigurationType const& configuration);

    //! \brief Make a dynamically-allocated copy.
    VectorFieldEvolver* clone() const;

    //! \brief Get the internal system.
    virtual const SystemType& system() const { return *_system; }

    //! \brief Make an enclosure from a computed box set.
    EnclosureType enclosure(ExactBoxType const&) const;

    //!@}

    //!@{
    //! \name Evolution using abstract sets.
    //! \brief Compute an approximation to the orbit set using upper semantics.
    OrbitType orbit(EnclosureType const& initial_set, TimeType const& time, Semantics semantics) const;
    OrbitType orbit(RealVariablesBox const& initial_set, TimeType const& time, Semantics semantics) const;
    OrbitType orbit(RealExpressionBoundedConstraintSet const& initial_set, TimeType const& time, Semantics semantics) const;

    //!@}

  protected:

    Void _process_timed_enclosure(WorkloadType::Access& workload, TimedEnclosureType const& current_timed_set,
                                  TimeType const& maximum_time, Semantics semantics, SharedPointer<SynchronisedOrbit> result) const;

    Void _process_timed_enclosure_step(WorkloadType::Access& workload, TimedEnclosureType const& current_timed_set,
                                       TimeType const& maximum_time, Semantics semantics, SharedPointer<SynchronisedOrbit> result) const;

    Void _append_initial_set(WorkloadType& workload, TimeStepType const& initial_time, EnclosureType const& current_set) const;

  private:
    SharedPointer<SystemType> _system;
};

} // namespace Ariadne

namespace ProNest {

using Ariadne::VectorFieldEvolver;
using Ariadne::ExactDouble;
using Ariadne::ApproximateDouble;
using Ariadne::Bool;
using Ariadne::IntegratorInterface;
using Ariadne::TaylorPicardIntegrator;
using ProNest::RangeConfigurationProperty;


//! \brief Configuration for a VectorFieldEvolver, essentially for controlling the accuracy of continuous evolution methods.
template<> struct Configuration<VectorFieldEvolver> final : public SearchableConfiguration {
    typedef Configuration<VectorFieldEvolver> C;
    typedef double RealType;
    typedef RangeConfigurationProperty<RealType> RealTypeProperty;
    typedef InterfaceListConfigurationProperty<IntegratorInterface> IntegratorProperty;

    Configuration() {
        add_property("enable_premature_termination",BooleanConfigurationProperty(false));
        add_property("enable_reconditioning",BooleanConfigurationProperty(true));
        add_property("enable_subdivisions",BooleanConfigurationProperty(false));
        add_property("enable_clobbering",BooleanConfigurationProperty(false));
        add_property("integrator", InterfaceListConfigurationProperty<IntegratorInterface>(TaylorPicardIntegrator(Configuration<TaylorPicardIntegrator>())));
        add_property("maximum_enclosure_radius",RealTypeProperty(std::numeric_limits<double>::max(),Log10SearchSpaceConverter<RealType>()));
        add_property("maximum_spacial_error",RealTypeProperty(std::numeric_limits<double>::max(),Log10SearchSpaceConverter<RealType>()));
        add_property("maximum_step_size",RealTypeProperty(std::numeric_limits<double>::max(),Log2SearchSpaceConverter<RealType>()));
    }

    //! \brief Enable premature termination of lower evolution
    Bool const& enable_premature_termination() const { return at<BooleanConfigurationProperty>("enable_premature_termination").get(); }
    C& set_enable_premature_termination(Bool const& value) { at<BooleanConfigurationProperty>("enable_premature_termination").set(value); return *this; }

    //! \brief Enable reconditioning of basic sets
    Bool const& enable_reconditioning() const { return at<BooleanConfigurationProperty>("enable_reconditioning").get(); }
    C& set_enable_reconditioning(Bool const& value) { at<BooleanConfigurationProperty>("enable_reconditioning").set(value); return *this; }
    C& set_both_enable_reconditioning() { at<BooleanConfigurationProperty>("enable_reconditioning").set_both(); return *this; }

    //! \brief Enable subdivisions for upper evolution
    Bool const& enable_subdivisions() const { return at<BooleanConfigurationProperty>("enable_subdivisions").get(); }
    C& set_enable_subdivisions(Bool const& value) { at<BooleanConfigurationProperty>("enable_subdivisions").set(value); return *this; }

    //! \brief Enable clobbering (i.e., removing model error) to evolve state function
    //! \details This makes the orbit result approximate
    Bool const& enable_clobbering() const { return at<BooleanConfigurationProperty>("enable_clobbering").get(); }
    C& set_enable_clobbering(Bool const& value) { at<BooleanConfigurationProperty>("enable_clobbering").set(value); return *this; }

    //! \brief The maximum allowable step size for integration.
    //! Decreasing this value increases the accuracy of the computation.
    RealType const& maximum_step_size() const { return at<RealTypeProperty>("maximum_step_size").get(); }
    C& set_maximum_step_size(RealType const& value) { at<RealTypeProperty>("maximum_step_size").set(value); return *this; }
    C& set_maximum_step_size(RealType const& lower, RealType const& upper) { at<RealTypeProperty>("maximum_step_size").set(lower,upper); return *this; }

    //! \brief The maximum allowable approximation error for reconditioning
    RealType const& maximum_spacial_error() const { return at<RealTypeProperty>("maximum_spacial_error").get(); }
    C& set_maximum_spacial_error(RealType const& value) { at<RealTypeProperty>("maximum_spacial_error").set(value); return *this; }
    C& set_maximum_spacial_error(RealType const& lower, RealType const& upper) { at<RealTypeProperty>("maximum_spacial_error").set(lower,upper); return *this; }

    //! \brief The maximum allowable radius of a basic set during integration.
    //! Decreasing this value increases the accuracy of the computation of an over-approximation.
    RealType const& maximum_enclosure_radius() const { return at<RealTypeProperty>("maximum_enclosure_radius").get(); }
    C& set_maximum_enclosure_radius(RealType const& value) { at<RealTypeProperty>("maximum_enclosure_radius").set(value); return *this; }

    //! \brief The integrator to be used.
    IntegratorInterface const& integrator() const { return at<IntegratorProperty>("integrator").get(); }
    C& set_integrator(IntegratorInterface const& integrator) { at<IntegratorProperty>("integrator").set(integrator); return *this; }
    C& set_integrator(SharedPointer<IntegratorInterface> const& integrator) { at<IntegratorProperty>("integrator").set(integrator); return *this; }
};

} // namespace ProNest

namespace pExplore {

using Ariadne::VectorFieldEvolver;
using Ariadne::LabelledEnclosure;
using Ariadne::IntegratorInterface;
using Ariadne::EffectiveVectorMultivariateFunction;
using Ariadne::Dyadic;

template<> struct TaskInput<VectorFieldEvolver> {
    TaskInput(EffectiveVectorMultivariateFunction const& dynamic_, LabelledEnclosure const& current_set_, Dyadic const& current_time_) :
        dynamic(dynamic_), current_set(current_set_), current_time(current_time_) { }
    EffectiveVectorMultivariateFunction const& dynamic;
    LabelledEnclosure const& current_set;
    Dyadic const& current_time;
    shared_ptr<IntegratorInterface> integrator_ptr;
};

template<> struct TaskOutput<VectorFieldEvolver> {
    TaskOutput(LabelledEnclosure const& evolve_, LabelledEnclosure const& reach_, Dyadic const& time_) :
        evolve(evolve_), reach(reach_), time(time_) { }
    LabelledEnclosure const evolve;
    LabelledEnclosure const reach;
    Dyadic const time;
};

template<> struct Task<VectorFieldEvolver> final: ParameterSearchTaskBase<VectorFieldEvolver> {
    TaskOutput<VectorFieldEvolver> run(TaskInput<VectorFieldEvolver> const& in, Configuration<VectorFieldEvolver> const& cfg) const override;
};

} // namespace pExplore

#endif // ARIADNE_VECTOR_FIELD_EVOLVER_HPP
