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


#include "../utility/tuple.hpp"

#include "../dynamics/vector_field.hpp"
#include "../function/function_interface.hpp"
#include "../solvers/configuration.hpp"
#include "../solvers/integrator_interface.hpp"
#include "../solvers/integrator.hpp"
#include "../dynamics/evolver_base.hpp"
#include "../concurrency/task_runner.tpl.hpp"
#include "../concurrency/task_appraisal_space.hpp"

#include "../output/logging.hpp"

namespace Ariadne {

class VectorField;
template<class ES> class Orbit;

class VectorFieldEvolver;
template<> class Configuration<VectorFieldEvolver>;
class VectorFieldFlowStepTask;

//! \brief A class for computing the evolution of a vector_field system.
//!
//! The actual evolution steps are performed by the Integrator class.
class VectorFieldEvolver
    : public EvolverBase<VectorField,LabelledEnclosure,typename VectorField::TimeType>,
      public Configurable<VectorFieldEvolver>,
      public TaskRunnable<VectorFieldEvolver>
{
  public:
    typedef VectorField SystemType;
    typedef typename VectorField::TimeType TimeType;
    typedef Dyadic TimeStepType;
    typedef TimeType TerminationType;
    typedef LabelledEnclosure EnclosureType;
    typedef Pair<TimeStepType, EnclosureType> TimedEnclosureType;
    typedef Orbit<EnclosureType> OrbitType;
    typedef ListSet<EnclosureType> EnclosureListType;
    typedef ValidatedFunctionModelDPFactory::Interface FunctionFactoryType;
  public:

    //! \brief Construct from parameters and an integrator to compute the flow.
    VectorFieldEvolver(
    		const SystemType& system,
            const IntegratorInterface& integrator);

    //! \brief Make a dynamically-allocated copy.
    VectorFieldEvolver* clone() const override { return new VectorFieldEvolver(*this); }

    //! \brief Get the internal system.
    virtual const SystemType& system() const override { return *_sys_ptr; }

    //! \brief Make an enclosure from a user set.
    EnclosureType enclosure(RealBox const&) const;
    EnclosureType enclosure(RealBox const&, EnclosureConfiguration const&) const;

    //! \brief Make an enclosure from a user set with variables.
    EnclosureType enclosure(RealVariablesBox const&) const;
    EnclosureType enclosure(RealVariablesBox const&, EnclosureConfiguration const&) const;

    //! \brief Make an enclosure from a computed box set.
    EnclosureType enclosure(ExactBoxType const&) const;
    EnclosureType enclosure(ExactBoxType const&, EnclosureConfiguration const&) const;

    //! \brief The class which constructs functions for the enclosures.
    const FunctionFactoryType& function_factory() const;

    //!@}

    //!@{
    //! \name Evolution using abstract sets.
    //! \brief Compute an approximation to the orbit set using upper semantics.
    Orbit<EnclosureType> orbit(const EnclosureType& initial_set, const TimeType& time, Semantics semantics=Semantics::UPPER) const override;

    using EvolverBase< VectorField, EnclosureType, TerminationType >::evolve;
    using EvolverBase< VectorField, EnclosureType, TerminationType >::reach;

    //! \brief Compute an approximation to the evolution set using upper semantics.
    EnclosureListType evolve(const EnclosureType& initial_set, const TimeType& time) const {
        EnclosureListType final; EnclosureListType reachable; EnclosureListType intermediate;
        this->_evolution(final,reachable,intermediate,initial_set,time,Semantics::UPPER,false);
        return final; }

    //! \brief Compute an approximation to the reachable set under upper semantics.
    EnclosureListType reach(const EnclosureType& initial_set, const TimeType& time) const {
        EnclosureListType final; EnclosureListType reachable; EnclosureListType intermediate;
        this->_evolution(final,reachable,intermediate,initial_set,time,Semantics::UPPER,true);
        return reachable; }
    //!@}

protected:

    virtual Void _evolution(EnclosureListType& final, EnclosureListType& reachable, EnclosureListType& intermediate,
                            const EnclosureType& initial, const TimeType& time,
                            Semantics semantics, Bool reach) const override;

    virtual Void _evolution_step(List< TimedEnclosureType >& working_sets,
                                 EnclosureListType& final, EnclosureListType& reachable, EnclosureListType& intermediate,
                                 const TimedEnclosureType& current_set, StepSizeType& last_step_size, const TimeType& time,
                                 Semantics semantics, Bool reach) const;

    virtual Void _append_initial_set(List<TimedEnclosureType>& working_sets, const TimeStepType& initial_time, const EnclosureType& current_set) const;

  private:
    SharedPointer<SystemType> _sys_ptr;
};

template class RangeConfigurationProperty<ExactDouble>;
template class ListConfigurationProperty<IntegratorInterface>;

//! \brief Configuration for a VectorFieldEvolver, essentially for controlling the accuracy of continuous evolution methods.
template<> class Configuration<VectorFieldEvolver> final : public SearchableConfiguration
{
  public:
    typedef ExactDouble RealType;
    typedef ApproximateDouble ApproximateRealType;
    typedef RangeConfigurationProperty<RealType> RealTypeProperty;
    typedef ListConfigurationProperty<IntegratorInterface> IntegratorProperty;

    //! \brief Default constructor gives reasonable values.
    Configuration() {
        add_property("enable_reconditioning",BooleanConfigurationProperty(false));
        add_property("integrator",ListConfigurationProperty<IntegratorInterface>(TaylorPicardIntegrator(1e-2)));
        add_property("maximum_spacial_error",RealTypeProperty(ExactDouble::infinity(),Log10SearchSpaceConverter<RealType>()));
        add_property("maximum_enclosure_radius",RealTypeProperty(ExactDouble::infinity(),Log10SearchSpaceConverter<RealType>()));
        add_property("maximum_step_size",RealTypeProperty(ExactDouble::infinity(),Log2SearchSpaceConverter<RealType>()));
    }

    //! \brief Enable reconditioning of basic sets
    Bool const& enable_reconditioning() const { return dynamic_cast<BooleanConfigurationProperty const&>(*properties().get("enable_reconditioning")).get(); }
    void set_enable_reconditioning() { dynamic_cast<BooleanConfigurationProperty&>(*properties().get("enable_reconditioning")).set(); }
    void set_enable_reconditioning(Bool const& value) { dynamic_cast<BooleanConfigurationProperty&>(*properties().get("enable_reconditioning")).set(value); }

    //! \brief The maximum allowable step size for integration.
    //! Decreasing this value increases the accuracy of the computation.
    RealType const& maximum_step_size() const { return dynamic_cast<RealTypeProperty const&>(*properties().get("maximum_step_size")).get(); }
    void set_maximum_step_size(ApproximateRealType const& value) { dynamic_cast<RealTypeProperty&>(*properties().get("maximum_step_size")).set(cast_exact(value)); }
    void set_maximum_step_size(ApproximateRealType const& lower, ApproximateRealType const& upper) { dynamic_cast<RealTypeProperty&>(*properties().get("maximum_step_size")).set(cast_exact(lower),cast_exact(upper)); }

    //! \brief The maximum allowable approximation error for reconditioning
    RealType const& maximum_spacial_error() const { return dynamic_cast<RealTypeProperty const&>(*properties().get("maximum_spacial_error")).get(); }
    void set_maximum_spacial_error(ApproximateRealType const& value) { dynamic_cast<RealTypeProperty&>(*properties().get("maximum_spacial_error")).set(cast_exact(value)); }
    void set_maximum_spacial_error(ApproximateRealType const& lower, ApproximateRealType const& upper) { dynamic_cast<RealTypeProperty&>(*properties().get("maximum_spacial_error")).set(cast_exact(lower),cast_exact(upper)); }

    //! \brief The maximum allowable radius of a basic set during integration.
    //! Decreasing this value increases the accuracy of the computation of an over-approximation.
    RealType const& maximum_enclosure_radius() const { return dynamic_cast<RealTypeProperty const&>(*properties().get("maximum_enclosure_radius")).get(); }
    void set_maximum_enclosure_radius(ApproximateRealType const& value) { dynamic_cast<RealTypeProperty&>(*properties().get("maximum_enclosure_radius")).set(cast_exact(value)); }
    void set_maximum_enclosure_radius(ApproximateRealType const& lower, ApproximateRealType const& upper) { dynamic_cast<RealTypeProperty&>(*properties().get("maximum_enclosure_radius")).set(cast_exact(lower),cast_exact(upper)); }

    //! \brief The integrator to be used.
    IntegratorInterface const& integrator() const { return dynamic_cast<IntegratorProperty const&>(*properties().get("integrator")).get(); }
    void set_integrator(IntegratorInterface const& integrator) { dynamic_cast<IntegratorProperty&>(*properties().get("integrator")).set(integrator); }
    void set_integrator(SharedPointer<IntegratorInterface> const& integrator) { dynamic_cast<IntegratorProperty&>(*properties().get("integrator")).set(integrator); }
};

template<> struct TaskInput<VectorFieldEvolver> {
    TaskInput(EffectiveVectorMultivariateFunction const& dynamic_, LabelledEnclosure const& current_set_,
                          FloatDPExactBox const& current_set_bounds_, Dyadic const& current_time_, Dyadic const& previous_step_size_) :
            dynamic(dynamic_), current_set(current_set_), current_set_bounds(current_set_bounds_),
            current_time(current_time_), previous_step_size(previous_step_size_) { }
    EffectiveVectorMultivariateFunction const& dynamic;
    LabelledEnclosure const& current_set;
    FloatDPExactBox const& current_set_bounds;
    Dyadic const& current_time;
    Dyadic const& previous_step_size;
};

template<> struct TaskOutput<VectorFieldEvolver> {
    TaskOutput(LabelledEnclosure const& evolve_, LabelledEnclosure const& reach_, Dyadic const& time_, Dyadic const& step_size_used_) :
            evolve(evolve_), reach(reach_), time(time_), step_size_used(step_size_used_) { }
    LabelledEnclosure const evolve;
    LabelledEnclosure const reach;
    Dyadic const time;
    Dyadic const step_size_used;
};

template<> class Task<VectorFieldEvolver> final: public ParameterSearchTaskBase<VectorFieldEvolver> {
    typedef VectorFieldEvolver R;
    typedef TaskInput<R> I;
    typedef TaskOutput<R> O;
    typedef Configuration<R> C;
public:
    Task() : ParameterSearchTaskBase<R>(
            "stp",
            TaskSearchSpace({
                                    MetricSearchParameter("starting_step_size_num_refinements", 2, 5,2),
                                    MetricSearchParameter("sweep_threshold", exp(-RealVariable("sweep_threshold") * log(RealConstant(10))), 8, 12,9),
                                    MetricSearchParameter("maximum_temporal_order", 9, 15,12)
                            }),
            TaskAppraisalSpaceBuilder<R>()
                    .add(execution_time_appraisal_parameter<R>)
                    .add(ScalarAppraisalParameter<R>("step_size_used",TaskAppraisalParameterOptimisation::MAXIMISE,[](I const& i,O const& o,DurationType const& d) { return o.step_size_used.get_d(); }),2)
                    /*.add(VectorAppraisalParameter<R>("final_set_width_increases",TaskAppraisalParameterOptimisation::MINIMISE,
                                                       [](I const& i,O const& o,DurationType const& d,SizeType const& idx) { return (o.evolve.euclidean_set().bounding_box()[idx].width() - i.current_set_bounds[idx].width()).get_d(); },
                                                       [](I const& i){ return i.current_set_bounds.dimension(); }))*/
                    .build()) { }

    C singleton_configuration(C const& cfg, TaskSearchPoint const& p) const override {
        auto result = cfg;
        TaylorPicardIntegrator const& old_integrator = static_cast<TaylorPicardIntegrator const&>(cfg.integrator());
        result.set_integrator(TaylorPicardIntegrator(
                MaximumError(old_integrator.maximum_error()),
                ThresholdSweeper<FloatDP>(DoublePrecision(),p.value("sweep_threshold")),
                LipschitzConstant(old_integrator.lipschitz_tolerance()),
                StartingStepSizeNumRefinements(p.value("starting_step_size_num_refinements").get_d()),
                StepMaximumError(old_integrator.step_maximum_error()),
                MinimumTemporalOrder(old_integrator.minimum_temporal_order()),
                MaximumTemporalOrder(p.value("maximum_temporal_order").get_d())
        ));
        return result;
    }

    O run_task(I const& in, C const& cfg) const override {
        LabelledEnclosure next_set = in.current_set;
        LabelledEnclosure reach_set = in.current_set;
        Dyadic next_time = in.current_time;
        Dyadic chosen_step_size = cfg.maximum_step_size();
        FlowStepModelType flow_model = cfg.integrator().flow_step(in.dynamic, in.current_set_bounds, in.previous_step_size,chosen_step_size);
        ARIADNE_LOG_PRINTLN_VAR_AT(1, chosen_step_size);
        ARIADNE_LOG_PRINTLN_VAR_AT(1, flow_model);
        next_time += chosen_step_size;
        ARIADNE_LOG_PRINTLN_VAR_AT(1, next_time);
        reach_set.apply_full_reach_step(flow_model);
        ARIADNE_LOG_PRINTLN_VAR_AT(1, reach_set);
        next_set.apply_fixed_evolve_step(flow_model, chosen_step_size);
        ARIADNE_LOG_PRINTLN_VAR_AT(1, next_set);
        return O(next_set, reach_set, next_time, chosen_step_size);
    }
};

} // namespace Ariadne

#endif // ARIADNE_VECTOR_FIELD_EVOLVER_HPP
