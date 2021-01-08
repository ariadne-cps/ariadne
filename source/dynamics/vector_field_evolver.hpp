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
#include "../solvers/configuration_interface.hpp"
#include "../solvers/integrator_interface.hpp"
#include "../solvers/integrator.hpp"
#include "../dynamics/evolver_base.hpp"

#include "../concurrency/task_runner_interface.hpp"

#include "../output/logging.hpp"

namespace Ariadne {

class VectorField;
template<class ES> class Orbit;

class VectorFieldEvolverConfiguration;

struct FlowStepRunnerInput {
    FlowStepRunnerInput(EffectiveVectorMultivariateFunction const& dynamic_, IntegratorInterface const& integrator_, LabelledEnclosure const& current_set_,
                        FloatDPExactBox const& current_set_bounds_, Dyadic const& current_time_, Dyadic const& previous_step_size_, Dyadic const& maximum_step_size_) :
            dynamic(dynamic_), integrator(integrator_), current_set(current_set_), current_set_bounds(current_set_bounds_),
            current_time(current_time_), previous_step_size(previous_step_size_), maximum_step_size(maximum_step_size_) { }
    EffectiveVectorMultivariateFunction const& dynamic;
    IntegratorInterface const& integrator;
    LabelledEnclosure const& current_set;
    FloatDPExactBox const& current_set_bounds;
    Dyadic const& current_time;
    Dyadic const& previous_step_size;
    Dyadic const& maximum_step_size;
};

struct FlowStepRunnerConfiguration {
    FlowStepRunnerConfiguration(SharedPointer<TaylorPicardIntegrator> const& integrator_) : integrator(integrator_){ }
    SharedPointer<TaylorPicardIntegrator> integrator;
};

struct FlowStepRunnerOutput {
    FlowStepRunnerOutput(LabelledEnclosure const& evolve_, LabelledEnclosure const& reach_, Dyadic const& time_, Dyadic const& step_size_used_) :
            evolve(evolve_), reach(reach_), time(time_), step_size_used(step_size_used_) { }
    LabelledEnclosure const evolve;
    LabelledEnclosure const reach;
    Dyadic const time;
    Dyadic const step_size_used;
};

inline TaskSearchSpace make_flow_step_runner_space() {
    RealVariable sssnr("starting_step_size_num_refinements"), st("sweep_threshold"), mto("maximum_temporal_order");
    return TaskSearchSpace({MetricSearchParameter(sssnr, 5, 2),
                            MetricSearchParameter(st, exp(-st * log(RealConstant(10))), 12, 9),
                            MetricSearchParameter(mto, 15, 12)
                             },(st*mto)/sssnr);
}

inline FlowStepRunnerConfiguration
to_configuration(FlowStepRunnerInput const& in, TaskSearchPoint const& p) {
    TaylorPicardIntegrator const& default_integrator = static_cast<TaylorPicardIntegrator const&>(in.integrator);
    SharedPointer<TaylorPicardIntegrator> integrator(new TaylorPicardIntegrator(
            MaximumError(default_integrator.maximum_error()),
            ThresholdSweeper<FloatDP>(DoublePrecision(),p.value("sweep_threshold")),
            LipschitzConstant(default_integrator.lipschitz_tolerance()),
            StartingStepSizeNumRefinements(p.value("starting_step_size_num_refinements").get_d()),
            StepMaximumError(default_integrator.step_maximum_error()),
            MinimumTemporalOrder(default_integrator.minimum_temporal_order()),
            MaximumTemporalOrder(p.value("maximum_temporal_order").get_d())
    ));
    return FlowStepRunnerConfiguration(integrator);
}

inline FlowStepRunnerOutput
run_task(FlowStepRunnerInput const& in, FlowStepRunnerConfiguration const& cfg) {
    LabelledEnclosure next_set = in.current_set;
    LabelledEnclosure reach_set = in.current_set;
    Dyadic next_time = in.current_time;
    Dyadic chosen_step_size = in.maximum_step_size;
    FlowStepModelType flow_model = cfg.integrator->flow_step(in.dynamic, in.current_set_bounds, in.previous_step_size,chosen_step_size);
    ARIADNE_LOG_PRINTLN("step_size = " << chosen_step_size);
    ARIADNE_LOG_PRINTLN_AT(1, "flow_model = " << flow_model);
    next_time += chosen_step_size;
    ARIADNE_LOG_PRINTLN_AT(1, "next_time = " << next_time)
    reach_set.apply_full_reach_step(flow_model);
    ARIADNE_LOG_PRINTLN_AT(1, "reach_set = " << reach_set);
    next_set.apply_fixed_evolve_step(flow_model, chosen_step_size);
    ARIADNE_LOG_PRINTLN_AT(1, "next_set = " << next_set);
    return FlowStepRunnerOutput(next_set, reach_set, next_time, chosen_step_size);
}

inline Set<TaskSearchPointCost>
evaluate(Map<TaskSearchPoint,TaskIOData<FlowStepRunnerInput,FlowStepRunnerOutput>> const& data) {
    Set<TaskSearchPointCost> result;

    Nat max_x = 0;
    for (auto entry : data) max_x = std::max(max_x,(Nat)entry.second.execution_time().count());
    for (auto entry : data) {
        CostType x = CostType(entry.second.execution_time().count())/max_x;
        CostType p = entry.second.output().step_size_used.get_d();
        result.insert(TaskSearchPointCost(entry.first, x/p));
    }
    return result;
}

struct VectorFieldEvolverFlowStepSerialRunner final : public SequentialRunnerBase<FlowStepRunnerInput,FlowStepRunnerOutput,FlowStepRunnerConfiguration> {
    VectorFieldEvolverFlowStepSerialRunner() : SequentialRunnerBase<FlowStepRunnerInput,FlowStepRunnerOutput,FlowStepRunnerConfiguration>(make_flow_step_runner_space()) { }
    FlowStepRunnerConfiguration to_configuration(FlowStepRunnerInput const& in, TaskSearchPoint const& p) const override { return Ariadne::to_configuration(in, p); }
    FlowStepRunnerOutput run_task(FlowStepRunnerInput const& in, FlowStepRunnerConfiguration const& cfg) const override { return Ariadne::run_task(in,cfg); }
    Set<TaskSearchPointCost> evaluate(Map<TaskSearchPoint,TaskIOData<FlowStepRunnerInput,FlowStepRunnerOutput>> const& data) const override { return Ariadne::evaluate(data); }
};

struct VectorFieldEvolverFlowStepDetachedRunner final : public DetachedRunnerBase<FlowStepRunnerInput,FlowStepRunnerOutput,FlowStepRunnerConfiguration> {
    VectorFieldEvolverFlowStepDetachedRunner() : DetachedRunnerBase<FlowStepRunnerInput,FlowStepRunnerOutput,FlowStepRunnerConfiguration>("step", make_flow_step_runner_space()) { }
    FlowStepRunnerConfiguration to_configuration(FlowStepRunnerInput const& in, TaskSearchPoint const& p) const override { return Ariadne::to_configuration(in, p); }
    FlowStepRunnerOutput run_task(FlowStepRunnerInput const& in, FlowStepRunnerConfiguration const& cfg) const override { return Ariadne::run_task(in,cfg); }
    Set<TaskSearchPointCost> evaluate(Map<TaskSearchPoint,TaskIOData<FlowStepRunnerInput,FlowStepRunnerOutput>> const& data) const override { return Ariadne::evaluate(data); }
};

struct VectorFieldEvolverFlowStepParameterSearchRunner final : public ParameterSearchRunnerBase<FlowStepRunnerInput,FlowStepRunnerOutput,FlowStepRunnerConfiguration> {
    VectorFieldEvolverFlowStepParameterSearchRunner(Nat concurrency) : ParameterSearchRunnerBase<FlowStepRunnerInput,FlowStepRunnerOutput,FlowStepRunnerConfiguration>("stp", make_flow_step_runner_space(),concurrency) { }
    FlowStepRunnerConfiguration to_configuration(FlowStepRunnerInput const& in, TaskSearchPoint const& p) const override { return Ariadne::to_configuration(in, p); }
    FlowStepRunnerOutput run_task(FlowStepRunnerInput const& in, FlowStepRunnerConfiguration const& cfg) const override { return Ariadne::run_task(in,cfg); }
    Set<TaskSearchPointCost> evaluate(Map<TaskSearchPoint,TaskIOData<FlowStepRunnerInput,FlowStepRunnerOutput>> const& data) const override { return Ariadne::evaluate(data); }
};

//! \brief A class for computing the evolution of a vector_field system.
//!
//! The actual evolution steps are performed by the Integrator class.
class VectorFieldEvolver
    : public EvolverBase< VectorField, LabelledEnclosure, typename VectorField::TimeType >,
      public TaskRunnableInterface<FlowStepRunnerInput,FlowStepRunnerOutput,FlowStepRunnerConfiguration>
{
  public:
    typedef TaskRunnerInterface<FlowStepRunnerInput,FlowStepRunnerOutput,FlowStepRunnerConfiguration> RunnerType;
    typedef VectorFieldEvolverConfiguration ConfigurationType;
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

    //! \brief Get the internal integrator.
    const IntegratorInterface* integrator() const { return _integrator.get(); }

    //! \brief Make an enclosure from a user set.
    EnclosureType enclosure(RealBox const&) const;
    EnclosureType enclosure(RealBox const&, EnclosureConfiguration const&) const;

    //! \brief Make an enclosure from a user set with variables.
    EnclosureType enclosure(RealVariablesBox const&) const;
    EnclosureType enclosure(RealVariablesBox const&, EnclosureConfiguration const&) const;

    //! \brief Make an enclosure from a computed box set.
    EnclosureType enclosure(ExactBoxType const&) const;
    EnclosureType enclosure(ExactBoxType const&, EnclosureConfiguration const&) const;

    //!@{
    //! \name Configuration for the class.
    //! \brief A reference to the configuration controlling the evolution.
    ConfigurationType& configuration() { return *this->_configuration; }
    const ConfigurationType& configuration() const { return *this->_configuration; }

    //! \brief The class which constructs functions for the enclosures.
    const FunctionFactoryType& function_factory() const;

    //!@}

    //! \brief Set the runner for the internal task
    void set_runner(SharedPointer<RunnerType> runner) override { this->_runner = runner; }

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
    SharedPointer<IntegratorInterface> _integrator;
    SharedPointer<ConfigurationType> _configuration;
    SharedPointer<RunnerType> _runner;
};


//! \brief Configuration for a VectorFieldEvolver, essentially for controlling the accuracy of continuous evolution methods.
class VectorFieldEvolverConfiguration : public ConfigurationInterface
{
  public:
    typedef ExactDouble RealType;
    typedef ApproximateDouble ApproximateRealType;

    //! \brief Default constructor gives reasonable values.
    VectorFieldEvolverConfiguration();

    virtual ~VectorFieldEvolverConfiguration() = default;

  private:

    //! \brief The maximum allowable step size for integration.
    //! Decreasing this value increases the accuracy of the computation.
    RealType _maximum_step_size;

    //! \brief The maximum allowable radius of a basic set during integration.
    //! Decreasing this value increases the accuracy of the computation of an over-approximation.
    RealType _maximum_enclosure_radius;

    //! \brief The maximum allowable approximation error in the parameter-to-space mapping of an enclosure set.
    //! Decreasing this value increases the accuracy of the computation of an over-approximation.
    RealType _maximum_spacial_error;

    //! \brief Enable reconditioning of basic sets (false by default).
    Bool _enable_reconditioning;

  public:

    const RealType& maximum_step_size() const { return _maximum_step_size; }
    Void set_maximum_step_size(const ApproximateRealType value) { _maximum_step_size = cast_exact(value); }

    const RealType& maximum_enclosure_radius() const { return _maximum_enclosure_radius; }
    Void set_maximum_enclosure_radius(const ApproximateRealType value) { _maximum_enclosure_radius = cast_exact(value); }

    const RealType& maximum_spacial_error() const { return _maximum_spacial_error; }
    Void set_maximum_spacial_error(const ApproximateRealType value) { _maximum_spacial_error = cast_exact(value); }

    const Bool& enable_reconditioning() const { return _enable_reconditioning; }
    Void set_enable_reconditioning(const Bool value) { _enable_reconditioning = value; }

  public:

    virtual OutputStream& _write(OutputStream& os) const;
};

} // namespace Ariadne

#endif // ARIADNE_VECTOR_FIELD_EVOLVER_HPP
