/***************************************************************************
 *            dynamics/vector_field_evolver.cpp
 *
 *  Copyright  2008-20  Alberto Casagrande, Pieter Collins
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

#include "function/functional.hpp"
#include "config.hpp"

#include "utility/macros.hpp"
#include "helper/array.hpp"
#include "utility/tuple.hpp"
#include "helper/stlio.hpp"
#include "helper/container.hpp"
#include "algebra/vector.hpp"
#include "function/function.hpp"
#include "function/constraint.hpp"
#include "dynamics/enclosure.hpp"
#include "dynamics/orbit.hpp"
#include "dynamics/inner_approximation.hpp"

#include "solvers/integrator.hpp"

#include "conclog/logging.hpp"

#include "dynamics/vector_field.hpp"
#include "dynamics/vector_field_evolver.hpp"

#include "symbolic/space.hpp"
#include "symbolic/assignment.hpp"
#include "symbolic/expression_set.hpp"

#include "pexplore/task_runner.tpl.hpp"

using namespace ConcLog;

namespace Ariadne {

// TODO: Move to Numeric module
inline PositiveValidatedUpperNumber abs(PositiveValidatedUpperNumber y) { return y; }

namespace {

template<class ES> List<ES> subdivide(const ES& enclosure) {
    List<ES> result;
    Pair<ES,ES> split=enclosure.split();
    result.append(split.first);
    result.append(split.second);
    return result;
}

} // namespace

VectorFieldEvolver::VectorFieldEvolver(SystemType const& system, ConfigurationType const& configuration) :
        TaskRunnable(configuration), _system(system.clone())
{ }

VectorFieldEvolver* VectorFieldEvolver::clone() const {
    return new VectorFieldEvolver(system(),configuration());
}

auto VectorFieldEvolver::enclosure(const ExactBoxType& box) const -> EnclosureType {
    return EnclosureType(box,this->system().state_space(),EnclosureConfiguration(TaylorFunctionFactory(ThresholdSweeper<FloatDP>(DoublePrecision(),Configuration<ThresholdSweeper<FloatDP>>().set_threshold(1e-7)))));
}

auto VectorFieldEvolver::orbit(RealVariablesBox const& initial_set, TimeType const& time, Semantics semantics) const -> Orbit<EnclosureType> {
    auto enclosure = EnclosureType(initial_set,this->system().state_space(),EnclosureConfiguration(TaylorFunctionFactory(ThresholdSweeper<FloatDP>(DoublePrecision(),Configuration<ThresholdSweeper<FloatDP>>().set_threshold(1e-7)))));
    enclosure.set_auxiliary(this->system().auxiliary_space(),this->system().auxiliary_mapping());
    return orbit(enclosure,time,semantics);
}

auto VectorFieldEvolver::orbit(RealExpressionBoundedConstraintSet const& initial_set, TimeType const& time, Semantics semantics) const -> Orbit<EnclosureType> {
    auto enclosure = EnclosureType(initial_set.euclidean_set(this->system().state_space()),this->system().state_space(),EnclosureConfiguration(TaylorFunctionFactory(ThresholdSweeper<FloatDP>(DoublePrecision(),Configuration<ThresholdSweeper<FloatDP>>().set_threshold(1e-7)))));
    enclosure.set_auxiliary(this->system().auxiliary_space(),this->system().auxiliary_mapping());
    return orbit(enclosure,time,semantics);
}

auto VectorFieldEvolver::orbit(EnclosureType const& initial_set, TimeType const& time, Semantics semantics) const -> Orbit<EnclosureType>
{
    CONCLOG_SCOPE_CREATE
    ARIADNE_PRECONDITION(this->system().state_auxiliary_space() == initial_set.state_auxiliary_space())
    auto result = std::make_shared<SynchronisedOrbit>(initial_set);

    try {
        WorkloadType workload([time](TimedEnclosureType const& timed_enclosure, SharedPointer<ProgressIndicator> indicator){
                                    indicator->update_current(timed_enclosure.first.get_d());
                                    indicator->update_final(time.get_d());
                                  },
                              std::bind_front(&VectorFieldEvolver::_process_timed_enclosure,this),time,semantics,result);
        _append_initial_set(workload,TimeStepType(0u),initial_set);
        workload.process();
    } catch (pExplore::NoActiveConstraintsException<VectorFieldEvolver>*) {
        CONCLOG_PRINTLN("Terminated early due to no active constraints, returning partial orbit.")
    }

    return std::move(*result);
}

Void VectorFieldEvolver::
_append_initial_set(WorkloadType& workload, TimeStepType const& initial_time, EnclosureType const& current_set) const
{
    if (decide(current_set.euclidean_set().bounding_box().radius() > configuration().maximum_enclosure_radius())) {
        CONCLOG_PRINTLN_AT(1,"set is too large, splitting")
        Pair<EnclosureType,EnclosureType> split_sets = current_set.split();
        if(!definitely(split_sets.first.is_empty())) { _append_initial_set(workload,initial_time,split_sets.first); }
        if(!definitely(split_sets.second.is_empty())) { _append_initial_set(workload,initial_time,split_sets.second); }
    } else {
        workload.append({initial_time,current_set});
    }
}

Void
VectorFieldEvolver::
_process_timed_enclosure(WorkloadType::Access& workload,
                         TimedEnclosureType const& current_timed_set,
                         TimeType const& maximum_time,
                         Semantics semantics,
                         SharedPointer<SynchronisedOrbit> result) const {
    CONCLOG_SCOPE_CREATE
    TimeStepType current_time=current_timed_set.first;
    EnclosureType current_set=current_timed_set.second;
    FloatDPUpperBound current_set_radius=current_set.euclidean_set().bounding_box().radius();

    CONCLOG_PRINTLN("#r="<<std::setw(5)<<std::left<<result->reach_size()
                             <<" t="<<std::setw(7)<<std::fixed<<current_time.get_d()
                             <<" p="<<std::setw(4)<<std::left<<current_set.number_of_parameters()
                             <<" r="<<std::setw(7)<<current_set.radius()
                             <<" c="<<current_set.centre())

    if (definitely(current_time>=maximum_time)) {
        result->adjoin_final(current_set);
    } else if (semantics == Semantics::UPPER and configuration().enable_subdivisions() and decide(current_set_radius>configuration().maximum_enclosure_radius())) {
        auto subdivisions=subdivide(current_set);
        for (auto const& sub : subdivisions)
            workload.append({current_time,sub});
    } else if(semantics == Semantics::LOWER && configuration().enable_premature_termination() && decide(current_set_radius>configuration().maximum_enclosure_radius())) {
        CONCLOG_PRINTLN("Terminating lower evolution at time " << current_time << " and set " << current_set << " due to maximum radius being exceeded.")
    } else {
        _process_timed_enclosure_step(workload,current_timed_set,maximum_time,semantics,result);
    }
}

Void
VectorFieldEvolver::
_process_timed_enclosure_step(WorkloadType::Access& workload,
                              TimedEnclosureType const& working_timed_set_model,
                              TimeType const& maximum_time,
                              Semantics semantics,
                              SharedPointer<SynchronisedOrbit> result) const
{
    CONCLOG_SCOPE_CREATE

    EnclosureType current_set;
    TimeStepType current_time;
    CONCLOG_PRINTLN_AT(1,"working_timed_set_model = "<<working_timed_set_model)
    make_lpair(current_time, current_set)=working_timed_set_model;

    CONCLOG_PRINTLN("current_time = "<<current_time)
    CONCLOG_PRINTLN("current_set = " << current_set)

    CONCLOG_PRINTLN("box = " << current_set.bounding_box())
    CONCLOG_PRINTLN("radius = " << current_set.euclidean_set().bounding_box().radius())

    /////////////// Main Evolution ////////////////////////////////

    runner()->push({system().dynamic_function(), current_set, current_time});
    auto out = runner()->pull();

    result->adjoin_reach(out.reach);

    if (this->configuration().enable_clobbering()) {
        auto state_function = out.evolve.euclidean_set().state_function();
        ValidatedVectorMultivariateTaylorFunctionModelDP tf = dynamic_cast<ValidatedVectorMultivariateTaylorFunctionModelDP const&>(state_function.reference());
        tf.clobber();
        Enclosure new_evolve(out.evolve.euclidean_set().domain(),tf, out.evolve.euclidean_set().time_function(), out.evolve.constraints(),out.evolve.euclidean_set().configuration());
        LabelledEnclosure new_labelled_evolve(new_evolve,out.evolve.state_space(),out.evolve.auxiliary_space());
        new_labelled_evolve.set_auxiliary_mapping(out.evolve.auxiliary_mapping());
        result->adjoin_intermediate(new_labelled_evolve);
        workload.append({out.time,new_labelled_evolve});
    } else {
        result->adjoin_intermediate(out.evolve);
        workload.append({out.time,out.evolve});
    }
}

}

namespace pExplore {

using Ariadne::cast_exact;
using Ariadne::NativeSimplex;
using Ariadne::ParallelLinearisationContractor;
using Ariadne::NonlinearCandidateValidationInnerApproximator;

auto Task<VectorFieldEvolver>::run(TaskInput<VectorFieldEvolver> const& in, Configuration<VectorFieldEvolver> const& cfg) const -> TaskOutput<VectorFieldEvolver> {
    LabelledEnclosure next_set = in.current_set;
    LabelledEnclosure reach_set = in.current_set;
    Dyadic next_time = in.current_time;
    Dyadic chosen_step_size = Ariadne::cast_exact(cfg.maximum_step_size());

    if(cfg.enable_reconditioning() and cast_exact(norm(next_set.state_function().errors())) > cast_exact(cfg.maximum_spacial_error())) {
        CONCLOG_PRINTLN_AT(1,"reconditioning from errors "<< next_set.state_function().errors());
        next_set.recondition();
    }

    auto set_bounds=cast_exact_box(next_set.euclidean_set().bounding_box());
    CONCLOG_PRINTLN_VAR_AT(1, set_bounds);

    auto flow_model = cfg.integrator().flow_step(in.dynamic, set_bounds,chosen_step_size);

    CONCLOG_PRINTLN_VAR_AT(1, chosen_step_size);
    CONCLOG_PRINTLN_VAR_AT(1, flow_model);
    next_time += chosen_step_size;
    CONCLOG_PRINTLN_VAR_AT(1, next_time);
    reach_set.apply_full_reach_step(flow_model);
    CONCLOG_PRINTLN_VAR_AT(1, reach_set);
    next_set.apply_fixed_evolve_step(flow_model, chosen_step_size);
    CONCLOG_PRINTLN_VAR_AT(1, next_set);

    Lazy<LabelledEnclosure> inner_next_set([next_set]{
        auto approximator = NonlinearCandidateValidationInnerApproximator(ParallelLinearisationContractor(NativeSimplex(),2,0));
        auto inner_evolve = approximator.compute_from(next_set);
        return inner_evolve.clone();
    });

    return {next_set, reach_set, inner_next_set, next_time};
}

} // namespace pExplore
