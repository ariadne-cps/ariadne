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

#include <array>

#include "function/functional.hpp"
#include "config.hpp"

#include "utility/macros.hpp"
#include "utility/array.hpp"
#include "utility/tuple.hpp"
#include "helper/stlio.hpp"
#include "utility/container.hpp"
#include "algebra/vector.hpp"
#include "function/function.hpp"
#include "function/taylor_function.hpp"
#include "function/constraint.hpp"
#include "dynamics/enclosure.hpp"
#include "dynamics/orbit.hpp"

#include "solvers/integrator.hpp"

#include "conclog/logging.hpp"

#include "dynamics/vector_field.hpp"
#include "dynamics/vector_field_evolver.hpp"

#include "symbolic/space.hpp"
#include "symbolic/assignment.hpp"
#include "symbolic/expression_set.hpp"

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

struct CarriedExpansionSnapshot {
    static constexpr std::size_t degree_slots=32u;
    std::array<unsigned long long,degree_slots> below_tight_count{};
    std::array<unsigned long long,degree_slots> bridge_count{};
    std::array<unsigned long long,degree_slots> above_loose_count{};
    std::array<double,degree_slots> below_tight_abs_mass{};
    std::array<double,degree_slots> bridge_abs_mass{};
    std::array<double,degree_slots> above_loose_abs_mass{};
};

CarriedExpansionSnapshot
carried_expansion_snapshot(ValidatedVectorMultivariateFunctionPatch const& mapping)
{
    CarriedExpansionSnapshot profile;
    for(SizeType component=0u; component!=mapping.result_size(); ++component) {
        auto const& generic_component=mapping.get(component);
        auto const* taylor_component=
            dynamic_cast<ValidatedScalarMultivariateTaylorFunctionModelDP const*>(
                generic_component.raw_pointer());
        ARIADNE_ASSERT(taylor_component!=nullptr);
        auto const& expansion=taylor_component->expansion();
        for(auto const& term : expansion) {
            auto degree=static_cast<std::size_t>(term.index().degree());
            if(degree>=CarriedExpansionSnapshot::degree_slots) {
                degree=CarriedExpansionSnapshot::degree_slots-1u;
            }
            const double magnitude=std::abs(term.coefficient().get_d());
            if(magnitude<3e-14) {
                ++profile.below_tight_count[degree];
                profile.below_tight_abs_mass[degree]+=magnitude;
            } else if(magnitude<1e-12) {
                ++profile.bridge_count[degree];
                profile.bridge_abs_mass[degree]+=magnitude;
            } else {
                ++profile.above_loose_count[degree];
                profile.above_loose_abs_mass[degree]+=magnitude;
            }
        }
    }
    return profile;
}

Void
print_carried_expansion_snapshot(
    String const& stage,
    SizeType step,
    ValidatedVectorMultivariateFunctionPatch const& mapping,
    Sweeper<FloatDP> const& sweeper)
{
    auto const profile=carried_expansion_snapshot(mapping);
    for(std::size_t degree=0u;
        degree!=CarriedExpansionSnapshot::degree_slots; ++degree) {
        const auto total_count=
            profile.below_tight_count[degree]
            +profile.bridge_count[degree]
            +profile.above_loose_count[degree];
        if(total_count==0u) { continue; }
        std::cerr << "[CarriedExpansionSnapshot]"
                  << " sweeper=" << sweeper
                  << " stage=" << stage
                  << " step=" << step
                  << " degree=" << degree
                  << " below_3e-14_count=" << profile.below_tight_count[degree]
                  << " bridge_3e-14_to_1e-12_count=" << profile.bridge_count[degree]
                  << " above_1e-12_count=" << profile.above_loose_count[degree]
                  << " below_3e-14_abs_mass=" << profile.below_tight_abs_mass[degree]
                  << " bridge_3e-14_to_1e-12_abs_mass=" << profile.bridge_abs_mass[degree]
                  << " above_1e-12_abs_mass=" << profile.above_loose_abs_mass[degree]
                  << std::endl;
    }
}

} // namespace

VectorFieldEvolver::VectorFieldEvolver(const SystemType& system, const IntegratorInterface& i)
    : _system(system.clone())
    , _integrator(i.clone())
    , _configuration(new ConfigurationType())
{
}

auto VectorFieldEvolver::enclosure(const ExactBoxType& box) const -> EnclosureType {
    return EnclosureType(box,this->system().state_space(),EnclosureConfiguration(this->function_factory()));
}

auto VectorFieldEvolver::function_factory() const -> FunctionFactoryType const& {
    return std::dynamic_pointer_cast<const IntegratorBase>(this->_integrator)->function_factory();
}

auto VectorFieldEvolver::orbit(RealVariablesBox const& initial_set, TimeType const& time, Semantics semantics) const -> Orbit<EnclosureType> {
    auto enclosure = EnclosureType(initial_set,this->system().state_space(),EnclosureConfiguration(this->function_factory()));
    enclosure.set_auxiliary(this->system().auxiliary_space(),this->system().auxiliary_mapping());
    return orbit(enclosure,time,semantics);
}

auto VectorFieldEvolver::orbit(RealExpressionBoundedConstraintSet const& initial_set, TimeType const& time, Semantics semantics) const -> Orbit<EnclosureType> {
    auto enclosure = EnclosureType(initial_set.euclidean_set(this->system().state_space()),this->system().state_space(),EnclosureConfiguration(this->function_factory()));
    enclosure.set_auxiliary(this->system().auxiliary_space(),this->system().auxiliary_mapping());
    return orbit(enclosure,time,semantics);
}

auto VectorFieldEvolver::orbit(EnclosureType const& initial_set, TimeType const& time, Semantics semantics) const -> Orbit<EnclosureType>
{
    CONCLOG_SCOPE_CREATE
    ARIADNE_PRECONDITION(this->system().state_auxiliary_space() == initial_set.state_auxiliary_space())

    auto result = std::make_shared<SynchronisedOrbit>(initial_set);
    WorkloadType workload([time](TimedEnclosureType const& timed_enclosure, SharedPointer<ProgressIndicator> indicator){
                                indicator->update_current(timed_enclosure.first.get_d());
                                indicator->update_final(time.get_d());
                              },
                          std::bind_front(&VectorFieldEvolver::_process_timed_enclosure,this),time,semantics,result);
    _append_initial_set(workload,TimeStepType(0u),initial_set);
    workload.process();

    return std::move(*result);
}

Void VectorFieldEvolver::
_append_initial_set(WorkloadType& workload, TimeStepType const& initial_time, EnclosureType const& current_set) const
{
    if (possibly(current_set.euclidean_set().bounding_box().radius() > this->_configuration->maximum_enclosure_radius())) {
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
    } else if (semantics == Semantics::UPPER and this->_configuration->enable_subdivisions() and decide(current_set_radius>this->_configuration->maximum_enclosure_radius())) {
        // Subdivide
        List< EnclosureType > subdivisions=subdivide(current_set);
        for(SizeType i=0; i!=subdivisions.size(); ++i) {
            EnclosureType const& subdivided_set_model=subdivisions[i];
            workload.append({current_time,subdivided_set_model});
        }
    } else if (semantics == Semantics::LOWER and decide(current_set_radius>this->_configuration->maximum_enclosure_radius())) {
        CONCLOG_PRINTLN("Terminating lower evolution at time " << current_time << " and set " << current_set << " due to maximum radius being exceeded.")
    } else {
        this->_process_timed_enclosure_step(workload,current_timed_set,maximum_time,semantics,result);
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
    typedef EffectiveVectorMultivariateFunction FunctionType;

    EnclosureType current_set=working_timed_set_model.second;
    TimeStepType current_time=working_timed_set_model.first;
    SharedPointer<PreconditionedTaylorSeriesState> carried_preconditioned_state=
        working_timed_set_model.preconditioned_state;
    CONCLOG_PRINTLN_AT(1,"working_timed_set_model time = "<<current_time)

    CONCLOG_PRINTLN("current_time = "<<current_time)
    CONCLOG_PRINTLN("current_set = " << current_set)

    CONCLOG_PRINTLN("box = " << current_set.bounding_box())
    CONCLOG_PRINTLN("radius = " << current_set.euclidean_set().bounding_box().radius())

    IntegratorInterface const* integrator=this->_integrator.operator->();
    auto const* preconditioned_integrator=
        dynamic_cast<PreconditionedGradedTaylorSeriesIntegrator const*>(integrator);

    // The preconditioned integrator carries the local-initial Taylor mapping
    // across steps. Standard reconditioning would deliberately discard part of
    // that dependence, so it is bypassed on this specialised path.
    if (preconditioned_integrator==nullptr
        && this->_configuration->enable_reconditioning()
        && possibly(norm(current_set.state_function().errors()) > this->_configuration->maximum_spacial_error())) {
        current_set.recondition();
        workload.append({current_time,current_set});
        return;
    }

    /////////////// Main Evolution ////////////////////////////////
    const FunctionType& dynamic=_system->function();

    // Set evolution parameters
    const StepSizeType maximum_step_size=this->_configuration->maximum_step_size();

    // Get bounding boxes for time and space bounding_box
    auto current_set_bounds=cast_exact_box(current_set.euclidean_set().bounding_box());
    CONCLOG_PRINTLN("current_set_bounds = "<<current_set_bounds)

    StepSizeType step_size;
    EnclosureType reach_set=current_set;
    EnclosureType next_set=current_set;
    SharedPointer<PreconditionedTaylorSeriesState> next_preconditioned_state;

    if(preconditioned_integrator!=nullptr) {
        if(preconditioned_integrator->diagnostics()
           && result->reach_size()==1u
           && preconditioned_integrator->preconditioning()==TaylorSeriesPreconditioning::QR) {
            PreconditionedGradedTaylorSeriesIntegrator identity_probe=
                *preconditioned_integrator;
            identity_probe.set_preconditioning(TaylorSeriesPreconditioning::IDENTITY);
            identity_probe.set_diagnostics(true);

            PreconditionedGradedTaylorSeriesIntegrator qr_probe=
                *preconditioned_integrator;
            qr_probe.set_preconditioning(TaylorSeriesPreconditioning::QR);
            qr_probe.set_diagnostics(true);

            auto identity_state=
                identity_probe.precondition(current_set.state_function());
            auto qr_state=
                qr_probe.precondition(current_set.state_function());

            std::cerr << "[PreconditionedSecondStepState]"
                      << " physical_error=" << current_set.state_function().error()
                      << " physical_range=" << current_set.state_function().range()
                      << " identity_domy=" << identity_state.local_domain()
                      << " identity_A=" << identity_state.linear_map()
                      << " identity_normalized_error=" << identity_state.normalised_mapping().error()
                      << " qr_domy=" << qr_state.local_domain()
                      << " qr_A=" << qr_state.linear_map()
                      << " qr_normalized_error=" << qr_state.normalised_mapping().error()
                      << std::endl;

            auto identity_step=
                identity_probe.step(dynamic,identity_state,suggest(maximum_step_size));
            auto qr_step=
                qr_probe.step(dynamic,qr_state,suggest(maximum_step_size));

            std::cerr << "[PreconditionedSecondStepResult]"
                      << " identity_h=" << identity_step.time_step()
                      << " identity_flow_error=" << identity_step.flowpipe_mapping().error()
                      << " qr_h=" << qr_step.time_step()
                      << " qr_flow_error=" << qr_step.flowpipe_mapping().error()
                      << std::endl;
        }

        PreconditionedTaylorSeriesState local_state=
            carried_preconditioned_state
                ? *carried_preconditioned_state
                : preconditioned_integrator->precondition(current_set.state_function());

        const SizeType diagnostic_step=result->reach_size()+1u;
        const bool carried_snapshot_step=
            diagnostic_step==500u || diagnostic_step==1000u
            || diagnostic_step==1500u || diagnostic_step==2000u;
        if(preconditioned_integrator->carried_expansion_diagnostics()
           && carried_snapshot_step) {
            print_carried_expansion_snapshot(
                "input",diagnostic_step,
                local_state.normalised_mapping(),
                preconditioned_integrator->sweeper());
        }

        PreconditionedTaylorSeriesStep local_step=
            preconditioned_integrator->step(
                dynamic,local_state,suggest(maximum_step_size));

        if(preconditioned_integrator->carried_expansion_diagnostics()
           && carried_snapshot_step) {
            print_carried_expansion_snapshot(
                "output",diagnostic_step,
                local_step.final_state().normalised_mapping(),
                preconditioned_integrator->sweeper());
        }

        step_size=local_step.time_step();

        reach_set.apply_parameterised_full_reach_step(
            local_step.flowpipe_mapping());

        ValidatedVectorMultivariateFunctionPatch const& final_mapping=
            local_step.evolved_mapping();
        next_set.apply_parameterised_fixed_evolve_step(
            final_mapping,step_size);
        next_preconditioned_state=
            std::make_shared<PreconditionedTaylorSeriesState>(
                local_step.final_state());

        // Temporary diagnostic: compare endpoint-first composition with the
        // old endpoint-after-full-flowpipe ordering.
        if(preconditioned_integrator->diagnostics()
           && result->reach_size()<20u) {
            auto const& final_state=local_step.final_state();
            ValidatedVectorMultivariateFunctionPatch old_order_endpoint=
                partial_evaluate(
                    local_step.flowpipe_mapping(),
                    local_step.flowpipe_mapping().argument_size()-1u,
                    step_size);
            std::cerr << "[PreconditionedStepDiagnostic]"
                      << " step=" << result->reach_size()
                      << " carried_state=" << (carried_preconditioned_state ? 1 : 0)
                      << " t=" << current_time
                      << " h=" << step_size
                      << " state_error=" << current_set.state_function().error()
                      << " normalized_error=" << local_state.normalised_mapping().error()
                      << " normalized_range=" << local_state.normalised_mapping().range()
                      << " A=" << local_state.linear_map()
                      << " flowpipe_error=" << local_step.flowpipe_mapping().error()
                      << " old_order_endpoint_error=" << old_order_endpoint.error()
                      << " final_mapping_error=" << final_mapping.error()
                      << " final_normalized_error=" << final_state.normalised_mapping().error()
                      << " final_normalized_range=" << final_state.normalised_mapping().range()
                      << " final_A=" << final_state.linear_map()
                      << std::endl;
        }
    } else {
        FlowStepModelType flow_model=
            integrator->flow_step(
                dynamic,current_set_bounds,suggest(maximum_step_size));

        step_size=static_cast<StepSizeType>(
            flow_model.domain()[flow_model.argument_size()-1u].upper_bound());

        reach_set.apply_full_reach_step(flow_model);
        next_set.apply_fixed_evolve_step(flow_model,step_size);
    }

    TimeStepType next_time=current_time+TimeStepType(step_size);

    result->adjoin_reach(reach_set);
    result->adjoin_intermediate(next_set);
    if(next_preconditioned_state) {
        workload.append({next_time,next_set,std::move(next_preconditioned_state)});
    } else {
        workload.append({next_time,next_set});
    }
}


VectorFieldEvolverConfiguration::VectorFieldEvolverConfiguration()
{
    set_maximum_step_size(1);
    set_maximum_enclosure_radius(100.0);
    set_maximum_spacial_error(1e-2);
    set_enable_reconditioning(true);
    set_enable_subdivisions(false);
}


OutputStream&
VectorFieldEvolverConfiguration::_write(OutputStream& os) const
{
    os << "VectorFieldEvolverConfiguration("
       << "\n maximum_step_size=" << maximum_step_size()
       << ",\n maximum_enclosure_radius=" << maximum_enclosure_radius()
       << ",\n maximum_spacial_error=" << maximum_spacial_error()
       << ",\n enable_reconditioning=" << enable_reconditioning()
       << ",\n enable_subdivisions=" << enable_subdivisions()
       << "\n)";
    return os;
}


}  // namespace Ariadne

