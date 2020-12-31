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

#include "../function/functional.hpp"
#include "../config.hpp"

#include "../utility/macros.hpp"
#include "../utility/array.hpp"
#include "../utility/tuple.hpp"
#include "../utility/stlio.hpp"
#include "../utility/container.hpp"
#include "../algebra/vector.hpp"
#include "../function/function.hpp"
#include "../function/constraint.hpp"
#include "../function/taylor_function.hpp"
#include "../dynamics/enclosure.hpp"
#include "../dynamics/orbit.hpp"

#include "../solvers/integrator.hpp"

#include "../output/logging.hpp"
#include "../output/progress_indicator.hpp"

#include "../concurrency/loggable_smart_thread.hpp"
#include "../concurrency/buffer.hpp"
#include "../concurrency/task_parameter_point.hpp"

#include "../dynamics/vector_field.hpp"
#include "../dynamics/vector_field_evolver.hpp"

#include "../symbolic/space.hpp"
#include "../symbolic/assignment.hpp"

namespace Ariadne {

namespace {

template<class ES> List<ES> subdivide(const ES& enclosure) {
    List<ES> result;
    Pair<ES,ES> split=enclosure.split();
    result.append(split.first);
    result.append(split.second);
    return result;
}

} // namespace

EffectiveVectorMultivariateFunction make_auxiliary_function(
    Space<Real> const& state_space,
    List<RealAssignment> const& algebraic);

EffectiveVectorMultivariateFunction make_dynamic_function(
    Space<Real> const& space,
    List<RealAssignment> const& algebraic,
    List<DottedRealAssignment> const& differential);


VectorField::VectorField(List<DottedRealAssignment> const& dynamics)
    : VectorField(dynamics, List<RealAssignment>())
{
}

VectorField::VectorField(List<DottedRealAssignment> const& dynamics, List<RealAssignment> const& auxiliary)
    : _dynamics(dynamics), _auxiliary(auxiliary)
    , _dynamic_function(make_dynamic_function(left_hand_sides(dynamics),auxiliary,dynamics))
    , _auxiliary_function(make_auxiliary_function(left_hand_sides(dynamics),auxiliary))
{
}

VectorField::VectorField(EffectiveVectorMultivariateFunction const& function)
    : _dynamic_function(function), _auxiliary_function(0u,function.domain())
{
    ARIADNE_PRECONDITION(function.result_size()==function.argument_size());
}

RealSpace VectorField::state_space() const {
    return RealSpace(left_hand_sides(this->_dynamics));
}

RealSpace VectorField::auxiliary_space() const {
    return RealSpace(left_hand_sides(this->_auxiliary));
}


OutputStream& operator<<(OutputStream& os, const VectorField& vf) {
    os << "VectorField( dynamic_function = " << vf.dynamic_function() << ", "
          "auxiliary_function = " << vf.auxiliary_function() << ", "
          "dynamics = " << vf._dynamics << ", "
          "auxiliary = " << vf._auxiliary << ")";
    return os;
}

// Allow subdivisions in upper evolution
const Bool ENABLE_SUBDIVISIONS = false;
// Allow premature termination of lower evolution
const Bool ENABLE_PREMATURE_TERMINATION = false;

using std::shared_ptr;

class DegenerateCrossingException { };



VectorFieldEvolver::VectorFieldEvolver(const SystemType& system, const IntegratorInterface& i)
    : _sys_ptr(system.clone())
    , _integrator(i.clone())
    , _configuration(new ConfigurationType())
{
}

typename VectorFieldEvolver::EnclosureType VectorFieldEvolver::enclosure(const ExactBoxType& box) const {
    return EnclosureType(box,this->system().state_space(),EnclosureConfiguration(this->function_factory()));
}

typename VectorFieldEvolver::EnclosureType VectorFieldEvolver::enclosure(const ExactBoxType& box, const EnclosureConfiguration& config) const {
    return EnclosureType(box,this->system().state_space(),config);
}

typename VectorFieldEvolver::EnclosureType VectorFieldEvolver::enclosure(const RealBox& box) const {
    return EnclosureType(box,this->system().state_space(),EnclosureConfiguration(this->function_factory()));
}

typename VectorFieldEvolver::EnclosureType VectorFieldEvolver::enclosure(const RealBox& box, const EnclosureConfiguration& config) const {
    return EnclosureType(box,this->system().state_space(),config);
}

typename VectorFieldEvolver::EnclosureType VectorFieldEvolver::enclosure(const RealVariablesBox& box) const {
    return EnclosureType(box,this->system().state_space(),EnclosureConfiguration(this->function_factory()));
}

typename VectorFieldEvolver::EnclosureType VectorFieldEvolver::enclosure(const RealVariablesBox& box, const EnclosureConfiguration& config) const {
    return EnclosureType(box,this->system().state_space(),config);
}

typename VectorFieldEvolver::FunctionFactoryType const& VectorFieldEvolver::function_factory() const {
    return std::dynamic_pointer_cast<const IntegratorBase>(this->_integrator)->function_factory();
}


Orbit<VectorFieldEvolver::EnclosureType>
VectorFieldEvolver::
orbit(const EnclosureType& initial_set,
      const TimeType& time,
      Semantics semantics) const
{
    Orbit<EnclosureType> orbit(initial_set);
    EnclosureListType final;
    EnclosureListType reachable;
    EnclosureListType intermediate;
    this->_evolution(final,reachable,intermediate,
                     initial_set,time,semantics,false);
    orbit.adjoin_intermediate(intermediate);
    orbit.adjoin_reach(reachable);
    orbit.adjoin_final(final);
    return orbit;
}

Void VectorFieldEvolver::
_append_initial_set(List<TimedEnclosureType>& working_sets,
                   const TimeStepType& initial_time,
                   const EnclosureType& current_set) const
{
    ARIADNE_LOG_SCOPE_CREATE;
    if (possibly(current_set.euclidean_set().bounding_box().radius() > this->_configuration->maximum_enclosure_radius())) {
        ARIADNE_LOG_PRINTLN("initial set too large, splitting");
        Pair<EnclosureType,EnclosureType> split_sets = current_set.split();
        if(!definitely(split_sets.first.is_empty())) { _append_initial_set(working_sets,initial_time,split_sets.first); }
        if(!definitely(split_sets.second.is_empty())) { _append_initial_set(working_sets,initial_time,split_sets.second); }
    } else {
        working_sets.push_back(make_pair(initial_time,current_set));
    }
}

class FlowStepInput {
  public:
    FlowStepInput(LabelledEnclosure const& current_set, FloatDPExactBox const& current_set_bounds, Dyadic const& current_time, Dyadic const& previous_step_size) :
                  _current_set(current_set), _current_set_bounds(current_set_bounds), _current_time(current_time), _previous_step_size(previous_step_size) { }
    LabelledEnclosure const& _current_set;
    FloatDPExactBox const& _current_set_bounds;
    Dyadic const& _current_time;
    Dyadic const& _previous_step_size;
};

class FlowStepPoint {

};

class FlowStepConfiguration {
  public:
    FlowStepConfiguration(TaylorPicardIntegrator const& integrator) : _integrator(integrator.clone()){ }
    TaylorPicardIntegrator* _integrator;
};

class FlowStepOutput {
  public:
    FlowStepOutput(LabelledEnclosure const& evolve, LabelledEnclosure const& reach, Dyadic const& time, Dyadic const& step_size_used) :
                   _evolve(evolve), _reach(reach), _time(time), _step_size_used(step_size_used) { }
    LabelledEnclosure const _evolve;
    LabelledEnclosure const _reach;
    Dyadic const _time;
    Dyadic const _step_size_used;
};

class ConcurrentRunner {
    typedef Buffer<Pair<FlowStepInput,FlowStepPoint>> InputBufferType;
    typedef Buffer<Pair<FlowStepOutput,FlowStepPoint>> OutputBufferType;
  public:
    ConcurrentRunner(FlowStepPoint const& initial_point, EffectiveVectorMultivariateFunction const& dynamic, TaylorPicardIntegrator const& integrator, Dyadic const& maximum_step_size);

    ~ConcurrentRunner();

    Void push(FlowStepInput const& input);

    FlowStepOutput pull();

  private:

    Void _loop() {
        ARIADNE_LOG_SCOPE_CREATE;
        while(true) {
            std::unique_lock<std::mutex> locker(_input_mutex);
            _input_availability.wait(locker, [this]() { return _input_buffer.size()>0 || _terminate; });
            if (_terminate) break;
            auto pkg = _input_buffer.pop();
            _output_buffer.push({_task(pkg.first,_to_configuration(pkg.second)),pkg.second});
            _output_availability.notify_all();
        }
    }

    FlowStepConfiguration _to_configuration(FlowStepPoint const& p) {
        FlowStepConfiguration result(_integrator);
        return result;
    }

    FlowStepOutput _task(FlowStepInput const& in, FlowStepConfiguration const& cfg) {

        LabelledEnclosure next_set = in._current_set;
        LabelledEnclosure reach_set = in._current_set;
        Dyadic next_time = in._current_time;
        Dyadic chosen_step_size = _maximum_step_size;
        FlowStepModelType flow_model = cfg._integrator->flow_step(_dynamic, in._current_set_bounds, in._previous_step_size,chosen_step_size);
        ARIADNE_LOG_PRINTLN("step_size = " << chosen_step_size);
        ARIADNE_LOG_PRINTLN_AT(1, "flow_model = " << flow_model);
        next_time += chosen_step_size;
        ARIADNE_LOG_PRINTLN_AT(1, "next_time = " << next_time)
        reach_set.apply_full_reach_step(flow_model);
        ARIADNE_LOG_PRINTLN_AT(1, "reach_set = " << reach_set);
        next_set.apply_fixed_evolve_step(flow_model, chosen_step_size);
        ARIADNE_LOG_PRINTLN_AT(1, "next_set = " << next_set);

        return FlowStepOutput(next_set,reach_set,next_time,chosen_step_size);
    }

  private:
    // Initial point
    FlowStepPoint const _initial_point;
    // Constants
    EffectiveVectorMultivariateFunction const _dynamic;
    TaylorPicardIntegrator const& _integrator;
    Dyadic const _maximum_step_size;
    // Synchronization
    LoggableSmartThread _thread;
    InputBufferType _input_buffer;
    OutputBufferType _output_buffer;
    std::atomic<bool> _terminate;
    std::mutex _input_mutex;
    std::condition_variable _input_availability;
    std::mutex _output_mutex;
    std::condition_variable _output_availability;
};

ConcurrentRunner::ConcurrentRunner(FlowStepPoint const& initial_point, EffectiveVectorMultivariateFunction const& dynamic, TaylorPicardIntegrator const& integrator, Dyadic const& maximum_step_size)
    :  _initial_point(initial_point), _dynamic(dynamic), _integrator(integrator), _maximum_step_size(maximum_step_size), _thread("step", [this]() { _loop(); }),
       _input_buffer(InputBufferType(1)),_output_buffer(OutputBufferType(1)),
       _terminate(false)
{
    _thread.activate();
}

ConcurrentRunner::~ConcurrentRunner() {
    _terminate = true;
    _input_availability.notify_all();
}

Void
ConcurrentRunner::push(FlowStepInput const& input)
{
    _input_buffer.push({input,_initial_point});
    _input_availability.notify_all();
}

FlowStepOutput
ConcurrentRunner::pull() {
    std::unique_lock<std::mutex> locker(_output_mutex);
    _output_availability.wait(locker, [this]() { return _output_buffer.size()>0; });
    auto result = _output_buffer.pop();
    return result.first;
}


Void
VectorFieldEvolver::
_evolution(EnclosureListType& final_sets,
           EnclosureListType& reach_sets,
           EnclosureListType& intermediate_sets,
           const EnclosureType& initial_set,
           const TimeType& maximum_time,
           Semantics semantics,
           Bool reach) const
{
    ARIADNE_LOG_SCOPE_CREATE;

    List< TimedEnclosureType > working_sets;

    {
        TimeStepType initial_time = 0u;
        // Append the initial set, possibly splitting it
        _append_initial_set(working_sets,initial_time,initial_set);
    }

    // Track the previous step size used to properly find the starting step size
    StepSizeType previous_step_size = 0;

    ProgressIndicator initials_indicator(working_sets.size());
    ProgressIndicator time_indicator(maximum_time.get_d());

    ConcurrentRunner runner(FlowStepPoint(),_sys_ptr->dynamic_function(),
                            *dynamic_cast<TaylorPicardIntegrator*>(this->_integrator.operator->()),
                            _configuration->maximum_step_size());

    while(!working_sets.empty()) {
        TimedEnclosureType current_timed_set=working_sets.back();
        working_sets.pop_back();
        TimeStepType current_time=current_timed_set.first;
        EnclosureType current_set_model=current_timed_set.second;
        FloatDPUpperBound current_set_radius=current_set_model.euclidean_set().bounding_box().radius();

        if(definitely(current_time>=maximum_time)) {
            final_sets.adjoin(current_set_model);
        } else if(semantics == Semantics::UPPER && ENABLE_SUBDIVISIONS
                  && decide(current_set_radius>this->_configuration->maximum_enclosure_radius())) {
            // Subdivide
            List< EnclosureType > subdivisions=subdivide(current_set_model);
            for(SizeType i=0; i!=subdivisions.size(); ++i) {
                EnclosureType const& subdivided_set_model=subdivisions[i];
                working_sets.push_back(make_pair(current_time,subdivided_set_model));
            }
        } else if(semantics == Semantics::LOWER && ENABLE_PREMATURE_TERMINATION && decide(current_set_radius>this->_configuration->maximum_enclosure_radius())) {
            ARIADNE_WARN("Terminating lower evolution at time " << current_time << " and set " << current_set_model << " due to maximum radius being exceeded.")
        } else {
            // Compute evolution
            this->_evolution_step(working_sets,
                                  final_sets,reach_sets,intermediate_sets,
                                  current_timed_set,previous_step_size,maximum_time,
                                  semantics,reach,runner);
        }

        ARIADNE_LOG_PRINTLN("#w="<<std::setw(4)<<working_sets.size()
                            <<"#r="<<std::setw(4)<<std::left<<reach_sets.size()
                            <<" t="<<std::setw(7)<<std::fixed<<current_time.get_d()
                            <<" p="<<std::setw(4)<<std::left<<current_set_model.number_of_parameters()
                            <<" r="<<std::setw(7)<<current_set_model.radius()
                            <<" c="<<current_set_model.centre());


        initials_indicator.update_current(final_sets.size());
        time_indicator.update_current(current_time.get_d());
        if (initials_indicator.final_value() > 1) { ARIADNE_LOG_SCOPE_PRINTHOLD("[" << time_indicator.symbol() << "] " << initials_indicator.percentage() << "% of sets, " << time_indicator.percentage() << "% of current set"); }
        else ARIADNE_LOG_SCOPE_PRINTHOLD("[" << time_indicator.symbol() << "] " << time_indicator.percentage() << "%");
    }
    ARIADNE_LOG_PRINTLN("Finished evolution");
}


Void
VectorFieldEvolver::
_evolution_step(List< TimedEnclosureType >& working_sets,
                EnclosureListType& final_sets,
                EnclosureListType& reach_sets,
                EnclosureListType& intermediate_sets,
                const TimedEnclosureType& working_timed_set,
                StepSizeType& previous_step_size,
                const TimeType& maximum_time,
                Semantics semantics,
                Bool reach,
                ConcurrentRunner& runner) const
{
    ARIADNE_LOG_SCOPE_CREATE;

    EnclosureType current_set;
    TimeStepType current_time;
    ARIADNE_LOG_PRINTLN_AT(1,"working_timed_set_model = "<<working_timed_set);
    make_lpair(current_time,current_set)=working_timed_set;

    ARIADNE_LOG_PRINTLN("current_time = "<<current_time);
    ARIADNE_LOG_PRINTLN("current_set_model = "<<current_set);

    ARIADNE_LOG_PRINTLN("box = "<<current_set.bounding_box());
    ARIADNE_LOG_PRINTLN("radius = "<<current_set.euclidean_set().bounding_box().radius());

    // Test to see if set requires reconditioning
    if(this->_configuration->enable_reconditioning() &&
       possibly(norm(current_set.state_function().errors()) > this->_configuration->maximum_spacial_error())) {
        ARIADNE_LOG_PRINTLN("reconditioning from errors "<<current_set.state_function().errors());
        current_set.recondition();
    }

    /////////////// Main Evolution ////////////////////////////////

    // Get bounding boxes for time and space bounding_box
    auto current_set_bounds=cast_exact_box(current_set.euclidean_set().bounding_box());
    ARIADNE_LOG_PRINTLN("current_set_bounds = "<<current_set_bounds);

    // Push inputs
    runner.push(FlowStepInput(current_set,current_set_bounds,current_time,previous_step_size));
    // Pull outputs
    auto out = runner.pull();
    // Save outputs
    reach_sets.adjoin(out._reach);
    intermediate_sets.adjoin(out._evolve);
    working_sets.push_back(make_pair(out._time,out._evolve));
    previous_step_size = out._step_size_used;
}


VectorFieldEvolverConfiguration::VectorFieldEvolverConfiguration()
{
    set_maximum_step_size(1);
    set_maximum_enclosure_radius(100.0);
    set_enable_reconditioning(true);
    set_maximum_spacial_error(1e-2);
}


OutputStream&
VectorFieldEvolverConfiguration::_write(OutputStream& os) const
{
    os << "VectorFieldEvolverConfiguration("
       << "\n  maximum_step_size=" << maximum_step_size()
       << ",\n  maximum_enclosure_radius=" << maximum_enclosure_radius()
       << ",\n  enable_reconditioning=" << enable_reconditioning()
       << ",\n  maximum_spacial_error=" << maximum_spacial_error()
       << "\n)";
    return os;
}


}  // namespace Ariadne

