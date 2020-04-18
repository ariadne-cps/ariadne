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
#include "../dynamics/enclosure.hpp"
#include "../dynamics/orbit.hpp"

#include "../solvers/integrator.hpp"

#include "../output/logging.hpp"

#include "../dynamics/vector_field.hpp"
#include "../dynamics/vector_field_evolver.hpp"

#include "../symbolic/space.hpp"
#include "../symbolic/assignment.hpp"

namespace {

using namespace Ariadne;

template<class ES> List<ES> subdivide(const ES& enclosure) {
    List<ES> result;
    Pair<ES,ES> split=enclosure.split();
    result.append(split.first);
    result.append(split.second);
    return result;
}


}



namespace Ariadne {

VectorField::VectorField(List<DottedRealAssignment> const& dynamics)
    : _variable_names(variable_names(left_hand_sides(dynamics)))
    , _function(make_function(left_hand_sides(dynamics),Vector<RealExpression>(right_hand_sides(dynamics))))
{
}

VectorField::VectorField(EffectiveVectorMultivariateFunction const& function) {
    List<Identifier> variable_names;
    for (auto i : range(0,function.result_size()))
        variable_names.append(Identifier("x"+std::to_string(i)));

    _variable_names = variable_names;
    _function = function;
}

RealSpace VectorField::state_space() const
{
    return real_space(this->_variable_names);
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
    this->charcode = "v";
}

typename VectorFieldEvolver::EnclosureType VectorFieldEvolver::enclosure(const ExactBoxType& box) const {
    return Enclosure(box,this->function_factory());
}

typename VectorFieldEvolver::EnclosureType VectorFieldEvolver::enclosure(const RealBox& box) const {
    return Enclosure(box,this->function_factory());
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
    ARIADNE_LOG(5,ARIADNE_PRETTY_FUNCTION<<"\n");

    List< TimedEnclosureType > working_sets;

    {
        // Set up initial timed set models
        ARIADNE_LOG(6,"initial_set = "<<initial_set<<"\n");
        TimeStepType initial_time = 0u;
        ARIADNE_LOG(6,"initial_time = "<<initial_time<<"\n");
        EnclosureType initial_set_model(initial_set);
        ARIADNE_LOG(6,"initial_set_model = "<<initial_set_model<<"\n");
        working_sets.push_back(make_pair(initial_time,initial_set_model));
    }


    while(!working_sets.empty()) {
        TimedEnclosureType current_timed_set=working_sets.back();
        working_sets.pop_back();
        TimeStepType current_time=current_timed_set.first;
        EnclosureType current_set_model=current_timed_set.second;
        FloatDPUpperBound current_set_radius=current_set_model.bounding_box().radius();
        if(definitely(current_time>=maximum_time)) {
            final_sets.adjoin(current_set_model);
        } else if(semantics == Semantics::UPPER && ENABLE_SUBDIVISIONS
                  && decide(current_set_radius>this->_configuration->maximum_enclosure_radius())) {
            // Subdivide
            List< EnclosureType > subdivisions=subdivide(current_set_model);
            for(Nat i=0; i!=subdivisions.size(); ++i) {
                EnclosureType const& subdivided_set_model=subdivisions[i];
                working_sets.push_back(make_pair(current_time,subdivided_set_model));
            }
        } else if(semantics == Semantics::LOWER && ENABLE_PREMATURE_TERMINATION && decide(current_set_radius>this->_configuration->maximum_enclosure_radius())) {
            ARIADNE_WARN("Terminating lower evolution at time " << current_time
                         << " and set " << current_set_model << " due to maximum radius being exceeded.");
        } else {
            // Compute evolution
            this->_evolution_step(working_sets,
                                  final_sets,reach_sets,intermediate_sets,
                                  current_timed_set,maximum_time,
                                  semantics,reach);
        }

        ARIADNE_LOG(1,"#w="<<std::setw(4)<<working_sets.size()
                        <<"#r="<<std::setw(4)<<std::left<<reach_sets.size()
                        <<" t="<<std::setw(7)<<std::fixed<<current_time.get_d()
                        <<" p="<<std::setw(4)<<std::left<<current_set_model.number_of_parameters()
                        <<" r="<<std::setw(7)<<current_set_model.radius()
                        <<" c="<<current_set_model.centre() << "\n");

    }
}





Void
VectorFieldEvolver::
_evolution_step(List< TimedEnclosureType >& working_sets,
                EnclosureListType& final_sets,
                EnclosureListType& reach_sets,
                EnclosureListType& intermediate_sets,
                const TimedEnclosureType& working_timed_set_model,
                const TimeType& maximum_time,
                Semantics semantics,
                Bool reach) const
{
    typedef EffectiveVectorMultivariateFunction FunctionType;

    EnclosureType current_set_model;
    TimeStepType current_time;
    ARIADNE_LOG(9,"working_timed_set_model = "<<working_timed_set_model<<"\n");
    make_lpair(current_time,current_set_model)=working_timed_set_model;

    ARIADNE_LOG(4,"current_time = "<<current_time<<"");
    ARIADNE_LOG(6,"current_set_model = "<<current_set_model<<"\n");

    ARIADNE_LOG(2,"box = "<<current_set_model.bounding_box()<<" ");
    ARIADNE_LOG(2,"radius = "<<current_set_model.bounding_box().radius()<<"\n\n");
    //const Nat nd=initial_set_model.result_size();
    //const Nat ng=initial_set_model.argument_size();


    // Test to see if set requires reconditioning
    if(this->_configuration->enable_reconditioning() &&
       possibly(norm(current_set_model.state_function().errors()) > this->_configuration->maximum_spacial_error())) {

        ARIADNE_LOG(4," reconditioning: errors "<<current_set_model.state_function().errors()<<"\n");
        current_set_model.recondition();
        working_sets.append(make_pair(current_time,current_set_model));
        return;
    }

    /////////////// Main Evolution ////////////////////////////////
    const FunctionType& dynamic=_sys_ptr->function();

    // Set evolution parameters
    const StepSizeType maximum_step_size=this->_configuration->maximum_step_size();
    //const FloatDP maximum_bounds_diameter=this->_parameters->maximum_enclosure_radius*2;
    //const FloatDP zero_time=0.0;

    // Get bounding boxes for time and space bounding_box
    ExactBoxType current_set_bounds=cast_exact_box(current_set_model.bounding_box());
    ARIADNE_LOG(4,"current_set_bounds = "<<current_set_bounds<<"\n");


    // Compute flow model
    // TODO: Modify this for general integrator interface
    //TaylorPicardIntegrator const* taylor_integrator=dynamic_cast<const TaylorPicardIntegrator*>(this->_integrator.operator->());
    IntegratorInterface const* integrator=this->_integrator.operator->();
    StepSizeType step_size=maximum_step_size;
    FlowStepModelType flow_model=integrator->flow_step(dynamic,current_set_bounds,step_size);
    ARIADNE_LOG(4,"step_size = "<<step_size<<"\n");
    ARIADNE_LOG(6,"flow_model = "<<flow_model<<"\n");
    FlowStepModelType flow_step_model=partial_evaluate(flow_model,flow_model.domain().size()-1u,step_size);
    ARIADNE_LOG(6,"flow_step_model = "<<flow_step_model<<"\n");

    // Compute the integration time model
    TimeStepType next_time=current_time+TimeStepType(step_size);
    ARIADNE_LOG(6,"next_time = "<<next_time<<"\n");
    // Compute the flow tube (reachable set) model and the final set
    ARIADNE_LOG(6,"product = "<<product(current_set_model,ExactIntervalType(0,step_size))<<"\n");
    EnclosureType reach_set_model=apply(flow_model,product(current_set_model,ExactIntervalType(0,step_size)));
    ARIADNE_LOG(6,"reach_set_model = "<<reach_set_model<<"\n");
    EnclosureType next_set_model=apply(flow_step_model,current_set_model);
    ARIADNE_LOG(6,"next_set_model = "<<next_set_model<<"\n");
    ARIADNE_LOG(4,"Done computing evolution\n");

    reach_sets.adjoin(reach_set_model);

    intermediate_sets.adjoin(EnclosureType(next_set_model));
    working_sets.push_back(make_pair(next_time,next_set_model));

    ARIADNE_LOG(2,"Done evolution_step.\n\n");

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
    os << "VectorFieldEvolverConfiguration"
       << ",\n  maximum_step_size=" << maximum_step_size()
       << ",\n  maximum_enclosure_radius=" << maximum_enclosure_radius()
       << ",\n  enable_reconditioning=" << enable_reconditioning()
       << ",\n  maximum_spacial_error=" << maximum_spacial_error()
       << "\n)\n";
    return os;
}


}  // namespace Ariadne

