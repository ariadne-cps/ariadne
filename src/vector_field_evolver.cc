/***************************************************************************
 *            vector_field_evolver.cc
 *
 *  Copyright  2008  Alberto Casagrande, Pieter Collins
 *
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include "macros.h"
#include "array.h"
#include "tuple.h"
#include "stlio.h"
#include "vector.h"
#include "function.h"
#include "taylor_function.h"
#include "taylor_set.h"
#include "orbit.h"
#include "evolution_parameters.h"

#include "integrator.h"

#include "logging.h"

#include "vector_field.h"
#include "vector_field_evolver.h"

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

// Allow subdivisions in upper evolution
const bool ENABLE_SUBDIVISIONS = false;
// Allow premature termination of lower evolution
const bool ENABLE_PREMATURE_TERMINATION = false;

static const int BLOCKING_EVENT = -2;
using boost::shared_ptr;

class DegenerateCrossingException { };



VectorFieldEvolver::VectorFieldEvolver(const EvolutionParametersType& p, const IntegratorInterface& i)
    : _parameters(new EvolutionParametersType(p))
    , _integrator(i.clone())
{
    ARIADNE_ASSERT_MSG(dynamic_cast<const TaylorPicardIntegrator*>(&i),"Only TaylorPicardIntegrator supported by VectorFieldEvolver\n");
}


Orbit<VectorFieldEvolver::EnclosureType>
VectorFieldEvolver::
orbit(const SystemType& system,
      const EnclosureType& initial_set,
      const TimeType& time,
      Semantics semantics) const
{
    Orbit<EnclosureType> orbit(initial_set);
    EnclosureListType final;
    EnclosureListType reachable;
    EnclosureListType intermediate;
    this->_evolution(final,reachable,intermediate,
                     system,initial_set,time,semantics,false);
    orbit.adjoin_intermediate(intermediate);
    orbit.adjoin_reach(reachable);
    orbit.adjoin_final(final);
    return orbit;
}




enum PredicateKind { INVARIANT, ACTIVATION, GUARD, TIME, MIXED };
enum CrossingKind { TRANSVERSE, TOUCHING, NONE, UNKNOWN };



void
VectorFieldEvolver::
_evolution(EnclosureListType& final_sets,
           EnclosureListType& reach_sets,
           EnclosureListType& intermediate_sets,
           const SystemType& system,
           const EnclosureType& initial_set,
           const TimeType& maximum_time,
           Semantics semantics,
           bool reach) const
{
    typedef RealVectorFunction FunctionType;
    typedef Vector<Interval> BoxType;
    typedef VectorTaylorFunction FunctionModelType;
    typedef VectorTaylorFunction FlowModelType;

    ARIADNE_LOG(5,ARIADNE_PRETTY_FUNCTION<<"\n");

    List< TimedEnclosureType > working_sets;

    {
        // Set up initial timed set models
        ARIADNE_LOG(6,"initial_set = "<<initial_set<<"\n");
        TimeType initial_time = 0.0;
        ARIADNE_LOG(6,"initial_time = "<<initial_time<<"\n");
        EnclosureType initial_set_model(initial_set);
        ARIADNE_LOG(6,"initial_set_model = "<<initial_set_model<<"\n");
        working_sets.push_back(make_pair(initial_time,initial_set_model));
    }


    while(!working_sets.empty()) {
        TimedEnclosureType current_timed_set=working_sets.back();
        working_sets.pop_back();
        TimeType current_time=current_timed_set.first;
        EnclosureType current_set_model=current_timed_set.second;
        Float current_set_radius=radius(current_set_model.bounding_box());
        if(current_time>=maximum_time) {
            final_sets.adjoin(current_set_model);
        } else if(UPPER_SEMANTICS && ENABLE_SUBDIVISIONS
                  && (current_set_radius>this->_parameters->maximum_enclosure_radius)) {
            // Subdivide
            List< EnclosureType > subdivisions=subdivide(current_set_model);
            for(uint i=0; i!=subdivisions.size(); ++i) {
                EnclosureType const& subdivided_set_model=subdivisions[i];
                working_sets.push_back(make_pair(current_time,subdivided_set_model));
            }
        } else if(LOWER_SEMANTICS && ENABLE_PREMATURE_TERMINATION && current_set_radius>this->_parameters->maximum_enclosure_radius) {
            ARIADNE_WARN("Terminating lower evolution at time " << current_time
                         << " and set " << current_set_model << " due to maximum radius being exceeded.");
        } else {
            // Compute evolution
            this->_evolution_step(working_sets,
                                  final_sets,reach_sets,intermediate_sets,
                                  system,current_timed_set,maximum_time,
                                  semantics,reach);
        }

        if(verbosity==1) {
            ARIADNE_LOG(1,"\r"
                        <<"#w="<<std::setw(4)<<working_sets.size()
                        <<"#r="<<std::setw(4)<<std::left<<reach_sets.size()
                        <<" t="<<std::setw(7)<<std::fixed<<current_time
                        <<" r="<<std::setw(7)<<current_set_model.radius()
                        <<" c="<<current_set_model.centre()
                        <<"                      ");
        }

    }

}





void
VectorFieldEvolver::
_evolution_step(List< TimedEnclosureType >& working_sets,
                EnclosureListType& final_sets,
                EnclosureListType& reach_sets,
                EnclosureListType& intermediate_sets,
                const SystemType& system,
                const TimedEnclosureType& working_timed_set_model,
                const TimeType& maximum_time,
                Semantics semantics,
                bool reach) const
{
    typedef RealVectorFunction FunctionType;
    typedef Vector<Interval> BoxType;
    typedef VectorTaylorFunction MapModelType;
    typedef VectorTaylorFunction FlowModelType;
    typedef TaylorConstrainedImageSet EnclosureType;

    EnclosureType current_set_model;
    TimeType current_time;
    ARIADNE_LOG(9,"working_timed_set_model = "<<working_timed_set_model<<"\n");
    make_lpair(current_time,current_set_model)=working_timed_set_model;

    ARIADNE_LOG(4,"current_time = "<<current_time<<"");
    ARIADNE_LOG(6,"current_set_model = "<<current_set_model<<"\n");

    ARIADNE_LOG(2,"box = "<<current_set_model.bounding_box()<<" ");
    ARIADNE_LOG(2,"radius = "<<radius(current_set_model.bounding_box())<<"\n\n");
    //const uint nd=initial_set_model.result_size();
    //const uint ng=initial_set_model.argument_size();


    /////////////// Main Evolution ////////////////////////////////
    const FunctionType& dynamic=system.function();

    // Set evolution parameters
    const Float maximum_step_size=this->_parameters->maximum_step_size;
    const Float maximum_bounds_diameter=this->_parameters->maximum_enclosure_radius*2;
    const Float zero_time=0.0;

    // Get bounding boxes for time and space bounding_box
    Vector<Interval> current_set_bounds=current_set_model.bounding_box();
    ARIADNE_LOG(4,"current_set_bounds = "<<current_set_bounds<<"\n");

    //ARIADNE_ASSERT(initial_time_bounding_box.width() <= maximum_step_size);


    // Compute flow bounds and find flow bounding box
    Vector<Interval> flow_bounds;
    Float step_size;
    make_lpair(step_size,flow_bounds)=this->_integrator->flow_bounds(dynamic,current_set_bounds,maximum_step_size);
    ARIADNE_LOG(4,"step_size = "<<step_size<<"\n");
    ARIADNE_LOG(4,"flow_bounds = "<<flow_bounds<<"\n");

    // Compute flow model
    // TODO: Modify this for general integrator interface
    TaylorPicardIntegrator const* taylor_integrator=dynamic_cast<const TaylorPicardIntegrator*>(this->_integrator.operator->());
    FlowModelType flow_model=taylor_integrator->flow_step(dynamic,current_set_bounds,step_size,flow_bounds);
    //FlowModelType flow_model=this->_integrator->flow_step(dynamic,current_set_bounds,step_size,flow_bounds);
    ARIADNE_LOG(6,"flow_model = "<<flow_model<<"\n");
    FlowModelType flow_step_model=partial_evaluate(flow_model,flow_model.domain().size()-1u,step_size);
    ARIADNE_LOG(6,"flow_step_model = "<<flow_step_model<<"\n");

    // Compute the integration time model
    TimeType next_time=current_time+step_size;
    ARIADNE_LOG(6,"next_time = "<<next_time<<"\n");
    // Compute the flow tube (reachable set) model and the final set
    //std::cerr<<"flow_model.argument_size()="<<flow_model.argument_size()<<"\n";
    //std::cerr<<"product(current_set_model,Interval(0,step_size))="<<product(current_set_model,Interval(0,step_size))<<"\n";
    EnclosureType reach_set_model=apply(flow_model,product(current_set_model,Interval(0.0,step_size)));
    ARIADNE_LOG(6,"reach_set_model = "<<reach_set_model<<"\n");
    EnclosureType next_set_model=apply(flow_step_model,current_set_model);
    ARIADNE_LOG(6,"next_set_model = "<<next_set_model<<"\n");
    ARIADNE_LOG(4,"Done computing evolution\n");

    reach_sets.adjoin(reach_set_model);

    intermediate_sets.adjoin(EnclosureType(next_set_model));
    working_sets.push_back(make_pair(next_time,next_set_model));





    ARIADNE_LOG(2,"Done evolution_step.\n\n");

}




}  // namespace Ariadne

