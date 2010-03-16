/***************************************************************************
 *            hybrid_evolver-constrained.cc
 *
 *  Copyright  2009  Pieter Collins
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

#include "numeric.h"
#include "vector.h"
#include "function.h"
#include "taylor_model.h"
#include "taylor_function.h"
#include "grid_set.h"
#include "hybrid_time.h"
#include "hybrid_automaton.h"
#include "hybrid_evolver-constrained.h"

#include "integrator.h"

namespace Ariadne {

/*
# Define extra functionality for TaylorFunctions
def taylor_set(f,j,g):
    """Set the value of the jth variable to g. The new domain is the old domain with the jth variable removed.
       Either g is a constant, or g is a scalar variable with domain equal to the new domain."""
    n=f.domain().size()
    subdomain=Box( [ f.domain()[i] for i in range(0,j)+range(j+1,n) ] )
    if not isinstance(g,ScalarTaylorFunction):
        g=ScalarTaylorFunction.constant(subdomain,g)
    assert subdomain==Box(g.domain())
    h=VectorTaylorFunction([ ScalarTaylorFunction.variable(subdomain,i) for i in range(0,j) ] + [ g ] + [ ScalarTaylorFunction.variable(subdomain,i) for i in range(j,n-1) ])
    return compose(f,h)

def taylor_substitute(f,j,g):
    """Substitute the value of the jth variable with g. 
       Either g is a constant, or g is a scalar variable with domain equal to the domain of f."""
    domain=Box(f.domain())
    if not isinstance(g,ScalarTaylorFunction):
        g = ScalarTaylorFunction.constant(domain,g)
    assert domain==Box(g.domain())
    h=VectorTaylorFunction.identity(domain)
    h[j]=g
    return compose(f,h)
*/



template<class T> std::string str(const T& t) { std::stringstream ss; ss<<t; return ss.str(); }
template<class T1> inline void print(const T1& t1) { std::cerr<<t1; }
template<class T1,class T2> inline void print(const T1& t1, const T2& t2) { std::cerr<<t1<<" "<<t2; }
template<class T1,class T2,class T3> inline void print(const T1& t1, const T2& t2, const T3& t3) { std::cerr<<t1<<" "<<t2<<" "<<t3; }

typedef Vector<Interval> IntervalVector;


struct ConstrainedImageSet::Data {
    IntervalVector _domain;
    Vector<TaylorModel> _models;
    Map<DiscreteEvent,TaylorModel> _positive;
    Map<DiscreteEvent,TaylorModel> _zero;
    Map<DiscreteEvent,TaylorModel> _negative;
    Map<DiscreteEvent,TaylorModel> _always_negative;
    Map<DiscreteEvent,TaylorModel> _always_negative_and_zero;

    VectorFunction function() const {
        return VectorTaylorFunction(this->_domain,this->_models).function(); }

    void _adjoin_outer_approximation_to(GridTreeSet& gts, const IntervalVector& subdomain, uint depth);
};

ConstrainedImageSet::ConstrainedImageSet() 
    : _data(new Data)
{
}

ConstrainedImageSet::ConstrainedImageSet(const ConstrainedImageSet& other) 
    : _data(new Data(other._data.operator*()))
{
}

ConstrainedImageSet& ConstrainedImageSet::operator=(const ConstrainedImageSet& other)
{
    this->_data.operator*() = other._data.operator*();
    return *this;
}

ConstrainedImageSet* ConstrainedImageSet::clone() const
{
    return new ConstrainedImageSet(*this);
}

ConstrainedImageSet::ConstrainedImageSet(IntervalVector domain) 
    : _data(new Data)
{
    _data->_domain=domain;
    _data->_models=VectorTaylorFunction::identity(domain).models();
}


ConstrainedImageSet::ConstrainedImageSet(IntervalVector domain, VectorFunction function) 
    : _data(new Data)
{
    _data->_domain=domain;
    _data->_models=VectorTaylorFunction(domain,function).models();
}


void ConstrainedImageSet::apply_map(VectorFunction map) {
    _data->_models=compose(map,_data->_models);
}

void ConstrainedImageSet::apply_flow(VectorFunction phi, Interval time_domain) {
    _data->_domain=join(_data->_domain,time_domain);
    _data->_models=compose(phi,combine(_data->_models,ScalarTaylorFunction::variable(time_domain).model()));
}

void ConstrainedImageSet::apply_flow(VectorTaylorFunction phi) {
    std::cerr<<"apply_flow(VectorTaylorFunction)\n";
    Interval time_domain=phi.domain()[phi.domain().size()-1];
    _data->_domain=join(_data->_domain,time_domain);
    _data->_models=unchecked_compose(phi.models(),phi.domain(),combine(_data->_models,ScalarTaylorFunction::variable(time_domain).model()));
}

void ConstrainedImageSet::new_invariant(DiscreteEvent event, ScalarFunction constraint, ScalarFunction derivative) {
    _data->_always_negative[event]=compose(constraint,_data->_models);
}

void ConstrainedImageSet::new_activation(DiscreteEvent event,ScalarFunction constraint,ScalarFunction derivative) {
    _data->_positive[event]=compose(constraint,_data->_models);
}

void ConstrainedImageSet::new_guard(DiscreteEvent event,ScalarFunction constraint,ScalarFunction derivative) {
    _data->_always_negative[event]=compose(constraint,_data->_models);
    _data->_positive[event]=compose(constraint,_data->_models);
}

void ConstrainedImageSet::new_equation(DiscreteEvent event,ScalarFunction constraint) {
    _data->_zero[event]=compose(constraint,_data->_models);
}

uint ConstrainedImageSet::dimension() const {
    return _data->_models.size();
}

uint ConstrainedImageSet::number_of_parameters() const {
    return _data->_domain.size();
}

Box ConstrainedImageSet::bounding_box() const {
    return Box(evaluate(_data->_models,_data->_domain)).bounding_box();
}

tribool ConstrainedImageSet::disjoint(Box bx) {
    ARIADNE_NOT_IMPLEMENTED;
}


GridTreeSet ConstrainedImageSet::outer_approximation(Grid grid, int depth) const {
    GridTreeSet gts(grid);
    this->_data->_adjoin_outer_approximation_to(gts,_data->_domain,depth);
    return gts;
}


void ConstrainedImageSet::draw(CanvasInterface& canvas) const {
    return this->outer_approximation(Grid(this->dimension()),4).draw(canvas);
};

std::ostream& ConstrainedImageSet::write(std::ostream& os) const {
    os << "\n  domain=" << _data->_domain << ","
       << "\n  function=" << _data->function() << ","
       << "\n  activations=" << "...,";
    return os;
}

void ConstrainedImageSet::Data::
_adjoin_outer_approximation_to(GridTreeSet& gts, const IntervalVector& subdomain, uint depth)
{
   typedef Map<DiscreteEvent,TaylorModel>::const_iterator iterator;
    double max_error=1.0/(1<<depth);
    const uint ds=this->_domain.size();
    for(iterator iter=this->_positive.begin(); iter!=this->_positive.end(); ++iter) {
        Interval constraint_range=iter->second.evaluate(subdomain);
        if(constraint_range.upper() < 0.0) { 
            return;
        }
    }
    for(iterator iter=this->_always_negative_and_zero.begin(); iter!=this->_always_negative_and_zero.end(); ++iter) {
        const TaylorModel& constraint=iter->second;
        Interval constraint_range=constraint.evaluate(subdomain);
        if(constraint_range.lower() > 0.0 or constraint_range.upper()<0.0) {
            return;
        }
        Vector<Interval> lowdomain=subdomain;
        while(lowdomain[ds-1].lower() >= this->_domain[ds-1].lower()) {
            if(constraint.evaluate(lowdomain).lower() > 0.0) {
                return;
            }
            lowdomain[ds-1]=lowdomain[ds-1]-lowdomain[ds-1].width();
        }
    }
    for(iterator iter=this->_always_negative.begin(); iter!=this->_always_negative.end(); ++iter) {
        const TaylorModel& constraint=iter->second;
        Vector<Interval> lowdomain=subdomain;
        while(lowdomain[ds-1].lower()>this->_domain[ds-1].lower()) {
            if(constraint.evaluate(lowdomain).lower() > 0.0) {
                return;
            }
            lowdomain[ds-1]=lowdomain[ds-1]-lowdomain[ds-1].width();
        }
    }
    Box range=evaluate(_models,subdomain);
    if(range.radius()<max_error) {
        gts.adjoin_outer_approximation(range,depth+2);
    } else {
        Vector<Interval> subdomain1,subdomain2;
        make_lpair(subdomain1,subdomain2)=split(subdomain);
        this->_adjoin_outer_approximation_to(gts,subdomain1,depth);
        this->_adjoin_outer_approximation_to(gts,subdomain2,depth);
    }
}





Orbit<HybridConstrainedImageSet>
ConstrainedImageSetHybridEvolver::
orbit(const HybridAutomaton& system,
      const HybridConstrainedImageSet& initial,
      const HybridTime& time,
      Semantics semantics) const
{
    EnclosureListType final, reachable, intermediate;
    this->_evolution(final,reachable,intermediate,system,initial,time,semantics,true);
    Orbit<HybridConstrainedImageSet> orbit(initial);
    orbit.adjoin_intermediate(intermediate);
    orbit.adjoin_reach(reachable);
    orbit.adjoin_final(final);
    return orbit;
}


void
ConstrainedImageSetHybridEvolver::
_evolution(EnclosureListType& final,
           EnclosureListType& reachable,
           EnclosureListType& intermediate,
           HybridAutomaton const& system,
           HybridConstrainedImageSet const& initial_set,
           HybridTime const& maximum_time,
           Semantics semantics,
           bool reach) const
{
    List<TimedHybridConstrainedImageSet> working;

    working.push_back(TimedHybridConstrainedImageSet(initial_set));

    while(!working.empty()) {
        this->_upper_evolution_step(working,final,reachable,intermediate,system,maximum_time);
    }
}

void
ConstrainedImageSetHybridEvolver::
_upper_evolution_step(List<TimedHybridConstrainedImageSet>& working_sets,
                      ListSet<HybridConstrainedImageSet>& reach_sets,
                      ListSet<HybridConstrainedImageSet>& evolve_sets,
                      ListSet<HybridConstrainedImageSet>& intermediate_sets,
                      HybridAutomaton const& system,
                      HybridTime const& maximum_hybrid_time) const
{
    typedef Map<DiscreteEvent,ScalarFunction>::const_iterator constraint_iterator;

    TaylorIntegrator integrator(5);
    Float maximum_time=maximum_hybrid_time.continuous_time();

    std::cerr << "UPPER EVOLUTION STEP";

    TimedHybridConstrainedImageSet current_data=working_sets.back();
    working_sets.pop_back();
    
    List<DiscreteEvent>& starting_events=current_data._events;
    DiscreteState starting_location=current_data._location;
    ConstrainedImageSet& starting_set=current_data._set;
    ScalarFunction& starting_time=current_data._time;

    print("\nstarting_events:",starting_events);
    print("\nstarting_location:",starting_location);
    print("\nstarting_set:",starting_set);

    // Set the dimension
    const uint n=starting_set.dimension();
    const uint m=starting_set.number_of_parameters();
    
    // Extract mode and transitions, dynamic and constraints
    DiscreteMode mode=system.mode(starting_location);
    std::list<DiscreteTransition> transitions=system.transitions(starting_location);

    VectorFunction dynamic=mode.dynamic();
    Map<DiscreteEvent,ScalarFunction> invariants(mode.scalar_invariants());

    Map<DiscreteEvent,ScalarFunction> guards;
    Map<DiscreteEvent,ScalarFunction> activations;
    Map<DiscreteEvent,VectorFunction> resets;
    Map<DiscreteEvent,DiscreteState> targets;
    
    for(std::list<DiscreteTransition>::const_iterator transition_iter=transitions.begin();
        transition_iter!=transitions.end(); ++transition_iter)
    {
        const DiscreteEvent& event=transition_iter->event();
        if(transition_iter->forced()) {
            guards[event]=transition_iter->scalar_activation();
        } else {
            activations[event]=transition_iter->scalar_activation();
        }
        resets[event]=transition_iter->reset();
        targets[event]=transition_iter->target();
    }

    print("\ninvariants:",invariants);
    print("\nguards:",guards);
    print("\nactivations:",activations);
    print("\nresets:",resets);
    print("\n\n");

    Map<DiscreteEvent,ScalarFunction> blocking=join(invariants,guards);
    Map<DiscreteEvent,ScalarFunction> enabling=join(activations,guards);
    
    // Compute initially active guards
    Box starting_bounding_box=Box(starting_set.bounding_box());
    Map<DiscreteEvent,Interval> starting_ranges;
    bool no_flow=false;
    for(constraint_iterator iter=blocking.begin(); iter!=blocking.end(); ++iter) {
        const DiscreteEvent& event=iter->first;
        const ScalarFunction& constraint=iter->second;
        Interval constraint_range=constraint(starting_bounding_box);
        print(str(event)+".starting_range():",constraint_range,"\n");
        starting_ranges[event]=constraint(starting_bounding_box);
        if(constraint_range.lower()>0.0) {
            no_flow=true;
        }
    }
    assert(no_flow==false);

    // Compute flow and actual time step size used
    VectorTaylorFunction flow_model=integrator.flow(dynamic,IntervalVector(starting_set.bounding_box()),maximum_time);
    print("flow_model:",flow_model);
    Box flow_domain=Box(flow_model.domain());
    Float step_size=flow_domain[n].upper();
    flow_domain[n]=Interval(0,step_size);
    flow_model=restrict(flow_model,flow_domain);
    print("flow_model:",flow_model);
    print("flow_model.domain():",flow_model.domain());
    print("flow_model.range():",flow_model.range());
    VectorFunction flow=flow_model.function();
    print("flow:",flow);
    ConstrainedImageSet flowed_set=starting_set;
    flowed_set.apply_flow(flow,flow_domain[n]);
    print("\nflowed_set:",flowed_set);

    // Compute restrictions on continuous evolution
    ConstrainedImageSet reached_set=flowed_set;
    for(constraint_iterator iter=blocking.begin(); iter!=blocking.end(); ++iter) {
        DiscreteEvent const& event=iter->first;
        ScalarFunction const& constraint=iter->second;
        ScalarFunction derivative=lie_derivative(constraint,dynamic);
        reached_set.new_invariant(event,constraint,derivative);
    }
    print(reached_set);
    reach_sets.adjoin(HybridConstrainedImageSet(starting_events,starting_location,reached_set));

    // Compute the set reached after one time step of the flow; corresponds to setting t=t0+delta
    // By "finishing" we mean completing one full time step without a jump; this corresponds to delta=h.
    // By "final", we mean completing the full evolution; this corresponds to t=t_max.
    print("\nstarting_time:",starting_time);
    ScalarFunction starting_time_function=embed(starting_time,1u);
    print("starting_time_function:",starting_time_function);
    ScalarFunction dwell_time_function=ScalarFunction::variable(n+1u,n);
    print("dwell_time_function:",dwell_time_function);
    ScalarFunction evolution_time_function=starting_time_function+dwell_time_function;
    print("evolution_time_function:",evolution_time_function);
    ScalarFunction finishing_constraint=dwell_time_function-step_size;
    print("finishing_constraint:",finishing_constraint);
    ScalarFunction final_constraint=evolution_time_function-maximum_time;
    print("final_constraint:",final_constraint);

    ConstrainedImageSet finishing_set(reached_set);
    finishing_set.new_equation(DiscreteEvent(0),finishing_constraint);
    ConstrainedImageSet final_set(reached_set);
    final_set.new_equation(DiscreteEvent(0),final_constraint);
    print("\nreached_set:",reached_set);
    print("\nfinishing_set:",finishing_set);
    print("\nfinal_set:",final_set);
    evolve_sets.adjoin(HybridConstrainedImageSet(starting_events,starting_location,final_set));
    working_sets.append(TimedHybridConstrainedImageSet(starting_events,starting_location,finishing_set,evolution_time_function));
    
    // Compute the reached set under a single event
    for(constraint_iterator iter=enabling.begin(); iter!=enabling.end(); ++iter) {
        DiscreteEvent const& event=iter->first;
        ScalarFunction const& constraint=iter->second;
        ScalarFunction derivative=lie_derivative(constraint,dynamic);
        DiscreteState target=targets[event];
        VectorFunction reset=resets[event];
        ConstrainedImageSet jumping_set(reached_set);
        jumping_set.new_activation(event,constraint,derivative);
        ConstrainedImageSet jumped_set(jumping_set);
        jumped_set.apply_map(reset);
        print("\njumped_set("+str(event)+"):",jumped_set);
        ScalarFunction dwell_time_function=ScalarFunction::variable(m+1u,m);
        ScalarFunction maximum_time_function=embed(starting_time,1u)+dwell_time_function;
        List<DiscreteEvent> jumped_events=starting_events;
        jumped_events.append(event);
        working_sets.append(TimedHybridConstrainedImageSet(jumped_events,target,jumped_set,maximum_time_function));
    }
            
    print("\nDONE STEP\n");

    return;
}

} // namespace Ariadne
