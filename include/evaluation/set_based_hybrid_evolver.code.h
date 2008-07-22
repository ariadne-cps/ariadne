/***************************************************************************
 *            set_based_hybrid_evolver.code.h
 *
 *  Copyright  2004-7  Alberto Casagrande,  Pieter Collins
 *  casagrande@dimi.uniud.it  Pieter.Collins@cwi.nl
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
 
#include "base/tuple.h"
#include "geometry/box_list_set.h"
#include "geometry/set_interface.h"
#include "geometry/hybrid_set.h"
#include "system/hybrid_automaton.h"
#include "evaluation/hybrid_time.h"
#include "evaluation/evolution_parameters.h"

#include "evaluation/applicator_interface.h"
#include "evaluation/integrator_interface.h"
#include "evaluation/satisfier_interface.h"

#include "evaluation/standard_applicator.h"
#include "evaluation/standard_integrator.h"
#include "evaluation/standard_approximator.h"
#include "evaluation/standard_satisfier.h"
#include "evaluation/standard_subdivider.h"
#include "evaluation/cascade_reducer.h"

#include "output/epsstream.h"
#include "output/logging.h"

#include "set_based_hybrid_evolver.h"


namespace Ariadne {
 
typedef Rational Q;
const uint BISECTION_STEPS=8;




template<class R> 
Evolver< HybridAutomaton<R>, Zonotope<R> >::
Evolver(const EvolutionParameters<R>& parameters, 
        const ApplicatorInterface<ES>& applicator, 
        const IntegratorInterface<ES>& integrator, 
        const SatisfierInterface<ES>& satisfier, 
        const SubdividerInterface<ES>& subdivider, 
        const ReducerInterface<ES>& reducer)
  : _parameters(parameters.clone()),
    _applicator(applicator.clone()),
    _integrator(integrator.clone()),
    _satisfier(satisfier.clone()),
    _subdivider(subdivider.clone()),
    _reducer(reducer.clone()),
    verbosity(parameters.verbosity())
{
}

template<class R> 
Evolver< HybridAutomaton<R>, Zonotope<R> >::
Evolver(const EvolutionParameters<R>& parameters, 
        const ApplicatorInterface<ES>& applicator, 
        const IntegratorInterface<ES>& integrator)
{
  StandardSatisfier< ES > satisfier;
  StandardSubdivider< ES > subdivider;
  CascadeReducer< ES > reducer(3);
  *this = Evolver(parameters,applicator,integrator,satisfier,subdivider,reducer);
}

template<class R> 
Evolver< HybridAutomaton<R>, Zonotope<R> >::
Evolver(const EvolutionParameters<R>& parameters)
{
  StandardApplicator< ES > applicator;
  StandardIntegrator< ES > integrator;
  StandardSatisfier< ES > satisfier;
  StandardSubdivider< ES > subdivider;
  CascadeReducer< ES > reducer(3);
  *this = Evolver(parameters,applicator,integrator,satisfier,subdivider,reducer);
}

template<class R> 
Evolver< HybridAutomaton<R>, Zonotope<R> >::
Evolver()
{
  EvolutionParameters<R> parameters;
  StandardApplicator< ES > applicator;
  StandardIntegrator< ES > integrator;
  StandardSatisfier< ES > satisfier;
  StandardSubdivider< ES > subdivider;
  CascadeReducer< ES > reducer(3);
  *this = Evolver(parameters,applicator,integrator,satisfier,subdivider,reducer);
}


template<class R>
bool
Evolver< HybridAutomaton<R>, Zonotope<R> >::
_satisfies(const ES& bs,
           const CS& inv, 
           const Semantics semantics) const
{
  ARIADNE_LOG(8,"  SetBasedHybridEvolver::_satisfies(...)\n");
  ARIADNE_LOG(9,"    bs="<<bs<<"\n    inv="<<inv<<"\n");
	ARIADNE_LOG(9,"    bs.bounding_box="<<bs.bounding_box()<<"\n");
  ARIADNE_LOG(9,(semantics==upper_semantics?"    upper_semantics\n":"    lower_semantics\n"));
  bool res;
  if(semantics==upper_semantics) {
    res = not definitely(disjoint(bs,inv));
		ARIADNE_LOG(9,"    result="<<res<<"\n");
    return res;
  } else { 
	  res = definitely(subset(bs,inv));
		ARIADNE_LOG(9,"    result="<<res<<"\n");
    return res;
  }
}



template<class R>
Rational
Evolver< HybridAutomaton<R>, Zonotope<R> >::
_initial_activation_time(const VF& vf, 
                         const CS& inv, 
                         const ES& bs,
                         const Q& maximum_time,
                         const Bx& bb,
                         const Semantics semantics) const 
{
  ARIADNE_LOG(8,"  SetBasedHybridEvolver<BS>::_initial_activation_time\n");
  // Compute the set of times for which inv is satisfied.
  // Assume that the set is connected and nonempty
  Q initial_time=0;
  Q initial_activation_time=0;
	Q minimum_step_size=this->minimum_step_size();

  // Compute initial activation time
  if(!_satisfies(bs,inv,semantics)) {
    Q lower_bound=initial_time;
    Q upper_bound=maximum_time;
    for(uint i=0; i!=BISECTION_STEPS && upper_bound-lower_bound > minimum_step_size; ++i) {
      Q approx=med(lower_bound,upper_bound);
      ES rs=continuous_evolution_step(vf,bs,lower_bound,approx,bb);
      if(_satisfies(rs,inv,semantics)) {
        upper_bound=approx;
      } else {
        lower_bound=approx;
      }
    }
    if(semantics==upper_semantics) {
      initial_activation_time=lower_bound;
    } else {
      initial_activation_time=upper_bound;
    }
  }
  return initial_activation_time;
}



template<class R>
Q
Evolver< HybridAutomaton<R>, Zonotope<R> >::
_final_activation_time(const VF& vf, const CS& inv, const ES& bs, 
                       const Q& maximum_time, const Bx& bb, const Semantics semantics) const 
{
  ARIADNE_LOG(8,"  SetBasedHybridEvolver<ES>::_final_activation_time\n");
  // Compute final activation time
  ES es=continuous_integration_step(vf,bs,maximum_time,bb);
  ARIADNE_LOG(9,"    es="<<es<<"\n");
  if(_satisfies(es,inv,semantics)) {
    return maximum_time;
  }

  Q minimum_step_size=this->minimum_step_size();
  Q final_activation_time=maximum_time;
  Q lower_bound=0;
  Q upper_bound=maximum_time;
  for(uint i=0; i!=BISECTION_STEPS && upper_bound-lower_bound > minimum_step_size; ++i) {
    Q approx=med(lower_bound,upper_bound);
		ARIADNE_LOG(9,"   lower_bound="<<lower_bound<<" upper_bound="<<upper_bound<<"\n");
    ES rs=continuous_evolution_step(vf,bs,approx,upper_bound,bb);
    if(_satisfies(rs,inv,semantics)) {
      lower_bound=approx;
    } else {
      upper_bound=approx;
    }
  }
  if(semantics==upper_semantics) {
    final_activation_time=upper_bound;
  } else {
    final_activation_time=lower_bound;
  }
  return final_activation_time;
}



template<class R>
tuple<Rational,Zonotope<R> >
Evolver< HybridAutomaton<R>, Zonotope<R> >::
_saltation_map(const VF& vf1, const VF& vf2, const Mp& rm, const CS& inv, 
               const ES& bs, const Q& h1, const Q& h2, const Bx& bb, const Semantics sem) const
{
  ARIADNE_LOG(8,"\nSetBasedHybridEvolver<BS>::_saltation_map\n");
  ARIADNE_LOG(9,(sem==upper_semantics?"  upper_semantics\n":"  lower_semantics\n"));
  ARIADNE_ASSERT(h1<=h2);
  if(sem==upper_semantics) {
    ES rs=continuous_evolution_step(vf1,bs,h1,h2,bb);
    ARIADNE_LOG(9,"  rs="<< rs << "\n");    
		return make_tuple(h1,apply(rm,rs));  } 
  else {
    Q hmed=med(h1,h2);
    ES rs=continuous_integration_step(vf1,bs,hmed,bb);
    ARIADNE_LOG(9,"  rs="<< rs << "\n");    
    return make_tuple(hmed,apply(rm,rs));  
  }
}



template<class R>
void
Evolver< HybridAutomaton<R>, Zonotope<R> >::
_step(HESL& evolve,
      HESL& reach,
      HESL& intermediate,
      THESL& working,
      const HA& automaton, 
      const Q& time,
      const Semantics semantics) const
{
  //const_cast<SetBasedHybridEvolver<ES>*>(this)->verbosity=9;
  ARIADNE_LOG(6,"\nSetBasedHybridEvolver<ES>::_step\n");
  ARIADNE_ASSERT(!working.empty());
  ARIADNE_LOG(9," evolution_time="<<time<<"\n");
  ARIADNE_LOG(9," working.size()="<<working.size()<<"\n");
  ARIADNE_LOG(9," evolve.size()="<<evolve.size()<<"\n");
  ARIADNE_LOG(9," reach.size()="<<reach.size()<<"\n");
  //ARIADNE_ASSERT(working.size()<32)
  //THES thbs=working.back();
  //ARIADNE_LOG(7,"  t="<<thbs.time()<<", n="<<thbs.steps()<<", ds="<<thbs.state()<<"\n");
  //ARIADNE_LOG(7,"  bs="<<thbs.set()<<"\n");
  

  Q t; Z n; DS ds(0); ES bs;
  make_ltuple(t,n,ds,bs)=working.pop();
  ARIADNE_LOG(7,"  t="<<t<<", n="<<n<<", ds="<<ds);
  ARIADNE_LOG(7,(semantics==upper_semantics?", upper_semantics\n":", lower_semantics\n"));
  ARIADNE_LOG(7,"  bs="<<bs<<"\n");
  ARIADNE_ASSERT(t<=time);

  if(t==time) {
    ARIADNE_LOG(6," reached end time\n");
    evolve.adjoin(HES(ds,bs));
  } else if(bs.radius()>this->maximum_enclosure_radius()) {
    ARIADNE_LOG(6," subdivide\n");
    ARIADNE_LOG(7,"  r="<<bs.radius()<<"; max_r="<<this->maximum_enclosure_radius()<<"\n");
    this->append_subdivision(working,THES(t,n,ds,bs));
  } else {
    ARIADNE_LOG(6," time step\n");
    ARIADNE_LOG(9," integrator="<<*this->_integrator<<"\n");
    const DM& mode=automaton.mode(ds);
    reference_vector<const DT> transitions=automaton.transitions(ds);
    const VF& vf=mode.dynamic();
    const CS& inv=mode.invariant();

    // Compute continous evolve and reach sets
    Q h; Bx bb;
    make_lpair(h,bb)=this->flow_bounds(vf,bs.bounding_box());
    if(Q(t+h)>time) { h=time-t; }
    ARIADNE_LOG(6,"  h="<<h<<"\n");
		ARIADNE_LOG(7," bb="<<bb<<"\n");
    ES ebs=this->continuous_integration_step(vf,bs,h,bb);
    ebs=this->reduce(ebs);
    ARIADNE_LOG(7,"  ebs="<<ebs<<"\n");
    Q rh=this->_final_activation_time(vf,inv,bs,h,bb,semantics);
    ES rbs=this->continuous_reachability_step(vf,bs,rh,bb);
    ARIADNE_LOG(6,"  rh="<<rh<< "\n");
		ARIADNE_LOG(7,"  rbs="<<rbs<<"\n");
    ESL rbsl(rbs);

    // Process evolved set
    ARIADNE_LOG(9,"  i="<<this->intersect(ebs,mode.invariant()));
    ARIADNE_LOG(9,", s="<<this->subset(ebs,mode.invariant()));
    ARIADNE_LOG(9,", inv="<<mode.invariant().bounding_box()<<"\n");
    if(rh==h) {
      ARIADNE_LOG(2,"Continuous step: push with t="<<Q(t+h)<<" n="<<n<<"\n");
      intermediate.adjoin(HES(ds,ebs));
      working.push(THES(Q(t+h),n,ds,ebs));
    }
  
    // Process reach set
    for(typename ESL::const_iterator rbs_iter=rbsl.begin(); rbs_iter!=rbsl.end(); ++rbs_iter) {
      reach.adjoin(HES(ds,*rbs_iter));
    }

    // Process transitions
    ARIADNE_LOG(6,"Checking for active transitions...\n");
    for(transitions_const_iterator transition_iter=transitions.begin();
        transition_iter!=transitions.end(); ++transition_iter)
    {
      Z imn=n+1;
      const DT& transition=*transition_iter;
      ARIADNE_LOG(7,"  transition = "<<transition<<"\n");
      DS imds=transition.destination().discrete_state();
      ARIADNE_LOG(9,"  i="<<this->intersect(ebs,transition.activation()));
      ARIADNE_LOG(9,", s="<<this->subset(ebs,transition.activation()));
      ARIADNE_LOG(9,", act="<<transition.activation().bounding_box()<<"\n");
      for(typename ESL::const_iterator rbs_iter=rbsl.begin(); rbs_iter!=rbsl.end(); ++rbs_iter) {
        const ES& rbs=*rbs_iter;
        const CS& act=transition.activation();
        // Always for activation by looking for possibility of intersection, which corresponds to upper semantics
        if(_satisfies(rbs,act,upper_semantics)) { 
          Q h1=_initial_activation_time(vf,act,bs,rh,bb,semantics);
          Q h2=_final_activation_time(vf,act,bs,rh,bb,semantics);
          ARIADNE_LOG(9,"    t=["<<h1.get_d()<<","<<h2.get_d()<<"]\n");
          if(h1>h2) {
            // Can't prove that a transition occurs
          } else {
            ARIADNE_LOG(7,"  Transition occurs at time t=["<<h1.get_d()<<","<<h2.get_d()<<"]\n");
            Q imh; DS imds=transition.destination().discrete_state(); ES imbs;
            make_ltuple(imh,imbs)=this->_saltation_map(vf,transition.destination().dynamic(),transition.reset(),transition.activation(),
                                                       bs,h1,h2,bb,semantics);
            Q imt=t+imh;
            imbs=this->reduce(imbs);
            ARIADNE_LOG(9,"    imbs="<<imbs<<"\n");
            ARIADNE_LOG(4,"Transition to "<<transition.destination().discrete_state()<<" at time="<<imt.get_d()<<"\n");
            ARIADNE_ASSERT(imt>=t);
            ARIADNE_ASSERT(imt<=Q(t+h));
            ARIADNE_LOG(2,"Transition: push with imt="<<imt<<" imn="<<imn<<"\n");
            intermediate.adjoin(HES(imds,imbs));
            working.push(THES(imt,imn,imds,imbs));
          }
        }
      }
    }

  }  // done step
  ARIADNE_LOG(9," done step\n")
  ARIADNE_LOG(6,"  working.size()="<<working.size()<<"\n");
  ARIADNE_LOG(9,"  evolve.size()="<<evolve.size()<<"\n");
  ARIADNE_LOG(9,"  reach.size()="<<reach.size()<<"\n\n");
}



template<class R> void
Evolver< HybridAutomaton<R>, Zonotope<R> >::
_evolution(HybridEnclosureSetList& final,
           HybridEnclosureSetList& reachable,
           HybridEnclosureSetList& intermediate,
           const Automaton& automaton,
           const HybridEnclosureSet& initial,
           const Time& time,
           Semantics semantics,
           bool reach) const
{
  THESL working=this->timed_enclosure_set_list(initial);

  ARIADNE_LOG(5,"  working_set="<<working<<"\n\n");	
  while(working.size()!=0) { 
    this->_step(final,reachable,intermediate,working,automaton,time,semantics);
    ARIADNE_LOG(3,"  working_set="<<working<<"\n\n");	
  }
  return;
}


}
