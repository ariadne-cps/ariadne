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
#include "geometry/rectangle_expression.h"
#include "geometry/box_list_set.h"
#include "geometry/set_interface.h"
#include "geometry/hybrid_set.h"
#include "system/hybrid_automaton.h"
#include "evaluation/hybrid_time.h"
#include "evaluation/evolution_parameters.h"
#include "evaluation/map_evolver.h"
#include "evaluation/vector_field_evolver.h"

#include "evaluation/applicator_interface.h"
#include "evaluation/integrator_interface.h"

#include "output/epsstream.h"
#include "output/logging.h"

#include "set_based_hybrid_evolver.h"

// For default applicator/integrator
#include "kuhn_applicator.h"
#include "kuhn_integrator.h"

namespace Ariadne {
 
typedef Numeric::Rational Q;


namespace Evaluation { 
const uint BISECTION_STEPS=8;
static int& verbosity = hybrid_evolver_verbosity; 
}


template<class BS> 
Evaluation::SetBasedHybridEvolver<BS>::SetBasedHybridEvolver(const EvolutionParameters<R>& parameters)
  : _parameters(parameters.clone()),
    _applicator(new KuhnApplicator<typename BS::real_type>(3)),
    _integrator(new KuhnIntegrator<typename BS::real_type>(4,3)), // 4th order, casacade size 3
    verbosity(parameters.verbosity())
{
}


template<class BS> 
Evaluation::SetBasedHybridEvolver<BS>::SetBasedHybridEvolver(const EvolutionParameters<R>& parameters, 
                                                             const ApplicatorInterface<BS>& applicator, 
                                                             const IntegratorInterface<BS>& integrator)
  : _parameters(parameters.clone()),
    _applicator(applicator.clone()),
    _integrator(integrator.clone()),
    verbosity(parameters.verbosity())
{
}


template<class BS>
bool
Evaluation::SetBasedHybridEvolver<BS>::_satisfies(const BS& bs,
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



template<class BS>
Q
Evaluation::SetBasedHybridEvolver<BS>::_initial_activation_time(const VF& vf, 
                                                                 const CS& inv, 
                                                                 const BS& bs,
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
      BS rs=continuous_evolution_step(vf,bs,lower_bound,approx,bb);
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



template<class BS>
Q
Evaluation::SetBasedHybridEvolver<BS>::
_final_activation_time(const VF& vf, const CS& inv, const BS& bs, 
                       const Q& maximum_time, const Bx& bb, const Semantics semantics) const 
{
  ARIADNE_LOG(8,"  SetBasedHybridEvolver<BS>::_final_activation_time\n");
  // Compute final activation time
  BS es=continuous_integration_step(vf,bs,maximum_time,bb);
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
    BS rs=continuous_evolution_step(vf,bs,approx,upper_bound,bb);
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



template<class BS>
tuple<Q,BS>
Evaluation::SetBasedHybridEvolver<BS>::
_saltation_map(const VF& vf1, const VF& vf2, const Mp& rm, const CS& inv, 
               const BS& bs, const Q& h1, const Q& h2, const Bx& bb, const Semantics sem) const
{
  ARIADNE_LOG(8,"\nEvaluation::SetBasedHybridEvolver<BS>::_saltation_map\n");
  ARIADNE_LOG(9,(sem==upper_semantics?"  upper_semantics\n":"  lower_semantics\n"));
  ARIADNE_ASSERT(h1<=h2);
  if(sem==upper_semantics) {
    BS rs=continuous_evolution_step(vf1,bs,h1,h2,bb);
    ARIADNE_LOG(9,"  rs="<< rs << "\n");    
		return make_tuple(h1,apply(rm,rs));  } 
  else {
    Q hmed=med(h1,h2);
    BS rs=continuous_integration_step(vf1,bs,hmed,bb);
    ARIADNE_LOG(9,"  rs="<< rs << "\n");    
    return make_tuple(hmed,apply(rm,rs));  
  }
}



template<class BS>
void
Evaluation::SetBasedHybridEvolver<BS>::_step(HBSL& evolve,
                                             HBSL& reach,
                                             THBSL& working,
                                             const HA& automaton, 
                                             const Q& time,
                                             const Semantics semantics) const
{
  //const_cast<SetBasedHybridEvolver<BS>*>(this)->verbosity=9;
  ARIADNE_LOG(6,"\nSetBasedHybridEvolver<BS>::_step\n");
  ARIADNE_ASSERT(!working.empty());
  ARIADNE_LOG(9," evolution_time="<<time<<"\n");
  ARIADNE_LOG(9," working.size()="<<working.size()<<"\n");
  ARIADNE_LOG(9," evolve.size()="<<evolve.size()<<"\n");
  ARIADNE_LOG(9," reach.size()="<<reach.size()<<"\n");
  //ARIADNE_ASSERT(working.size()<32)
  //THBS thbs=working.back();
  //ARIADNE_LOG(7,"  t="<<thbs.time()<<", n="<<thbs.steps()<<", ds="<<thbs.state()<<"\n");
  //ARIADNE_LOG(7,"  bs="<<thbs.set()<<"\n");
  

  Q t; Z n; DS ds(0); BS bs;
  make_ltuple(t,n,ds,bs)=working.pop();
  ARIADNE_LOG(7,"  t="<<t<<", n="<<n<<", ds="<<ds);
  ARIADNE_LOG(7,(semantics==upper_semantics?", upper_semantics\n":", lower_semantics\n"));
  ARIADNE_LOG(7,"  bs="<<bs<<"\n");
  ARIADNE_ASSERT(t<=time);

  if(t==time) {
    ARIADNE_LOG(6," reached end time\n");
    evolve.adjoin(HBS(ds,bs));
  } else if(bs.radius()>this->maximum_basic_set_radius()) {
    ARIADNE_LOG(6," subdivide\n");
    ARIADNE_LOG(7,"  r="<<bs.radius()<<"; max_r="<<this->maximum_basic_set_radius()<<"\n");
    this->append_subdivision(working,THBS(t,n,ds,bs));
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
    BS ebs=this->continuous_integration_step(vf,bs,h,bb);
    ARIADNE_LOG(7,"  ebs="<<ebs<<"\n");
    Q rh=this->_final_activation_time(vf,inv,bs,h,bb,semantics);
    BS rbs=this->continuous_reachability_step(vf,bs,rh,bb);
    ARIADNE_LOG(6,"  rh="<<rh<< "\n");
		ARIADNE_LOG(7,"  rbs="<<rbs<<"\n");
    BSL rbsl(rbs);

    // Process evolved set
    ARIADNE_LOG(9,"  i="<<this->intersect(ebs,mode.invariant()));
    ARIADNE_LOG(9,", s="<<this->subset(ebs,mode.invariant()));
    ARIADNE_LOG(9,", inv="<<mode.invariant().bounding_box()<<"\n");
    if(rh==h) {
			ARIADNE_LOG(3,"Continuous step: push with t="<<Q(t+h)<<" n="<<n<<"\n");
      working.push(THBS(Q(t+h),n,ds,ebs));
    }
  
    // Process reach set
    for(typename BSL::const_iterator rbs_iter=rbsl.begin(); rbs_iter!=rbsl.end(); ++rbs_iter) {
      reach.adjoin(make_pair(ds,*rbs_iter));
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
      for(typename BSL::const_iterator rbs_iter=rbsl.begin(); rbs_iter!=rbsl.end(); ++rbs_iter) {
        const BS& rbs=*rbs_iter;
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
            Q imh; DS imds=transition.destination().discrete_state(); BS imbs;
            make_ltuple(imh,imbs)=this->_saltation_map(vf,transition.destination().dynamic(),transition.reset(),transition.activation(),
                                                       bs,h1,h2,bb,semantics);
            Q imt=t+imh;
						ARIADNE_LOG(9,"    imbs="<<imbs<<"\n");
						ARIADNE_LOG(6,"Transition to "<<transition.destination().discrete_state()<<" at time="<<imt.get_d()<<"\n");
            ARIADNE_ASSERT(imt>=t);
						ARIADNE_ASSERT(imt<=Q(t+h));
            ARIADNE_LOG(3,"Transition: push with imt="<<imt<<" imn="<<imn<<"\n");
            working.push(THBS(imt,imn,imds,imbs));
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



template<class BS>
Geometry::HybridGridMaskSet<typename BS::real_type> 
Evaluation::SetBasedHybridEvolver<BS>::lower_evolve(const System::HybridAutomaton<R>& automaton, 
                                                    const Geometry::HybridSet<R>& initial_set,
                                                    const Numeric::Rational& time) const
{
  ARIADNE_LOG(2,"SetBasedHybridEvolver::lower_evolve(automaton,set,time)\n");
  ARIADNE_LOG(3,"  initial_set="<<initial_set<<"\n\n");
  HBSL reach(automaton.locations());
	HBSL evolve(automaton.locations());
	THBSL working=this->timed_basic_set_list(this->lower_approximation(initial_set,this->grid(automaton.locations())));

  ARIADNE_LOG(5,"  working_set="<<working<<"\n\n");	
	while(working.size()!=0) { 
    this->_step(evolve,reach,working,automaton,time,lower_semantics);
		ARIADNE_LOG(5,"  working_set="<<working<<"\n\n");	
  }
  HGCLS result=this->outer_approximation(evolve,this->grid(automaton.locations()));
  return HGMS(result);
}

template<class BS>
Geometry::HybridGridMaskSet<typename BS::real_type> 
Evaluation::SetBasedHybridEvolver<BS>::lower_reach(const System::HybridAutomaton<R>& automaton, 
                                                   const Geometry::HybridSet<R>& initial_set,
                                                   const Numeric::Rational& time) const
{
  ARIADNE_LOG(2,"SetBasedHybridEvolver::lower_reach(automaton,set,time)\n");
  ARIADNE_LOG(3,"  initial_set="<<initial_set<<"\n\n");
  HBSL reach(automaton.locations()), evolve(automaton.locations());
  THBSL working=this->timed_basic_set_list(this->lower_approximation(initial_set,this->grid(automaton.locations())));
  
  ARIADNE_LOG(5,"  working_set="<<working<<"\n\n");	
  while(working.size()!=0) { 
    this->_step(evolve,reach,working,automaton,time,lower_semantics);
		ARIADNE_LOG(5,"  working_set="<<working<<"\n\n");			
  }
  HGCLS result=this->outer_approximation(reach,this->grid(automaton.locations()));
  return HGMS(result);
}

template<class BS>
Geometry::HybridGridMaskSet<typename BS::real_type> 
Evaluation::SetBasedHybridEvolver<BS>::upper_evolve(const System::HybridAutomaton<R>& automaton, 
                                                    const Geometry::HybridSet<R>& initial_set,
                                                    const Numeric::Rational& time) const
{
  ARIADNE_LOG(2,"SetBasedHybridEvolver::upper_evolve(automaton,set,time)\n");
  ARIADNE_LOG(3,"  initial_set="<<initial_set<<"\n\n");
  HBSL reach(automaton.locations()), evolve(automaton.locations());
  THBSL working=this->timed_basic_set_list(this->outer_approximation(initial_set,this->grid(automaton.locations())));
  ARIADNE_LOG(5,"  working_set="<<working<<"\n\n");	
  while(working.size()!=0) { 
    this->_step(evolve,reach,working,automaton,time,upper_semantics);
		ARIADNE_LOG(5,"  working_set="<<working<<"\n\n");			
  }
  HGCLS result=this->outer_approximation(evolve,this->grid(automaton.locations()));
  return HGMS(result);
}

template<class BS>
Geometry::HybridGridMaskSet<typename BS::real_type> 
Evaluation::SetBasedHybridEvolver<BS>::upper_reach(const System::HybridAutomaton<R>& automaton, 
                                                   const Geometry::HybridSet<R>& initial_set,
                                                   const Numeric::Rational& time) const
{
  ARIADNE_LOG(2,"SetBasedHybridEvolver::upper_reach(automaton,set,time)\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<"\n\n");
  HBSL reach(automaton.locations()), evolve(automaton.locations());
  THBSL working=this->timed_basic_set_list(this->outer_approximation(initial_set,this->grid(automaton.locations())));
  ARIADNE_LOG(5,"  working_set="<<working<<"\n\n");	
  while(working.size()!=0) { 
    this->_step(evolve,reach,working,automaton,time,upper_semantics);
		ARIADNE_LOG(5,"  working_set="<<working<<"\n\n");			
  }
  HGCLS result=this->outer_approximation(reach,this->grid(automaton.locations()));
  return HGMS(result);
}


template<class BS>
Geometry::HybridGridMaskSet<typename BS::real_type> 
Evaluation::SetBasedHybridEvolver<BS>::chainreach(const System::HybridAutomaton<R>& automaton, 
                                                  const Geometry::HybridSet<R>& initial_set) const
{
  ARIADNE_LOG(2,"SetBasedHybridEvolver::chainreach(autonaton, initial_set)\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<"\n\n");	
  HGr grid=this->grid(automaton.locations());
	HGMS domain=this->outer_approximation(this->domain(automaton.locations()),grid);
	ARIADNE_LOG(3,"domain="<<this->domain(automaton.locations())<<"\n");	
	ARIADNE_LOG(3,"domain.size()="<<domain.size()<<"\n");	
  HGCLS result(grid);
  Numeric::Rational time=this->lock_to_grid_time();
  HGCLS found=this->outer_approximation(initial_set,grid);
	found.restrict(domain);
	ARIADNE_LOG(3,"initial_set.size()="<<found.size()<<"\n");
  while(!found.empty()) {
		ARIADNE_LOG(3,"found.size()="<<found.size()<<"\n");
		THBSL working=this->timed_basic_set_list(found);
		ARIADNE_LOG(3,"  working_set="<<working<<"\n\n");	
		HBSL reach(automaton.locations()), evolve(automaton.locations());
		while(working.size()!=0) { 
			this->_step(evolve,reach,working,automaton,time,upper_semantics);
			ARIADNE_LOG(3,"  working_set="<<working<<"\n\n");			
		}
		found.adjoin(this->outer_approximation(evolve,this->grid(automaton.locations())));
		found.remove(result);
		found.restrict(domain);
		ARIADNE_LOG(3,"found "<<found.size()<<" new cells\n");
		result.adjoin(this->outer_approximation(reach,this->grid(automaton.locations())));
		ARIADNE_LOG(3,"result.size()="<<result.size()<<"\n");
  }
	result.restrict(domain);
	ARIADNE_LOG(3,"Final result.size()="<<result.size()<<"\n");
  return HGMS(result);
}

}
