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
#include "evaluation/evolution_profiler.h"
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
    _profiler(new EvolutionProfiler),
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
    _profiler(new EvolutionProfiler),
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
  if(semantics==upper_semantics) {
    return not definitely(disjoint(bs,inv));
  } else { 
    return definitely(subset(bs,inv));
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
  // Compute the set of times for which inv is satisfied.
  // Assume that the set is connected and nonempty
  Q initial_time=0;
  Q initial_activation_time=0;

  // Compute initial activation time
  if(!_satisfies(bs,inv,semantics)) {
    Q lower_bound=initial_time;
    Q upper_bound=maximum_time;
    for(uint i=0; i!=BISECTION_STEPS; ++i) {
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

  Q final_activation_time=maximum_time;
  Q lower_bound=0;
  Q upper_bound=maximum_time;
  for(uint i=0; i!=BISECTION_STEPS; ++i) {
    Q approx=med(lower_bound,upper_bound);
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
tuple<Q, Q>
Evaluation::SetBasedHybridEvolver<BS>::_activation_times(const VF& vf, 
                                                         const CS& inv, 
                                                         const BS& bs,
                                                         const Q& maxt,
                                                         const Bx& bb,
                                                         const Semantics sem) const 
{
  return make_tuple(_initial_activation_time(vf,inv,bs,maxt,bb,sem),
                    _final_activation_time(vf,inv,bs,maxt,bb,sem));
}




template<class BS>
tuple<Q,BS>
Evaluation::SetBasedHybridEvolver<BS>::_saltation_map(const VF& vf1, 
                                                      const VF& vf2, 
                                                      const Mp& rm, 
                                                      const CS& inv, 
                                                      const BS& bs,
                                                      const Q& maxt,
                                                      const Bx& bb,
                                                      const Semantics sem) const
{
  // Compute activation times up to the given time and apply the saltation map 
  Q h1=_initial_activation_time(vf1,inv,bs,maxt,bb,sem);
  Q h2=_final_activation_time(vf1,inv,bs,maxt,bb,sem);
  return _saltation_map(vf1,vf2,rm,inv,bs,h1,h2,bb,sem);
}

template<class BS>
tuple<Q,BS>
Evaluation::SetBasedHybridEvolver<BS>::
_saltation_map(const VF& vf1, const VF& vf2, const Mp& rm, const CS& inv, 
               const BS& bs, const Q& h1, const Q& h2, const Bx& bb, const Semantics sem) const
{
  ARIADNE_ASSERT(h1<=h2);
  if(sem==upper_semantics) {
    BS rs=continuous_evolution_step(vf1,bs,h1,h2,bb);
    return make_tuple(h1,apply(rm,rs));  } 
  else {
    Q hmed=med(h1,h2);
    BS rs=continuous_integration_step(vf1,bs,hmed,bb);
    return make_tuple(hmed,apply(rm,rs));  
  }
}




template<class BS>
BS
Evaluation::SetBasedHybridEvolver<BS>::_continuous_reachability_step(const VF& vf, 
                                                                     const CS& inv, 
                                                                     const BS& bs,
                                                                     const Q& maximum_time,
                                                                     const Bx& bb,
                                                                     const Semantics semantics) const
{
  Q h=_final_activation_time(vf,inv,bs,maximum_time,bb,semantics);
  return this->continuous_reachability_step(vf,bs,h,bb);
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
    ARIADNE_LOG(7," reached end time\n");
    evolve.adjoin(HBS(ds,bs));
  } else if(bs.radius()>this->maximum_basic_set_radius()) {
    ARIADNE_LOG(7," subdivide\n");
    ++this->_profiler->subdivisions;
    ARIADNE_LOG(7,"  r="<<bs.radius()<<"; max_r="<<this->maximum_basic_set_radius()<<"\n");
    this->append_subdivision(working,THBS(t,n,ds,bs));
  } else {
    ARIADNE_LOG(7," time step\n");
    ++this->_profiler->time_steps;
    ARIADNE_LOG(9," integrator="<<*this->_integrator<<"\n");
    const DM& mode=automaton.mode(ds);
    reference_vector<const DT> transitions=automaton.transitions(ds);
    const VF& vf=mode.dynamic();
    const CS& inv=mode.invariant();

    // Compute continous evolve and reach sets
    Q h; Bx bb;
    make_lpair(h,bb)=this->flow_bounds(vf,this->bounding_box(bs));
    this->_profiler->minimum_time_step=std::min(h,this->_profiler->minimum_time_step);
    this->_profiler->total_stepping_time+=h;
    if(Q(t+h)>time) { h=time-t; }
    ARIADNE_LOG(7,"  h="<<h<<", bb="<<bb<<"\n");
    BS ebs=this->continuous_integration_step(vf,bs,h,bb);
    ARIADNE_LOG(7,"  ebs="<<ebs<<"\n");
    Q rh=this->_final_activation_time(vf,inv,bs,h,bb,semantics);
    BS rbs=this->continuous_reachability_step(vf,bs,rh,bb);
    ARIADNE_LOG(7,"  rh="<<rh<<", rbs="<<rbs<<"\n");
    BSL rbsl(rbs);

    // Process evolved set
    ARIADNE_LOG(9,"  i="<<this->intersect(ebs,mode.invariant()));
    ARIADNE_LOG(9,", s="<<this->subset(ebs,mode.invariant()));
    ARIADNE_LOG(9,", inv="<<mode.invariant().bounding_box()<<"\n");
    if(rh==h) {
      working.push(THBS(Q(t+h),n,ds,ebs));
    }
  
    // Process reach set
    for(typename BSL::const_iterator rbs_iter=rbsl.begin(); rbs_iter!=rbsl.end(); ++rbs_iter) {
      reach.adjoin(make_pair(ds,*rbs_iter));
    }

    // Process transitions
    for(transitions_const_iterator transition_iter=transitions.begin();
        transition_iter!=transitions.end(); ++transition_iter)
    {
      Z imn=n+1;
      const DT& transition=*transition_iter;
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
            ++this->_profiler->transitions;
            Q imt; DS imds=transition.destination().discrete_state(); BS imbs;
            make_ltuple(imt,imbs)=this->_saltation_map(vf,transition.destination().dynamic(),transition.reset(),transition.activation(),
                                                       bs,h1,h2,bb,semantics);
            ARIADNE_LOG(9,"    imbs="<<imbs<<"\n");
            working.push(THBS(imt,imn,imds,imbs));
          }
        }
      }
    }

  }  // done step
  ARIADNE_LOG(9," done step\n")
  ARIADNE_LOG(9,"  working.size()="<<working.size()<<"\n");
  ARIADNE_LOG(9,"  evolve.size()="<<evolve.size()<<"\n");
  ARIADNE_LOG(9,"  reach.size()="<<reach.size()<<"\n\n");
}




template<class BS>
void
Evaluation::SetBasedHybridEvolver<BS>::_evolve(HBSL& evolve, 
                                               HBSL& reach, 
                                               const HBSL& initial, 
                                               const HA& ha,
                                               const Q& time, 
                                               const Semantics semantics) const
{
  //const_cast<SetBasedHybridEvolver<BS>*>(this)->verbosity=6;
  ARIADNE_LOG(5,"SetBasedHybridEvolver::_evolve(...)\n");
  THBSL working=initial;
  ARIADNE_LOG(5,"initial="<<initial<<"\n");
  ARIADNE_LOG(5,"working="<<working<<"\n");
  while(working.size()!=0) {
    _step(evolve,reach,working,ha,time,semantics);
  }
}

template<class BS>
typename Evaluation::SetBasedHybridEvolver<BS>::HBSL
Evaluation::SetBasedHybridEvolver<BS>::basic_set_evolve(const HA& automaton, 
                                                        const HBSL& sets,
                                                        const Q& time,
                                                        const Semantics semantics) const
{
  HBSL reach(automaton.locations());
  HBSL evolve(automaton.locations());
  THBSL working(sets);
  while(!working.empty()) {
    this->_step(evolve,reach,working,automaton,time,semantics);
  }
  return evolve;
}


template<class BS>
typename Evaluation::SetBasedHybridEvolver<BS>::HBSL
Evaluation::SetBasedHybridEvolver<BS>::basic_set_reach(const HA& automaton, 
                                                       const HBSL& sets,
                                                       const Q& time,
                                                       const Semantics semantics) const
{
  ARIADNE_LOG(1,"SetBasedHybridEvolver::basic_set_reach(...)\n");
  ARIADNE_LOG(1,"sets="<<sets<<"\n")

  HBSL reach(automaton.locations());
  HBSL evolve(automaton.locations());
  THBSL working(sets);
  while(!working.empty()) {
    this->_step(evolve,reach,working,automaton,time,semantics);
  }
  return reach;
}


template<class BS>
void
Evaluation::SetBasedHybridEvolver<BS>::_upper_evolve(HGCLS& result, const HA& ha, const HGCLS& initial_set, const T& t) const
{
  ARIADNE_LOG(4,"SetBasedHybridEvolver::_upper_evolve(result,automaton,grid_set,time)\n");
  HBSL reach(ha.locations());
  HBSL evolve(ha.locations());
  HBSL initial=this->basic_set_list(initial_set);
  this->_evolve(evolve,reach,initial,ha,t.time(),upper_semantics);
  result=this->outer_approximation(evolve,initial_set.grid());
}

template<class BS>
void
Evaluation::SetBasedHybridEvolver<BS>::_upper_reach(HGCLS& result, const HA& ha, const HGCLS& initial_set, const T& t) const
{
  ARIADNE_LOG(4,"SetBasedHybridEvolver::_upper_reach(result,automaton,grid_set,time)\n");
  HBSL reach(ha.locations());
  HBSL evolve(ha.locations());
  HBSL initial=this->basic_set_list(initial_set);
  ARIADNE_LOG(4,"initial_working="<<initial<<")\n");
  this->_evolve(evolve,reach,initial,ha,t.time(),upper_semantics);
  ARIADNE_LOG(4,"evolve="<<evolve<<")\n");
  ARIADNE_LOG(4,"reach="<<reach<<")\n");
  result=this->outer_approximation(reach,initial_set.grid());
}


template<class BS>
typename Evaluation::SetBasedHybridEvolver<BS>::HGCLS
Evaluation::SetBasedHybridEvolver<BS>::_upper_evolve(const HA& ha, const HGCLS& initial_set, const T& t) const
{
  ARIADNE_LOG(4,"SetBasedHybridEvolver::_upper_evolve(automaton,grid_set,time)\n");
  HBSL reach(ha.locations());
  HBSL evolve(ha.locations());
  HBSL initial=this->basic_set_list(initial_set);
  this->_evolve(evolve,reach,initial,ha,t.time(),upper_semantics);
  return this->outer_approximation(evolve,initial_set.grid());
}

template<class BS>
typename Evaluation::SetBasedHybridEvolver<BS>::HGCLS
Evaluation::SetBasedHybridEvolver<BS>::_upper_reach(const HA& ha, const HGCLS& initial_set, const T& t) const
{
  ARIADNE_LOG(4,"SetBasedHybridEvolver::_upper_reach(automaton,grid_set,time)\n");
  ARIADNE_LOG(5,"initial_set="<<initial_set<<"\n");
  HBSL reach(ha.locations());
  HBSL evolve(ha.locations());
  ARIADNE_LOG(5,"empty="<<evolve<<"\n");
  HBSL initial=this->basic_set_list(initial_set);
  ARIADNE_LOG(5,"initial_basic_sets="<<initial<<"\n");
  this->_evolve(evolve,reach,initial,ha,t.time(),upper_semantics);
  return this->outer_approximation(reach,initial_set.grid());
}


template<class BS>
Geometry::HybridGridMaskSet<typename BS::real_type> 
Evaluation::SetBasedHybridEvolver<BS>::lower_evolve(const System::HybridAutomaton<R>& automaton, 
                                                    const Geometry::HybridSet<R>& initial_set,
                                                    const Numeric::Rational& time) const
{
  this->_profiler->reset();
  HBSL reach, evolve;
  HBSL initial;
  THBSL working=this->timed_basic_set_list(this->lower_approximation(initial_set,this->grid(automaton.locations())));
  while(working.size()!=0) { 
    this->_step(evolve,reach,working,automaton,time,lower_semantics);
  }
  HGCLS result=this->outer_approximation(evolve,this->grid(automaton.locations()));
  ARIADNE_LOG(2,*this->_profiler);
  ARIADNE_LOG(2,"initial.size()="<<initial.size()<<" final.size()="<<evolve.size()<<"\n");
  return HGMS(result);
}

template<class BS>
Geometry::HybridGridMaskSet<typename BS::real_type> 
Evaluation::SetBasedHybridEvolver<BS>::lower_reach(const System::HybridAutomaton<R>& automaton, 
                                                   const Geometry::HybridSet<R>& initial_set,
                                                   const Numeric::Rational& time) const
{
  this->_profiler->reset();
  HBSL reach(automaton.locations()), evolve(automaton.locations());
  HBxLS initial(this->lower_approximation(initial_set,this->grid(automaton.locations())));
  THBSL working=this->timed_basic_set_list(initial);
  //THBSL working=this->timed_basic_set_list(this->lower_approximation(initial_set,this->grid(automaton.locations())));
  while(working.size()!=0) { 
    this->_step(evolve,reach,working,automaton,time,lower_semantics);
  }
  HGCLS result=this->outer_approximation(reach,this->grid(automaton.locations()));
  ARIADNE_LOG(2,*this->_profiler);
  ARIADNE_LOG(2,"initial.size()="<<initial.size()<<" final.size()="<<evolve.size()<<" reach.size()="<<reach.size()<<"\n");
  return HGMS(result);
}

template<class BS>
Geometry::HybridGridMaskSet<typename BS::real_type> 
Evaluation::SetBasedHybridEvolver<BS>::upper_evolve(const System::HybridAutomaton<R>& automaton, 
                                                    const Geometry::HybridSet<R>& initial_set,
                                                    const Numeric::Rational& time) const
{
  ARIADNE_LOG(2,"SetBasedHybridEvolver::upper_evolve(automaton,set,time)\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<"\n\n");
  this->_profiler->reset();
  HGr grid=this->grid(automaton.locations());
  ARIADNE_LOG(3,"grid="<<grid<<"\n");
  HGCLS initial=this->outer_approximation(initial_set,grid);
  ARIADNE_LOG(3,"working_set="<<initial<<"\n\n");
  HGCLS evolve=initial;
  HGCLS reach=evolve;
  Q lock_time=this->lock_to_grid_time();
  Z steps=floor(Q(time/lock_time));
  Q last_time=time-steps*lock_time;
  for(size_type i=0; i!=steps; ++i) {
    ARIADNE_LOG(5,"step i="<<i<<",  time="<<lock_time<<"\n");
    evolve=this->_upper_evolve(automaton,evolve,lock_time);
    ARIADNE_LOG(5,"working_set="<<evolve<<"\n\n");
  }
  evolve=this->_upper_evolve(automaton,evolve,lock_time);
  ARIADNE_LOG(2,*this->_profiler);
  ARIADNE_LOG(2,"initial.size()="<<initial.size()<<" final.size()="<<evolve.size()<<"\n");
  return HGMS(evolve);
}

template<class BS>
Geometry::HybridGridMaskSet<typename BS::real_type> 
Evaluation::SetBasedHybridEvolver<BS>::upper_reach(const System::HybridAutomaton<R>& automaton, 
                                                   const Geometry::HybridSet<R>& initial_set,
                                                   const Numeric::Rational& time) const
{
  ARIADNE_LOG(2,"SetBasedHybridEvolver::upper_reach(automaton,set,time)\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<"\n\n");
  this->_profiler->reset();
  HGr grid=this->grid(automaton.locations());
  ARIADNE_LOG(3,"grid="<<grid<<"\n");
  HGCLS initial=this->outer_approximation(initial_set,grid);
  ARIADNE_LOG(3,"working_set="<<initial<<"\n\n");
  HGCLS evolve=initial;
  HGCLS reach=evolve;
  Q lock_time=this->lock_to_grid_time();
  Z steps=floor(Q(time/lock_time));
  Q last_time=time-steps*lock_time;
  for(size_type i=0; i!=steps; ++i) {
    ARIADNE_LOG(5,"step i="<<i<<",  time="<<lock_time<<"\n");
    reach.adjoin(this->_upper_reach(automaton,evolve,lock_time));
    evolve=this->_upper_evolve(automaton,evolve,lock_time);
    ARIADNE_LOG(5,"working_set="<<evolve<<"\n\n");
  }
  reach.adjoin(this->_upper_evolve(automaton,evolve,lock_time));
  const_cast<int&>(verbosity)=3;
  ARIADNE_LOG(2,*this->_profiler);
  ARIADNE_LOG(2,"initial.size()="<<initial.size()<<" final.size()="<<evolve.size()<<" reach.size()="<<reach.size()<<"\n");
  return HGMS(reach);
}


template<class BS>
Geometry::HybridGridMaskSet<typename BS::real_type> 
Evaluation::SetBasedHybridEvolver<BS>::chainreach(const System::HybridAutomaton<R>& automaton, 
                                                  const Geometry::HybridSet<R>& initial_set) const
{
  //const_cast<SetBasedHybridEvolver<BS>*>(this)->verbosity=5;
  ARIADNE_LOG(2,"SetBasedHybridEvolver::chainreach(...)\n");
  this->_profiler->reset();
  HGr grid=this->grid(automaton.locations());
  HGCLS result(grid);
  T time=this->lock_to_grid_time();
  HGCLS found=this->outer_approximation(initial_set,grid);
  ARIADNE_LOG(3,"found.size()="<<found.size()<<"\n");
  this->_upper_reach(found,automaton,found,time);
  ARIADNE_LOG(3,"found.size()="<<found.size()<<"\n");
  while(!found.empty()) {
    result.adjoin(found);
    this->_upper_evolve(found,automaton,found,time);
    found.remove(result);
  }
  ARIADNE_LOG(2,*this->_profiler);
  return HGMS(result);
}






}
