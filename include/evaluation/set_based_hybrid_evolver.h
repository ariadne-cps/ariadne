/***************************************************************************
 *            set_based_hybrid_evolver.h
 *
 *  Copyright  2006-7  Alberto Casagrande,  Pieter Collins
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
 
/*! \file set_based_hybrid_evolver.h
 *  \brief Evolver for hybrid systems defined using invariant and activation sets.
 */

#ifndef ARIADNE_SET_BASED_HYBRID_EVOLVER_H
#define ARIADNE_SET_BASED_HYBRID_EVOLVER_H

#include <string>
#include <vector>
#include <list>
#include <iostream>

#include "base/types.h"
#include "base/tuple.h"
#include "base/stack.h"
#include "base/tribool.h"

#include "geometry/declarations.h"
#include "system/declarations.h"
#include "evaluation/declarations.h"

#include "geometry/grid_approximation_scheme.h"
#include "geometry/grid_approximation.h"
#include "geometry/list_set.h"
#include "geometry/hybrid_space.h"
#include "geometry/hybrid_set.h"
#include "geometry/hybrid_timed_set.h"
#include "system/hybrid_automaton.h"
#include "evaluation/hybrid_time.h"
#include "evaluation/evolver_interface.h"



#include "output/logging.h"

namespace Ariadne {  
  
  
    class HybridTime;
    template<class Sys, class ES> class Evolver;

    /*! \ingroup Evolvers
     *  \brief A class for computing the evolution of a hybrid system.
     */
    template<class R>
    class Evolver< HybridAutomaton<R>, Zonotope<R> >
      : public EvolverBase< HybridAutomaton<R>, HybridBasicSet< Zonotope<R> > >
    {
      typedef Integer Z;
      typedef Rational Q;
      typedef Zonotope<R> ES;
      typedef ListSet<ES> ESL;
      typedef HybridBasicSet<ES> HES;
      typedef ListSet<HES> HESL;
      
      typedef HybridAutomaton<R> Automaton;
      typedef HybridAutomaton<R> Sys;
      typedef Rational Time;

      typedef HES HybridEnclosureSet;
      typedef HESL HybridEnclosureSetList;
     public:
      //@{
      //! \name Constructors and destructors

      /*! \brief Construct from evolution parameters, an applicator and an integrator. */
      Evolver(const EvolutionParameters<R>& parameters, 
              const ApplicatorInterface<ES>& applicator, 
              const IntegratorInterface<ES>& integrator, 
              const SatisfierInterface<ES>& satisfier, 
              const SubdividerInterface<ES>& reducer, 
              const ReducerInterface<ES>& reducer);

      /*! \brief Make a dynamically-allocated copy. */
      Evolver<Sys,ES>* clone() const { return new Evolver<Sys,ES>(*this); }

      //@}

      //@{
      //! \name Evolution using abstract sets.
     
      /*! \brief Compute a lower approximation to the evolution up to time \a time using lower semantics. */
      virtual void evolution(HybridEnclosureSetList& final,
                             HybridEnclosureSetList& intermediate,
                             const Automaton& automaton,
                             const HybridEnclosureSet& initial,
                             const Time& time,
                             Semantics semantics,
                             bool reach) const;
      //@}
      
     private:
      // Simplifying typedefs
      typedef DiscreteState DS;
      typedef Box<R> Bx;
      typedef SetInterface<R> SI;
      typedef ConstraintSet<R> CS;
      typedef Map<R> Mp;
      typedef VectorField<R> VF;
      typedef HybridSpace HSp;
      typedef HybridSet<R> HS;
      typedef HybridBox<R> HBx;
      typedef HybridTimedSet<ES> THES;
      typedef stack<THES> THESL;
      typedef HybridAutomaton<R> HA;
      typedef DiscreteMode<R> DM;
      typedef DiscreteTransition<R> DT;
     private:
      // Services provided by other classes
      std::pair<Q,Bx> flow_bounds(const VF& vf, const Bx& bx) const {
        return this->_integrator->flow_bounds(vf,bx,this->maximum_step_size()); }
      ES apply(const Mp& ha, const ES& es) const {
        return this->_applicator->apply(ha,es); }
      ES continuous_integration_step(const VF& vf, const ES& es, const Q& h, const Bx& bb) const {
        return this->_integrator->integration_step(vf,es,h,bb); }
      ES continuous_reachability_step(const VF& vf, const ES& es, const Q& h, const Bx& bb) const {
        return this->_integrator->reachability_step(vf,es,h,bb); }
      ES continuous_evolution_step(const VF& vf, const ES& es, const Q& h1, const Q& h2, const Bx& bb) const {
        return this->_integrator->evolution_step(vf,es,h1,h2,bb); }
      tribool disjoint(const ES& es, const CS& cs) const {
        return Ariadne::disjoint(es,cs); }
      tribool subset(const ES& es, const CS& cs) const {
        return Ariadne::subset(es,cs); }
      tribool intersect(const ES& es, const CS& cs) const {
        return Ariadne::intersects(es,cs); }
      ES reduce(const ES& es) const {
        return this->_reducer->over_approximate(es); }
      ESL subdivide(const ES& es) const {
        return Ariadne::subdivide(es,this->maximum_basic_set_radius()); }
     private:
      // Helper functions for timed sets
      THES integration_step(const VF& vf, const THES& thes, const Q& h, const Bx& bb) const {
        return THES(Q(thes.time()+h),thes.steps(),thes.state(),this->continuous_integration_step(vf,thes.set(),h,bb)); }
      THESL timed_enclosure_set_list(const HES& hes) const {
        return THESL(HESL(hes)); }
      THESL timed_enclosure_set_list(const HESL& hesl) const {
        return THESL(hesl); }
      R radius(const THES& thes) const {
        return thes.set().radius(); }
      void append_subdivision(THESL& working, const THES& thes) const {
        ESL sets=this->subdivide(orthogonal_over_approximation(thes.set()));
        for(typename ESL::const_iterator iter=sets.begin(); iter!=sets.end(); ++iter) {
          working.push(THES(thes.time(),thes.steps(),thes.state(),*iter)); } }
     private:
      typedef typename reference_vector<const DT>::const_iterator transitions_const_iterator;
     private:
      // Helper functions for computing orbits sets
      bool _satisfies(const ES&, const CS&, const Semantics) const;
      Q _initial_activation_time(const VF&, const CS&, const ES&, const Q&, const Bx&, const Semantics) const;
      Q _final_activation_time(const VF&, const CS&, const ES&, const Q&, const Bx&, const Semantics) const;
      tuple<Q,ES> _saltation_map(const VF&, const VF&, const Mp&, const CS&, const ES&, const Q&, const Q&, const Bx&, const Semantics) const;
      void _step(HESL& evolve, HESL& reach, THESL& working, const HA& ha, const Q& time, const Semantics semantics) const;
     private:
      // Helper functions for accessing parameters
      HS domain(const HSp& s) const { return this->_parameters->hybrid_bounding_domain(s); }
      Q maximum_step_size() const { return this->_parameters->maximum_step_size(); }
      Q minimum_step_size() const { return this->_parameters->minimum_step_size(); }
      R maximum_basic_set_radius() const { return this->_parameters->maximum_basic_set_radius(); }
     private:
      boost::shared_ptr< EvolutionParameters<R> > _parameters;
      boost::shared_ptr< ApplicatorInterface<ES> > _applicator;
      boost::shared_ptr< IntegratorInterface<ES> > _integrator;
      boost::shared_ptr< SatisfierInterface<ES> > _satisfier;
      boost::shared_ptr< SubdividerInterface<ES> > _subdivider;
      boost::shared_ptr< ReducerInterface<ES> > _reducer;
      boost::shared_ptr< EvolutionProfiler > _profiler;
      int verbosity;
    };

    inline std::ostream& operator<<(std::ostream& os, const Semantics& semantics) {
      return os << (semantics==lower_semantics ? "lower_semantics" : "upper_semantics"); 
    }

  
} // namespace Ariadne

#endif /* ARIADNE_SET_BASED_HYBRID_EVOLVER_H */
