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

#include "geometry/grid_approximation.h"
#include "geometry/list_set.h"
#include "geometry/hybrid_space.h"
#include "geometry/hybrid_set.h"
#include "geometry/hybrid_timed_set.h"
#include "system/hybrid_automaton.h"
#include "evaluation/hybrid_time.h"


#include "evaluation/approximator_interface.h"

namespace Ariadne {  
  namespace Evaluation {
  
    class HybridTime;

    /*! \ingroup Evolve
     *  \brief A class for computing the evolution of a hybrid system.
     */
    template<class BS>
    class SetBasedHybridEvolver
    {
      typedef Numeric::Integer Z;
      typedef Numeric::Rational Q;
      typedef typename BS::real_type R;
     public:
      //@{
      //! \name Constructors and destructors

      /*! \brief Construct from evolution parameters, using a default applicator and integrator. */
      SetBasedHybridEvolver(const EvolutionParameters<R>& parameters);

      /*! \brief Construct from evolution parameters, an applicator and an integrator. */
      SetBasedHybridEvolver(const EvolutionParameters<R>& parameters, 
                            const ApplicatorInterface<BS>& applicator, 
                            const IntegratorInterface<BS>& integrator);

      //@}

      //@{
      //! \name Evolution using abstract sets.
     
      /*! \brief Compute a lower approximation to the evolution up to time \a time using lower semantics. */
      Geometry::HybridGridMaskSet<R> lower_evolve(const System::HybridAutomaton<R>& automaton, 
                                                  const Geometry::HybridSet<R>& initial_set,
                                                  const Numeric::Rational& time) const;

      /*! \brief Compute a lower approximation to the reachable set up to \a time using lower semantics. */
      Geometry::HybridGridMaskSet<R> lower_reach(const System::HybridAutomaton<R>& automaton, 
                                                 const Geometry::HybridSet<R>& initial_set,
                                                 const Numeric::Rational& time) const;

      /*! \brief Compute an over approximation to the evolution up to time \a time using lower semantics. */
      Geometry::HybridGridMaskSet<R> upper_evolve(const System::HybridAutomaton<R>& automaton, 
                                                  const Geometry::HybridSet<R>& initial_set,
                                                  const Numeric::Rational& time) const;

      /*! \brief Compute an over approximation to the reachable set up to \a time using lower semantics. */
      Geometry::HybridGridMaskSet<R> upper_reach(const System::HybridAutomaton<R>& automaton, 
                                                 const Geometry::HybridSet<R>& initial_set,
                                                 const Numeric::Rational& time) const;

      /*! \brief Compute an over approximation to the chain-reachable set using upper semantics. */
      Geometry::HybridGridMaskSet<R> chainreach(const System::HybridAutomaton<R>& automaton, 
                                                const Geometry::HybridSet<R>& initial_set) const;


      //@}
      
     private:
      // Simplifying typedefs
      typedef Geometry::DiscreteState DS;
      typedef Geometry::ListSet<BS> BSL;
      typedef Geometry::Box<R> Bx;
      typedef Geometry::BoxListSet<R> BxLS;
      typedef Geometry::Grid<R> Gr;
      typedef Geometry::GridCellListSet<R> GCLS;
      typedef Geometry::GridMaskSet<R> GMS;
      typedef Geometry::SetInterface<R> SI;
      typedef Geometry::ConstraintSet<R> CS;
      typedef System::Map<R> Mp;
      typedef System::VectorField<R> VF;
      typedef Geometry::HybridSpace HSp;
      typedef Geometry::HybridSet<R> HS;
      typedef Geometry::HybridBasicSet<BS> HBS;
      typedef Geometry::HybridListSet<BS> HBSL;
      typedef Geometry::HybridBox<R> HBx;
      typedef Geometry::HybridBoxListSet<R> HBxLS;
      typedef Geometry::HybridGrid<R> HGr;
      typedef Geometry::HybridGridCell<R> HGC;
      typedef Geometry::HybridGridCellListSet<R> HGCLS;
      typedef Geometry::HybridGridMaskSet<R> HGMS;
      typedef Geometry::HybridTimedSet<BS> THBS;
      typedef Base::stack<THBS> THBSL;
      typedef System::HybridAutomaton<R> HA;
      typedef System::DiscreteMode<R> DM;
      typedef System::DiscreteTransition<R> DT;
     private:
      // Services provided by other classes
      std::pair<Q,Bx> flow_bounds(const VF& vf, const Bx& bx) const {
        return this->_integrator->flow_bounds(vf,bx,this->maximum_step_size()); }
      BS apply(const Mp& ha, const BS& bs) const {
        return this->_applicator->apply(ha,bs); }
      BS continuous_integration_step(const VF& vf, const BS& bs, const Q& h, const Bx& bb) const {
        return this->_integrator->integration_step(vf,bs,h,bb); }
      BS continuous_reachability_step(const VF& vf, const BS& bs, const Q& h, const Bx& bb) const {
        return this->_integrator->reachability_step(vf,bs,h,bb); }
      BS continuous_evolution_step(const VF& vf, const BS& bs, const Q& h1, const Q& h2, const Bx& bb) const {
        return this->_integrator->evolution_step(vf,bs,h1,h2,bb); }
      tribool disjoint(const BS& bs, const CS& cs) const {
        return Geometry::disjoint(bs,cs); }
      tribool subset(const BS& bs, const CS& cs) const {
        return Geometry::subset(bs,cs); }
      tribool intersect(const BS& bs, const CS& cs) const {
        return Geometry::intersects(bs,cs); }
      BS basic_set(const Bx& bx) const {
        return BS(bx); }
      BSL subdivide(const BS& bs) const {
        return Geometry::subdivide(bs,this->maximum_basic_set_radius()); }
      HGCLS outer_approximation(const HBSL& bsl, const HGr& g) const {
        return Geometry::outer_approximation(bsl,g); }
      HGCLS outer_approximation(const HS& s, const HGr& g) const { 
        HGCLS hgcls(g); hgcls.adjoin_outer_approximation(s); return hgcls; }
      HBxLS lower_approximation(const HS& s, const HGr& g) const { 
        HBxLS hbxls(s.locations()); 
        for(typename HS::locations_const_iterator loc_iter=s.locations_begin(); loc_iter!=s.locations_end(); ++loc_iter) {
          DS ds=loc_iter->first; 
					try{
						hbxls[ds]=Geometry::lower_approximation(s[ds],g[ds]);
					}
					catch(Geometry::EmptyInterior& e) {
						// basic sets with empty interior are not added to lower approximation
					}
				} 
        return hbxls; }
      HBxLS point_approximation(const HS& s, const HGr& g) const { 
        HBxLS hbxls(s.locations()); 
        for(typename HS::locations_const_iterator loc_iter=s.locations_begin(); loc_iter!=s.locations_end(); ++loc_iter) {
          DS ds=loc_iter->first; hbxls[ds]=Geometry::point_approximation(s[ds],g[ds]); }
        return hbxls; }
     private:
      // Helper functions for timed sets
      THBS integration_step(const VF& vf, const THBS& thbs, const Q& h, const Bx& bb) const {
        return THBS(Q(thbs.time()+h),thbs.steps(),thbs.state(),this->continuous_integration_step(vf,thbs.set(),h,bb)); }
      THBSL timed_basic_set_list(const HGCLS& hgcls) const {
        return THBSL(hgcls); }
      THBSL timed_basic_set_list(const HBxLS& hbxls) const {
        return THBSL(hbxls); }
      R radius(const THBS& thbs) const {
        return this->_approximator->radius(thbs.set()); }
      void append_subdivision(THBSL& working, const THBS& thbs) const {
        BSL sets=this->subdivide(Geometry::orthogonal_over_approximation(thbs.set()));
        for(typename BSL::const_iterator iter=sets.begin(); iter!=sets.end(); ++iter) {
          working.push(THBS(thbs.time(),thbs.steps(),thbs.state(),*iter)); } }
     private:
      typedef typename reference_vector<const DT>::const_iterator transitions_const_iterator;
     private:
      // Helper functions for computing orbits sets
      bool _satisfies(const BS&, const CS&, const Semantics) const;
      Q _initial_activation_time(const VF&, const CS&, const BS&, const Q&, const Bx&, const Semantics) const;
      Q _final_activation_time(const VF&, const CS&, const BS&, const Q&, const Bx&, const Semantics) const;
      tuple<Q,BS> _saltation_map(const VF&, const VF&, const Mp&, const CS&, const BS&, const Q&, const Q&, const Bx&, const Semantics) const;
      void _step(HBSL& evolve, HBSL& reach, THBSL& working, const HA& ha, const Q& time, const Semantics semantics) const;
     private:
      // Helper functions for accessing parameters
			HGr grid(const HSp& s) const { return this->_parameters->hybrid_grid(s); }
      HS domain(const HSp& s) const { return this->_parameters->hybrid_bounding_domain(s); }
      Q lock_to_grid_time() const { return this->_parameters->lock_to_grid_time(); }
      Q maximum_step_size() const { return this->_parameters->maximum_step_size(); }
			Q minimum_step_size() const { return this->_parameters->minimum_step_size(); }
      R maximum_basic_set_radius() const { return this->_parameters->maximum_basic_set_radius(); }
     private:
      boost::shared_ptr< EvolutionParameters<R> > _parameters;
      boost::shared_ptr< SatisfierInterface<BS> > _satisfier;
      boost::shared_ptr< ApproximatorInterface<BS> > _approximator;
      boost::shared_ptr< ApplicatorInterface<BS> > _applicator;
      boost::shared_ptr< IntegratorInterface<BS> > _integrator;
      boost::shared_ptr< EvolutionProfiler > _profiler;
      int verbosity;
    };

    std::ostream& operator<<(std::ostream& os, const Semantics& semantics) {
      return os << (semantics==lower_semantics ? "lower_semantics" : "upper_semantics"); 
    }

  }
}

#endif /* ARIADNE_SET_BASED_HYBRID_EVOLVER_H */
