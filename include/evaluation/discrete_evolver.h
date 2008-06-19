/***************************************************************************
 *            discrete_evolver.h
 *
 *  Copyright  2006-8  Alberto Casagrande, Pieter Collins
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
 
/*! \file discrete_evolver.h
 *  \brief Methods for computing the evolution of systems on grids/pavings.
 */

#ifndef ARIADNE_DISCRETE_EVOLVER_H
#define ARIADNE_DISCRETE_EVOLVER_H

#include <boost/smart_ptr.hpp>

#include "base/types.h"
#include "base/declarations.h"
#include "geometry/declarations.h"
#include "system/declarations.h"
#include "evaluation/declarations.h"

#include "evaluation/evolver_interface.h"
#include "evaluation/approximator_interface.h"
#include "evaluation/discrete_evolver_interface.h"

// For approximation
#include "geometry/box_list_set.h"
#include "geometry/grid_cell_list_set.h"
#include "geometry/grid_approximation.h"
#include "geometry/timed_list_set.h"


namespace Ariadne {
  
  

 
    /*! \ingroup Evolvers
     *  \brief A class for computing the evolution of a discrete-time autonomous system.
     */
    template<class Sys, class Aprx, class ES>
    class DiscreteEvolver
      : public DiscreteEvolverInterface<Sys,Aprx>
    {
      typedef typename Sys::time_type T;
      typedef typename Sys::real_type R;
      typedef ListSet<ES> ESL;

      typedef typename Aprx::Cover Cover;
      typedef typename Aprx::Paving Paving;
      typedef typename Aprx::BasicSet BasicSet;
      typedef typename Aprx::CoverListSet CoverListSet;
      typedef typename Aprx::PartitionListSet PartitionListSet;
      typedef Sys System;
      typedef T Time;
     private:
      typedef Ariadne::NumericalSystemInterface<T,ES> NSysI;
      typedef Ariadne::NumericalSystem<Sys,ES> NSys;
     private:
      boost::shared_ptr< EvolverInterface<Sys,ES>  > _evolver;
      boost::shared_ptr< ApproximatorInterface<Aprx,ES>  > _approximator;
     private:
      // Functions for converting to numerical system discretiser
      NSys _numerical_system(const Sys& system) const { return NSys(system,*this->_evolver); }
      ESL _enclosure_set_list(const BasicSet& bs) const { return ESL(ES(bs)); }
      template<class LS> ESL _enclosure_set_list(const LS& ls) const { 
        ESL result; for(typename LS::const_iterator iter=ls.begin(); iter!=ls.end(); ++iter) {
          result.adjoin(ES(*iter)); } return result; }
     public:
      //@{
      //! \name Constructors and destructors

      /*! \brief Construct from evolution parameters and a method for evolving basic sets, and a scheme for approximating sets. */
      DiscreteEvolver(const EvolverInterface<Sys,ES>& evolver, 
                      const ApproximatorInterface<Aprx,ES>& approximator)
        : _evolver(evolver.clone()), _approximator(approximator.clone()) { }
      
      DiscreteEvolver<Sys, Aprx, ES>* clone() const { return new DiscreteEvolver<Sys,Aprx,ES>(*this); }
 
      //@}

       /*! \brief Compute a lower approximation to the evolved set of \a system starting in \a initial_set after \a time. */
      //@{
      //! \name Evaluation on basic sets.
    
      /*! \brief Compute a lower approximation to the evolved set of \a system starting in \a initial_set after \a time. */
      virtual void lower_evolution(CoverListSet& evolve, CoverListSet& reach, const System& system, const BasicSet& initial_set, const Time& time) const { 
        ESL initial_enclosures=this->_enclosure_set_list(initial_set); ESL evolved_enclosures, reached_enclosures;
        this->_evolver->evolution(evolved_enclosures,reached_enclosures,system,initial_enclosures,time,lower_semantics);
        this->_approximator->adjoin_over_approximations(evolve,evolved_enclosures); 
        this->_approximator->adjoin_over_approximations(reach,reached_enclosures); 
      }
    
      /*! \brief Compute a lower approximation to the evolved set of \a system starting in \a initial_set after \a time. */
      virtual void lower_evolution(CoverListSet& evolve, const System& system, const BasicSet& initial_set, const Time& time) const { 
        ESL initial_enclosures=this->_enclosure_set_list(initial_set); ESL evolved_enclosures;
        this->_evolver->evolution(evolved_enclosures,system,initial_enclosures,time,lower_semantics);
        this->_approximator->adjoin_over_approximations(evolve,evolved_enclosures); 
      }
    
      /*! \brief Compute a lower approximation to the evolved set of \a system starting in \a initial_set after \a time. */
      virtual void upper_evolution(PartitionListSet& evolve, PartitionListSet& reach, const System& system, const BasicSet& initial_set, const Time& time) const {
        ESL initial_enclosures=this->_enclosure_set_list(initial_set); ESL evolved_enclosures; ESL reached_enclosures;
        this->_evolver->evolution(evolved_enclosures,reached_enclosures,system,initial_set,time,lower_semantics);
        this->_approximator->adjoin_outer_approximation(evolve,evolved_enclosures); 
        this->_approximator->adjoin_outer_approximation(reach,reached_enclosures); 
      }

      /*! \brief Compute a lower approximation to the evolved set of \a system starting in \a initial_set after \a time. */
      virtual void upper_evolution(PartitionListSet& evolve, const System& system, const BasicSet& initial_set, const Time& time) const {
        ESL initial_enclosures=this->_enclosure_set_list(initial_set); ESL evolved_enclosures;
        this->_evolver->evolution(evolved_enclosures,system,initial_enclosures,time,lower_semantics);
        this->_approximator->adjoin_outer_approximation(evolve,evolved_enclosures); 
      }

      //@{
      //! \name Evaluation on list sets.
    
      /*! \brief Compute a lower approximation to the evolved set of \a system starting in \a initial_set after \a time. */
      virtual void lower_evolution(CoverListSet& evolve, CoverListSet& reach, const System& system, const CoverListSet& initial_set, const Time& time) const { 
        ESL initial_enclosures=this->_enclosure_set_list(initial_set); ESL evolved_enclosures, reached_enclosures;
        this->_evolver->evolution(evolved_enclosures,reached_enclosures,system,initial_enclosures,time,lower_semantics);
        this->_approximator->adjoin_over_approximations(evolve,evolved_enclosures); 
        this->_approximator->adjoin_over_approximations(reach,reached_enclosures); 
      }
    
      /*! \brief Compute a lower approximation to the evolved set of \a system starting in \a initial_set after \a time. */
      virtual void lower_evolution(CoverListSet& evolve, const System& system, const CoverListSet& initial_set, const Time& time) const { 
        ESL initial_enclosures=this->_enclosure_set_list(initial_set); ESL evolved_enclosures;
        this->_evolver->evolution(evolved_enclosures,system,initial_enclosures,time,lower_semantics);
        this->_approximator->adjoin_over_approximations(evolve,evolved_enclosures); 
      }
    
      /*! \brief Compute a lower approximation to the evolved set of \a system starting in \a initial_set after \a time. */
      virtual void upper_evolution(PartitionListSet& evolve, PartitionListSet& reach, const System& system, const PartitionListSet& initial_set, const Time& time) const {
        ESL initial_enclosures=this->_enclosure_set_list(initial_set); ESL evolved_enclosures; ESL reached_enclosures;
        this->_evolver->evolution(evolved_enclosures,reached_enclosures,system,initial_enclosures,time,lower_semantics);
        this->_approximator->adjoin_outer_approximation(evolve,evolved_enclosures); 
        this->_approximator->adjoin_outer_approximation(reach,reached_enclosures); 
      }

      /*! \brief Compute a lower approximation to the evolved set of \a system starting in \a initial_set after \a time. */
      virtual void upper_evolution(PartitionListSet& evolve, const System& system, const PartitionListSet& initial_set, const Time& time) const {
        ESL initial_enclosures=this->_enclosure_set_list(initial_set); ESL evolved_enclosures;
        this->_evolver->evolution(evolved_enclosures,system,initial_enclosures,time,lower_semantics);
        this->_approximator->adjoin_outer_approximation(evolve,evolved_enclosures); 
      }


    public:
      /*! \brief Compute a lower approximation to the evolved set of \a system starting in \a initial_set after \a time, giving cells in \a cover which are guarenteed to contain an evolved point. */
     virtual CoverListSet lower_evolve(const System& system, const CoverListSet& initial_set, const Time& time, const Cover& cover) const {
        CoverListSet evolve; this->lower_evolution(evolve,system,initial_set,time); return evolve; }
      /*! \brief Compute a lower approximation to the reach set of \a system starting in \a initial_set after \a time, giving cells in \a cover which are guarenteed to contain an evolved point. */
      virtual CoverListSet lower_reach(const System& system, const CoverListSet& initial_set, const Time& time, const Cover& cover) const {
        CoverListSet evolve; CoverListSet reach; this->lower_evolution(evolve,reach,system,initial_set,time); return reach; }
      /*! \brief Compute a lower approximation to the reach set of \a system starting in \a initial_set after \a time, giving cells in \a paving which are guarenteed to contain all evolved points. */
      virtual PartitionListSet upper_evolve(const System& system, const PartitionListSet& initial_set, const Time& time, const Paving& paving) const {
        PartitionListSet evolve(paving); this->upper_evolution(evolve,system,initial_set,time); return evolve; }
      /*! \brief Compute a lower approximation to the reach set of \a system starting in \a initial_set after \a time, giving cells in \a paving which are guarenteed to contain all reached points. */
      virtual PartitionListSet upper_reach(const System& system, const PartitionListSet& initial_set, const Time& time, const Paving& paving) const {
         PartitionListSet evolve(paving); PartitionListSet reach(paving); this->upper_evolution(evolve,reach,system,initial_set,time); return reach; }
    };


  
} // namespace Ariadne



#endif /* ARIADNE_DISCRETE_EVOLVER_H */
