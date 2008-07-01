/***************************************************************************
 *            map_evolver.h
 *
 *  Copyright  2008  Pieter Collins
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
 
/*! \file map_evolver.h
 *  \brief Methods for computing a single time step of the evolution of a map.
 */

#ifndef ARIADNE_MAP_EVOLVER_H
#define ARIADNE_MAP_EVOLVER_H

#include <boost/smart_ptr.hpp>

#include "base/types.h"
#include "base/declarations.h"
#include "geometry/declarations.h"
#include "system/declarations.h"
#include "evaluation/declarations.h"


#include "evaluation/evolver_interface.h"
#include "evaluation/applicator_interface.h"
#include "evaluation/subdivider_interface.h"
#include "evaluation/reducer_interface.h"

#include "evaluation/evolver_base.h"

// For approximation
#include "geometry/grid_approximation.h"
#include "geometry/timed_list_set.h"

#include "evaluation/evolution_parameters.h"

namespace Ariadne {
  
   template<class Sys, class ES> class Evolver;

    /*! \ingroup Evolvers
     *  \brief A class for evolving a discrete-time dynamical system.
     */
    template<class ES> 
    class Evolver< Map<typename ES::real_type>, ES >
      : public EvolverBase< Map<typename ES::real_type>, ES>
    {
      typedef typename ES::real_type R;
      typedef Map<R> Sys;
      typedef typename Sys::time_type T;

      typedef ListSet<ES> ESL;
      typedef TimedSet<T,ES> TES;
      typedef ListSet<TES> TESL;
      typedef Sys Mp;
     public:
      Evolver(const EvolutionParameters<R>&,const ApplicatorInterface<ES>&, const SubdividerInterface<ES>&, const ReducerInterface<ES>&);
      virtual Evolver<Sys,ES>* clone() const { return new Evolver<Sys,ES>(*this); }
     public:
      /*! \brief Compute an approximation to the evolution set. */
      ESL evolve(const Sys& system, const ES& initial_set, const T& time) const {
        ESL final; ESL intermediate; this->_evolution(final,intermediate,system,initial_set,time,upper_semantics,false); return final; }
      /*! \brief Compute an approximation to the evolution set under the given semantics. */
      ESL reach(const Sys& system, const ES& initial_set, const T& time) const {
        ESL final; ESL intermediate; this->_evolution(final,intermediate,system,initial_set,time,upper_semantics,true); 
        return intermediate; }
     protected:
      virtual void _evolution(ESL& final, ESL& intermediate, 
                              const Sys& system, const ES& initial, const T& time, 
                              Semantics semantics, bool reach) const;
     private:
      uint verbosity() const { 
        return this->_parameters->verbosity(); }
      R maximum_enclosure_radius() const { 
        return this->_parameters->maximum_enclosure_radius(); }
      R radius(const TES& tes) const { return tes.set().radius(); }
      void adjoin_subdivision(TESL& tls, const TES& ts) const { 
        T t=ts.time();
        ESL subdivisions=this->_subdivider->subdivide(ts.set(),this->maximum_enclosure_radius());
        for(size_type i=0; i!=subdivisions.size(); ++i) {
          tls.adjoin(TES(t,subdivisions[i]));
        }
      }
      TES apply(const Sys& sys, const TES& ts) const { 
        return TES(ts.time()+1,this->_applicator->apply(sys,ts.set()));
      }
      TES reduce(const TES& ts) const { 
        return TES(ts.time(),this->_reducer->over_approximate(ts.set()));
      }
     private:
      boost::shared_ptr< EvolutionParameters<R> > _parameters;
      boost::shared_ptr< ApplicatorInterface<ES> > _applicator;
      boost::shared_ptr< SubdividerInterface<ES> > _subdivider;
      boost::shared_ptr< ReducerInterface<ES> > _reducer;

    };



  
} // namespace Ariadne



#endif /* ARIADNE_MAP_EVOLVER_H */
