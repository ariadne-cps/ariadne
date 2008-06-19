/***************************************************************************
 *            vector_field_evolver.h
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
 
/*! \file vector_field_evolver.h
 *  \brief Methods for computing a single time step of the evolution of a vector_field.
 */

#ifndef ARIADNE_VECTOR_FIELD_EVOLVER_H
#define ARIADNE_VECTOR_FIELD_EVOLVER_H

#include <boost/smart_ptr.hpp>

#include "base/types.h"
#include "base/declarations.h"
#include "geometry/declarations.h"
#include "system/declarations.h"
#include "evaluation/declarations.h"

#include "base/tuple.h"

#include "evaluation/integrator_interface.h"
#include "evaluation/subdivider_interface.h"
#include "evaluation/reducer_interface.h"
#include "evaluation/evolver_interface.h"

// For approximation
#include "geometry/grid_approximation.h"
#include "geometry/timed_list_set.h"

#include "evaluation/evolution_parameters.h"



namespace Ariadne {
  
    template<class Sys, class ES> class Evolver;
 
    /*! \ingroup Evolvers
     *  \brief A class for evolving a continuous-time dynamical system.
     */
    template<class ES> 
    class Evolver< VectorField<typename ES::real_type>, ES>
      : public EvolverBase< VectorField<typename ES::real_type>, ES>
    {
      typedef typename ES::real_type R;
      typedef VectorField<R> Sys;
      typedef typename Sys::time_type T;
      typedef Box<R> Bx;
      typedef ListSet<ES> ESL;
      typedef TimedSet<T,ES> TES;
      typedef ListSet<TES> TESL;
      typedef Sys VF;
     public:
      Evolver(const EvolutionParameters<R>&,const IntegratorInterface<ES>&, const SubdividerInterface<ES>&, const ReducerInterface<ES>&);
      virtual Evolver<Sys,ES>* clone() const { return new Evolver<Sys,ES>(*this); }
      virtual void evolution(ESL& final, ESL& intermediate, const Sys& system, const ES& initial, const T& time, Semantics semantics, bool reach) const;
     public:
      /*! \brief Compute an approximation to the evolution set. */
      ESL evolve(const Sys& system, const ES& initial_set, const T& time) const {
        ESL final; ESL intermediate; this->evolution(final,intermediate,system,initial_set,time,upper_semantics,false); return final; }
      /*! \brief Compute an approximation to the evolution set under the given semantics. */
      ESL reach(const Sys& system, const ES& initial_set, const T& time) const {
        ESL final; ESL intermediate; this->evolution(final,intermediate,system,initial_set,time,upper_semantics,true); return intermediate; }
     private:
      // Helper functions for accessing parameters
      uint verbosity() const { 
        return this->_parameters->verbosity(); }
      T maximum_step_size() const { 
        return this->_parameters->maximum_step_size(); }
      R maximum_basic_set_radius() const { 
        return this->_parameters->maximum_basic_set_radius(); }
     private:
      // Services provided by other classes
      std::pair<T,Bx> flow_bounds(const VF& vf, const Bx& bx) const {
        return this->_integrator->flow_bounds(vf,bx,this->maximum_step_size()); }
      std::pair<T,Bx> flow_bounds(const VF& vf, const Bx& bx, const T& h) const {
        return this->_integrator->flow_bounds(vf,bx,h); }
      TES integration_step(const VF& vf, const TES& ts, const T& h, const Bx& bb) const {
        return TES(ts.time()+h,this->_integrator->integration_step(vf,ts.set(),h,bb)); }
      ES reachability_step(const VF& vf, const ES& s, const T& h, const Bx& bb) const {
        return this->_integrator->reachability_step(vf,s,h,bb); }
      R radius(const ES& s) const {
        return s.radius(); }
      Bx bounding_box(const ES& s) const {
        return s.bounding_box(); }
      ESL subdivide(const ES& s) const {
        return this->_subdivider->subdivide(s,this->maximum_basic_set_radius()); }
     private:
      // Helper functions for timed sets
      R radius(const TES& tes) const { return tes.set().radius(); }
      void adjoin_subdivision(TESL& tls, const TES& ts) const { 
        T t=ts.time();
        ESL subdivisions=this->_subdivider->subdivide(ts.set(),this->maximum_basic_set_radius());
        for(size_type i=0; i!=subdivisions.size(); ++i) {
          tls.adjoin(TES(t,subdivisions[i]));
        }
      }
      TES reduce(const TES& ts) const {
        return TES(ts.time(),this->_reducer->over_approximate(ts.set())); }
     private:
      boost::shared_ptr< EvolutionParameters<R> > _parameters;
      boost::shared_ptr< IntegratorInterface<ES> >  _integrator;
      boost::shared_ptr< SubdividerInterface<ES> > _subdivider;
      boost::shared_ptr< ReducerInterface<ES> > _reducer;

      boost::shared_ptr< EvolutionProfiler >  _profiler;
 
    };



  
} // namespace Ariadne



#endif /* ARIADNE_VECTOR_FIELD_EVOLVER_H */
