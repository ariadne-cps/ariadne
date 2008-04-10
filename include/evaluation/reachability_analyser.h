/***************************************************************************
 *            reachability_analyser.h
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
 
/*! \file reachability_analyser.h
 *  \brief Methods for computing abstract reachable sets.
 */

#ifndef ARIADNE_REACHABILITY_ANALYSER_H
#define ARIADNE_REACHABILITY_ANALYSER_H

#include <boost/smart_ptr.hpp>

#include "base/types.h"
#include "base/declarations.h"
#include "geometry/declarations.h"
#include "system/declarations.h"
#include "evaluation/declarations.h"


#include "system/transition_system_interface.h"
#include "system/transition_system.h"
#include "evaluation/evolver_interface.h"
#include "evaluation/discretiser_interface.h"
#include "evaluation/reachability_analyser_interface.h"

// For approximation
#include "geometry/grid_approximation.h"
#include "geometry/timed_list_set.h"


namespace Ariadne {
  namespace Evaluation {

    template<class Sys,class Aprx> class ReachabilityAnalyser;

    /*! \brief A class for computing the evolution of a discrete-time autonomous system.
     *  \ingroup Analysers
     */
    template<class T,class Aprx>
    class ReachabilityAnalyser< System::TransitionSystemInterface<T,Aprx>, Aprx >
      : public ReachabilityAnalyserInterface< System::TransitionSystemInterface<T,Aprx>, typename Aprx::basic_set_type >
    {
      typedef System::TransitionSystemInterface<T,Aprx> Sys;
      typedef typename Aprx::real_type R;
      typedef typename Aprx::basic_set_type BS;
      typedef T Time;

      typedef System::TransitionSystemInterface<T,Aprx> System;
      typedef Geometry::SetInterface<BS> Set;
      typedef Geometry::SetInterface<BS>* SetPointer;

     private:
      boost::shared_ptr< EvolutionParameters<R> > _parameters;
     public:
       /*! \brief The underlying state space of the system. */
      typedef typename Sys::state_space_type state_space_type;

      //@{
      //! \name Constructors and destructors
      /*! \brief Construct from evolution parameters and a method for iterating basic sets. */
      ReachabilityAnalyser(const EvolutionParameters<R>& parameters);
      //@}

      //@{ 
      //! \name Methods to set and get the parameters controlling the accuracy.
      /*! \brief The parameters controlling the accuracy. */
      virtual const EvolutionParameters<R>& parameters() const;
      /*! \brief A reference to the parameters controlling the accuracy. */
      virtual EvolutionParameters<R>& parameters();
      //@}


      //@{
      //! \name Evaluation of systems on abstract sets
      /*! \brief Compute a lower-approximation to the set obtained by evolving \a system for \a time starting in \a initial_set. */
      virtual SetPointer lower_evolve(const System& system, const Set& initial_set, const Time& time) const;
    
      /*! \brief Compute a lower-approximation to the reachable set of \a system starting in \a initial_set up to \a time \a. */
      virtual SetPointer lower_reach(const System& system, const Set& initial_set, const Time& time) const;
    
      /*! \brief Compute an approximation to the set obtained by iterating \a time times \a system starting in \a initial_set. */
      virtual SetPointer upper_evolve(const System& system, const Set& initial_set, const Time& time) const;
    
      /*! \brief Compute an approximation to the reachable set of \a system starting in \a initial_set iterating at most \a time times. */
      virtual SetPointer upper_reach(const System& system, const Set& initial_set, const Time& time) const;
    
      /*! \brief Compute an outer-approximation to the chain-reachable set of \a system starting in \a initial_set. */
      virtual SetPointer chain_reach(const System& system, const Set& initial_set) const;
    
      /*! \brief Compute an outer-approximation to the viability kernel of \a system within \a bounding_set. */
      virtual SetPointer viable(const System& system, const Set& bounding_set) const;
    
      /*! \brief Attempt to verify that the reachable set of \a system starting in \a initial_set remains in \a safe_set. */
      virtual tribool verify(const System& system, const Set& initial_set, const Set& safe_set) const;
     //@}


     private:
      // Simplifying typedefs
      typedef typename Aprx::Paving Pv;
      typedef typename Aprx::CoverListSet CLS;
      typedef typename Aprx::PartitionListSet PLS;
      typedef typename Aprx::PartitionTreeSet PTS;
      typedef Numeric::Interval<R> I;
      typedef Numeric::Integer Integer;
      typedef Geometry::BoxListSet<R> BxLS;
      typedef Geometry::Grid<R> Gr;
      typedef Geometry::GridCell<R> GC;
      typedef Geometry::GridBlock<R> GB;
      typedef Geometry::GridCellListSet<R> GCLS;
      typedef Geometry::GridMaskSet<R> GMS;
      typedef Geometry::SetInterface<BS> SI;
      typedef Geometry::Box<R> Bx;
      typedef Geometry::ListSet<BS> BSL;
      typedef Geometry::TimedSet<T,BS> TBS;
      typedef Geometry::ListSet<TBS> TBSL;
     private:
      // Services provided by other classes
      CLS _lower_evolve(const Sys& sys, const BS& bs, const T& t) const {
        return sys.lower_evolve(bs,t); }
      CLS _lower_reach(const Sys& sys, const BS& bs, const T& t) const {
        return sys.lower_reach(bs,t); }
      std::pair<CLS,CLS> _lower_reach_evolve(const Sys& sys, const BS& bs, const T& t) const {
        return sys.lower_reach_evolve(bs,t); }
      PLS _upper_evolve(const Sys& sys, const BS& bs, const T& t) const {
        return sys.upper_evolve(bs,t); }
      PLS _upper_reach(const Sys& sys, const BS& bs, const T& t) const {
        return sys.upper_reach(bs,t); }
      std::pair<PLS,PLS> _upper_reach_evolve(const Sys& sys, const BS& bs, const T& t) const {
        return sys.upper_reach_evolve(bs,t); }
     public:
      // Helper functions for operators on lists of sets.
      PLS _upper_reach(const Sys& sys, const PLS& set, const T& time) const {
        PLS result(set.grid()); 
        for(typename PLS::const_iterator bs=set.begin(); bs!=set.end(); ++bs) {
          result.adjoin(this->_upper_reach(sys,*bs,time)); }
        return result; }
      PLS _upper_evolve(const Sys& sys, const PLS& set, const T& time) const {
        PLS result(set.grid()); 
        for(typename PLS::const_iterator bs=set.begin(); bs!=set.end(); ++bs) {
          result.adjoin(this->_upper_evolve(sys,*bs,time)); }
        return result; }
     private:
      // Helper functions for approximating sets
      Geometry::Grid<R> grid(Geometry::EuclideanSpace espc) const { 
        return Geometry::Grid<R>(espc.dimension(),this->_parameters->grid_length()); }
      Geometry::HybridGrid<R> grid(Geometry::HybridSpace hspc) const { 
        return Geometry::HybridGrid<R>(hspc,this->_parameters->grid_length()); }
      PLS _outer_approximation(const SI& s) const {
        return Geometry::outer_approximation(s,this->grid(s.space())); }
      CLS _lower_approximation(const SI& s) const {
        return Geometry::lower_approximation(s,this->grid(s.space())); }
      PLS _inner_approximation(const SI& s) const {
        return Geometry::inner_approximation(s,this->grid(s.space())); }
      PLS _outer_approximation(const BxLS& s) const {
        return GMS(Geometry::outer_approximation(s,this->grid(s.space()))); };
     private:
      // Helper functions for accessing parameters
      T lock_to_grid_time() const;
      R maximum_basic_set_radius() const { return this->_parameters->maximum_basic_set_radius(); }
      Bx bounding_domain(const Sys& system) const { return this->_parameters->bounding_domain(system.state_space().dimension()); }
      uint verbosity() const { return this->_parameters->verbosity(); }
    };






    /*! \brief A class for computing the evolution of a discrete-time autonomous system.
     *  \ingroup Analysers
     */
    template<class Sys,class Aprx>
    class ReachabilityAnalyser 
      : public ReachabilityAnalyserInterface<Sys,typename Aprx::basic_set_type>
    {
      typedef typename Sys::time_type T;
      typedef typename Sys::real_type R;
      typedef typename Aprx::basic_set_type BS;

      typedef System::TransitionSystemInterface<T,Aprx> TSI;
      typedef System::TransitionSystem<Sys,Aprx> TS;
      typedef ReachabilityAnalyser<TSI,Aprx> TSRA;

      typedef Sys System;
      typedef Geometry::SetInterface<BS> Set;
      typedef Geometry::SetInterface<BS>* SetPointer;
      typedef T Time;
     private:
      boost::shared_ptr< DiscretiserInterface<Sys,Aprx> > _discretiser;
      boost::shared_ptr< ReachabilityAnalyser<TSI,Aprx> > _analyser;
     private:
      TS _transition_system(const Sys& system) const { return TS(system,*this->_discretiser); }
     public:
      //@{
      //! \name Constructors and destructors
      /*! \brief Construct from evolution parameters and a method for iterating basic sets. */
      ReachabilityAnalyser(const EvolutionParameters<R>& parameters,
                           const DiscretiserInterface<Sys,Aprx>& discretiser)
        : _discretiser(discretiser.clone()), _analyser(new ReachabilityAnalyser<TSI,Aprx>(parameters)) { }
      //@}


      //@{ 
      //! \name Methods to set and get the parameters controlling the accuracy.
      /*! \brief The parameters controlling the accuracy. */
      virtual const EvolutionParameters<R>& parameters() const { return this->_analyser->parameters(); }

      /*! \brief A reference to the parameters controlling the accuracy. */
      virtual EvolutionParameters<R>& parameters() { return this->_analyser->parameters(); }
      //@}


      //@{
      //! \name Evaluation of systems on abstract sets

      /*! \brief Compute a lower-approximation to the set obtained by evolving \a system for \a time starting in \a initial_set. */
      virtual SetPointer lower_evolve(const System& system, const Set& initial_set, const Time& time) const {
        return this->_analyser->lower_evolve(this->_transition_system(system),initial_set,time); }
    
      /*! \brief Compute a lower-approximation to the reachable set of \a system starting in \a initial_set up to \a time \a. */
      virtual SetPointer lower_reach(const System& system, const Set& initial_set, const Time& time) const {
        return this->_analyser->lower_reach(this->_transition_system(system),initial_set,time); }
    
      /*! \brief Compute an approximation to the set obtained by iterating \a time times \a system starting in \a initial_set. */
      virtual SetPointer upper_evolve(const System& system, const Set& initial_set, const Time& time) const {
        return this->_analyser->upper_evolve(this->_transition_system(system),initial_set,time); }
    
      /*! \brief Compute an approximation to the reachable set of \a system starting in \a initial_set iterating at most \a time times. */
      virtual SetPointer upper_reach(const System& system, const Set& initial_set, const Time& time) const {
        return this->_analyser->upper_reach(this->_transition_system(system),initial_set,time); }
    
      /*! \brief Compute an outer-approximation to the chain-reachable set of \a system starting in \a initial_set. */
      virtual SetPointer chain_reach(const System& system, const Set& initial_set) const {
        return this->_analyser->chain_reach(this->_transition_system(system),initial_set); }
    
      /*! \brief Compute an outer-approximation to the viability kernel of \a system within \a bounding_set. */
      virtual SetPointer viable(const System& system, const Set& bounding_set) const {
        return this->_analyser->viable(this->_transition_system(system),bounding_set); }
    
      /*! \brief Compute an outer-approximation to the viability kernel of \a system within \a bounding_set. */
      virtual tribool verify(const System& system, const Set& initial_set, const Set& safe_set) const {
        return this->_analyser->verify(this->_transition_system(system),initial_set,safe_set); }
    
      //@}
    };








  }
}



#endif /* ARIADNE_REACHABILITY_ANALYSER_H */
