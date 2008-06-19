/***************************************************************************
 *            discretiser.h
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
 
/*! \file discretiser.h
 *  \brief Methods for computing the evolution of systems on grids/pavings.
 */

#ifndef ARIADNE_DISCRETISER_H
#define ARIADNE_DISCRETISER_H

#include <boost/smart_ptr.hpp>

#include "base/types.h"
#include "base/declarations.h"
#include "geometry/declarations.h"
#include "system/declarations.h"
#include "evaluation/declarations.h"

#include "evaluation/evolution_parameters.h"
#include "evaluation/approximator_interface.h"
#include "evaluation/discretiser_interface.h"

// For approximation
#include "geometry/grid_approximation.h"
#include "geometry/timed_list_set.h"


namespace Ariadne {
  
  

    template<class Sys, class Aprx, class ES> class Discretiser;


    /*! \ingroup Evolvers
     *  \brief A class for discretising the evolution of a numerically-integrated system.
     */
    template<class T, class Aprx, class ES>
    class Discretiser<NumericalSystemInterface<T,ES>,Aprx,ES>
      : public DiscretiserInterface<NumericalSystemInterface<T,ES>,Aprx>
    {
      typedef NumericalSystemInterface<T,ES> Sys;
      typedef typename Sys::real_type R;

      typedef typename Aprx::BasicSet BasicSet;
      typedef typename Aprx::Paving Paving;
      typedef typename Aprx::CoverListSet CoverListSet;
      typedef typename Aprx::PartitionListSet PartitionListSet;
      typedef Sys System;
      typedef T Time;
     private:
      boost::shared_ptr< EvolutionParameters<R> > _parameters;
      boost::shared_ptr< ApproximatorInterface<Aprx,ES> > _approximator;
     public:
  
      //@{
      //! \name Constructors and destructors

      /*! \brief Construct from evolution parameters and a method for evolving basic sets, and a scheme for approximating sets. */
      Discretiser(const EvolutionParameters<R>& parameters, 
                  const ApproximatorInterface<Aprx,ES>& approximator);
      Discretiser<Sys,Aprx,ES>* clone() const { return new Discretiser<Sys,Aprx,ES>(*this); }
      //@}

      //@{
      //! \name Evaluation on basic sets.
    
       /*! \brief Compute a lower approximation to the evolved set of \a system after \a time starting in \a initial_set. */
      virtual CoverListSet lower_evolve(const System& system, const BasicSet& initial_set, const Time& time) const;
    
      /*! \brief Compute an approximation to the reachable set of \a system starting in \a initial_set iterating at most \a time times. */
      virtual CoverListSet lower_reach(const System& system, const BasicSet& initial_set, const Time& time) const;
    
      /*! \brief Compute an approximation to the reachable and evolved sets of \a system starting in \a initial_set iterating at most \a time times. */
      virtual std::pair<CoverListSet,CoverListSet> lower_reach_evolve(const System& system, const BasicSet& initial_set, const Time& time) const;
    
      /*! \brief Compute an upper approximation to the evolved set of \a system after \a time starting in \a initial_set. */
      virtual PartitionListSet upper_evolve(const System& system, const BasicSet& initial_set, const Time& time) const;
    
      /*! \brief Compute an approximation to the reachable set of \a system starting in \a initial_set iterating at most \a time times. */
      virtual PartitionListSet upper_reach(const System& system, const BasicSet& initial_set, const Time& time) const;

      /*! \brief Compute an approximation to the reachable and evolved sets of \a system starting in \a initial_set iterating at most \a time times. */
      virtual std::pair<PartitionListSet,PartitionListSet> upper_reach_evolve(const System& system, const BasicSet& initial_set, const Time& time) const;
     private:
      typedef ListSet<ES> ESL;
      typedef BasicSet BS;
      typedef CoverListSet CLS;
      typedef PartitionListSet PLS;

      ES _over_approximation(const BS& bs) const {
        return this->_approximator->enclosure_set(bs); }
      CLS _lower_approximation(const ESL& esl) const {
         CLS result; for(size_type i=0; i!=esl.size(); ++i) { result.adjoin(this->_approximator->bounding_box(esl[i])); } return result; }
      PLS _outer_approximation(const ESL& esl) const {
         return this->_approximator->outer_approximation(esl); }
      std::pair<CLS,CLS> _lower_approximation(const std::pair<ESL,ESL>& esl) const {
        CLS first=this->_lower_approximation(esl.first); CLS second=this->_lower_approximation(esl.second); return std::make_pair(first,second); }
      std::pair<PLS,PLS> _outer_approximation(const std::pair<ESL,ESL>& esl) const {
        PLS first=this->_outer_approximation(esl.first); PLS second=this->_outer_approximation(esl.second); return std::make_pair(first,second); }
    };





    /*! \ingroup Evolvers
     *  \brief A class for computing the evolution of a discrete-time autonomous system.
     */
    template<class Sys, class Aprx, class ES>
    class Discretiser
      : public DiscretiserInterface<Sys,Aprx>
    {
      typedef typename Sys::time_type T;
      typedef typename Sys::real_type R;

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
      boost::shared_ptr< Discretiser<NSysI,Aprx,ES>  > _discretiser;
     private:
      // Functions for converting to numerical system discretiser
      NSys _numerical_system(const Sys& system) const { return NSys(system,*this->_evolver); }
     public:
      //@{
      //! \name Constructors and destructors

      /*! \brief Construct from evolution parameters and a method for evolving basic sets, and a scheme for approximating sets. */
      Discretiser(const EvolutionParameters<R>& parameters, 
                  const EvolverInterface<Sys,ES>& evolver, 
                  const ApproximatorInterface<Aprx,ES>& approximator)
        : _evolver(evolver.clone()), _discretiser(new Discretiser<NSysI,Aprx,ES>(parameters,approximator)) { }
      
      Discretiser<Sys, Aprx, ES>* clone() const { return new Discretiser<Sys,Aprx,ES>(*this); }
 
      //@}

      //@{
      //! \name Evaluation on basic sets.
    
      /*! \brief Compute a lower approximation to the evolved set of \a system starting in \a initial_set after \a time. */
      virtual CoverListSet lower_evolve(const System& system, const BasicSet& initial_set, const Time& time) const { 
        NSys numerical_system = this->_numerical_system(system);
        NSysI& numerical_system_interface = numerical_system;
        return this->_discretiser->lower_evolve(numerical_system_interface,initial_set,time); }
    
      /*! \brief Compute a lower approximation to the reachable set of \a system starting in \a initial_set over \a time. */
      virtual CoverListSet lower_reach(const System& system, const BasicSet& initial_set, const Time& time) const { 
        return this->_discretiser->lower_reach(this->_numerical_system(system),initial_set,time); }
    
      /*! \brief Compute lower approximations to the reachable and evolved sets of \a system starting in \a initial_set over \a time. */
      virtual std::pair<CoverListSet,CoverListSet> lower_reach_evolve(const System& system, const BasicSet& initial_set, const Time& time) const { 
        return this->_discretiser->lower_reach_evolve(this->_numerical_system(system),initial_set,time); }
    
      /*! \brief Compute an upper approximation to the evolved set of \a system starting in \a initial_set after \a time. */
      virtual PartitionListSet upper_evolve(const System& system, const BasicSet& initial_set, const Time& time) const { 
        return this->_discretiser->upper_evolve(this->_numerical_system(system),initial_set,time); }
    
      /*! \brief Compute an upper approximation to the reachable set of \a system starting in \a initial_set over \a time. */
      virtual PartitionListSet upper_reach(const System& system, const BasicSet& initial_set, const Time& time) const { 
        return this->_discretiser->upper_reach(this->_numerical_system(system),initial_set,time); }

      /*! \brief Compute upper approximations to the reachable and evolved sets of \a system starting in \a initial_set over \a time. */
      virtual std::pair<PartitionListSet,PartitionListSet> upper_reach_evolve(const System& system, const BasicSet& initial_set, const Time& time) const { 
        return this->_discretiser->upper_reach_evolve(this->_numerical_system(system),initial_set,time); }
    };


  
} // namespace Ariadne



#endif /* ARIADNE_DISCRETISER_H */
