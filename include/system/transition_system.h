/***************************************************************************
 *            transition_system.h
 *
 *  Copyright  2007  Pieter Collins
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
 
/*! \file transition_system.h
 *  \brief A class under defined by orbits of systems.
 */

#ifndef ARIADNE_TRANSITION_SYSTEM_H
#define ARIADNE_TRANSITION_SYSTEM_H

#include <boost/smart_ptr.hpp>

#include "base/types.h"

#include "evaluation/declarations.h"

#include "system/numerical_system_interface.h"
#include "system/transition_system_interface.h"
#include "evaluation/evolver_interface.h"
#include "evaluation/discretiser_interface.h"
#include "evaluation/approximator_interface.h"

#include "system/numerical_system.h"
#include "evaluation/discretiser.h"

namespace Ariadne {
  namespace System {

    /*! \ingroup NumericalSystem
     * \brief A numerical discretisation of a dynamic system, exposing reach and evolve interfaces to an approximation scheme.
     *
     * This class is explicitly templated on the system type, since it uses the partition evolver interface to compute the evolution of the system. 
     * It is not templated on the enclosure set type, since this abstraction is handled by the partition evolver.
     */
    template<class Sys, class Aprx>
    class TransitionSystem
      : public TransitionSystemInterface<typename Sys::time_type,Aprx>
    {
      typedef typename Sys::time_type T;
      typedef typename Aprx::Real R;
      typedef typename Aprx::BasicSet BasicSet;
      typedef typename Aprx::CoverListSet CoverListSet;
      typedef typename Aprx::PartitionListSet PartitionListSet;
      typedef Sys System;
     public:
      /*! \brief The type used to denote time. */
      typedef T time_type;
      /*! \brief The type of the system being discretised. */
      typedef Sys system_type;
      /*! \brief The type used to denote time. */
      typedef typename Aprx::space_type state_space_type;
      /*! \brief The type used to represent basic sets. */
      typedef BasicSet basic_set_type;
      /*! \brief The type used to represent lists of basic open sets. */
      typedef CoverListSet cover_list_set_type;
      /*! \brief The type used to represent lists of basic compact sets. */
      typedef PartitionListSet partition_list_set_type;

     public:
      /*! \brief Construct from a system, an evolver and an approximator. */
      template<class ES> TransitionSystem(const Sys& system, 
                                          const Evaluation::EvolutionParameters<R>& parameters, 
                                          const Evaluation::EvolverInterface<Sys,ES>& evolver, 
                                          const Evaluation::ApproximatorInterface<Aprx,ES>& approximator)
        : _system(system.clone()), _discretiser(new Evaluation::Discretiser<Sys,Aprx,ES>(parameters,evolver,approximator)) { }

      /*! \brief Construct from a system and an evolution algorithm. */
      TransitionSystem(const Sys& system, const Evaluation::DiscretiserInterface<Sys,Aprx>& discretiser)
        : _system(system.clone()), _discretiser(discretiser.clone()) { }

      /*! \brief Construct from a system and an evolution algorithm. */
      TransitionSystem(shared_ptr<Sys> system_ptr, const Evaluation::DiscretiserInterface<Sys,Aprx>& discretiser)
        : _system(system_ptr), _discretiser(discretiser.clone()) { }

      /*! \brief Cloning operator. */
      virtual TransitionSystem<Sys,Aprx>* clone() const { 
        return new TransitionSystem<Sys,Aprx>(*this); }
      /*! \brief Return the discretised system. */
      const Sys& system() const { return *this->_system; }
      /*! \brief The state space of the system. */
      virtual state_space_type state_space() const { 
        return this->_system->state_space(); }

      /*! \brief Compute a lower-approximation to the evolved set under the system evolution. */
      virtual CoverListSet lower_evolve(const BasicSet& s, const T& t) const { return _discretiser->lower_evolve(this->system(),s,t); }
      /*! \brief Compute a lower-approximation to the reachable set under the system evolution. */
      virtual CoverListSet lower_reach(const BasicSet& s, const T& t) const { return _discretiser->lower_reach(this->system(),s,t); }
      /*! \brief Compute a lower-approximation to the reachable and evolved sets under the system evolution. */
      virtual std::pair<CoverListSet,CoverListSet> lower_reach_evolve(const BasicSet& s, const T& t) const { return _discretiser->lower_reach_evolve(this->system(),s,t); }
      /*! \brief Compute an upper-approximation to the evolved set under the system evolution. */
      virtual PartitionListSet upper_evolve(const BasicSet& s, const T& t) const { return _discretiser->upper_evolve(this->system(),s,t); }
      /*! \brief Compute an upper-approximation to the reachable set under the system evolution. */
      virtual PartitionListSet upper_reach(const BasicSet& s, const T& t) const { return _discretiser->upper_reach(this->system(),s,t); }
      /*! \brief Compute an upper-approximation to the reachable set under the system evolution. */
      virtual std::pair<PartitionListSet,PartitionListSet> upper_reach_evolve(const BasicSet& s, const T& t) const { return _discretiser->upper_reach_evolve(this->system(),s,t); }
     private:
      boost::shared_ptr< Sys > _system;
      boost::shared_ptr< Evaluation::DiscretiserInterface<Sys,Aprx> > _discretiser;
    };


    /*! \ingroup NumericalSystem
     *  \brief A discretisation of a numerical dynamic system, exposing reach and evolve interfaces to an approximation scheme. 
     *
     *  This class is not explicitly templated on the system type, since it uses the numerical evolution directly, but is templated on the enclosure set.
     */
    template<class T, class Aprx, class ES>
    class TransitionSystem<NumericalSystemInterface<T,ES>,Aprx>
      : public TransitionSystemInterface<T,Aprx>
    {
      typedef NumericalSystemInterface<T,ES> NSys;

      typedef typename Aprx::Real R;

      typedef typename Aprx::BasicSet BasicSet;
      typedef typename Aprx::CoverListSet CoverListSet;
      typedef typename Aprx::PartitionListSet PartitionListSet;
      typedef T Time;
     public:
      /*! \brief The type used to denote time. */
      typedef Time time_type;
      /*! \brief The type of the system being discretised. */
      typedef NSys system_type;
      /*! \brief The type used to represent basic sets. */
      typedef BasicSet basic_set_type;
      /*! \brief The type used to represent lists of basic open sets. */
      typedef CoverListSet cover_list_set_type;
      /*! \brief The type used to represent lists of basic compact sets. */
      typedef PartitionListSet partition_list_set_type;
     public:

      /*! \brief Construct from a system, an evolution algorithm and an approximation scheme. */
      template<class Sys> TransitionSystem(const Sys& system, 
                                           const Evaluation::EvolutionParameters<R>& parameters, 
                                           const Evaluation::EvolverInterface<Sys,ES>& evolver,
                                           const Evaluation::ApproximatorInterface<Aprx,ES>& approximator)
        : _system(new NumericalSystem<T,ES>(system,evolver)), 
          _discretiser(new Evaluation::Discretiser<NSys,Aprx,ES>(parameters,approximator)) { }

      /*! \brief Construct from a numerical system and an approximation scheme. */
      TransitionSystem(const NumericalSystem<T,ES>& system, 
                       const Evaluation::EvolutionParameters<R>& parameters, 
                       const Evaluation::ApproximatorInterface<Aprx,ES>& approximator)
        : _system(system.clone()), _discretiser(new Evaluation::Discretiser<NSys,Aprx,ES>(parameters,approximator)) { }

      /*! \brief Construct from a numerical system and a discretisation algorithm. */
      TransitionSystem(const NumericalSystem<T,ES>& system, 
                       const Evaluation::Discretiser<NSys,Aprx,ES>& discretiser)
        : _system(system.clone()), _discretiser(discretiser.clone()) { }

      /*! \brief Cloning operator. */
      virtual TransitionSystem<NSys,Aprx>* clone() const { 
        return new TransitionSystem<NSys,Aprx>(*this); }

      /*! \brief Return the discretised system. */
      const NumericalSystemInterface<T,ES>& system() const { return *this->_system; }


      /*! \brief Compute a lower-approximation to the evolved set under the system evolution. */
      virtual CoverListSet lower_evolve(const BasicSet& s, const Time& t) const {
        return this->_discretiser->lower_evolve(*this->_system,s,t); }
      /*! \brief Compute a lower-approximation to the reachable set under the system evolution. */
      virtual CoverListSet lower_reach(const BasicSet& s, const Time& t) const {
        return this->_discretiser->lower_reach(*this->_system,s,t); }
      /*! \brief Compute a lower-approximation to the evolved and reachable sets under the system evolution. */
      virtual std::pair<CoverListSet,CoverListSet> lower_reach_evolve(const BasicSet& s, const Time& t) const {
        return this->_discretiser->lower_reach_evolve(*this->_system,s,t); }
      /*! \brief Compute an upper-approximation to the evolved set under the system evolution. */
      virtual PartitionListSet upper_evolve(const BasicSet& s, const Time& t) const {
        return this->_discretiser->upper_evolve(*this->_system,s,t); }
      /*! \brief Compute an upper-approximation to the reachable set under the system evolution. */
      virtual PartitionListSet upper_reach(const BasicSet& s, const Time& t) const {
        return this->_discretiser->upper_reach(*this->_system,s,t); }
      /*! \brief Compute an upper-approximation to the evolved and reachable sets under the system evolution. */
      virtual std::pair<PartitionListSet,PartitionListSet> upper_reach_evolve(const BasicSet& s, const Time& t) const {
        return this->_discretiser->upper_reach_evolve(*this->_system,s,t); }
     private:
      boost::shared_ptr< NSys > _system;
      boost::shared_ptr< Evaluation::DiscretiserInterface<NSys,Aprx> > _discretiser;
    };


  }
}



#endif /* ARIADNE_TRANSITION_SYSTEM_H */
