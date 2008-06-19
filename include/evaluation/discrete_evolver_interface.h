/***************************************************************************
 *            discrete_evolver_interface.h
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
 
/*! \file discrete_evolver_interface.h
 *  \brief Interface for computing a single time step of the evolution of a system.
 */

#ifndef ARIADNE_DISCRETE_EVOLVER_INTERFACE_H
#define ARIADNE_DISCRETE_EVOLVER_INTERFACE_H

#include "evaluation/declarations.h"

namespace Ariadne {
  
  
    /*! \ingroup EvaluatorInterfaces \ingroup Evolvers
     *  \brief Interface for evolving a dynamic system and discretising the result on a grid.
     */
    template<class Sys, class Aprx>
    class DiscreteEvolverInterface 
    {
      typedef typename Sys::time_type T;
      typedef typename Aprx::BasicSet BasicSet;
      typedef typename Aprx::Cover Cover;
      typedef typename Aprx::Paving Paving;
      typedef typename Aprx::CoverListSet CoverListSet;
      typedef typename Aprx::PartitionListSet PartitionListSet;
      typedef Sys System;
      typedef T Time;
     public:
      /*! \brief The type used to denote time. */
      typedef T time_type;
      /*! \brief The type of the system. */
      typedef Sys system_type;
      /*! \brief The type used to represent basic sets. */
      typedef BasicSet basic_set_type;
      /*! \brief The type used to represent lists of basic open sets. */
      typedef CoverListSet cover_list_set_type;
      /*! \brief The type used to represent lists of basic compact sets. */
      typedef PartitionListSet partition_list_set_type;
     public:
      /*! \brief Destructor. */
      virtual ~DiscreteEvolverInterface() { };

      /*! \brief Make a dynamically-allocated copy. */
      virtual DiscreteEvolverInterface<Sys,Aprx>* clone() const = 0;

      /*! \brief Compute a lower-approximation to the evolved and reachable sets under the system evolution. */
      virtual CoverListSet lower_evolve(const System& f, const CoverListSet& s, const Time& t, const Cover& cv) const = 0;
      /*! \brief Compute a lower-approximation to the evolved sets under the system evolution. */
      virtual CoverListSet lower_reach(const System& f, const CoverListSet& s, const Time& t, const Cover& cv) const = 0;
      /*! \brief Compute an outer-approximation to the evolved and reachable sets under the system evolution. */
      virtual PartitionListSet upper_evolve(const System& f, const PartitionListSet& s, const Time& t, const Paving& pv) const = 0;
      /*! \brief Compute an outer-approximation to the evolved sets under the system evolution. */
      virtual PartitionListSet upper_reach(const System& f, const PartitionListSet& s, const Time& t, const Paving& pv) const = 0;
     public:
      /*! \brief Compute a lower-approximation to the evolved and reachable sets under the system evolution. */
      virtual void lower_evolution(CoverListSet& e, CoverListSet& r, const System& f, const CoverListSet& s, const Time& t) const = 0;
      /*! \brief Compute a lower-approximation to the evolved sets under the system evolution. */
      virtual void lower_evolution(CoverListSet& e, const System& f, const CoverListSet& s, const Time& t) const = 0;
      /*! \brief Compute an outer-approximation to the evolved and reachable sets under the system evolution. */
      virtual void upper_evolution(PartitionListSet& e, PartitionListSet& r, const System& f, const PartitionListSet& s, const Time& t) const = 0;
      /*! \brief Compute an outer-approximation to the evolved sets under the system evolution. */
      virtual void upper_evolution(PartitionListSet& e, const System& f, const PartitionListSet& s, const Time& t) const = 0;
     public:
      /*! \brief Compute a lower-approximation to the evolved and reachable sets under the system evolution. */
      virtual void lower_evolution(CoverListSet& e, CoverListSet& r, const System& f, const BasicSet& s, const Time& t) const = 0;
      /*! \brief Compute a lower-approximation to the evolved sets under the system evolution. */
      virtual void lower_evolution(CoverListSet& e, const System& f, const BasicSet& s, const Time& t) const = 0;
      /*! \brief Compute an outer-approximation to the evolved and reachable sets under the system evolution. */
      virtual void upper_evolution(PartitionListSet& e, PartitionListSet& r, const System& f, const BasicSet& s, const Time& t) const = 0;
      /*! \brief Compute an outer-approximation to the evolved sets under the system evolution. */
      virtual void upper_evolution(PartitionListSet& e, const System& f, const BasicSet& s, const Time& t) const = 0;
    };


  
} // namespace Ariadne



#endif /* ARIADNE_DISCRETE_EVOLVER_INTERFACE_H */
