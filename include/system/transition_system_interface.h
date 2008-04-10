/***************************************************************************
 *            transition_system_interface.h
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
 
/*! \file transition_system_interface.h
 *  \brief Interface for discrete transition systems with evolve and reach steps.
 */

#ifndef ARIADNE_TRANSITION_SYSTEM_INTERFACE_H
#define ARIADNE_TRANSITION_SYSTEM_INTERFACE_H

#include <boost/smart_ptr.hpp>

#include "base/types.h"
#include "base/declarations.h"
#include "geometry/declarations.h"
#include "system/declarations.h"
#include "evaluation/declarations.h"

namespace Ariadne {

  namespace System {

    /*! \ingroup System
     * \brief Interface for transition systems defined by reach and evolve operators. 
     */
    template<class T, class Aprx>
    class TransitionSystemInterface 
    {
      typedef typename Aprx::space_type StateSpace;
      typedef typename Aprx::Real Real;
      typedef typename Aprx::BasicSet BasicSet;
      typedef typename Aprx::CoverListSet CoverListSet;
      typedef typename Aprx::PartitionListSet PartitionListSet;
      typedef T Time;
     public:
      /*! \brief The type used to describe the state space. */
      typedef typename Aprx::space_type state_space_type;
      /*! \brief The type used to represent real numbers. */
      typedef Real real_type;
      /*! \brief The type used for a single basic set. */
      typedef BasicSet basic_set_type;
      /*! \brief The type used for a list of basic open sets. */
      typedef CoverListSet cover_list_set_type; 
      /*! \brief The type used for a list of basic compact sets. */
      typedef PartitionListSet partition_list_set_type; 
      /*! \brief The type used to represent time. */
      typedef Time time_type;

      /*! \brief Virtual destructor. */
      virtual ~TransitionSystemInterface() { }
      /*! \brief Cloning operation. */
      virtual TransitionSystemInterface<T,Aprx>* clone() const = 0;

      /*! \brief The underlying state space of the system. */
      virtual StateSpace state_space() const = 0;

      /*! \brief Compute a lower approximation to the evolution for time \a t starting in set \a s. */
      virtual CoverListSet lower_evolve(const BasicSet& s, const Time& t) const = 0;
      /*! \brief Compute a lower approximation to the evolution for times up to \a t starting in set \a s. */
      virtual CoverListSet lower_reach(const BasicSet&, const Time&) const = 0;
      /*! \brief Compute a lower approximation to the reachable and evolved sets for time \a t starting in set \a s. */
      virtual std::pair<CoverListSet,CoverListSet> lower_reach_evolve(const BasicSet& s, const Time& t) const = 0;
      /*! \brief Compute an upper approximation to the evolution for time \a t starting in set \a s. */
      virtual PartitionListSet upper_evolve(const BasicSet&, const Time&) const = 0;
      /*! \brief Compute a upper approximation to the evolution for times up to \a t starting in set \a s. */
      virtual PartitionListSet upper_reach(const BasicSet&, const Time&) const = 0;
      /*! \brief Compute a upper approximation to the reachable and evolved sets for time \a t starting in set \a s. */
      virtual std::pair<PartitionListSet,PartitionListSet> upper_reach_evolve(const BasicSet& s, const Time& t) const = 0;
    };

  }
}

#endif /* ARIADNE_TRANSITION_SYSTEM_INTERFACE_H */
