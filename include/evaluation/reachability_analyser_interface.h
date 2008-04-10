/***************************************************************************
 *            reachability_analyser_interface.h
 *
 *  Copyright  2006-7  Alberto Casagrande, Pieter Collins
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
 
/*! \file reachability_analyser_interface.h
 *  \brief Interface for computing the evolution of sets.
 */

#ifndef ARIADNE_REACHABILITY_ANALYSER_INTERFACE_H
#define ARIADNE_REACHABILITY_ANALYSER_INTERFACE_H

#include <boost/smart_ptr.hpp>

#include "base/types.h"
#include "base/declarations.h"
#include "geometry/declarations.h"
#include "system/declarations.h"
#include "evaluation/declarations.h"


namespace Ariadne {

  namespace Evaluation {

    /*! \brief Interface for computing (chain) reachable sets of a dynamic system.
     *  \ingroup EvaluatorInterfaces \ingroup Analysers
     */
    template<class Sys, class BS>
    class ReachabilityAnalyserInterface {
      typedef typename Sys::time_type T;
      typedef typename Sys::real_type R;

      typedef Sys System;
      typedef T Time;
      typedef Geometry::SetInterface<BS> Set;
      typedef Geometry::SetInterface<BS>* SetPointer;

     public:
      /*! \brief Virtual destructor. */
      virtual ~ReachabilityAnalyserInterface() { }

      //@{
      //! \name Evaluation of maps on abstract sets

      /*! \brief Compute an approximation to the set obtained by iterating \a steps times \a system starting in \a initial_set. */
      virtual SetPointer lower_evolve(const System& system, const Set& initial_set, const Time& steps) const = 0;
    
      /*! \brief Compute an approximation to the reachable set of \a system starting in \a initial_set iterating at most \a steps times. */
      virtual SetPointer lower_reach(const System& system, const Set& initial_set, const Time& steps) const = 0;
    
      /*! \brief Compute an approximation to the set obtained by iterating \a steps times \a system starting in \a initial_set. */
      virtual SetPointer upper_evolve(const System& system, const Set& initial_set, const Time& steps) const = 0;
    
      /*! \brief Compute an approximation to the reachable set of \a system starting in \a initial_set iterating at most \a steps times. */
      virtual SetPointer upper_reach(const System& system, const Set& initial_set, const Time& steps) const = 0;
    
      /*! \brief Compute an outer-approximation to the chain-reachable set of \a system starting in \a initial_set. */
      virtual SetPointer chain_reach(const System& system, const Set& initial_set) const = 0;
    
      /*! \brief Compute an outer-approximation to the viability kernel of \a system within \a bounding_set. */
      virtual SetPointer viable(const System& system, const Set& bounding_set) const = 0;
    
      /*! \brief Attempt to verify that the reachable set of \a system starting in \a initial_set remains in \a safe_set. */
      virtual tribool verify(const System& system, const Set& initial_set, const Set& safe_set) const = 0;
      //@}

    };



  }
}

 


#endif /* ARIADNE_REACHABILITY_ANALYSER_INTERFACE_H */
