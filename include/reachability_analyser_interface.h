/***************************************************************************
 *      reachability_analyser_interface.h
 *
 *  Copyright  2006-8  Pieter Collins
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
 *  \brief Interface for performing reachability analysis.
 */

#ifndef ARIADNE_REACHABILITY_ANALYSER_INTERFACE_H
#define ARIADNE_REACHABILITY_ANALYSER_INTERFACE_H

#include <boost/smart_ptr.hpp>

#include "set_interface.h"
#include "hybrid_set_interface.h"

namespace Ariadne {



template<class SYS> class ReachabilityAnalyserInterface;

//! \ingroup EvaluationModule
//! \brief Interface for computing (chain) reachable sets of a dynamic system.
//!
//! \sa \link Ariadne::EvolverInterface \c EvolverInterface<SYS,ES> \endlink
template<class SYS> class ReachabilityAnalyserInterface {
  public:
    //! \brief The type of the system.
    typedef SYS SystemType;
    //! \brief The type used to define the elapsed evolution time for the system type.
    typedef typename SystemType::TimeType TimeType;
    //! \brief The type used to describe the state space of the system evolution.
    typedef typename SystemType::StateSpaceType StateSpaceType;
    //! \brief The type used to describe the interface for overt sets in the state space, which are used as initial sets for lower reachability analysis.
    typedef typename StateSpaceType::OvertSetInterfaceType OvertSetInterfaceType;
    //! \brief The type used to describe the interface for compact sets in the state space, which are used as initial sets for upper reachability analysis.
    typedef typename StateSpaceType::CompactSetInterfaceType CompactSetInterfaceType;
    //! \brief The type used to describe the interface for located sets in the state space, which are used as initial sets for verification.
    typedef typename StateSpaceType::LocatedSetInterfaceType LocatedSetInterfaceType;
    //! \brief The type used to describe the interface for regular sets in the state space, which are used as safe sets for verification.
    typedef typename StateSpaceType::RegularSetInterfaceType RegularSetInterfaceType;
    //! \brief The type used to describe the interface for bounded sets in the state space.
    typedef typename StateSpaceType::BoundedSetInterfaceType BoundedSetInterfaceType;
    //! \brief The type used to describe the type used for concrete approximations to sets.
    typedef typename StateSpaceType::SetApproximationType SetApproximationType;
  public:
    //! \brief Virtual destructor.
    virtual ~ReachabilityAnalyserInterface() { }

    //@{
    //! \name Evaluation of maps on abstract sets

    //! \brief Compute an approximation to the set obtained by iterating \a steps times \a system starting in \a initial_set.
    virtual SetApproximationType
    lower_evolve(const SystemType& system,
     const OvertSetInterfaceType& initial_set,
     const TimeType& steps) const = 0;

    //! \brief Compute an approximation to the reachable set of \a system starting in \a initial_set iterating at most \a steps times.
    virtual SetApproximationType
    lower_reach(const SystemType& system,
    const OvertSetInterfaceType& initial_set,
    const TimeType& steps) const = 0;

    //! \brief Compute an approximation to the set obtained by iterating \a steps times \a system starting in \a initial_set.
    virtual SetApproximationType
    upper_evolve(const SystemType& system,
     const CompactSetInterfaceType& initial_set,
     const TimeType& steps) const = 0;

    //! \brief Compute an approximation to the reachable set
    //! of \a system starting in \a initial_set iterating at most \a steps times.
    virtual SetApproximationType
    upper_reach(const SystemType& system,
    const CompactSetInterfaceType& initial_set,
    const TimeType& steps) const = 0;

    //! \brief Compute an outer-approximation to the chain-reachable set
    //! of \a system starting in \a initial_set.
    virtual SetApproximationType
    chain_reach(const SystemType& system,
    const CompactSetInterfaceType& initial_set) const = 0;

    //! \brief Compute an outer-approximation to the viability kernel
    //! of \a system within \a bounding_set.
    virtual SetApproximationType
    viable(const SystemType& system,
     const CompactSetInterfaceType& bounding_set) const = 0;

    //! \brief Attempt to verify that the reachable set
    //! of \a system starting in \a initial_set remains in \a safe_set.
    virtual tribool
    verify(const SystemType& system,
     const LocatedSetInterfaceType& initial_set,
     const RegularSetInterfaceType& safe_set) const = 0;
    //@}

};


} // namespace Ariadne




#endif // ARIADNE_ANALYSER_INTERFACE_H
