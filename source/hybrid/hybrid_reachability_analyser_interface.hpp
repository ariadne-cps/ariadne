/***************************************************************************
 *            hybrid_reachability_analyser_interface.hpp
 *
 *  Copyright  2011  Luca Geretti
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

/*! \file hybrid_reachability_analyser_interface.hpp
 *  \brief Interface for analysis in the hybrid case.
 */

#ifndef ARIADNE_HYBRID_REACHABILITY_ANALYSER_INTERFACE_HPP
#define ARIADNE_HYBRID_REACHABILITY_ANALYSER_INTERFACE_HPP


#include "../hybrid/hybrid_set_interface.hpp"
#include "../dynamics/evolver_interface.hpp"
#include "../dynamics/reachability_analyser_interface.hpp"

#include "../utility/logging.hpp"

namespace Ariadne {

class HybridAutomatonInterface;
class HybridGridTreeSet;

/*! \brief A class for performing reachability analysis on a hybrid system.
 */
class HybridReachabilityAnalyserInterface
    : public virtual ReachabilityAnalyserInterface<HybridAutomatonInterface>
    , public Loggable
{
  public:
    typedef HybridAutomatonInterface SystemType;
    typedef SystemType::StateSpaceType StateSpaceType;
    typedef SystemType::TimeType TimeType;
    typedef HybridOvertSetInterface OvertSetInterfaceType;
    typedef HybridCompactSetInterface CompactSetInterfaceType;
    typedef HybridLocatedSetInterface LocatedSetInterfaceType;
    typedef HybridRegularSetInterface RegularSetInterfaceType;
    typedef HybridGridTreeSet SetApproximationType;

  public:
    //@{
    //! \name Constructors and destructors
    /*! \brief Virtual destructor */
    virtual ~HybridReachabilityAnalyserInterface() { }

    /*! \brief Make a dynamically-allocated copy. */
    virtual HybridReachabilityAnalyserInterface* clone() const = 0;
    //@}

};


//! \brief Factory for hybrid reachability analyser interface classes.
class HybridReachabilityAnalyserFactoryInterface
    : public ReachabilityAnalyserFactoryInterface<HybridAutomatonInterface>
{
  public:

    //! \brief Create a copy of the factory.
    virtual HybridReachabilityAnalyserFactoryInterface* clone() const = 0;
    //! \brief Create a hybrid reachability analyser interface object around a \a system.
    virtual HybridReachabilityAnalyserInterface* create(const HybridAutomatonInterface& system) const = 0;
};


} // namespace Ariadne


#endif // ARIADNE_HYBRID_REACHABILITY_ANALYSER_INTERFACE_HPP
