/***************************************************************************
 *            evolver_interface.h
 *
 *  Copyright  2008-10  Pieter Collins, Luca Geretti
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
 
/*! \file evolver_interface.h
 *  \brief Interface for computing a single time step of the evolution of a system.
 */

#ifndef ARIADNE_EVOLVER_INTERFACE_H
#define ARIADNE_EVOLVER_INTERFACE_H

#include "evolution_parameters.h"
#include "evolution_statistics.h"
#include "disprove_data.h"

namespace Ariadne {

using std::pair;


template<class ES> class ListSet;
template<class ES> class Orbit;

enum Semantics { LOWER_SEMANTICS, UPPER_SEMANTICS }; 

  
/*! \brief Interface for evolving a dynamic system.
 *
 * \sa \link Ariadne::CalculusInterface \c CalculusInterface<S,M,F> \endlink
 *   , \link Ariadne::ReachabilityAnalyserInterface \c ReachabilityAnalyserInterface<SYS> \endlink
 */
template<class SYS, class ES> 
class EvolverInterface 
{
  public:
    typedef SYS SystemType;
    typedef ES EnclosureType;
    typedef typename SystemType::TimeType TimeType;
    typedef ListSet<EnclosureType> EnclosureListType;
	typedef ContinuousEvolutionParameters EvolutionParametersType;


    //! \brief Virtual destructor. 
    virtual ~EvolverInterface() {};

    //! \brief Cloning operator.
    virtual EvolverInterface<SYS,ES>* clone() const = 0;

    //! \brief Write to an output stream. 
    virtual std::ostream& write(std::ostream& os) const = 0;

	//! \brief Set the parameters of the evolution.
	virtual EvolutionParametersType& parameters() = 0;

	//! \brief Get the parameters of the evolution.
    virtual const EvolutionParametersType& parameters() const = 0;


  public:
    //! \brief Compute an approximation to the evolved set under the given semantics. 
    virtual 
    Orbit<EnclosureType>
    orbit(const SystemType& system, 
          const EnclosureType& initial_set, 
          const TimeType& time, 
          Semantics semantics) const = 0;

    //! \brief Compute an approximation to the evolved set for upper semantics, with continuous evolution only. 
    virtual 
    Orbit<EnclosureType>
    upper_orbit_continuous(const SystemType& system, 
          				   const EnclosureType& initial_set, 
          				   const TimeType& time) const = 0;

    //! \brief Compute an approximation to the evolved set under the given semantics. 
    virtual 
    EnclosureListType 
    evolve(const SystemType& system, 
           const EnclosureType& initial_set, 
           const TimeType& time, 
           Semantics semantics) const = 0;

    //! \brief Compute an approximation to the reachable set under the given semantics. 
    virtual 
    EnclosureListType 
    reach(const SystemType& system, 
          const EnclosureType& initial_set, 
          const TimeType& time, 
          Semantics semantics) const = 0;

    //! \brief Compute an approximation to the evolved and reachable sets under the given semantics. 
    virtual 
    pair<EnclosureListType,EnclosureListType> 
    reach_evolve(const SystemType& system, 
                 const EnclosureType& initial_set, 
                 const TimeType& time, 
                 Semantics semantics) const = 0;

    //! \brief Compute an approximation to the evolved and reachable sets under lower semantics for chain reachability,
    //! where the disproving result is provided.
    virtual
    tuple<EnclosureListType,EnclosureListType,DisproveData>
    lower_chain_reach_evolve_disprove(const SystemType& system,
									  const EnclosureType& initial_set,
									  const TimeType& time, const HybridBoxes& disprove_bounds,
									  const bool& skip_if_disproved) const = 0;
  
    //! \brief Compute an approximation to the evolved set under the given semantics. 
    virtual 
    void 
    evolution(EnclosureListType& final, 
              const SystemType& system, 
              const EnclosureType& initial, 
              const TimeType& time, 
              Semantics semantics) const = 0;

    //! \brief Compute an approximation to the evolved and reachable sets 
    //! under the given semantics. 
    virtual void evolution(EnclosureListType& final, 
                           EnclosureListType& intermediate, 
                           const SystemType& system, 
                           const EnclosureType& initial, 
                           const TimeType& time, 
                           Semantics semantics) const = 0;
  

    //! \brief Compute an approximation to the evolved set under the given semantics, 
    //! starting from a list of enclosure sets. 
    virtual 
    void 
    evolution(EnclosureListType& final, 
              const SystemType& system, 
              const EnclosureListType& initial, 
              const TimeType& time, 
              Semantics semantics) const = 0;

    //! \brief Compute an approximation to the evolved and reachable sets 
    //! under the given semantics starting from a list of enclosure sets. 
    virtual 
    void 
    evolution(EnclosureListType& final, 
              EnclosureListType& intermediate, 
              const SystemType& system, 
              const EnclosureListType& initial, 
              const TimeType& time, 
              Semantics semantics) const = 0;
  

};


template<class SYS, class ES> inline
std::ostream& 
operator<<(std::ostream& os, const EvolverInterface<SYS,ES>& e) {
    return e.write(os); 
}


} // namespace Ariadne



#endif // ARIADNE_EVOLVER_INTERFACE_H
