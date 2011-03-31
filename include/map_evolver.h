/***************************************************************************
 *            map_evolver.h
 *
 *  Copyright  2007-8  Alberto Casagrande, Pieter Collins
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
 
/*! \file map_evolver.h
 *  \brief Evolver for map systems.
 */

#ifndef ARIADNE_MAP_EVOLVER_H
#define ARIADNE_MAP_EVOLVER_H

#include <string>
#include <vector>
#include <list>
#include <iostream>

#include <boost/smart_ptr.hpp>

#include "tuple.h"

#include "map.h"
#include "function_interface.h"
#include "evolver_base.h"
#include "settings.h"

#include "logging.h"

namespace Ariadne {  
  
template<class Sys, class BS> class Evolver;

class TaylorSet;
class IteratedMap;
template<class ES> class Orbit;

class EvolutionSettings;
class TaylorModel;
template<class Var> class CalculusInterface;

class EvolutionProfiler;


/*! \brief A class for computing the evolution of a map system. 
 *
 * The actual evolution steps are performed by the MapEvolver class.
 */
class MapEvolver
    : public EvolverBase< IteratedMap, TaylorSet>
    , public Loggable
{
    typedef TaylorModel VariableType;
  public:
    typedef ContinuousEvolutionSettings EvolutionSettingsType;
    typedef IteratedMap::TimeType TimeType;
    typedef IteratedMap SystemType;
    typedef TaylorSet SetModelType;
    typedef TaylorSet EnclosureType;
    typedef Orbit<EnclosureType> OrbitType;
    typedef ListSet<EnclosureType> EnclosureListType;
  public:
    
    //! \brief Default constructor.
    MapEvolver();
  
    //! \brief Construct from parameters using a default integrator.
    MapEvolver(const EvolutionSettingsType& parameters);
  
    /*! \brief Make a dynamically-allocated copy. */
    MapEvolver* clone() const { return new MapEvolver(*this); }

    //@{
    //! \name Settings controlling the evolution.

    EvolutionSettingsType& settings() { return *this->_settings; }
    const EvolutionSettingsType& settings() const { return *this->_settings; }

    //@}
  

    //@{
    //! \name Evolution using abstract sets.
    //! \brief Compute an approximation to the orbit set using upper semantics. 
    Orbit<EnclosureType> orbit(const SystemType& system, const EnclosureType& initial_set, const TimeType& time, Semantics semantics=UPPER_SEMANTICS) const;

    //! \name Evolution using abstract sets.
    //! \brief Compute an approximation to the orbit set using upper semantics, with only the continuous part. 
    Orbit<EnclosureType> upper_orbit_continuous(const SystemType& system, const EnclosureType& initial_set, const TimeType& time) {
		ARIADNE_NOT_IMPLEMENTED;
	}

    //! \brief Compute an approximation to the evolution set using upper semantics. 
    EnclosureListType evolve(const SystemType& system, const EnclosureType& initial_set, const TimeType& time) const {
        EnclosureListType final; EnclosureListType reachable; EnclosureListType intermediate; 
        this->_evolution(final,reachable,intermediate,system,initial_set,time,UPPER_SEMANTICS);
        return final; }

    //! \brief Compute an approximation to the evolution set under upper semantics. 
    EnclosureListType reach(const SystemType& system, const EnclosureType& initial_set, const TimeType& time) const {
        EnclosureListType final; EnclosureListType reachable; EnclosureListType intermediate; 
        this->_evolution(final,reachable,intermediate,system,initial_set,time,UPPER_SEMANTICS);
        return intermediate; }

    //! \brief Compute an approximation to the evolution set under the given semantics, returning the reached and final sets, and the information
    //! on having disproved.
    tuple<EnclosureListType,EnclosureListType,DisproveData> lower_chain_reach_evolve_disprove(const SystemType& system, const EnclosureType& initial_set,
																					  const TimeType& time, const HybridBoxes& disprove_bounds,
																					  const bool& skip_if_disproved) const {
        ARIADNE_NOT_IMPLEMENTED; }

  protected:
    virtual void _evolution(EnclosureListType& final, EnclosureListType& reachable, EnclosureListType& intermediate, 
                            const SystemType& system, const EnclosureType& initial, const TimeType& time, 
                            Semantics semantics) const;

    typedef tuple<TimeType, SetModelType> TimedSetType;
    virtual void _evolution_step(std::vector< TimedSetType >& working_sets, 
                                 EnclosureListType& final, EnclosureListType& reachable, EnclosureListType& intermediate,  
                                 const SystemType& system, const TimedSetType& current_set, const TimeType& time, 
                                 Semantics semantics) const;

  private:
    boost::shared_ptr< EvolutionSettingsType > _settings;
    boost::shared_ptr< CalculusInterface<TaylorModel> > _toolbox;
    //boost::shared_ptr< EvolutionProfiler >  _profiler;
};


  
} // namespace Ariadne

#endif // ARIADNE_MAP_EVOLVER_H
