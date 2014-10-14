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

#include "utility/tuple.h"

#include "dynamics/map.h"
#include "function/function_interface.h"
#include "solvers/configuration_interface.h"
#include "dynamics/evolver_base.h"

#include "utility/logging.h"

namespace Ariadne {

template<class Sys, class BS, class TRM> class Evolver;

class Enclosure;
class IteratedMap;
template<class ES> class Orbit;

class MapEvolverConfiguration;
class EvolutionProfiler;


/*! \brief A class for computing the evolution of a map system.
 *
 * The actual evolution steps are performed by the MapEvolver class.
 */
class MapEvolver
    : public EvolverBase< IteratedMap, Enclosure, Integer>
    , public Loggable
{
    typedef TaylorModel<ValidatedNumber> VariableType;
  public:
    typedef MapEvolverConfiguration ConfigurationType;
    typedef IteratedMap SystemType;
    typedef Integer TimeType;
    typedef Integer TerminationType;
    typedef Enclosure EnclosureType;
    typedef Pair<TerminationType, EnclosureType> TimedEnclosureType;
    typedef Orbit<EnclosureType> OrbitType;
    typedef ListSet<EnclosureType> EnclosureListType;
  public:

    //! \brief Default constructor.
    MapEvolver(const SystemType& system);

    /*! \brief Make a dynamically-allocated copy. */
    MapEvolver* clone() const { return new MapEvolver(*this); }

    /* \brief Get the internal system. */
    virtual const SystemType& system() const { return *_sys_ptr; }

    //@{
    //! \name Configuration for the class.

    //! \brief A reference to the configuration.
    ConfigurationType& configuration() { return *this->_configuration; }
    const ConfigurationType& configuration() const { return *this->_configuration; }

    //@}


    //@{
    //! \name Evolution using abstract sets.
    //! \brief Compute an approximation to the orbit set using upper semantics.
    Orbit<EnclosureType> orbit(const EnclosureType& initial_set, const TerminationType& termination, Semantics semantics=UPPER_SEMANTICS) const;


    //! \brief Compute an approximation to the evolution set using upper semantics.
    EnclosureListType evolve(const EnclosureType& initial_set, const TerminationType& termination) const {
        EnclosureListType final; EnclosureListType reachable; EnclosureListType intermediate;
        this->_evolution(final,reachable,intermediate,initial_set,termination,UPPER_SEMANTICS,false);
        return final; }

    //! \brief Compute an approximation to the evolution set under upper semantics.
    EnclosureListType reach(const EnclosureType& initial_set, const TerminationType& termination) const {
        EnclosureListType final; EnclosureListType reachable; EnclosureListType intermediate;
        this->_evolution(final,reachable,intermediate,initial_set,termination,UPPER_SEMANTICS,true);
        return intermediate; }

  protected:
    virtual void _evolution(EnclosureListType& final, EnclosureListType& reachable, EnclosureListType& intermediate,
                            const EnclosureType& initial, const TerminationType& termination,
                            Semantics semantics, bool reach) const;

    virtual void _evolution_step(List< TimedEnclosureType >& working_sets,
                                 EnclosureListType& final, EnclosureListType& reachable, EnclosureListType& intermediate,
                                 const TimedEnclosureType& current_set, const TerminationType& termination,
                                 Semantics semantics, bool reach) const;

  private:
    std::shared_ptr< SystemType > _sys_ptr;
    //std::shared_ptr< EvolutionProfiler >  _profiler;
    std::shared_ptr< ConfigurationType > _configuration;
};


//! \brief Configuration for a MapEvolver, essentially to control accuracy of evolution.
class MapEvolverConfiguration : public ConfigurationInterface
{
  public:
    typedef double RealType;

    //! \brief Default constructor gives reasonable values.
    MapEvolverConfiguration();

  private:

    //! \brief The maximum allowable radius of a basic set.
    //! Decreasing this value increases the accuracy of the computation of an over-approximation.
    RealType _maximum_enclosure_radius;

  public:

    const RealType& maximum_enclosure_radius() const { return _maximum_enclosure_radius; }
    void maximum_enclosure_radius(const RealType value) { _maximum_enclosure_radius = value; }

  public:

    virtual std::ostream& write(std::ostream& os) const;

};


} // namespace Ariadne

#endif // ARIADNE_MAP_EVOLVER_H
