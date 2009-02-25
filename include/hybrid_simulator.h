/***************************************************************************
 *            hybrid_simulator.h
 *
 *  Copyright  2009  Pieter Collins
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
 
/*! \file hybrid_simulator.h
 *  \brief Simulator for hybrid systems.
 */

#ifndef ARIADNE_HYBRID_SIMULATOR_H
#define ARIADNE_HYBRID_SIMULATOR_H

#include "evolver_base.h"
#include "logging.h"

namespace Ariadne {

class SimulationToolboxInterface;
typedef int DiscreteState;
class Point;
class HybridPoint;
class HybridAutomaton;

/*! \brief A class for computing the evolution of a hybrid system. 
 *
 * The actual evolution steps are performed by the HybridEvolver class.
 */
class HybridSimulator
    : public EvolverBase<HybridAutomaton, HybridPoint >
    , public Loggable
{
    typedef HybridPoint EnclosureType;
  public:
    
    //! \brief Default constructor.
    HybridSimulator();
  
    //! \brief Construct from parameters using a default integrator.
    HybridSimulator(const EvolutionParameters& parameters);
  
    /*! \brief Make a dynamically-allocated copy. */
    HybridSimulator* clone() const;

    //@{
    //! \name Parameters controlling the evolution.
    //! \brief A reference to the parameters controlling the evolution.
    EvolutionParameters& parameters() { return *this->_parameters; }
    const EvolutionParameters& parameters() const { return *this->_parameters; }

    //@}
  

    //@{
    //! \name Evolution using abstract sets.
    //! \brief Compute an approximation to the orbit set using upper semantics. 
    Orbit<EnclosureType> orbit(const SystemType& system, const EnclosureType& initial_point, const TimeType& time, Semantics semantics=UPPER_SEMANTICS) const;

  protected:
    virtual void _evolution(EnclosureListType&,EnclosureListType&,EnclosureListType&,const SystemType& system, const EnclosureType& initial_point, const TimeType& time, Semantics semantics, bool) const;

  private:
    boost::shared_ptr< EvolutionParameters > _parameters;
    boost::shared_ptr< SimulationToolboxInterface > _toolbox;
    //boost::shared_ptr< EvolutionProfiler >  _profiler;
};


  
} // namespace Ariadne

#endif // ARIADNE_HYBRID_EVOLVER_H
