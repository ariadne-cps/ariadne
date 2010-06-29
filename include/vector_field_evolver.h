/***************************************************************************
 *            vector_field_evolver.h
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
 
/*! \file vector_field_evolver.h
 *  \brief Evolver for vector_field systems.
 */

#ifndef ARIADNE_VECTOR_FIELD_EVOLVER_H
#define ARIADNE_VECTOR_FIELD_EVOLVER_H

#include <string>
#include <vector>
#include <list>
#include <iostream>

#include <boost/smart_ptr.hpp>

#include "tuple.h"

#include "vector_field.h"
#include "function_interface.h"
#include "evolver_base.h"
#include "evolution_parameters.h"

#include "logging.h"

namespace Ariadne {  
  
template<class Sys, class BS> class Evolver;

class VectorField;
template<class ES> class Orbit;

class EvolutionParameters;
class TaylorModel;
class TaylorSet;
template<class Var> class CalculusInterface;

class EvolutionProfiler;


/*! \brief A class for computing the evolution of a vector_field system. 
 *
 * The actual evolution steps are performed by the VectorFieldEvolver class.
 */
class VectorFieldEvolver
    : public EvolverBase< VectorField, TaylorSet >
    , public Loggable
{
  public:
    typedef ContinuousEvolutionParameters EvolutionParametersType;
    typedef VectorField::TimeType TimeType;
    typedef VectorField SystemType;
    typedef TaylorSet EnclosureType;
    typedef Orbit<EnclosureType> OrbitType;
    typedef ListSet<EnclosureType> EnclosureListType;
  public:
    
    //! \brief Default constructor.
    VectorFieldEvolver();
  
    //! \brief Construct from parameters using a default integrator.
    VectorFieldEvolver(const EvolutionParametersType& parameters);
  
    /*! \brief Make a dynamically-allocated copy. */
    VectorFieldEvolver* clone() const { return new VectorFieldEvolver(*this); }

    //@{
    //! \name Parameters controlling the evolution.
    //! \brief A reference to the parameters controlling the evolution.
    EvolutionParametersType& parameters() { return *this->_parameters; }
    const EvolutionParametersType& parameters() const { return *this->_parameters; }

    //@}
  

    //@{
    //! \name Evolution using abstract sets.
    //! \brief Compute an approximation to the orbit set using upper semantics. 
    Orbit<EnclosureType> orbit(const SystemType& system, const EnclosureType& initial_set, const TimeType& time, Semantics semantics=UPPER_SEMANTICS) const;

    //! \name Evolution using abstract sets.
    //! \brief Compute an approximation to the orbit set using upper semantics, with only the continuous part. 
    Orbit<EnclosureType> upper_orbit_continuous(const SystemType& system, const EnclosureType& initial_set, const TimeType& time) const {
		ARIADNE_NOT_IMPLEMENTED;
	}

    using EvolverBase< VectorField, TaylorSet >::evolve;
    using EvolverBase< VectorField, TaylorSet >::reach;

    //! \brief Compute an approximation to the evolution set using upper semantics.
    EnclosureListType evolve(const SystemType& system, const EnclosureType& initial_set, const TimeType& time) const {
        EnclosureListType final; EnclosureListType reachable; EnclosureListType intermediate; 
        this->_evolution(final,reachable,intermediate,system,initial_set,time,UPPER_SEMANTICS,false); 
        return final; }

    //! \brief Compute an approximation to the reachable set under upper semantics.
    EnclosureListType reach(const SystemType& system, const EnclosureType& initial_set, const TimeType& time) const {
        EnclosureListType final; EnclosureListType reachable; EnclosureListType intermediate; 
        this->_evolution(final,reachable,intermediate,system,initial_set,time,UPPER_SEMANTICS,true); 
        return reachable; }

    //! \brief Compute an approximation to the evolution set under the given semantics, returning the reached and final sets, and the information
    //! on having disproved.
    tuple<EnclosureListType,EnclosureListType,bool> lower_chain_reach_evolve_disprove(const SystemType& system, const EnclosureType& initial_set,
    																					  const TimeType& time, const HybridBoxes& disprove_bounds,
    																					  const bool& skip_if_disproved) const {
            ARIADNE_NOT_IMPLEMENTED; }

  protected:
    virtual void _evolution(EnclosureListType& final, EnclosureListType& reachable, EnclosureListType& intermediate, 
                            const SystemType& system, const EnclosureType& initial, const TimeType& time, 
                            Semantics semantics, bool reach) const;

    typedef tuple<TimeType, EnclosureType> TimedSetType;
    virtual void _evolution_step(std::vector< TimedSetType >& working_sets, 
                                 EnclosureListType& final, EnclosureListType& reachable, EnclosureListType& intermediate,  
                                 const SystemType& system, const TimedSetType& current_set, const TimeType& time, 
                                 Semantics semantics, bool reach) const;

  private:
    boost::shared_ptr< EvolutionParametersType > _parameters;
    boost::shared_ptr< CalculusInterface<TaylorModel> > _toolbox;
    //boost::shared_ptr< EvolutionProfiler >  _profiler;
};


  
} // namespace Ariadne

#endif // ARIADNE_VECTOR_FIELD_EVOLVER_H
