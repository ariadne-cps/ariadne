/***************************************************************************
 *            vector_field_evolver.hpp
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

/*! \file vector_field_evolver.hpp
 *  \brief Evolver for vector_field systems.
 */

#ifndef ARIADNE_VECTOR_FIELD_EVOLVER_HPP
#define ARIADNE_VECTOR_FIELD_EVOLVER_HPP

#include <string>
#include <vector>
#include <list>
#include <iostream>

#include <boost/smart_ptr.hpp>

#include "utility/tuple.hpp"

#include "dynamics/vector_field.hpp"
#include "function/function_interface.hpp"
#include "solvers/configuration_interface.hpp"
#include "solvers/integrator_interface.hpp"
#include "dynamics/evolver_base.hpp"

#include "utility/logging.hpp"

namespace Ariadne {

template<class Sys, class BS, class TRM> class Evolver;

class VectorField;
template<class ES> class Orbit;

class VectorFieldEvolverConfiguration;
class Enclosure;

class EvolutionProfiler;


/*! \brief A class for computing the evolution of a vector_field system.
 *
 * The actual evolution steps are performed by the VectorFieldEvolver class.
 */
class VectorFieldEvolver
    : public EvolverBase< VectorField, Enclosure, Float64 >
    , public Loggable
{
  public:
    typedef VectorFieldEvolverConfiguration ConfigurationType;
    typedef VectorField SystemType;
    typedef Float64 TimeType;
    typedef Float64 TerminationType;
    typedef Enclosure EnclosureType;
    typedef Pair<TimeType, EnclosureType> TimedEnclosureType;
    typedef Orbit<EnclosureType> OrbitType;
    typedef ListSet<EnclosureType> EnclosureListType;
  public:

    //! \brief Construct from parameters and an integrator to compute the flow.
    VectorFieldEvolver(
    		const SystemType& system,
            const IntegratorInterface& integrator);

    /*! \brief Make a dynamically-allocated copy. */
    VectorFieldEvolver* clone() const { return new VectorFieldEvolver(*this); }

    /* \brief Get the internal system. */
    virtual const SystemType& system() const { return *_sys_ptr; }

    //@{
    //! \name Configuration for the class.
    //! \brief A reference to the configuration controlling the evolution.
    ConfigurationType& configuration() { return *this->_configuration; }
    const ConfigurationType& configuration() const { return *this->_configuration; }

    //@}


    //@{
    //! \name Evolution using abstract sets.
    //! \brief Compute an approximation to the orbit set using upper semantics.
    Orbit<EnclosureType> orbit(const EnclosureType& initial_set, const TimeType& time, Semantics semantics=UPPER_SEMANTICS) const;

    using EvolverBase< VectorField, EnclosureType, TerminationType >::evolve;
    using EvolverBase< VectorField, EnclosureType, TerminationType >::reach;

    //! \brief Compute an approximation to the evolution set using upper semantics.
    EnclosureListType evolve(const EnclosureType& initial_set, const TimeType& time) const {
        EnclosureListType final; EnclosureListType reachable; EnclosureListType intermediate;
        this->_evolution(final,reachable,intermediate,initial_set,time,UPPER_SEMANTICS,false);
        return final; }

    //! \brief Compute an approximation to the reachable set under upper semantics.
    EnclosureListType reach(const EnclosureType& initial_set, const TimeType& time) const {
        EnclosureListType final; EnclosureListType reachable; EnclosureListType intermediate;
        this->_evolution(final,reachable,intermediate,initial_set,time,UPPER_SEMANTICS,true);
        return reachable; }

  protected:
    virtual Void _evolution(EnclosureListType& final, EnclosureListType& reachable, EnclosureListType& intermediate,
                            const EnclosureType& initial, const TimeType& time,
                            Semantics semantics, Bool reach) const;

    virtual Void _evolution_step(List< TimedEnclosureType >& working_sets,
                                 EnclosureListType& final, EnclosureListType& reachable, EnclosureListType& intermediate,
                                 const TimedEnclosureType& current_set, const TimeType& time,
                                 Semantics semantics, Bool reach) const;

  private:
    std::shared_ptr< SystemType > _sys_ptr;
    std::shared_ptr< IntegratorInterface > _integrator;
    //std::shared_ptr< EvolutionProfiler >  _profiler;
    std::shared_ptr< ConfigurationType > _configuration;
};


//! \brief Configuration for a VectorFieldEvolver, essentially for controlling the accuracy of continuous evolution methods.
class VectorFieldEvolverConfiguration : public ConfigurationInterface
{
  public:
    typedef double RealType;

    //! \brief Default constructor gives reasonable values.
    VectorFieldEvolverConfiguration();

  private:

    //! \brief The maximum allowable step size for integration.
    //! Decreasing this value increases the accuracy of the computation.
    RealType _maximum_step_size;

    //! \brief The maximum allowable radius of a basic set during integration.
    //! Decreasing this value increases the accuracy of the computation of an over-approximation.
    RealType _maximum_enclosure_radius;

  public:

    const RealType& maximum_step_size() const { return _maximum_step_size; }
    Void maximum_step_size(const RealType value) { _maximum_step_size = value; }

    const RealType& maximum_enclosure_radius() const { return _maximum_enclosure_radius; }
    Void maximum_enclosure_radius(const RealType value) { _maximum_enclosure_radius = value; }

  public:

    virtual OutputStream& write(OutputStream& os) const;
};

} // namespace Ariadne

#endif // ARIADNE_VECTOR_FIELD_EVOLVER_HPP
