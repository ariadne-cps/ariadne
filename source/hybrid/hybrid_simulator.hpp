/***************************************************************************
 *            hybrid/hybrid_simulator.hpp
 *
 *  Copyright  2009-20  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

/*! \file hybrid/hybrid_simulator.hpp
 *  \brief Simulator for hybrid systems.
 */

#ifndef ARIADNE_HYBRID_SIMULATOR_HPP
#define ARIADNE_HYBRID_SIMULATOR_HPP

#include "conclog/include/logging.hpp"
#include "solvers/configuration_interface.hpp"
#include "hybrid/hybrid_set.decl.hpp"

using namespace ConcLog;

namespace Ariadne {

class HybridTerminationCriterion;
class HybridAutomatonInterface;
class DiscreteEvent;
class DiscreteLocation;

class HybridSimulatorConfiguration;

template<class T> class Orbit;

/*! \brief A class for computing the evolution of a hybrid system.
 */
class HybridSimulator
{
  public:
    typedef HybridPoint<FloatDPApproximation> HybridApproximatePointType;
    typedef Point<FloatDPApproximation> ApproximatePointType;
    typedef HybridSimulatorConfiguration ConfigurationType;
    typedef HybridAutomatonInterface SystemType;
    typedef HybridApproximatePointType EnclosureType;
    typedef Orbit<HybridApproximatePointType> OrbitType;
    typedef HybridTerminationCriterion TerminationType;
  private:
    SharedPointer<SystemType> _sys_ptr;
    SharedPointer<ConfigurationType> _configuration;
  public:

    //! \brief Default constructor.
    HybridSimulator(const SystemType& system);

    //!@{
    //! \name Configuration for the class.
    //! \brief A reference to the configuration controlling the evolution.
    ConfigurationType& configuration() { return *this->_configuration; }
    const ConfigurationType& configuration() const { return *this->_configuration; }

    //!@{
    //! \name Evolution using abstract sets.
    //! \brief Compute an approximation to the orbit set using upper semantics.
    Orbit<HybridApproximatePointType> orbit(const HybridApproximatePointType& initial_point, const TerminationType& termination) const;
    Orbit<HybridApproximatePointType> orbit(const HybridRealPoint& initial_point, const TerminationType& termination) const;
    //! \brief Currently simulates from the midpoint of the set in the first location only.
    Orbit<HybridApproximatePointType> orbit(const HybridBoundedConstraintSet& initial_set, const TerminationType& termination) const;

  private:
    Map<DiscreteEvent,EffectiveScalarMultivariateFunction> _guard_functions(const DiscreteLocation& location) const;
    Bool _satisfies_invariants(const DiscreteLocation& location, const Point<FloatDPApproximation>& point) const;
};

class HybridSimulatorConfiguration : public ConfigurationInterface {
public:
    //! \brief Default constructor gives reasonable values.
    HybridSimulatorConfiguration();

    virtual ~HybridSimulatorConfiguration() = default;

private:

    //! \brief The step size for integration.
    //! Decreasing this value increases the accuracy of the computation.
    FloatDPApproximation _step_size;

public:

    const FloatDPApproximation& step_size() const { return _step_size; }
    Void set_step_size(ApproximateDouble value) { _step_size = FloatDPApproximation(value,double_precision); }

    virtual OutputStream& _write(OutputStream& os) const;
};

} // namespace Ariadne

#endif // ARIADNE_HYBRID_SIMULATOR_HPP
