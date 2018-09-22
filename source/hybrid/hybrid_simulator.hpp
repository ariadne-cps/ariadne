/***************************************************************************
 *            hybrid_simulator.hpp
 *
 *  Copyright  2009  Pieter Collins
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

/*! \file hybrid_simulator.hpp
 *  \brief Simulator for hybrid systems.
 */

#ifndef ARIADNE_HYBRID_SIMULATOR_HPP
#define ARIADNE_HYBRID_SIMULATOR_HPP

#include "../output/logging.hpp"
#include "../hybrid/hybrid_set.decl.hpp"

namespace Ariadne {

class HybridTime;
class HybridAutomatonInterface;

template<class T> class Orbit;



/*! \brief A class for computing the evolution of a hybrid system.
 */
class HybridSimulator
    : public Loggable
{
    typedef HybridApproximatePoint EnclosureType;
    double _step_size;
  public:

    //! \brief Default constructor.
    HybridSimulator();
    Void set_step_size(double h);

    //@{
    //! \name Evolution using abstract sets.
    //! \brief Compute an approximation to the orbit set using upper semantics.
    Orbit<HybridApproximatePoint> orbit(const HybridAutomatonInterface& system, const HybridApproximatePoint& initial_point, const HybridTime& time) const;
    Orbit<HybridApproximatePoint> orbit(const HybridAutomatonInterface& system, const HybridRealPoint& initial_point, const HybridTime& time) const;
};



} // namespace Ariadne

#endif // ARIADNE_HYBRID_SIMULATOR_HPP
