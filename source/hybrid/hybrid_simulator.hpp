/***************************************************************************
 *            hybrid_simulator.hpp
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

/*! \file hybrid_simulator.hpp
 *  \brief Simulator for hybrid systems.
 */

#ifndef ARIADNE_HYBRID_SIMULATOR_HPP
#define ARIADNE_HYBRID_SIMULATOR_HPP

#include "utility/logging.hpp"
#include "hybrid/hybrid_set.decl.hpp"

namespace Ariadne {

class HybridTime;
class HybridAutomatonInterface;

template<class T> class Orbit;



/*! \brief A class for computing the evolution of a hybrid system.
 */
class HybridSimulator
    : public Loggable
{
    typedef HybridExactPoint EnclosureType;
    double _step_size;
  public:

    //! \brief Default constructor.
    HybridSimulator();
    Void set_step_size(double h);

    //@{
    //! \name Evolution using abstract sets.
    //! \brief Compute an approximation to the orbit set using upper semantics.
    Orbit<HybridExactPoint> orbit(const HybridAutomatonInterface& system, const HybridExactPoint& initial_point, const HybridTime& time) const;
};



} // namespace Ariadne

#endif // ARIADNE_HYBRID_SIMULATOR_HPP
