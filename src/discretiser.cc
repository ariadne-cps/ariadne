/***************************************************************************
 *            discretiser.cc
 *
 *  Copyright  2006-8  Alberto Casagrande, Pieter Collins
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
#include "config.h"

#include "discretiser.h"

#include "stlio.h"
#include "function.h"
#include "enclosure.h"
#include "orbit.h"
#include "taylor_function.h"
#include "grid_set.h"
#include "map.h"
#include "vector_field.h"

namespace Ariadne {

inline Sweeper default_sweeper() { return Sweeper(); }


template<class Sys,class ES>
Orbit<typename Discretiser<Sys,ES>::BasicSetType>
Discretiser<Sys,ES>::
evolution(const SystemType& system,
          const BasicSetType& initial_set,
          const TimeType& time,
          const Grid& grid,
          const AccuracyType accuracy,
          const Semantics semantics) const
{
    ARIADNE_LOG(3,ARIADNE_PRETTY_FUNCTION<<"\n");
    EnclosureType enclosure=this->_enclosure(initial_set);
    ARIADNE_LOG(4,"enclosure="<<enclosure<<"\n");
    Orbit<EnclosureType> continuous_orbit=this->_evolver->orbit(system,enclosure,time,semantics);
    ARIADNE_LOG(5,"continuous_orbit="<<continuous_orbit<<"\nOK\n");
    Orbit<BasicSetType> discrete_orbit=this->_discretise(continuous_orbit,initial_set,grid,accuracy);
    ARIADNE_LOG(5,"discrete_orbit="<<discrete_orbit<<"\n");
    return discrete_orbit;
}

template<class Sys,class ES>
Orbit<typename Discretiser<Sys,ES>::BasicSetType>
Discretiser<Sys,ES>::
lower_evolution(const SystemType& system,
                const BasicSetType& initial_set,
                const TimeType& time,
                const Grid& grid,
                const AccuracyType accuracy) const
{
    return this->evolution(system, initial_set, time, grid, accuracy, LOWER_SEMANTICS);
}

template<class Sys,class ES>
Orbit<typename Discretiser<Sys,ES>::BasicSetType>
Discretiser<Sys,ES>::
upper_evolution(const SystemType& system,
                const BasicSetType& initial_set,
                const TimeType& time,
                const Grid& grid,
                const AccuracyType accuracy) const
{
    return this->evolution(system, initial_set, time, grid, accuracy, UPPER_SEMANTICS);
}

template<class Sys, class ES>
typename Discretiser<Sys,ES>::EnclosureType
Discretiser<Sys,ES>::
_enclosure(const BasicSetType& initial_set) const
{
    return EnclosureType(initial_set.box(),TaylorFunctionFactory(Sweeper()));
}

template<class Sys, class ES>
Orbit<typename Discretiser<Sys,ES>::BasicSetType>
Discretiser<Sys, ES>::
_discretise(const Orbit<EnclosureType>& continuous_orbit,
            const BasicSetType& initial_set,
            const Grid& grid,
            const int accuracy) const
{
    ARIADNE_LOG(3,ARIADNE_PRETTY_FUNCTION<<"\n");
    ARIADNE_LOG(6,"continuous_reach_set="<<continuous_orbit.reach().bounding_boxes()<<"\n");
    DenotableSetType reach_set
        = outer_approximation(continuous_orbit.reach(),
                              grid,
                              accuracy);
    ARIADNE_LOG(4,"reach_set="<<reach_set<<"\n");
    DenotableSetType intermediate_set
        = outer_approximation(continuous_orbit.intermediate(),
                              grid,
                              accuracy);
    ARIADNE_LOG(4,"intermediate_set="<<intermediate_set<<"\n");
    DenotableSetType final_set
        = outer_approximation(continuous_orbit.final(),
                              grid,
                              accuracy);
    ARIADNE_LOG(4,"final_set="<<final_set<<"\n");
    Orbit<BasicSetType> orbit(initial_set,reach_set,intermediate_set,final_set);
    ARIADNE_LOG(4,"orbit="<<orbit<<"\n");
    return orbit;

}

template class Discretiser<VectorField,Enclosure>;
template class Discretiser<IteratedMap,Enclosure>;



} // namespace Ariadne


