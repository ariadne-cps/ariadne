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
 
#include "discretiser.h"

#include "approximate_taylor_model.h"
#include "orbit.h"
#include "grid_set.h"
#include "hybrid_set.h"
#include "hybrid_automaton.h"

namespace Ariadne {

typedef ApproximateTaylorModel DefaultModelType;
typedef ApproximateTaylorModel DefaultEnclosureType;
typedef std::pair<DiscreteState,DefaultEnclosureType> DefaultHybridEnclosureType;

template<class ES>
Orbit<typename HybridDiscretiser<ES>::BasicSetType> 
HybridDiscretiser<ES>::
lower_evolve(const SystemType& system, 
             const BasicSetType& initial_set, 
             const TimeType& time) const
{
  return this->_discretise(this->_evolver->orbit(system,this->_enclosure(initial_set),time,LOWER_SEMANTICS));
}

template<class ES>
Orbit<typename HybridDiscretiser<ES>::BasicSetType> 
HybridDiscretiser<ES>::
upper_evolve(const SystemType& system, 
             const BasicSetType& initial_set, 
             const TimeType& time) const
{
  return this->_discretise(this->_evolver->orbit(system,this->_enclosure(initial_set),time,UPPER_SEMANTICS));
}

template<class ES>
typename HybridDiscretiser<ES>::EnclosureType 
HybridDiscretiser<ES>::
_enclosure(const BasicSetType& initial_set) const
{
  return EnclosureType(initial_set.first,ES(initial_set.second.box()));
}

template<class ES>
Orbit<typename HybridDiscretiser<ES>::BasicSetType> 
HybridDiscretiser<ES>::
_discretise(const Orbit<EnclosureType>& orbit) const
{
  ARIADNE_NOT_IMPLEMENTED;
}

template class HybridDiscretiser<DefaultEnclosureType>;

} // namespace Ariadne
