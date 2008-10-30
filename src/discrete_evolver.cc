/***************************************************************************
 *            discrete_evolver.cc
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
 
#include "discrete_evolver.h"

#include "approximate_taylor_model.h"
#include "orbit.h"
#include "function_set.h"
#include "grid_set.h"
#include "hybrid_set.h"
#include "hybrid_automaton.h"

namespace Ariadne {

template<class ES>
GridTreeSet 
outer_approximation(const ListSet<ES>& ls,
                    const Grid& gr,
                    const uint accuracy)
{
  GridTreeSet result(gr);
  for(typename ListSet<ES>::const_iterator 
        iter=ls.begin(); iter!=ls.end(); ++iter)
  {
    result.adjoin_outer_approximation(*iter,accuracy);
  }
  return result;
}


template<class ES>
HybridGridTreeSet 
outer_approximation(const HybridListSet<ES>& hls,
                    const HybridGrid& hgr,
                    const int accuracy)
{
  HybridGridTreeSet result;
  for(typename HybridListSet<ES>::const_iterator 
        iter=hls.begin(); iter!=hls.end(); ++iter)
  {
    DiscreteState loc=iter->first;
    const ES& es=iter->second;
    if(result.find(loc)==result.locations_end()) {
      result.insert(make_pair(loc,GridTreeSet(hgr[loc])));
    }
    GridTreeSet& gts=result[loc];
    gts.adjoin_outer_approximation(ImageSet(es.range()),accuracy);
    //gts.adjoin_outer_approximation(ModelSet<ES>(es),accuracy);
  }
  return result;
}





typedef ApproximateTaylorModel DefaultModelType;
typedef ApproximateTaylorModel DefaultEnclosureType;
typedef std::pair<DiscreteState,DefaultEnclosureType> DefaultHybridEnclosureType;

template<class ES>
Orbit<typename HybridDiscreteEvolver<ES>::BasicSetType> 
HybridDiscreteEvolver<ES>::
lower_evolution(const SystemType& system, 
                const BasicSetType& initial_set, 
                const TimeType& time,
                const int accuracy) const
{
  return this->_discretise(this->_evolver->orbit(system,this->_enclosure(initial_set),time,LOWER_SEMANTICS),initial_set,accuracy);
}

template<class ES>
Orbit<typename HybridDiscreteEvolver<ES>::BasicSetType> 
HybridDiscreteEvolver<ES>::
upper_evolution(const SystemType& system, 
                const BasicSetType& initial_set, 
                const TimeType& time,
                const int accuracy) const
{
  ARIADNE_LOG(3,ARIADNE_PRETTY_FUNCTION);
  EnclosureType enclosure=this->_enclosure(initial_set);
  ARIADNE_LOG(4,"enclosure="<<enclosure<<"\n");
  Orbit<EnclosureType> continuous_orbit=this->_evolver->orbit(system,enclosure,time,UPPER_SEMANTICS);
  ARIADNE_LOG(5,"continuous_orbit="<<continuous_orbit<<"\nOK\n");
  Orbit<BasicSetType> discrete_orbit=this->_discretise(continuous_orbit,initial_set,accuracy);
  ARIADNE_LOG(5,"discrete_orbit="<<discrete_orbit<<"\n");
  return discrete_orbit;
}

template<class ES>
typename HybridDiscreteEvolver<ES>::EnclosureType 
HybridDiscreteEvolver<ES>::
_enclosure(const BasicSetType& initial_set) const
{
  return EnclosureType(initial_set.first,ES(initial_set.second.box()));
}

template<class ES>
Orbit<typename HybridDiscreteEvolver<ES>::BasicSetType> 
HybridDiscreteEvolver<ES>::
_discretise(const Orbit<EnclosureType>& continuous_orbit,
            const BasicSetType& initial_set,
            const int accuracy) const
{
  ARIADNE_LOG(3,ARIADNE_PRETTY_FUNCTION<<"\n");
  ARIADNE_LOG(6,"continuous_reach_set="<<continuous_orbit.reach()<<"\n");
  DenotableSetType reach_set
    = outer_approximation(continuous_orbit.reach(),
                          HybridGrid(continuous_orbit.reach().space(),Float(1)),
                          accuracy);
  ARIADNE_LOG(4,"reach_set="<<reach_set<<"\n");
  DenotableSetType intermediate_set
    = outer_approximation(continuous_orbit.intermediate(),
                          HybridGrid(continuous_orbit.intermediate().space(),Float(1)),
                          accuracy);
  ARIADNE_LOG(4,"intermediate_set="<<intermediate_set<<"\n");
  DenotableSetType final_set
    = outer_approximation(continuous_orbit.final(),
                          HybridGrid(continuous_orbit.final().space(),Float(1)),
                          accuracy);
  ARIADNE_LOG(4,"final_set="<<final_set<<"\n");
  return Orbit<BasicSetType>(initial_set,reach_set,intermediate_set,final_set);
 
}

template class HybridDiscreteEvolver<DefaultEnclosureType>;

} // namespace Ariadne
