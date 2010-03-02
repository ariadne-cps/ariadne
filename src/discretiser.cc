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

#include "taylor_set.h"
#include "orbit.h"
#include "function_set.h"
#include "polytope.h"
#include "grid_set.h"
#include "hybrid_set.h"
#include "map.h"
#include "hybrid_automaton.h"
#include "hybrid_evolver-constrained.h"

namespace Ariadne {


template<class ES>
HybridGridTreeSet 
outer_approximation(const ListSet<HybridBasicSet<ES> >& hls,
                    const HybridGrid& hgr,
                    const int accuracy)
{
    HybridGridTreeSet result(hgr);
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

template<class ES>
HybridGridTreeSet 
outer_approximation(const HybridBasicSet<ES>& hs,
                    const HybridGrid& hgr,
                    const int accuracy)
{
    HybridGridTreeSet result(hgr);
    DiscreteState loc=hs.location();
    const ES& es=hs.continuous_state_set();
    if(result.find(loc)==result.locations_end()) {
        result.insert(make_pair(loc,GridTreeSet(hgr[loc])));
    }
    GridTreeSet& gts=result[loc];
    gts.adjoin_outer_approximation(ImageSet(es.range()),accuracy);
    return result;
}

template<>
HybridGridTreeSet 
outer_approximation(const HybridBasicSet<ConstrainedImageSet>& hs,
                    const HybridGrid& hgr,
                    const int accuracy)
{
    HybridGridTreeSet result(hgr);
    DiscreteState loc=hs.location();
    const ConstrainedImageSet& es=hs.continuous_state_set();
    if(result.find(loc)==result.locations_end()) {
        result.insert(make_pair(loc,GridTreeSet(hgr[loc])));
    }
    GridTreeSet& gts=result[loc];
    gts.adjoin(es.outer_approximation(gts.grid(),accuracy));
    return result;
}

template<>
HybridGridTreeSet 
outer_approximation(const ListSet< HybridBasicSet<ConstrainedImageSet> >& hls,
                    const HybridGrid& hgr,
                    const int accuracy)
{
    typedef ConstrainedImageSet ES;
    HybridGridTreeSet result;
    //for(typename HybridListSet<ES>::const_iterator 
    for(ListSet< HybridBasicSet<ES> >::const_iterator
            iter=hls.begin(); iter!=hls.end(); ++iter)
        {
            HybridBasicSet<ConstrainedImageSet> hes(*iter);
            DiscreteState loc=hes.location();
            const ES& es=hes.continuous_state_set();
            if(result.find(loc)==result.locations_end()) {
                Grid grid=hgr[loc];
                result.insert(make_pair(loc,GridTreeSet(grid)));
            }
            GridTreeSet& gts=result[loc];
            gts.adjoin(es.outer_approximation(gts.grid(),accuracy));
            //gts.adjoin_outer_approximation(ModelSet<ES>(es),accuracy);
        }
    return result;
}


//typedef ApproximateTaylorModel DefaultModelType;
//typedef ApproximateTaylorModel DefaultEnclosureType;
//typedef std::pair<DiscreteState,DefaultEnclosureType> DefaultHybridEnclosureType;

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
    return EnclosureType(initial_set.box());
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

template class Discretiser<VectorField,TaylorSet>;
template class Discretiser<IteratedMap,TaylorSet>;





template<class ES>
std::pair<HybridGridTreeSet,HybridGridTreeSet>
HybridDiscretiser<ES>::
evolution(const SystemType& system, 
          const EnclosureType& initial_set, 
          const TimeType& time,
          const AccuracyType accuracy,
          const Semantics semantics) const
{
    ARIADNE_LOG(3,ARIADNE_PRETTY_FUNCTION<<"\n");
    Orbit<EnclosureType> continuous_orbit=this->_evolver->orbit(system,initial_set,time,semantics);
    ARIADNE_LOG(5,"continuous_orbit reach size="<<continuous_orbit.reach().size()<<"\n");
    ARIADNE_LOG(5,"continuous_orbit final size="<<continuous_orbit.final().size()<<"\nOK\n");
    HybridGridTreeSet reach=this->_discretise(continuous_orbit.reach(),system.grid(),accuracy);
    HybridGridTreeSet final=this->_discretise(continuous_orbit.final(),system.grid(),accuracy);
	reach.adjoin(final); // Always adjoin the reached region with the final region (preferable for consistency with _evolver.reach() )
    ARIADNE_LOG(5,"discretised reach size="<<reach.size()<<"\n");
    ARIADNE_LOG(5,"discretised final size="<<final.size()<<"\n");

	ARIADNE_ASSERT_MSG(possibly(final.subset(reach.bounding_box())), "The final region is not included into the reached region!"); // TO BE REMOVED as soon as the "THE FINAL always included in REACH" property is verified

    return make_pair(reach,final);
}


template<class ES>
HybridGridTreeSet 
HybridDiscretiser<ES>::
reach(const SystemType& system, 
            const EnclosureType& initial_set, 
            const TimeType& time,
            const AccuracyType accuracy,
            const Semantics semantics) const
{
    return this->_discretise(this->_evolver->reach(system,initial_set,time,semantics),
                             system.grid(),accuracy);
}

template<class ES>
HybridGridTreeSet 
HybridDiscretiser<ES>::
evolve(const SystemType& system, 
             const EnclosureType& initial_set, 
             const TimeType& time,
             const AccuracyType accuracy,
             const Semantics semantics) const
{
    ListSet<EnclosureType> final_enclosures=this->_evolver->evolve(system,initial_set,time,semantics);
    HybridGrid grid=system.grid();
    return this->_discretise(final_enclosures,grid,accuracy);
}

template<class ES>
std::pair<HybridGridTreeSet,HybridGridTreeSet>
HybridDiscretiser<ES>::
lower_evolution(const SystemType& system, 
                const EnclosureType& initial_set, 
                const TimeType& time,
                const AccuracyType accuracy) const 
{ 
    return this->evolution(system, initial_set, time, accuracy, LOWER_SEMANTICS);
}

template<class ES>
std::pair<HybridGridTreeSet,HybridGridTreeSet>
HybridDiscretiser<ES>::
upper_evolution(const SystemType& system, 
                const EnclosureType& initial_set, 
                const TimeType& time,
                const AccuracyType accuracy) const 
{ 
    return this->evolution(system, initial_set, time, accuracy, UPPER_SEMANTICS);
}


template<class ES>
typename HybridDiscretiser<ES>::EnclosureType 
HybridDiscretiser<ES>::
enclosure(const BasicSetType& initial_set) const
{
    return EnclosureType(initial_set.first,ES(initial_set.second.box()));
}


template<class ES>
HybridGridTreeSet
HybridDiscretiser<ES>::
_discretise(const ListSet<EnclosureType>& enclosure_list_set,
            const HybridGrid& hgrid,
            const int accuracy) const
{
    ARIADNE_LOG(3,ARIADNE_PRETTY_FUNCTION<<"\n");
    ARIADNE_LOG(6,"enclosure_list_set="<<enclosure_list_set<<"\n");

    DenotableSetType discretised_set
        = outer_approximation(enclosure_list_set,hgrid,accuracy);
    ARIADNE_LOG(4,"discretised_set="<<discretised_set<<"\n");
    return discretised_set; 
}


template class HybridDiscretiser<TaylorSet>;
template class HybridDiscretiser<ConstrainedImageSet>;

} // namespace Ariadne


