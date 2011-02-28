/***************************************************************************
 *            hybrid_discretiser.cc
 *
 *  Copyright  2006-11  Alberto Casagrande, Pieter Collins
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

#include "hybrid_discretiser.h"

#include "stlio.h"
#include "taylor_set.h"
#include "function_set.h"
#include "grid_set.h"
#include "hybrid_evolver.h"
#include "hybrid_set.h"
#include "hybrid_orbit.h"
#include "hybrid_space.h"

#include "hybrid_automaton_interface.h"

namespace Ariadne {

inline Sweeper default_sweeper() { return Sweeper(); }

template<class ES>
HybridGridTreeSet
//outer_approximation(const HybridListSet<ES>& hls,
outer_approximation(const ListSet<HybridBasicSet<ES> >& hls,
                    const HybridGrid& hgr,
                    const int accuracy)
{
    HybridGridTreeSet result(hgr);
    for(typename HybridListSet<ES>::const_iterator
            iter=hls.begin(); iter!=hls.end(); ++iter)
        {
            DiscreteLocation loc=iter->first;
            const ES& es=iter->second;
            if(!result.has_location(loc)) {
                result.insert(make_pair(loc,GridTreeSet(hgr[loc])));
            }
            GridTreeSet& gts=result[loc];
            gts.adjoin_outer_approximation(ImageSet(es.range()),accuracy);
            //gts.adjoin_outer_approximation(ModelSet<ES>(es),accuracy);
        }
    return result;
}




template<>
HybridGridTreeSet
outer_approximation(const ListSet< HybridBasicSet<TaylorConstrainedImageSet> >& hls,
//outer_approximation(const HybridListSet< TaylorConstrainedImageSet >& hls,
                    const HybridGrid& grid,
                    const int accuracy)
{
    typedef TaylorConstrainedImageSet ES;
    HybridGridTreeSet result(grid);
    //for(typename HybridListSet<ES>::const_iterator
    for(ListSet< HybridBasicSet<ES> >::const_iterator
            iter=hls.begin(); iter!=hls.end(); ++iter)
        {
            HybridBasicSet<ES> hes(*iter);
            DiscreteLocation loc=hes.location();
            const ES& es=hes.continuous_state_set();
            GridTreeSet& gts=result[loc];
            gts.adjoin(es.outer_approximation(gts.grid(),accuracy));
            //gts.adjoin_outer_approximation(ModelSet<ES>(es),accuracy);
        }
    return result;
}


template<class HES>
HybridDiscretiser<HES>::
HybridDiscretiser(const HybridEvolverInterface& evolver)
    : _evolver_ptr(evolver.clone()), _scaling_ptr(new HybridScaling)
{
}

template<class HES>
void
HybridDiscretiser<HES>::
set_scaling(const HybridScalingInterface& scaling)
{
    this->_scaling_ptr=shared_ptr<HybridScalingInterface>(scaling.clone());
}


template<class HES>
Orbit<HybridGridCell>
HybridDiscretiser<HES>::
evolution(const SystemType& system,
          const BasicSetType& initial_set,
          const TimeType& time,
          const AccuracyType accuracy,
          const Semantics semantics) const
{
    ARIADNE_LOG(3,ARIADNE_PRETTY_FUNCTION<<"\n");
    EnclosureType enclosure=this->_enclosure(initial_set);
    ARIADNE_LOG(4,"enclosure"<<enclosure<<"\n");
    Orbit<EnclosureType> continuous_orbit=this->_evolver_ptr->orbit(system,enclosure,time,semantics);
    ARIADNE_LOG(5,"continuous_orbit reach size="<<continuous_orbit.reach().size()<<"\n");
    ARIADNE_LOG(5,"continuous_orbit final size="<<continuous_orbit.final().size()<<"\n");
    HybridGrid hgrid=this->_hybrid_grid(system);
    ARIADNE_LOG(5,"hybrid_grid="<<hgrid<<"\n");
    Orbit<BasicSetType> discrete_orbit=this->_discretise(continuous_orbit,initial_set,hgrid,accuracy);
    ARIADNE_LOG(5,"discrete_orbit reach size="<<discrete_orbit.reach().size()<<"\n");
    ARIADNE_LOG(5,"discrete_orbit final size="<<discrete_orbit.final().size()<<"\n");
    return discrete_orbit;
}


template<class HES>
HybridGridTreeSet
HybridDiscretiser<HES>::
reach(const SystemType& system,
            const BasicSetType& initial_set,
            const TimeType& time,
            const AccuracyType accuracy,
            const Semantics semantics) const
{
    return this->_discretise(this->_evolver_ptr->reach(system,this->_enclosure(initial_set),time,semantics),
                             initial_set,this->_hybrid_grid(system),accuracy);
}

template<class HES>
HybridGridTreeSet
HybridDiscretiser<HES>::
evolve(const SystemType& system,
             const BasicSetType& initial_set,
             const TimeType& time,
             const AccuracyType accuracy,
             const Semantics semantics) const
{
    EnclosureType initial_enclosure=this->_enclosure(initial_set);
    ListSet<EnclosureType> final_enclosures=this->_evolver_ptr->evolve(system,initial_enclosure,time,semantics);
    HybridGrid grid=this->_hybrid_grid(system);
    return this->_discretise(final_enclosures,initial_set,grid,accuracy);
}

template<class HES>
Orbit<HybridGridCell>
HybridDiscretiser<HES>::
lower_evolution(const SystemType& system,
                const BasicSetType& initial_set,
                const TimeType& time,
                const AccuracyType accuracy) const
{
    return this->evolution(system, initial_set, time, accuracy, LOWER_SEMANTICS);
}

template<class HES>
Orbit<HybridGridCell>
HybridDiscretiser<HES>::
upper_evolution(const SystemType& system,
                const BasicSetType& initial_set,
                const TimeType& time,
                const AccuracyType accuracy) const
{
    return this->evolution(system, initial_set, time, accuracy, UPPER_SEMANTICS);
}


template<class HES>
typename HybridDiscretiser<HES>::EnclosureType
HybridDiscretiser<HES>::
_enclosure(const BasicSetType& initial_set) const
{
    return this->_evolver_ptr->enclosure(initial_set.box());
    ContinuousEnclosureType continuous_enclosure(initial_set.second.box(),default_sweeper());
    return EnclosureType(initial_set.first,continuous_enclosure);
}

template<class HES>
Orbit<typename HybridDiscretiser<HES>::BasicSetType>
HybridDiscretiser<HES>::
_discretise(const Orbit<EnclosureType>& continuous_orbit,
            const BasicSetType& initial,
            const HybridGrid& grid,
            const int accuracy) const
{
    ARIADNE_LOG(3,"HybridDiscretiser<HES>::_discretise(...)"<<"\n");
    ARIADNE_LOG(6,"continuous_orbit="<<continuous_orbit<<"\n");
    ARIADNE_LOG(6,"initial="<<initial<<"\n");
    ARIADNE_LOG(6,"grid="<<grid<<"\n");
    DenotableSetType initial_set(grid);
    initial_set.adjoin(initial);
    ARIADNE_LOG(4,"initial_set size="<<initial_set.size()<<"\n");
    ARIADNE_LOG(6,"initial_set="<<initial_set<<"\n");
    DenotableSetType reach_set
        = outer_approximation(continuous_orbit.reach(),grid,accuracy);
    ARIADNE_LOG(4,"reach_set size="<<reach_set.size()<<"\n");
    ARIADNE_LOG(6,"reach_set="<<reach_set<<"\n");
    DenotableSetType intermediate_set
        = outer_approximation(continuous_orbit.intermediate(),grid,accuracy);
    ARIADNE_LOG(4,"intermediate_set size="<<intermediate_set.size()<<"\n");
    ARIADNE_LOG(6,"intermediate_set="<<intermediate_set<<"\n");
    DenotableSetType final_set
        = outer_approximation(continuous_orbit.final(),grid,accuracy);
    ARIADNE_LOG(4,"final_set size="<<final_set.size()<<"\n");
    ARIADNE_LOG(6,"final_set="<<final_set<<"\n");
    return Orbit<BasicSetType>(initial_set,reach_set,intermediate_set,final_set);

}

template<class HES>
HybridGridTreeSet
HybridDiscretiser<HES>::
_discretise(const ListSet<EnclosureType>& enclosure_list_set,
            const BasicSetType& initial_set,
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


template<class HES>
HybridGrid
HybridDiscretiser<HES>::
_hybrid_grid(const SystemType& system) const
{
    return HybridGrid(system.state_space(),*this->_scaling_ptr);
}



template class HybridDiscretiser<HybridEnclosure>;

} // namespace Ariadne


