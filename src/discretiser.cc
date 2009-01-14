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
            //            result.adjoin_outer_approximation(*iter,accuracy);
            result.adjoin_outer_approximation(ImageSet(iter->range()),accuracy);
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
Orbit<typename ContinuousDiscretiser<ES>::BasicSetType> 
ContinuousDiscretiser<ES>::
lower_evolution(const SystemType& system, 
                const BasicSetType& initial_set, 
                const TimeType& time,
                const int accuracy) const
{
    return this->_discretise(this->_evolver->orbit(system,this->_enclosure(initial_set),time,LOWER_SEMANTICS),initial_set,accuracy);
}

template<class ES>
Orbit<typename ContinuousDiscretiser<ES>::BasicSetType> 
ContinuousDiscretiser<ES>::
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
typename ContinuousDiscretiser<ES>::EnclosureType 
ContinuousDiscretiser<ES>::
_enclosure(const BasicSetType& initial_set) const
{
    return EnclosureType(initial_set.box());
}

template<class ES>
Orbit<typename ContinuousDiscretiser<ES>::BasicSetType> 
ContinuousDiscretiser<ES>::
_discretise(const Orbit<EnclosureType>& continuous_orbit,
            const BasicSetType& initial_set,
            const int accuracy) const
{
    ARIADNE_LOG(3,ARIADNE_PRETTY_FUNCTION<<"\n");
    ARIADNE_LOG(6,"continuous_reach_set="<<continuous_orbit.reach()<<"\n");
    DenotableSetType reach_set
        = outer_approximation(continuous_orbit.reach(),
                              Grid(continuous_orbit.reach().dimension(),Float(1)),
                              accuracy);
    ARIADNE_LOG(4,"reach_set="<<reach_set<<"\n");
    DenotableSetType intermediate_set
        = outer_approximation(continuous_orbit.intermediate(),
                              Grid(continuous_orbit.intermediate().dimension(),Float(1)),
                              accuracy);
    ARIADNE_LOG(4,"intermediate_set="<<intermediate_set<<"\n");
    DenotableSetType final_set
        = outer_approximation(continuous_orbit.final(),
                              Grid(continuous_orbit.final().dimension(),Float(1)),
                              accuracy);
    ARIADNE_LOG(4,"final_set="<<final_set<<"\n");
    return Orbit<BasicSetType>(initial_set,reach_set,intermediate_set,final_set);
 
}

template class ContinuousDiscretiser<DefaultEnclosureType>;





template<class ES>
Orbit<typename HybridDiscretiser<ES>::BasicSetType> 
HybridDiscretiser<ES>::
lower_evolution(const SystemType& system, 
                const BasicSetType& initial_set, 
                const TimeType& time,
                const int accuracy) const
{
    return this->_discretise(this->_evolver->orbit(system,this->_enclosure(initial_set),time,LOWER_SEMANTICS),initial_set,accuracy);
}

template<class ES>
Orbit<typename HybridDiscretiser<ES>::BasicSetType> 
HybridDiscretiser<ES>::
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
typename HybridDiscretiser<ES>::EnclosureType 
HybridDiscretiser<ES>::
_enclosure(const BasicSetType& initial_set) const
{
    return EnclosureType(initial_set.first,ES(initial_set.second.box()));
}

template<class ES>
Orbit<typename HybridDiscretiser<ES>::BasicSetType> 
HybridDiscretiser<ES>::
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

template class HybridDiscretiser<DefaultEnclosureType>;

} // namespace Ariadne


typedef unsigned char uchar;

namespace Ariadne {

Cell make_cell(const GridCell& gc) {
    const uint max_height=112;
    const int max_depth=112;
    const uint max_size=112;
    assert(gc.height()<max_height);
    assert(gc.depth()<max_depth);
    assert(gc.height()+gc.depth()<max_size);
    
    std::bitset<112> wd;
    for(uchar i=0; i!=gc.word().size(); ++i) {
        wd[i]=gc.word()[i];
    }
    Cell c= { gc.height(),gc.depth(),wd };
    return c;
}

std::ostream& operator<<(std::ostream& os, const Cell& c) {
    os<<"("<<int(c.height)<<","<<int(c.depth)<<",";
    for(uint i=0; i!=c.height+c.depth; ++i) { os<<c.word[i]; }
    return os<<")";
}

GridCell make_grid_cell(const Cell& c, const Grid& g) {
    BinaryWord wd;
    for(uchar i=0; i!=wd.size(); ++i) {
        wd.push_back(c.word[i]);
    }
    return GridCell(g,c.height,wd);
}

std::vector<Cell> 
successor(const DiscretiserInterface<VectorField,GridCell>& discretiser, 
          const Grid& grid, const VectorField& system, const Cell& cell, Float time, char reachevolve)
{
    GridCell grid_cell=make_grid_cell(cell,grid);
    //std::cerr<<"initial_grid_cell="<<grid_cell<<"\n";
    GridTreeSet grid_tree_set;
    ARIADNE_ASSERT(reachevolve=='r' || reachevolve=='e');
    Orbit<GridCell> orbit=discretiser.upper_evolution(system,grid_cell,time,grid_cell.depth());
    //std::cerr<<"orbit="<<orbit<<"\n\n";
    if(reachevolve=='r') {
        grid_tree_set=orbit.reach();
    } else if(reachevolve=='e') {
        grid_tree_set=orbit.final();
    }
    grid_tree_set.mince(grid_cell.depth());
    std::vector<Cell> cells;
    for(GridTreeSet::const_iterator iter=grid_tree_set.begin();
        iter!=grid_tree_set.end(); ++iter)
    {
        cells.push_back(make_cell(*iter));
    }
    return cells;
}

} // namespace Ariadne


