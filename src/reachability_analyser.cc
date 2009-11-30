/***************************************************************************
 *            reachability_analyser.cc
 *
 *  Copyright  2006-9  Alberto Casagrande, Pieter Collins, Davide Bresolin
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

#include <string>
#include <sstream>
#include <algorithm>

#include <list>
#include <set>
#include <vector>
#include <valarray>

#include "exceptions.h"

#include "numeric.h"

#include "vector.h"
#include "matrix.h"

#include "box.h"
#include "list_set.h"
#include "grid_set.h"

#include "orbit.h"

#include "hybrid_time.h"
#include "hybrid_automaton.h"

#include "evolution_parameters.h"
#include "evolver_interface.h"

#include "discretiser.h"
#include "reachability_analyser.h"
#include "logging.h"

#include "graphics.h"


namespace Ariadne {


static const double DEFAULT_MAXIMUM_ENCLOSURE_RADIUS=0.25;
static const double DEFAULT_GRID_LENGTH=0.125;

HybridReachabilityAnalyser::
~HybridReachabilityAnalyser()
{
}


HybridReachabilityAnalyser::
HybridReachabilityAnalyser(const HybridDiscretiser<HybridEvolver::ContinuousEnclosureType>& discretiser)
    : _parameters(new EvolutionParametersType())
    , _discretiser(discretiser.clone())
{
}






// Helper functions for operators on lists of sets.
HybridGridTreeSet
HybridReachabilityAnalyser::_upper_reach(const HybridAutomaton& sys,
                                         const HybridGridTreeSet& set,
                                         const HybridTime& time,
                                         const int accuracy) const
{
    Gr grid=sys.grid();
    HybridGridTreeSet result(grid);
    HybridGridTreeSet cells=set;
    cells.mince(accuracy);
    for(HybridGridTreeSet::const_iterator iter=cells.begin(); iter!=cells.end(); ++iter) {
        EnclosureType enclosure=this->_discretiser->enclosure(*iter);
        result.adjoin(this->_discretiser->reach(sys,enclosure,time,accuracy,UPPER_SEMANTICS));
    }
    return result;
}


HybridGridTreeSet
HybridReachabilityAnalyser::_upper_evolve(const HybridAutomaton& sys,
                                          const HybridGridTreeSet& set,
                                          const HybridTime& time,
                                          const int accuracy) const
{
    Gr grid=sys.grid();
    GTS result(grid); GTS cells=set; cells.mince(accuracy);
    for(HybridGridTreeSet::const_iterator iter=cells.begin(); iter!=cells.end(); ++iter) {
        ARIADNE_LOG(5,"Evolving cell = "<<*iter<<"\n");
        EnclosureType enclosure=this->_discretiser->enclosure(*iter);
        result.adjoin(this->_discretiser->evolve(sys,enclosure,time,accuracy,UPPER_SEMANTICS));
    }
    ARIADNE_LOG(4,"_upper_evolve result size = "<<result.size()<<"\n");
    return result;
}


std::pair<HybridGridTreeSet,HybridGridTreeSet>
HybridReachabilityAnalyser::_upper_reach_evolve(const HybridAutomaton& sys,
                                                const HybridGridTreeSet& set,
                                                const HybridTime& time,
                                                const int accuracy) const
{
    ARIADNE_LOG(2,"HybridReachabilityAnalyser::_upper_reach_evolve(...)\n");
    Gr grid=sys.grid();
    std::pair<GTS,GTS> result=make_pair(GTS(grid),GTS(grid));
    GTS& reach=result.first; GTS& evolve=result.second;
    GTS cells=set; cells.mince(accuracy);

    for(HybridGridTreeSet::const_iterator iter=cells.begin(); iter!=cells.end(); ++iter) {
        EnclosureType enclosure=this->_discretiser->enclosure(*iter);
        GTS cell_reach, cell_final;
        make_lpair(cell_reach,cell_final)=this->_discretiser->evolution(sys,enclosure,time,accuracy,UPPER_SEMANTICS);
        ARIADNE_LOG(4,"  evolution reach size= "<<cell_reach.size()<<"\n");
        ARIADNE_LOG(4,"  evolution final size= "<<cell_final.size()<<"\n");
        reach.adjoin(cell_reach);
        evolve.adjoin(cell_final);
    }
    ARIADNE_LOG(3,"  final reach size = "<<reach.size()<<"\n");
    ARIADNE_LOG(3,"  final evolve size = "<<evolve.size()<<"\n");
    ARIADNE_LOG(2,"Done.\n");
    return result;
}


HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
lower_evolve(const SystemType& system,
             const HybridImageSet& initial_set,
             const TimeType& time) const
{
    ARIADNE_LOG(2,"HybridReachabilityAnalyser::lower_evolve(...)\n");
    int grid_depth = this->_parameters->maximum_grid_depth;
    int grid_height = this->_parameters->maximum_grid_height;
    Gr grid(system.grid());
    GTS initial; GTS final;

    // For each location, test if the radius of the set is smaller than the grid cell
    for(HybridImageSet::locations_const_iterator loc_iter=initial_set.locations_begin();
        loc_iter!=initial_set.locations_end(); ++loc_iter) 
    {
        Vector<Float> cell = grid[loc_iter->first].lengths();
        Float cell_radius = (min(cell))/(1 << (grid_depth+4));
        Vector<Float> origin = grid[loc_iter->first].origin();
        Box bbox = loc_iter->second.bounding_box();
        if (radius(bbox) > cell_radius) {
            // if bigger, map to the grid
            // First of all, test if the bounding box lies on cell boundaries or not
            for(uint i = 0 ; i < bbox.dimension() ; ++i) {
                // test if the i-th dimension is a singleton interval AND
                // if it lies on a cell boundary
                cell_radius = cell[i]/(1 << (grid_depth+4));
                Float intpart, fractpart;
                fractpart = modf((bbox[i].lower()-origin[i])/cell_radius,&intpart);
                if(bbox[i].singleton() && (fractpart == 0.0)) {
                    // the set lies on the boundary, shift the grid center by half cell size
                    origin[i]+=cell_radius/2.0;
                    grid[loc_iter->first].set_origin(origin);
                }
            }        
            // if bigger, map to the grid
            ARIADNE_LOG(3,"Adjoining initial set for location "<<loc_iter->first<<" to the grid...\n");
            initial.insert(make_pair(loc_iter->first,grid[loc_iter->first]));
            initial[loc_iter->first].adjoin_lower_approximation(loc_iter->second,grid_height,grid_depth+4);
        } else {
            ARIADNE_LOG(3,"Computing evolution for initial set in location "<<loc_iter->first<<" directly...\n");
            // if smaller, compute the evolution directly
            EnclosureType initial_enclosure(loc_iter->first,ContinuousEnclosureType(loc_iter->second));
            GTS cell_final=this->_discretiser->evolve(system,initial_enclosure,time,grid_depth,LOWER_SEMANTICS);
            final.adjoin(cell_final);            
            initial.insert(make_pair(loc_iter->first,grid[loc_iter->first]));
        }
    }

    ARIADNE_LOG(4,"grid="<<grid<<"\n");

    ARIADNE_LOG(3,"initial.size()="<<initial.size()<<"\n");
    if(!initial.empty()) {
        ARIADNE_LOG(3,"computing lower evolution from the grid.");
        for(GTS::const_iterator bs_iter=initial.begin(); bs_iter!=initial.end(); ++bs_iter) {
            ARIADNE_LOG(3,".");
            EnclosureType enclosure=this->_discretiser->enclosure(*bs_iter);
            GTS cell_final=this->_discretiser->evolve(system,enclosure,time,grid_depth,LOWER_SEMANTICS);
            final.adjoin(cell_final);
        }
    }
    ARIADNE_LOG(3,"\n");
    return final;
}


HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
lower_reach(const SystemType& system,
            const HybridImageSet& initial_set,
            const TimeType& time) const
{
    ARIADNE_LOG(2,"HybridReachabilityAnalyser::lower_reach(...)\n");
    int grid_depth = this->_parameters->maximum_grid_depth;
    int grid_height = this->_parameters->maximum_grid_height;
    Gr grid(system.grid());
    GTS initial; GTS reach;

    // For each location, test if the radius of the set is smaller than the grid cell
    for(HybridImageSet::locations_const_iterator loc_iter=initial_set.locations_begin();
        loc_iter!=initial_set.locations_end(); ++loc_iter) 
    {
        Vector<Float> cell = grid[loc_iter->first].lengths();
        Float cell_radius = (min(cell))/(1 << (grid_depth+4));
        Vector<Float> origin = grid[loc_iter->first].origin();
        Box bbox = loc_iter->second.bounding_box();
        if (radius(bbox) > cell_radius) {
            // if bigger, map to the grid
            // First of all, test if the bounding box lies on cell boundaries or not
            for(uint i = 0 ; i < bbox.dimension() ; ++i) {
                // test if the i-th dimension is a singleton interval AND
                // if it lies on a cell boundary
                cell_radius = cell[i]/(1 << (grid_depth+4));
                Float intpart, fractpart;
                fractpart = modf((bbox[i].lower()-origin[i])/cell_radius,&intpart);
                if(bbox[i].singleton() && (fractpart == 0.0)) {
                    // the set lies on the boundary, shift the grid center by half cell size
                    origin[i]+=cell_radius/2.0;
                    grid[loc_iter->first].set_origin(origin);
                }
            }        
            // if bigger, map to the grid
            ARIADNE_LOG(3,"Adjoining initial set for location "<<loc_iter->first<<" to the grid...\n");
            initial.insert(make_pair(loc_iter->first,grid[loc_iter->first]));
            initial[loc_iter->first].adjoin_lower_approximation(loc_iter->second,grid_height,grid_depth+4);
        } else {
            ARIADNE_LOG(3,"Computing evolution for initial set in location "<<loc_iter->first<<" directly...\n");
            // if smaller, compute the evolution directly
            EnclosureType initial_enclosure(loc_iter->first,ContinuousEnclosureType(loc_iter->second));
            GTS cell_reach=this->_discretiser->reach(system,initial_enclosure,time,grid_depth,LOWER_SEMANTICS);
            reach.adjoin(cell_reach);            
            initial.insert(make_pair(loc_iter->first,grid[loc_iter->first]));
        }
    }

    ARIADNE_LOG(4,"grid="<<grid<<"\n");

    ARIADNE_LOG(3,"initial.size()="<<initial.size()<<"\n");
    if(!initial.empty()) {
        ARIADNE_LOG(3,"Computing lower reach set from the grid...");
        for(GTS::const_iterator bs_iter=initial.begin(); bs_iter!=initial.end(); ++bs_iter) {
            ARIADNE_LOG(3,".");
            EnclosureType enclosure=this->_discretiser->enclosure(*bs_iter);
            GTS cell_reach = this->_discretiser->reach(system,enclosure,time,grid_depth,LOWER_SEMANTICS);
            reach.adjoin(cell_reach);
        }
    }
    ARIADNE_LOG(3,"\n");
    return reach;
}



std::pair<HybridReachabilityAnalyser::SetApproximationType,HybridReachabilityAnalyser::SetApproximationType>
HybridReachabilityAnalyser::
lower_reach_evolve(const SystemType& system,
                   const HybridImageSet& initial_set,
                   const TimeType& time) const
{
    ARIADNE_LOG(2,"HybridReachabilityAnalyser::lower_reach_evolve(...)\n");
    int grid_depth = this->_parameters->maximum_grid_depth;
    int grid_height = this->_parameters->maximum_grid_height;

    Gr grid=system.grid();
    GTS initial;
    GTS reach; GTS evolve;

    // For each location, test if the radius of the set is smaller than the grid cell
    for(HybridImageSet::locations_const_iterator loc_iter=initial_set.locations_begin();
        loc_iter!=initial_set.locations_end(); ++loc_iter) 
    {
        Vector<Float> cell = grid[loc_iter->first].lengths();
        Float cell_radius = (min(cell))/(1 << (grid_depth+4));
        Vector<Float> origin = grid[loc_iter->first].origin();
        Box bbox = loc_iter->second.bounding_box();
        if (radius(bbox) > cell_radius) {
            // if bigger, map to the grid
            // First of all, test if the bounding box lies on cell boundaries or not
            for(uint i = 0 ; i < bbox.dimension() ; ++i) {
                // test if the i-th dimension is a singleton interval AND
                // if it lies on a cell boundary
                cell_radius = cell[i]/(1 << (grid_depth+4));
                Float intpart, fractpart;
                fractpart = modf((bbox[i].lower()-origin[i])/cell_radius,&intpart);
                if(bbox[i].singleton() && (fractpart == 0.0)) {
                    // the set lies on the boundary, shift the grid center by half cell size
                    origin[i]+=cell_radius/2.0;
                    grid[loc_iter->first].set_origin(origin);
                }
            }        
            // if bigger, map to the grid
            ARIADNE_LOG(3,"Adjoining initial set for location "<<loc_iter->first<<" to the grid...\n");
            initial.insert(make_pair(loc_iter->first,grid[loc_iter->first]));
            initial[loc_iter->first].adjoin_lower_approximation(loc_iter->second,grid_height,grid_depth+4);
        } else {
            ARIADNE_LOG(3,"Computing evolution for initial set in location "<<loc_iter->first<<" directly...\n");
            // if smaller, compute the evolution directly
            EnclosureType initial_enclosure(loc_iter->first,ContinuousEnclosureType(loc_iter->second));
            GTS cell_reach, cell_final;
            make_lpair(cell_reach,cell_final)=this->_discretiser->evolution(system,initial_enclosure,time,grid_depth,LOWER_SEMANTICS);
            reach.adjoin(cell_reach);
            evolve.adjoin(cell_final);
            initial.insert(make_pair(loc_iter->first,grid[loc_iter->first]));           
        }
    }

    ARIADNE_LOG(4,"grid="<<grid<<"\n");

    ARIADNE_LOG(3,"initial.size()="<<initial.size()<<"\n");
    if(!initial.empty()) {
        ARIADNE_LOG(3,"computing lower evolution from the grid.");
        for(GTS::const_iterator bs_iter=initial.begin(); bs_iter!=initial.end(); ++bs_iter) {
            ARIADNE_LOG(3,".");
            EnclosureType enclosure=this->_discretiser->enclosure(*bs_iter);
            GTS cell_reach,cell_final;
            make_lpair(cell_reach,cell_final) = this->_discretiser->evolution(system,enclosure,time,grid_depth,LOWER_SEMANTICS);
            reach.adjoin(cell_reach);
            evolve.adjoin(cell_final);
        }
    }
    ARIADNE_LOG(3,"\n");
    return make_pair(reach,evolve);
}


HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
upper_evolve(const SystemType& system,
             const HybridImageSet& initial_set,
             const TimeType& time) const
{
    ARIADNE_LOG(2,"HybridReachabilityAnalyser::upper_evolve(...)\n");
    Gr grid=system.grid();
    GTS evolve(grid), initial(grid);
    int grid_depth = this->_parameters->maximum_grid_depth;
    Float real_time=time.continuous_time();
    uint discrete_steps=time.discrete_time();
    Float lock_to_grid_time=this->_parameters->lock_to_grid_time;
    uint time_steps=uint(real_time/lock_to_grid_time);
    Float remainder_time=real_time-time_steps*lock_to_grid_time;
    if(time_steps == 0) {
        time_steps=1;
        remainder_time=0.0;
        lock_to_grid_time=real_time;
    }
    HybridTime hybrid_lock_to_grid_time(lock_to_grid_time,discrete_steps);
    HybridTime hybrid_remainder_time(remainder_time,discrete_steps);
    ARIADNE_LOG(3,"real_time="<<real_time<<"\n");
    ARIADNE_LOG(3,"time_steps="<<time_steps<<"  lock_to_grid_time="<<lock_to_grid_time<<"\n");
    ARIADNE_LOG(3,"computing first reachability step...\n");
    // For each location, test if the radius of the set is smaller than the grid cell
    for(HybridImageSet::locations_const_iterator loc_iter=initial_set.locations_begin();
        loc_iter!=initial_set.locations_end(); ++loc_iter) 
    {
        Vector<Float> cell = grid[loc_iter->first].lengths();
        Float cell_radius = (min(cell))/(1 << (grid_depth));
        if (radius(loc_iter->second.bounding_box()) > cell_radius) {
            // if bigger, map to the grid
            ARIADNE_LOG(3,"Adjoining initial set for location "<<loc_iter->first<<" to the grid...\n");
            initial[loc_iter->first].adjoin_outer_approximation(loc_iter->second,grid_depth);
        } else {
            ARIADNE_LOG(3,"Computing evolution for initial set in location "<<loc_iter->first<<" directly...\n");
            // if smaller, compute the evolution directly
            EnclosureType initial_enclosure(loc_iter->first,ContinuousEnclosureType(loc_iter->second));
            GTS cell_final=this->_discretiser->evolve(system,initial_enclosure,time,grid_depth,UPPER_SEMANTICS);
            evolve.adjoin(cell_final);
        }
    }
    if(!initial.empty()) {
        ARIADNE_LOG(3,"computing evolution from the grid...\n");
        ARIADNE_LOG(4,"initial_evolve.size()="<<initial.size()<<"\n");
        initial=this->_upper_evolve(system,initial,hybrid_lock_to_grid_time,grid_depth);
        evolve.adjoin(initial);
    }
    // time steps evolution loop
    for(uint i=1; i<time_steps; ++i) {
        ARIADNE_LOG(3,"computing "<<i+1<<"-th reachability step...\n");
        evolve=this->_upper_evolve(system,evolve,hybrid_lock_to_grid_time,grid_depth);
    }
    ARIADNE_LOG(3,"remainder_time="<<remainder_time<<"\n");
    if(!evolve.empty() && remainder_time > 0) {
        ARIADNE_LOG(3,"computing evolution for remainder time...\n");
        evolve=this->_upper_evolve(system,evolve,hybrid_remainder_time,grid_depth);
    }
    evolve.recombine();
    ARIADNE_LOG(4,"final_evolve.size()="<<evolve.size()<<"\n");
    return evolve;
}



HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
upper_reach(const SystemType& system,
            const HybridImageSet& initial_set,
            const TimeType& time) const
{
    ARIADNE_LOG(2,"HybridReachabilityAnalyser::upper_reach(system,set,time)\n");
    ARIADNE_LOG(4,"initial_set="<<initial_set<<"\n");
    Gr grid=system.grid();
    GTS found(grid),evolve(grid),reach(grid),initial(grid);
    int grid_depth = this->_parameters->maximum_grid_depth;
    Float real_time=time.continuous_time();
    uint discrete_steps=time.discrete_time();
    Float lock_to_grid_time=this->_parameters->lock_to_grid_time;
    uint time_steps=uint(real_time/lock_to_grid_time);
    Float remainder_time=real_time-time_steps*lock_to_grid_time;
    if(time_steps == 0) {
        time_steps=1;
        remainder_time=0.0;
        lock_to_grid_time=real_time;
    }
    HybridTime hybrid_lock_to_grid_time(lock_to_grid_time,discrete_steps);
    HybridTime hybrid_remainder_time(remainder_time,discrete_steps);
    ARIADNE_LOG(3,"real_time="<<real_time<<"\n");
    ARIADNE_LOG(3,"time_steps="<<time_steps<<"  lock_to_grid_time="<<lock_to_grid_time<<"\n");
    ARIADNE_LOG(3,"computing first reachability step...\n");
    // For each location, test if the radius of the set is smaller than the grid cell
    for(HybridImageSet::locations_const_iterator loc_iter=initial_set.locations_begin();
        loc_iter!=initial_set.locations_end(); ++loc_iter) 
    {
        Vector<Float> cell = grid[loc_iter->first].lengths();
        Float cell_radius = (min(cell))/(1 << (grid_depth));
        if (radius(loc_iter->second.bounding_box()) > cell_radius) {
            // if bigger, map to the grid
            ARIADNE_LOG(3,"Adjoining initial set for location "<<loc_iter->first<<" to the grid...\n");
            initial[loc_iter->first].adjoin_outer_approximation(loc_iter->second,grid_depth);
        } else {
            ARIADNE_LOG(3,"Computing evolution for initial set in location "<<loc_iter->first<<" directly...\n");
            // if smaller, compute the evolution directly
            EnclosureType initial_enclosure(loc_iter->first,ContinuousEnclosureType(loc_iter->second));
            GTS cell_reach,cell_final;
            make_lpair(cell_reach,cell_final)=this->_discretiser->evolution(system,initial_enclosure,time,grid_depth,UPPER_SEMANTICS);
            reach.adjoin(cell_reach);
            evolve.adjoin(cell_final);
        }
    }
    if(!initial.empty()) {
        ARIADNE_LOG(3,"computing evolution from the grid...\n");
        ARIADNE_LOG(4,"initial_evolve.size()="<<initial.size()<<"\n");
        make_lpair(found,initial)=this->_upper_reach_evolve(system,initial,hybrid_lock_to_grid_time,grid_depth);
        evolve.adjoin(initial);
        reach.adjoin(found);
    }
    // time steps evolution loop        
    for(uint i=1; i<time_steps; ++i) {
        ARIADNE_LOG(3,"computing "<<i+1<<"-th reachability step...\n");
        make_lpair(found,evolve)=this->_upper_reach_evolve(system,evolve,hybrid_lock_to_grid_time,grid_depth);
        ARIADNE_LOG(5,"found.size()="<<found.size()<<"\n");
        ARIADNE_LOG(5,"evolve.size()="<<evolve.size()<<"\n");
        evolve.remove(reach);
        reach.adjoin(found);
        ARIADNE_LOG(3,"  found "<<found.size()<<" cells, of which "<<evolve.size()<<" are new.\n");
        if(evolve.empty()) break;
        ARIADNE_LOG(6,"evolve="<<evolve<<"\n");
    }
    ARIADNE_LOG(3,"remainder_time="<<remainder_time<<"\n");
    if(!evolve.empty() && remainder_time > 0) {
        ARIADNE_LOG(3,"computing evolution for remainder time...\n");
        make_lpair(found,evolve)=this->_upper_reach_evolve(system,evolve,hybrid_remainder_time,grid_depth);
        reach.adjoin(found);
    }
    // This last step is necessary to add the final set to the result.
    reach.adjoin(evolve);
    reach.recombine();
    ARIADNE_LOG(4,"final_reach size = "<<reach.size()<<"\n");
    return reach;
}



std::pair<HybridReachabilityAnalyser::SetApproximationType,HybridReachabilityAnalyser::SetApproximationType>
HybridReachabilityAnalyser::
upper_reach_evolve(const SystemType& system,
                   const HybridImageSet& initial_set,
                   const TimeType& time) const
{
    ARIADNE_LOG(2,"HybridReachabilityAnalyser::upper_reach_evolve(system,set,time)\n");
    ARIADNE_LOG(4,"initial_set="<<initial_set<<"\n");
    Gr grid=system.grid();
    GTS found(grid),evolve(grid),reach(grid),initial(grid);
    int grid_depth = this->_parameters->maximum_grid_depth;
    Float real_time=time.continuous_time();
    uint discrete_steps=time.discrete_time();
    Float lock_to_grid_time=this->_parameters->lock_to_grid_time;
    uint time_steps=uint(real_time/lock_to_grid_time);
    Float remainder_time=real_time-time_steps*lock_to_grid_time;
    if(time_steps == 0) {
        time_steps=1;
        remainder_time=0.0;
        lock_to_grid_time=real_time;
    }
    HybridTime hybrid_lock_to_grid_time(lock_to_grid_time,discrete_steps);
    HybridTime hybrid_remainder_time(remainder_time,discrete_steps);
    ARIADNE_LOG(3,"real_time="<<real_time<<"\n");
    ARIADNE_LOG(3,"time_steps="<<time_steps<<"  lock_to_grid_time="<<lock_to_grid_time<<"\n");
    ARIADNE_LOG(3,"computing first reachability step...\n");
    // For each location, test if the radius of the set is smaller than the grid cell
    for(HybridImageSet::locations_const_iterator loc_iter=initial_set.locations_begin();
        loc_iter!=initial_set.locations_end(); ++loc_iter) 
    {
        Vector<Float> cell = grid[loc_iter->first].lengths();
        Float cell_radius = (min(cell))/(1 << (grid_depth));
        if (radius(loc_iter->second.bounding_box()) > cell_radius) {
            // if bigger, map to the grid
            ARIADNE_LOG(3,"Adjoining initial set for location "<<loc_iter->first<<" to the grid...\n");
            initial[loc_iter->first].adjoin_outer_approximation(loc_iter->second,grid_depth);
        } else {
            ARIADNE_LOG(3,"Computing evolution for initial set in location "<<loc_iter->first<<" directly...\n");
            // if smaller, compute the evolution directly
            EnclosureType initial_enclosure(loc_iter->first,ContinuousEnclosureType(loc_iter->second));
            GTS cell_reach,cell_final;
            make_lpair(cell_reach,cell_final)=this->_discretiser->evolution(system,initial_enclosure,time,grid_depth,UPPER_SEMANTICS);
            reach.adjoin(cell_reach);
            evolve.adjoin(cell_final);
        }
    }
    if(!initial.empty()) {
        ARIADNE_LOG(3,"computing evolution from the grid...\n");
        ARIADNE_LOG(4,"initial_evolve.size()="<<initial.size()<<"\n");
        make_lpair(found,initial)=this->_upper_reach_evolve(system,initial,hybrid_lock_to_grid_time,grid_depth);
        evolve.adjoin(initial);
        reach.adjoin(found);
    }
    // time steps evolution loop        
    for(uint i=1; i<time_steps; ++i) {
        ARIADNE_LOG(3,"computing "<<i+1<<"-th reachability step...\n");
        make_lpair(found,evolve) = this->_upper_reach_evolve(system,evolve,hybrid_lock_to_grid_time,grid_depth);
        ARIADNE_LOG(5,"found.size()="<<found.size()<<"\n");
        ARIADNE_LOG(5,"evolve.size()="<<evolve.size()<<"\n");
        reach.adjoin(found);
        ARIADNE_LOG(3,"  found "<<found.size()<<" cells.\n");
    }
    ARIADNE_LOG(3,"remainder_time="<<remainder_time<<"\n");
    if(!evolve.empty() && remainder_time > 0) {
        ARIADNE_LOG(3,"computing evolution for remainder time...\n");
        make_lpair(found,evolve) = this->_upper_reach_evolve(system,evolve,hybrid_remainder_time,grid_depth);
        reach.adjoin(found);
    }
    reach.recombine();
    ARIADNE_LOG(4,"reach="<<reach<<"\n");
    evolve.recombine();
    ARIADNE_LOG(4,"evolve="<<evolve<<"\n");
    return std::make_pair(reach,evolve);
}






HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
chain_reach(const SystemType& system,
            const HybridImageSet& initial_set) const
{
    ARIADNE_LOG(2,"HybridReachabilityAnalyser::chain_reach(system,set)\n");
    Float transient_time = this->_parameters->transient_time;
    int transient_steps = this->_parameters->transient_steps;
    Float lock_to_grid_time=this->_parameters->lock_to_grid_time;
    int lock_to_grid_steps=this->_parameters->lock_to_grid_steps;
    HybridTime hybrid_lock_to_grid_time(lock_to_grid_time,lock_to_grid_steps);
    int maximum_grid_depth = this->_parameters->maximum_grid_depth;
    int maximum_grid_height = this->_parameters->maximum_grid_height;
    ARIADNE_LOG(3,"transient_time=("<<transient_time<<","<<transient_steps<<")\n");
    ARIADNE_LOG(3,"lock_to_grid_time=("<<lock_to_grid_time<<","<<lock_to_grid_steps<<")\n");

    ARIADNE_LOG(5,"initial_set="<<initial_set<<"\n");

    Gr grid=system.grid();
    GTS evolve(grid), reach(grid), initial(grid), found(grid);

    if(transient_time <= 0.0 || transient_steps <= 0) {
        transient_time = lock_to_grid_time;
        transient_steps = lock_to_grid_steps;
    }
    HybridTime hybrid_transient_time(transient_time, transient_steps);
    ARIADNE_LOG(3,"Computing first evolution step...\n");
    // For each location, test if the radius of the set is smaller than the grid cell
    for(HybridImageSet::locations_const_iterator loc_iter=initial_set.locations_begin();
        loc_iter!=initial_set.locations_end(); ++loc_iter) 
    {
        Vector<Float> cell = grid[loc_iter->first].lengths();
        Float cell_radius = (min(cell))/(1 << (maximum_grid_depth));
        if (radius(loc_iter->second.bounding_box()) > cell_radius) {
            // if bigger, map to the grid
            ARIADNE_LOG(3,"Adjoining initial set for location "<<loc_iter->first<<" to the grid...\n");
            initial[loc_iter->first].adjoin_outer_approximation(loc_iter->second,maximum_grid_depth);
        } else {
            ARIADNE_LOG(3,"Computing evolution for initial set in location "<<loc_iter->first<<" directly...\n");
            // if smaller, compute the evolution directly
            EnclosureType initial_enclosure(loc_iter->first,ContinuousEnclosureType(loc_iter->second));
            GTS cell_reach,cell_final;
            make_lpair(cell_reach,cell_final)=this->_discretiser->evolution(system,initial_enclosure,hybrid_transient_time,maximum_grid_depth,UPPER_SEMANTICS);
            reach.adjoin(cell_reach);
            evolve.adjoin(cell_final);
        }
    }
    if(!initial.empty()) {
        ARIADNE_LOG(3,"Computing evolution on the grid...\n")
        make_lpair(found,initial)=this->_upper_reach_evolve(system,initial,hybrid_transient_time,maximum_grid_depth);
        ARIADNE_LOG(3,"  found "<<found.size()<<" cells.\n");
        reach.adjoin(found);
        evolve.adjoin(initial);
    }
    ARIADNE_LOG(3,"Computing recurrent evolution...\n");
    while(!evolve.empty()) {
        make_lpair(found,evolve)=this->_upper_reach_evolve(system,evolve,hybrid_lock_to_grid_time,maximum_grid_depth);
        ARIADNE_LOG(5,"found.size()="<<found.size()<<"\n");
        ARIADNE_LOG(5,"evolve.size()="<<evolve.size()<<"\n");
        evolve.remove(reach);
        found.restrict_to_height(maximum_grid_height);
        reach.adjoin(found);
        ARIADNE_LOG(3,"  found "<<found.size()<<" cells, of which "<<evolve.size()<<" are new.\n");
        // ARIADNE_LOG(5,"bounded_new_size="<<found.size()<<"\n");
    }
    reach.recombine();
    return reach;
}


HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
chain_reach(const SystemType& system,
            const HybridImageSet& initial_set,
            const HybridBoxes& bounding_set) const
{
    // FIXME: Use tree sets throughout

    HybridBoxes bounding_domain=bounding_set;
    ARIADNE_LOG(2,"HybridReachabilityAnalyser::chain_reach(system,set,bounding_set)\n");
    Float transient_time = this->_parameters->transient_time;
    int transient_steps = this->_parameters->transient_steps;
    Float lock_to_grid_time=this->_parameters->lock_to_grid_time;
    int lock_to_grid_steps=this->_parameters->lock_to_grid_steps;
    int maximum_grid_depth = this->_parameters->maximum_grid_depth;
    ARIADNE_LOG(3,"transient_time=("<<transient_time<<","<<transient_steps<<")\n");
    ARIADNE_LOG(3,"lock_to_grid_time=("<<lock_to_grid_time<<","<<lock_to_grid_steps<<")\n");

    ARIADNE_LOG(5,"bounding_domain="<<bounding_domain<<"\n");
    ARIADNE_LOG(5,"initial_set="<<initial_set<<"\n");

    Gr grid=system.grid();
    //GTS bounding; bounding.adjoin_inner_approximation(bounding_domain,maximum_grid_depth);
    GTS bounding(grid);
    bounding.adjoin_outer_approximation(bounding_domain,maximum_grid_depth); bounding.recombine();
    ARIADNE_LOG(5,"bounding_size="<<bounding.size()<<"\n");

    GTS evolve(grid), reach(grid), initial(grid), found(grid);

    if(transient_time <= 0.0 || transient_steps <= 0) {
        transient_time = lock_to_grid_time;
        transient_steps = lock_to_grid_steps;
    }
    HybridTime hybrid_transient_time(transient_time, transient_steps);
    ARIADNE_LOG(3,"Computing first evolution step...\n");
    // For each location, test if the radius of the set is smaller than the grid cell
    for(HybridImageSet::locations_const_iterator loc_iter=initial_set.locations_begin();
        loc_iter!=initial_set.locations_end(); ++loc_iter) 
    {
        Vector<Float> cell = grid[loc_iter->first].lengths();
        Float cell_radius = (min(cell))/(1 << (maximum_grid_depth));
        if (radius(loc_iter->second.bounding_box()) > cell_radius) {
            // if bigger, map to the grid
            ARIADNE_LOG(3,"Adjoining initial set for location "<<loc_iter->first<<" to the grid...\n");
            initial[loc_iter->first].adjoin_outer_approximation(loc_iter->second,maximum_grid_depth);
        } else {
            ARIADNE_LOG(3,"Computing evolution for initial set in location "<<loc_iter->first<<" directly...\n");
            // if smaller, compute the evolution directly
            EnclosureType initial_enclosure(loc_iter->first,ContinuousEnclosureType(loc_iter->second));
            GTS cell_reach,cell_final;
            make_lpair(cell_reach,cell_final)=this->_discretiser->evolution(system,initial_enclosure,hybrid_transient_time,maximum_grid_depth,UPPER_SEMANTICS);
            reach.adjoin(cell_reach);
            evolve.adjoin(cell_final);
        }
    }
    if(!initial.empty()) {
        ARIADNE_LOG(3,"Computing evolution on the grid...\n")
        make_lpair(found,initial)=this->_upper_reach_evolve(system,initial,hybrid_transient_time,maximum_grid_depth);
        ARIADNE_LOG(3,"  found "<<found.size()<<" cells.\n");
        reach.adjoin(found);
        evolve.adjoin(initial);
    }
    evolve.restrict(bounding);
    ARIADNE_LOG(3,"  found "<<reach.size()<<" cells, of which "<<evolve.size()<<" are new.\n");   
    
    ARIADNE_LOG(3,"Computing recurrent evolution...\n");
    HybridTime hybrid_lock_to_grid_time(lock_to_grid_time,lock_to_grid_steps);
    while(!evolve.empty()) {
        make_lpair(found,evolve)=this->_upper_reach_evolve(system,evolve,hybrid_lock_to_grid_time,maximum_grid_depth);
        ARIADNE_LOG(5,"found.size()="<<found.size()<<"\n");
        ARIADNE_LOG(5,"evolve.size()="<<evolve.size()<<"\n");
        evolve.remove(reach);
        evolve.restrict(bounding);
        reach.adjoin(found);
        ARIADNE_LOG(3,"  found "<<found.size()<<" cells, of which "<<evolve.size()<<" are new.\n");
        // ARIADNE_LOG(5,"bounded_new_size="<<found.size()<<"\n");
    }
    reach.recombine();
    reach.restrict(bounding);
    return reach;
}





HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
viable(const SystemType& system,
       const HybridImageSet& bounding_set) const
{
    ARIADNE_NOT_IMPLEMENTED;
}



tribool
HybridReachabilityAnalyser::
verify(const SystemType& system,
       const HybridImageSet& initial_set,
       const HybridImageSet& safe_set) const
{
    ARIADNE_NOT_IMPLEMENTED;
}












} // namespace Ariadne
