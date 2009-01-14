/***************************************************************************
 *            successor.cc
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
 
#include "successor.h"

#include "discretiser_interface.h"
#include "orbit.h"
#include "vector_field.h"
#include "grid_set.h"

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
    for(uchar i=0; i!=c.height+c.depth; ++i) {
        wd.push_back(c.word[i]);
    }
    return GridCell(g,c.height,wd);
}

std::vector<Cell> 
successor(const DiscretiserInterface<VectorField,GridCell>& discretiser, 
          const Grid& grid, const VectorField& system, const Cell& cell, Float time, EvolutionKind reachevolve)
{
    GridCell grid_cell=make_grid_cell(cell,grid);
    GridTreeSet grid_tree_set;

    Orbit<GridCell> orbit=discretiser.upper_evolution(system,grid_cell,time,grid_cell.depth());
    if(reachevolve==REACH) {
        grid_tree_set=orbit.reach();
    } else if(reachevolve==EVOLVE) {
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


