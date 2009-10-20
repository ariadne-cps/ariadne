/***************************************************************************
 *            hybrid_set.cc
 *
 *  Copyright 2008  Pieter Collins
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
 
#include "hybrid_set.h"

namespace Ariadne {

HybridGrid::HybridGrid(uint n) : _component_grids(n) { }

HybridGrid::HybridGrid(const HybridSpace& hs, double l) {
    ARIADNE_DEPRECATED("HybridGrid(HybridSpace,double)","No current replacement");
    ARIADNE_NOT_IMPLEMENTED;
}

void HybridGrid::insert(uint i, AtomicDiscreteLocation q, const Grid& g) {
    this->_component_grids[i].insert(q,g); }

Grid HybridGrid::operator[](const DiscreteLocation& loc) const {
    assert(loc.size()==this->_component_grids.size());
    Vector<Float> lengths;
    for(uint i=0; i!=this->_component_grids.size(); ++i) {
        lengths=join(lengths,this->_component_grids[i][loc[i]].lengths());
    }
    return Grid(lengths);
}

Grid& HybridGrid::operator[](const DiscreteLocation& loc) {
    ARIADNE_DEPRECATED("Grid& HybridGrid::operator[](DiscreteLocation)",
                       "Cannot directly set grid in a location.");
    ARIADNE_NOT_IMPLEMENTED;
}

bool HybridGrid::has_location(const AtomicDiscreteLocation& q) const {
    return this->_component_grids[0].find(q) != this->_component_grids[0].end();
}

bool HybridGrid::has_location(const DiscreteLocation& loc) const {
    ARIADNE_NOT_IMPLEMENTED;
}

} // namespace Ariadne
