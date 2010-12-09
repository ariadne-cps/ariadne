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
#include "hybrid_automaton_interface.h"

namespace Ariadne {

HybridGrid::HybridGrid() : _grids(), _system_ptr(0) { }

HybridGrid::HybridGrid(const HybridAutomatonInterface& ha) : _grids(),  _system_ptr(&ha) { }

HybridGrid::HybridGrid(const HybridSpace& hs) : _grids(), _system_ptr(0)  {
    for(HybridSpace::const_iterator iter=hs.begin(); iter!=hs.end(); ++iter) {
        this->_grids.insert(iter->first,Grid(iter->second));
    }
}

void HybridGrid::insert(DiscreteLocation q, const Grid& g) {
    this->_grids.insert(q,g);
}

void HybridGrid::insert(const std::pair<DiscreteLocation,Grid>& qg) {
    this->_grids.insert(qg);
}

Grid HybridGrid::operator[](const DiscreteLocation& q) const {
    if(this->_grids.has_key(q)) { return this->_grids[q]; }
    if(_system_ptr!=0) { this->_grids.insert(q,_system_ptr->grid(q)); return this->_grids[q]; }

    ARIADNE_THROW(std::runtime_error,"HybridGrid::operator[](DiscreteLocation)",
                  "No location "<<q<<" in grid.");
}

bool HybridGrid::has_location(const DiscreteLocation& q) const {
    return this->_grids.has_key(q);
}


} // namespace Ariadne
