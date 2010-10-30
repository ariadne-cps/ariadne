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

HybridGrid::HybridGrid() : _grids(), _system_ptr() { }

HybridGrid::HybridGrid(const HybridAutomatonInterface& ha) : _grids(),  _system_ptr(&ha) { }

HybridGrid::HybridGrid(const HybridSpace& hs, double l) {
    for(HybridSpace::const_iterator iter=hs.begin(); iter!=hs.end(); ++iter) {
        this->_grids.insert(iter->first,Grid(iter->second,l));
    }
}

void HybridGrid::insert(DiscreteLocation q, const Grid& g) {
    this->_grids.insert(q,g);
}

Grid HybridGrid::operator[](const DiscreteLocation& loc) const {
    if(_system_ptr) { return _system_ptr->grid(loc); }
    else { return this->_grids[loc]; }
}

Grid& HybridGrid::operator[](const DiscreteLocation& loc) {
    if(_system_ptr) { this->_grids[loc]=this->_system_ptr->grid(loc); return this->_grids[loc]; }
    else { return this->_grids[loc]; }
}

bool HybridGrid::has_location(const DiscreteLocation& q) const {
    return this->_grids.find(q) != this->_grids.end();
}


} // namespace Ariadne
