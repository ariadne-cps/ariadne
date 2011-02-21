/***************************************************************************
 *            hybrid_grid.h
 *
 *  Copyright 2008-11  Pieter Collins
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

/*! \file hybrid_grid.h
 *  \brief Grids in hybrid spaces.
 */

#ifndef ARIADNE_HYBRID_GRID_H
#define ARIADNE_HYBRID_GRID_H

#include <map>

#include "container.h"
#include "stlio.h"
#include "space.h"
#include "hybrid_space.h"

#include "grid.h"

namespace Ariadne {

class HybridGridTreeSet;


//! \ingroup HybridModule
//! \brief A class which defines the state space grid to use in location \a loc given the continuous state variables \a spc.
class HybridScalingInterface
{
  public:
    //!
    virtual Grid grid(const DiscreteLocation& loc, const RealSpace& spc) const = 0;
    virtual std::ostream& write(std::ostream& os) const = 0;
    friend std::ostream& operator<<(std::ostream& os, const HybridScalingInterface& hsc) { return hsc.write(os); }
};

class HybridScaling
    : public HybridScalingInterface
{
    Map<String,Float> _scalings;
  public:
    HybridScaling() : _scalings() { }
    virtual Grid grid(const DiscreteLocation& loc, const RealSpace& spc) const;
    virtual std::ostream& write(std::ostream& os) const { return os << "HybridScaling( " << this->_scalings << " )"; }
};

inline Grid
HybridScaling::grid(const DiscreteLocation& loc, const RealSpace& spc) const
{
    FloatVector lengths(spc.dimension(),1.0);
    for(uint i=0; i!=lengths.size(); ++i) {
        if(this->_scalings.has_key(spc[i].name())) {
            lengths[i] = this->_scalings[spc[i].name()];
        }
    }
    return Grid(lengths);
}




//! \ingroup HybridModule
//! \brief A grid in a hybrid space
class HybridGrid
{
    // NOTE: The use of the system is to allow the "Grid" of a compositional
    // hybrid automaton to be computed on-the-fly since it might not be
    // feasible to compute the reachable states.
    // TODO: There should be a better way of doing this, but this might
    // involve changing the HybridReachabilityAnalyser or HybridDiscretiser code.
    HybridSpace _space;
    shared_ptr<const HybridScalingInterface> _scaling_ptr;
  public:
    //! Test whether the grid has location \a q.
    bool has_location(const DiscreteLocation& q) const {
        return this->_space.has_location(q); }
    //! The grid in location \a loc.
    Grid operator[](const DiscreteLocation& loc) const {
        return this->_scaling_ptr->grid(loc,this->_space.operator[](loc)); }
    const HybridSpace& space() const { return this->_space; }
  public:
    //! The grid corresponding to unit resolution on each dimension.
    HybridGrid(const HybridSpace& hsp) : _space(hsp), _scaling_ptr(new HybridScaling()) { }
    //!
    HybridGrid(const HybridSpace& hsp, const HybridScaling& hsc)
        : _space(hsp)
        , _scaling_ptr(new HybridScaling(hsc)) { }

    //!
    friend inline std::ostream&
    operator<<(std::ostream& os, const HybridGrid& hgrid) {
        return os << "HybridGrid( space="<<hgrid._space << ", scalings=" << *hgrid._scaling_ptr << " )"; }
};




}

#endif // ARIADNE_HYBRID_GRID_H

