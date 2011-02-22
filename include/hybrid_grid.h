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
    virtual HybridScalingInterface* clone() const = 0;
    virtual Float scaling(const DiscreteLocation& loc, const RealVariable& var) const = 0;
    virtual std::ostream& write(std::ostream& os) const = 0;
    friend std::ostream& operator<<(std::ostream& os, const HybridScalingInterface& hsc) { return hsc.write(os); }
};

class HybridScaling
    : public HybridScalingInterface
{
    Map<String,Float> _scalings;
  public:
    HybridScaling() : _scalings() { }
    HybridScaling(const HybridScalingInterface& hsc) : _scalings(dynamic_cast<const HybridScaling&>(hsc)._scalings) { }
    void set_scaling(const RealVariable& var, Float res) { ARIADNE_ASSERT(res>0.0); _scalings[var.name()]=res; }
    virtual HybridScaling* clone() const { return new HybridScaling(*this); }
    virtual Float scaling(const DiscreteLocation& loc, const RealVariable& var) const {
        return (this->_scalings.has_key(var.name())) ? this->_scalings[var.name()] : Float(1.0); }
    virtual std::ostream& write(std::ostream& os) const { return os << "HybridScaling( " << this->_scalings << " )"; }
};



//! \ingroup HybridModule
//! \brief A grid in a hybrid space
class HybridGrid
{
    // NOTE: The use of the system is to allow the "Grid" of a compositional
    // hybrid automaton to be computed on-the-fly since it might not be
    // feasible to compute the reachable states.
    HybridSpace _space;
    HybridScaling _scaling;
  public:
    //! Test whether the grid has location \a q.
    bool has_location(const DiscreteLocation& q) const {
        return this->_space.has_location(q); }
    //! The grid in location \a loc.
    Grid operator[](const DiscreteLocation& loc) const;
    //! The underlying hybrid space.
    const HybridSpace& space() const { return this->_space; }
  public:
    //! The grid corresponding to unit resolution on each dimension.
    HybridGrid(const HybridSpace& hsp) : _space(hsp), _scaling(HybridScaling()) { }
    //! The grid corresponding to scaling each variable in each location according to the policy \a hsc.
    HybridGrid(const HybridSpace& hsp, const HybridScaling& hsc)
        : _space(hsp), _scaling(hsc) { }

    //!
    friend inline std::ostream&
    operator<<(std::ostream& os, const HybridGrid& hgrid) {
        return os << "HybridGrid( space="<<hgrid._space << ", scalings=" << hgrid._scaling << " )"; }
};

inline
Grid HybridGrid::operator[](const DiscreteLocation& loc) const
{
    RealSpace continuous_space = this->_space[loc];
    FloatVector lengths(continuous_space.size());
    for(uint i=0; i!=continuous_space.size(); ++i) {
        lengths[i] = this->_scaling.scaling(loc,continuous_space.variable(i));
    }
    return Grid(lengths);
}


}

#endif // ARIADNE_HYBRID_GRID_H

