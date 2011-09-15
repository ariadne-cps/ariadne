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

#include "real.h"

#include "expression_set.h"
#include "function_set.h"

#include "hybrid_space.h"
#include "hybrid_time.h"
#include "hybrid_orbit.h"
#include "hybrid_automaton_interface.h"
#include <boost/concept_check.hpp>
#include <include/rounding.h>

namespace Ariadne {



Orbit<HybridPoint>::Orbit(const HybridPoint& hpt)
    : _curves(new std::vector<HybridInterpolatedCurve>(1u,HybridInterpolatedCurve(hpt.location(),hpt.space(),InterpolatedCurve(hpt.point()))))
{ }

uint
Orbit<HybridPoint>::size() const
{
    return this->_curves->size();
}

const InterpolatedCurve&
Orbit<HybridPoint>::curve(uint m) const
{
    return (*this->_curves)[m].third;
}

void
Orbit<HybridPoint>::insert(HybridTime ht, const HybridPoint& hpt)
{
    ARIADNE_ASSERT((uint)ht.discrete_time()<=this->size());
    if(this->size()==(uint)ht.discrete_time()) {
        this->_curves->push_back(HybridInterpolatedCurve(hpt.location(),hpt.space(),InterpolatedCurve(hpt.point())));
    } else {
        (*this->_curves)[ht.discrete_time()].third.insert(ht.continuous_time(),hpt.point());
    }
}


template<>
std::ostream&
operator<<(std::ostream& os, const Orbit< HybridPoint >& orb)
{
    return os << orb.curves();
}


struct Orbit<HybridGridCell>::Data {
    Data(const HybridGrid& grid)
        : initial(grid), reach(grid), intermediate(grid), final(grid) { }
    HybridGridTreeSet initial;
    HybridGridTreeSet reach;
    HybridGridTreeSet intermediate;
    HybridGridTreeSet final;
};

Orbit<HybridGridCell>::
Orbit(const HybridGridTreeSet& initial_set)
    : _data(new Data(initial_set.grid()))
{
    this->_data->initial=initial_set;
}

Orbit<HybridGridCell>::
Orbit(const HybridGridTreeSet& initial_set,
      const HybridGridTreeSet& reach_set,
      const HybridGridTreeSet& intermediate_set,
      const HybridGridTreeSet& final_set)
    : _data(new Data(initial_set.grid()))
{
    this->_data->initial=initial_set;
    this->_data->reach=reach_set;
    this->_data->intermediate=intermediate_set;
    this->_data->final=final_set;
}

HybridGridTreeSet const&
Orbit<HybridGridCell>::
initial() const
{
    return this->_data->initial;
}

HybridGridTreeSet const&
Orbit<HybridGridCell>::
reach() const
{
    return this->_data->reach;
}

HybridGridTreeSet const&
Orbit<HybridGridCell>::
intermediate() const
{
    return this->_data->intermediate;
}

HybridGridTreeSet const&
Orbit<HybridGridCell>::
final() const
{
    return this->_data->final;
}

void draw(CanvasInterface& graphic, const Orbit<HybridPoint>& orbit)
{
    for(uint i=0; i<=orbit.size(); ++i) {
        orbit.curve(i).draw(graphic);
    }
}


Map<RealVariable,RealInterval> make_map(const List<RealVariableInterval>& b) {
    Map<RealVariable,RealInterval> res;
    for(uint i=0; i!=b.size(); ++i) {
        res.insert(b[i].variable(),RealInterval(b[i].lower(),b[i].upper()));
    }
    return res;
}

HybridSet::HybridSet(const DiscreteLocation& q, const List<RealVariableInterval>& b, const List<ContinuousPredicate>& c)
    : _location(q), _bounds(make_map(b)), _constraints(c)
{
}

HybridSet::HybridSet(const DiscreteLocation& q, const RealExpressionSet& s)
    : _location(q), _bounds(s.bounds()), _constraints(s.constraints())
{
}

OutputStream& operator<<(OutputStream& os, const HybridSet& hs) {
    return os << "HybridSet( " << hs.location() << ", " << hs.bounds() << ", " << hs.constraints() << ")";
}

// Map<DiscreteLocation,ConstrainedImageSet> HybridBoundedConstraintSet::_sets;
// HybridSpace HybridBoundedConstraintSet::_space;

HybridBasicSet<Box>::HybridBasicSet(const DiscreteLocation& q, const List<RealVariableInterval>& b)
    : _location(q), _space(), _set(b.size())
{
    List<Identifier> variables;
    for(uint i=0; i!=b.size(); ++i) {
        Interval bl=b[i].lower();
        Interval bu=b[i].upper();
        ARIADNE_ASSERT_MSG(bl.lower()==bl.upper(),"Cannot convert "<<Real(b[i].lower())<<" exactly to a Float.");
        ARIADNE_ASSERT_MSG(bu.lower()==bu.upper(),"Cannot convert "<<b[i].upper()<<" exactly to a Float.");
        _set[i]=Interval(numeric_cast<Float>(b[i].lower()),numeric_cast<Float>(b[i].upper()));
        variables.append(b[i].variable().name());
    }
    _space=RealSpace(variables);
}

HybridBoundedConstraintSet::HybridBoundedConstraintSet()
    : _sets(), _spaces()
{
}

HybridBoundedConstraintSet::HybridBoundedConstraintSet(const HybridBox& hbx)
    : _sets(), _spaces()
{
    _sets.insert(hbx.location(),BoundedConstraintSet(hbx.continuous_state_set()));
    _spaces.insert(hbx.location(),hbx.space());
}

HybridBoundedConstraintSet* HybridBoundedConstraintSet::clone() const {
    return new HybridBoundedConstraintSet(*this);
}

HybridSpace HybridBoundedConstraintSet::space() const {
    return MonolithicHybridSpace(this->_spaces);
}

tribool HybridBoundedConstraintSet::overlaps(const HybridBox& bx) const {
    if(this->_sets.has_key(bx.location())) {
        return this->_sets[bx.location()].overlaps(bx.continuous_state_set());
    } else {
        return false;
    }
}

tribool HybridBoundedConstraintSet::covers(const HybridBox& bx) const {
    if(this->_sets.has_key(bx.location())) {
        return this->_sets[bx.location()].covers(bx.continuous_state_set());
    } else {
        return bx.continuous_state_set().empty();
    }
}

tribool HybridBoundedConstraintSet::disjoint(const HybridBox& bx) const {
    if(this->_sets.has_key(bx.location())) {
        return this->_sets[bx.location()].disjoint(bx.continuous_state_set());
    } else {
        return true;
    }
}


tribool HybridBoundedConstraintSet::inside(const HybridBoxes& bxs) const {
    tribool result=true;
    for(Map<DiscreteLocation,BoundedConstraintSet>::const_iterator iter=this->_sets.begin(); iter!=this->_sets.end(); ++iter) {
        result = result && iter->second.inside(bxs[iter->first]);
    }
    return result;
}

BoundedConstraintSet const& HybridBoundedConstraintSet::operator[](DiscreteLocation loc) const {
    return this->_sets[loc];
}

Set<DiscreteLocation> HybridBoundedConstraintSet::locations() const {
    return this->_sets.keys();
}

HybridBoxes HybridBoundedConstraintSet::bounding_box() const {
    HybridBoxes result;
    for(Map<DiscreteLocation,BoundedConstraintSet>::const_iterator iter=this->_sets.begin(); iter!=this->_sets.end(); ++iter) {
        result.insert(iter->first,iter->second.bounding_box());
    }
    return result;
}

std::ostream& HybridBoundedConstraintSet::write(std::ostream& os) const {
    return os << "HybridBoundedConstraintSet( "<< this->_spaces << ", " << this->_sets << " )";
}

} // namespace Ariadne
