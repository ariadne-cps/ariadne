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

#include "functional.h"
#include "config.h"

#include "hybrid_set.h"

#include "real.h"

#include "expression_set.h"
#include "function_set.h"

#include "hybrid_space.h"
#include "hybrid_time.h"
#include "hybrid_orbit.h"
#include "hybrid_automaton_interface.h"
#include "graphics.h"
#include "hybrid_graphics.h"
#include <boost/concept_check.hpp>
#include "rounding.h"
#include "assignment.h"
#include "graphics_interface.h"
#include "function_set.h"

namespace Ariadne {



Orbit<HybridPoint>::Orbit(const HybridPoint& hpt)
    : _curves_ptr(new std::vector<HybridInterpolatedCurve>(1u,HybridInterpolatedCurve(hpt.location(),hpt.space(),InterpolatedCurve(hpt.point()))))
{ }

uint
Orbit<HybridPoint>::size() const
{
    return this->_curves_ptr->size();
}

const InterpolatedCurve&
Orbit<HybridPoint>::curve(uint m) const
{
    return (*this->_curves_ptr)[m].continuous_set();
}

void
Orbit<HybridPoint>::insert(HybridTime ht, const HybridPoint& hpt)
{
    ARIADNE_ASSERT((uint)ht.discrete_time()<=this->size());
    if(this->size()==(uint)ht.discrete_time()) {
        this->_curves_ptr->push_back(HybridInterpolatedCurve(hpt.location(),hpt.space(),InterpolatedCurve(Float(ht.continuous_time()),hpt.point())));
    } else {
        (*this->_curves_ptr)[ht.discrete_time()].continuous_set().insert(static_cast<Float>(ht.continuous_time()),hpt.point());
    }
}

void Orbit<HybridPoint>::draw(CanvasInterface& canvas, const Set<DiscreteLocation>& locations, const Variables2d& axes) const {
    const Orbit<HybridPoint>& orbit=*this;
    for(uint i=0; i!=orbit._curves_ptr->size(); ++i) {
        HybridInterpolatedCurve const& hcurve=this->_curves_ptr->at(i);
        if(locations.empty() || locations.contains(hcurve.location())) {
            RealSpace const& space=hcurve.space();
            if(valid_axis_variables(space,axes)) {
                hcurve.continuous_set().draw(canvas,projection(space,axes));
            }
        }
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



void Orbit<HybridEnclosure>::draw(CanvasInterface& c, const Set<DiscreteLocation>& l, const Variables2d& v) const {
    this->reach().draw(c,l,v);
}

template<>
std::ostream&
operator<<(std::ostream& os, const Orbit< HybridEnclosure >& orb)
{
    os << "Orbit(\n  initial=" << orb.initial()
       << "\n  intermediate=" << orb.intermediate()
       << "\n  reach=" << orb.reach()
       << "\n  final=" << orb.final()
       << ")\n";
    return os;
}







Map<RealVariable,IntervalSet> make_map(const List<RealVariableInterval>& b) {
    Map<RealVariable,IntervalSet> res;
    for(uint i=0; i!=b.size(); ++i) {
        res.insert(b[i].variable(),IntervalSet(b[i].lower(),b[i].upper()));
    }
    return res;
}



HybridPoint::HybridPoint(const DiscreteLocation& q, const Map<Identifier,Real>& x)
    : HybridBasicSet<Point>(q,make_list(x.keys()),Point(x.size()))
{
    uint i=0;
    for(Map<Identifier,Real>::const_iterator iter=x.begin(); iter!=x.end(); ++iter, ++i) {
        this->point()[i]=numeric_cast<Float>(iter->second);
    }
}

HybridPoint::HybridPoint(const DiscreteLocation& q, const Map<Identifier,Float>& x)
    : HybridBasicSet<Point>(q,make_list(x.keys()),Point(x.size()))
{
    uint i=0;
    for(Map<Identifier,Float>::const_iterator iter=x.begin(); iter!=x.end(); ++iter, ++i) {
        this->point()[i]=iter->second;
    }
}

HybridPoint::HybridPoint(const DiscreteLocation& q, const List<RealConstantAssignment>& x)
    : HybridBasicSet<Point>(q,left_hand_sides(x),Point(x.size()))
{
    uint i=0;
    for(List<RealConstantAssignment>::const_iterator iter=x.begin(); iter!=x.end(); ++iter, ++i) {
        this->point()[i]=static_cast<Float>(iter->right_hand_side());
    }
}

Map<RealVariable,Float> HybridPoint::values() const {
    Map<RealVariable,Float> r;
    for(uint i=0; i!=this->space().dimension(); ++i) {
        r.insert(this->space()[i],this->point()[i]);
    }
    return r;
}




template<class BS> void draw(CanvasInterface& canvas, const DiscreteLocation& location, const Variables2d& axes, const HybridBasicSet<BS>& set)
{
    if(set.location()==location) {
        Projection2d projection(set.continuous_set().dimension(),set.space().index(axes.x_variable()),set.space().index(axes.y_variable()));
        set.continuous_set().draw(canvas,projection);
    }
}

void HybridBox::draw(CanvasInterface& c, const Set<DiscreteLocation>& q, const Variables2d& p) const {
    if(q.empty() || q.contains(this->location())) {
        this->continuous_set().draw(c,projection(this->space(),p));
    }
}

void HybridGridTreeSet::draw(CanvasInterface& canvas, const Set<DiscreteLocation>& locations, const Variables2d& axis_variables) const {
    for(locations_const_iterator loc_iter=this->locations_begin(); loc_iter!=this->locations_end(); ++loc_iter) {
        if(locations.empty() || locations.contains(loc_iter->first)) {
            RealSpace const& space=this->space(loc_iter->first);
            Projection2d projection(space.dimension(),space.index(axis_variables.x_variable()),space.index(axis_variables.y_variable()));
            loc_iter->second.draw(canvas,projection);
        }
    }
}



HybridConstraintSet::HybridConstraintSet()
    : _sets()
{
}

HybridConstraintSet::HybridConstraintSet(const DiscreteLocation& loc,
                                                 const List<ContinuousPredicate>& cnstr)
    : _sets()
{
    _sets.insert(loc,RealExpressionConstraintSet(cnstr));
}

HybridConstraintSet* HybridConstraintSet::clone() const {
    return new HybridConstraintSet(*this);
}

Set<RealVariable> HybridConstraintSet::variables(DiscreteLocation loc) const {
    ARIADNE_ASSERT(this->_sets.has_key(loc));
    return _sets[loc].variables();
}

ConstraintSet const HybridConstraintSet::euclidean_set(DiscreteLocation loc, RealSpace spc) const {
    ARIADNE_ASSERT(this->_sets.has_key(loc));
    return ConstraintSet(this->_sets[loc].euclidean_set(spc));
}

RegularSetInterface* HybridConstraintSet::_euclidean_set(DiscreteLocation loc, RealSpace spc) const {
    return new ConstraintSet(this->euclidean_set(loc,spc));
}

tribool HybridConstraintSet::overlaps(const HybridBox& bx) const {
    if(this->_sets.has_key(bx.location())) {
        return this->_sets[bx.location()].euclidean_set(bx.space()).overlaps(bx.continuous_set());
    } else {
        return false;
    }
}

tribool HybridConstraintSet::covers(const HybridBox& bx) const {
    if(this->_sets.has_key(bx.location())) {
        return this->_sets[bx.location()].euclidean_set(bx.space()).covers(bx.continuous_set());
    } else {
        return bx.continuous_set().empty();
    }
}

tribool HybridConstraintSet::separated(const HybridBox& bx) const {
    if(this->_sets.has_key(bx.location())) {
        return this->_sets[bx.location()].euclidean_set(bx.space()).separated(bx.continuous_set());
    } else {
        return true;
    }
}

std::ostream& HybridConstraintSet::write(std::ostream& os) const {
    return os << "HybridConstraintSet( "<< this->_sets << " )";
}



HybridBoundedConstraintSet::HybridBoundedConstraintSet()
    : _sets()
{
}

HybridBoundedConstraintSet::HybridBoundedConstraintSet(const DiscreteLocation& loc,
                                                               const List<RealVariableInterval>& bnd)
    : _sets()
{
    _sets.insert(loc,RealExpressionBoundedConstraintSet(bnd));
}

HybridBoundedConstraintSet::HybridBoundedConstraintSet(const DiscreteLocation& loc,
                                                               const List<RealVariableInterval>& bnd,
                                                               const List<ContinuousPredicate>& cnstr)
    : _sets()
{
    _sets.insert(loc,RealExpressionBoundedConstraintSet(bnd,cnstr));
}

HybridBoundedConstraintSet::HybridBoundedConstraintSet(const DiscreteLocation& loc,
                                                               const RealVariablesBox& bx)
    : _sets()
{
    _sets.insert(loc,RealExpressionBoundedConstraintSet(bx));
}

HybridBoundedConstraintSet* HybridBoundedConstraintSet::clone() const {
    return new HybridBoundedConstraintSet(*this);
}

Set<RealVariable> HybridBoundedConstraintSet::variables(DiscreteLocation loc) const {
    ARIADNE_ASSERT(this->_sets.has_key(loc));
    return _sets[loc].variables();
}

BoundedConstraintSet const HybridBoundedConstraintSet::euclidean_set(DiscreteLocation loc, RealSpace spc) const {
    // FIXME: Should be no need to cache Euclidean sets.
    ARIADNE_ASSERT(this->_sets.has_key(loc));
    return BoundedConstraintSet(this->_sets[loc].euclidean_set(spc));
}

BoundedConstraintSet* HybridBoundedConstraintSet::_euclidean_set(DiscreteLocation loc, RealSpace spc) const {
    return new BoundedConstraintSet(this->euclidean_set(loc,spc));
}

tribool HybridBoundedConstraintSet::overlaps(const HybridBox& bx) const {
    if(this->_sets.has_key(bx.location())) {
        return this->_sets[bx.location()].euclidean_set(bx.space()).overlaps(bx.continuous_set());
    } else {
        return false;
    }
}

tribool HybridBoundedConstraintSet::covers(const HybridBox& bx) const {
    if(this->_sets.has_key(bx.location())) {
        return this->_sets[bx.location()].euclidean_set(bx.space()).covers(bx.continuous_set());
    } else {
        return bx.continuous_set().empty();
    }
}

tribool HybridBoundedConstraintSet::separated(const HybridBox& bx) const {
    if(this->_sets.has_key(bx.location())) {
        return this->_sets[bx.location()].euclidean_set(bx.space()).separated(bx.continuous_set());
    } else {
        return true;
    }
}

tribool HybridBoundedConstraintSet::inside(const HybridBoxes& bxs) const {
    tribool result=true;
    for(Map<DiscreteLocation,RealExpressionBoundedConstraintSet>::const_iterator iter=this->_sets.begin(); iter!=this->_sets.end(); ++iter) {
        DiscreteLocation const& loc=iter->first;
        RealExpressionBoundedConstraintSet const& set = iter->second;
        Set<RealVariable> vars=set.variables();
        RealSpace const& spc = bxs[loc].space();
        Box const& bx = bxs[loc].continuous_set();
        result = result && set.euclidean_set(spc).inside(bx);
    }
    return result;
}

DiscreteLocation HybridBoundedConstraintSet::location() const {
    ARIADNE_ASSERT(this->_sets.size()==1);
    return this->_sets.begin()->first;
}

Set<DiscreteLocation> HybridBoundedConstraintSet::locations() const {
    return this->_sets.keys();
}

HybridBoxes HybridBoundedConstraintSet::bounding_box() const {
    HybridBoxes result;
    for(Map<DiscreteLocation,RealExpressionBoundedConstraintSet>::const_iterator iter=this->_sets.begin(); iter!=this->_sets.end(); ++iter) {
        result.insert(iter->first,over_approximation(RealVariablesBox(iter->second.bounds())));
    }
    return result;
}

std::ostream& HybridBoundedConstraintSet::write(std::ostream& os) const {
    return os << "HybridBoundedConstraintSet( "<< this->_sets << " )";
}

void HybridBoundedConstraintSet::draw(CanvasInterface& c, const Set<DiscreteLocation>& q, const Variables2d& p) const {
    if(q.empty() || q.contains(this->location())) {
        Set<RealVariable> variables=this->variables(this->location());
        RealSpace space(List<RealVariable>(variables.begin(),variables.end()));
        this->euclidean_set(this->location(),space).draw(c,projection(space,p));
    }
}

} // namespace Ariadne
