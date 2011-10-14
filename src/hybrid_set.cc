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
#include "graphics.h"
#include "hybrid_graphics.h"
#include <boost/concept_check.hpp>
#include <include/rounding.h>
#include <include/assignment.h>
#include <include/graphics_interface.h>
#include <include/function_set.h>

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
    return (*this->_curves_ptr)[m].third;
}

void
Orbit<HybridPoint>::insert(HybridTime ht, const HybridPoint& hpt)
{
    ARIADNE_ASSERT((uint)ht.discrete_time()<=this->size());
    if(this->size()==(uint)ht.discrete_time()) {
        this->_curves_ptr->push_back(HybridInterpolatedCurve(hpt.location(),hpt.space(),InterpolatedCurve(Float(ht.continuous_time()),hpt.point())));
    } else {
        (*this->_curves_ptr)[ht.discrete_time()].third.insert(ht.continuous_time(),hpt.point());
    }
}

void Orbit<HybridPoint>::draw(CanvasInterface& canvas, const Set<DiscreteLocation>& locations, const Variables2d& axes) const {
    const Orbit<HybridPoint>& orbit=*this;
    for(uint i=0; i!=orbit._curves_ptr->size(); ++i) {
        HybridInterpolatedCurve const& hcurve=this->_curves_ptr->at(i);
        if(locations.empty() || locations.contains(hcurve.location())) {
            RealSpace const& space=hcurve.space();
            if(valid_axis_variables(space,axes)) {
                hcurve.continuous_state_set().draw(canvas,projection(space,axes));
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

RealBoundedConstraintSet HybridSet::continuous_state_set(const RealSpace& space) const {
    ARIADNE_ASSERT_MSG(this->_bounds.size()==space.dimension()," set="<<*this<<", space="<<space<<"\n");
    RealBox domain(this->_bounds.size());
    for(uint i=0; i!=domain.size(); ++i) {
        domain[i]=_bounds[space[i]];
    }
    List< RealNonlinearConstraint> constraints;
    for(uint i=0; i!=this->_constraints.size(); ++i) {
        constraints.append( make_function(indicator(this->_constraints[i],POSITIVE),space) <= 0 );
    }
    return RealBoundedConstraintSet(domain,constraints);
}

void HybridSet::draw(CanvasInterface& c, const Set<DiscreteLocation>& q, const Variables2d& p) const {
    if(q.empty() || q.contains(this->location())) {
        Set<RealVariable> variables=this->variables();
        RealSpace space(List<RealVariable>(variables.begin(),variables.end()));
        this->continuous_state_set(space).draw(c,projection(space,p));
    }
}

OutputStream& operator<<(OutputStream& os, const HybridSet& hs) {
    return os << "HybridSet( " << hs.location() << ", " << hs.bounds() << ", " << hs.constraints() << ")";
}

// Map<DiscreteLocation,ConstrainedImageSet> HybridBoundedConstraintSet::_sets;
// HybridSpace HybridBoundedConstraintSet::_space;

HybridPoint::HybridPoint(const DiscreteLocation& q, const Map<Identifier,Real>& x)
    : Tuple<DiscreteLocation,RealSpace,Point>(q,RealSpace(),Point(x.size()))
{
    uint i=0;
    for(Map<Identifier,Real>::const_iterator iter=x.begin(); iter!=x.end(); ++iter, ++i) {
        this->second.append(RealVariable(iter->first));
        this->third[i]=numeric_cast<Float>(iter->second);
    }
}

HybridPoint::HybridPoint(const DiscreteLocation& q, const Map<Identifier,Float>& x)
    : Tuple<DiscreteLocation,RealSpace,Point>(q,RealSpace(),Point(x.size()))
{
    uint i=0;
    for(Map<Identifier,Float>::const_iterator iter=x.begin(); iter!=x.end(); ++iter, ++i) {
        this->second.append(RealVariable(iter->first));
        this->third[i]=iter->second;
    }
}

Map<RealVariable,Float> HybridPoint::values() const {
    Map<RealVariable,Float> r;
    for(uint i=0; i!=this->second.dimension(); ++i) {
        r.insert(this->second[i],this->third[i]);
    }
    return r;
}



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

template<class BS> void draw(CanvasInterface& canvas, const DiscreteLocation& location, const Variables2d& axes, const HybridBasicSet<BS>& set)
{
    if(set.location()==location) {
        Projection2d projection(set.continuous_state_set().dimension(),set.space().index(axes.x_variable()),set.space().index(axes.y_variable()));
        set.continuous_state_set().draw(canvas,projection);
    }
}

void HybridBasicSet<Box>::draw(CanvasInterface& c, const Set<DiscreteLocation>& q, const Variables2d& p) const {
    if(q.empty() || q.contains(this->location())) {
        this->continuous_state_set().draw(c,projection(this->space(),p));
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

tribool HybridBoundedConstraintSet::separated(const HybridBox& bx) const {
    if(this->_sets.has_key(bx.location())) {
        return this->_sets[bx.location()].separated(bx.continuous_state_set());
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
