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

#include "function/functional.h"
#include "config.h"

#include "hybrid/hybrid_set.h"

#include "numeric/real.h"

#include "expression/expression_set.h"
#include "geometry/function_set.h"

#include "hybrid/hybrid_space.h"
#include "hybrid/hybrid_time.h"
#include "hybrid/hybrid_orbit.h"
#include "hybrid/hybrid_automaton_interface.h"
#include "output/graphics.h"
#include "hybrid/hybrid_graphics.h"
#include <boost/concept_check.hpp>
#include "numeric/rounding.h"
#include "expression/assignment.h"
#include "output/graphics_interface.h"
#include "geometry/function_set.h"

namespace Ariadne {

template<> inline Float64Value numeric_cast<Float64Value>(Real const& r) {
    return cast_exact(Float64Approximation(r));
}

Orbit<HybridPoint>::Orbit(const HybridPoint& hpt)
    : _curves_ptr(new std::vector<HybridInterpolatedCurve>(1u,HybridInterpolatedCurve(hpt.location(),hpt.space(),InterpolatedCurve(0,hpt.point()))))
{ }

Nat
Orbit<HybridPoint>::size() const
{
    return this->_curves_ptr->size();
}

const InterpolatedCurve&
Orbit<HybridPoint>::curve(Nat m) const
{
    return (*this->_curves_ptr)[m].continuous_set();
}

Void
Orbit<HybridPoint>::insert(HybridTime ht, const HybridPoint& hpt)
{
    ARIADNE_ASSERT(ht.discrete_time()<=this->size());
    // FIXME: Should allow non-exact times
    Real time=ht.continuous_time();
    Float64Value flt_time=cast_exact(time);
    ARIADNE_ASSERT(decide(Real(flt_time)==time));
    if(this->size()==ht.discrete_time()) {
        this->_curves_ptr->push_back(HybridInterpolatedCurve(hpt.location(),hpt.space(),InterpolatedCurve(flt_time,hpt.point())));
    } else {
        (*this->_curves_ptr)[ht.discrete_time().get_si()].continuous_set().insert(flt_time,hpt.point());
    }
}

Void Orbit<HybridPoint>::draw(CanvasInterface& canvas, const Set<DiscreteLocation>& locations, const Variables2d& axes) const {
    const Orbit<HybridPoint>& orbit=*this;
    for(Nat i=0; i!=orbit._curves_ptr->size(); ++i) {
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
OutputStream&
operator<<(OutputStream& os, const Orbit< HybridPoint >& orb)
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



Void Orbit<HybridEnclosure>::draw(CanvasInterface& c, const Set<DiscreteLocation>& l, const Variables2d& v) const {
    this->reach().draw(c,l,v);
}

template<>
OutputStream&
operator<<(OutputStream& os, const Orbit< HybridEnclosure >& orb)
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
    for(Nat i=0; i!=b.size(); ++i) {
        res.insert(b[i].variable(),RealInterval(b[i].lower(),b[i].upper()));
    }
    return res;
}

HybridPoint::HybridPoint(const DiscreteLocation& q, const Map<Identifier,Real>& x)
    : HybridBasicSet<ExactPoint>(q,make_list(x.keys()),ExactPoint(x.size()))
{
    Nat i=0;
    for(Map<Identifier,Real>::ConstIterator iter=x.begin(); iter!=x.end(); ++iter, ++i) {
        this->point()[i]=numeric_cast<Float64Value>(iter->second);
    }
}

HybridPoint::HybridPoint(const DiscreteLocation& q, const Map<Identifier,Float64Value>& x)
    : HybridBasicSet<ExactPoint>(q,make_list(x.keys()),ExactPoint(x.size()))
{
    Nat i=0;
    for(Map<Identifier,Float64Value>::ConstIterator iter=x.begin(); iter!=x.end(); ++iter, ++i) {
        this->point()[i]=iter->second;
    }
}

HybridPoint::HybridPoint(const DiscreteLocation& q, const List<RealConstantAssignment>& x)
    : HybridBasicSet<ExactPoint>(q,left_hand_sides(x),ExactPoint(x.size()))
{
    Nat i=0;
    for(List<RealConstantAssignment>::ConstIterator iter=x.begin(); iter!=x.end(); ++iter, ++i) {
        this->point()[i]=numeric_cast<Float64Value>(iter->right_hand_side());
    }
}

HybridPoint::HybridPoint(const DiscreteLocation& q, const InitializerList<RealConstantAssignment>& x)
    : HybridPoint(q,List<RealConstantAssignment>(x))
{
}

Map<RealVariable,Float64Value> HybridPoint::values() const {
    Map<RealVariable,Float64Value> r;
    for(Nat i=0; i!=this->space().dimension(); ++i) {
        r.insert(this->space()[i],this->point()[i]);
    }
    return r;
}




template<class BS> Void draw(CanvasInterface& canvas, const DiscreteLocation& location, const Variables2d& axes, const HybridBasicSet<BS>& set)
{
    if(set.location()==location) {
        Projection2d projection(set.continuous_set().dimension(),set.space().index(axes.x_variable()),set.space().index(axes.y_variable()));
        set.continuous_set().draw(canvas,projection);
    }
}

Void HybridBoxType::draw(CanvasInterface& c, const Set<DiscreteLocation>& q, const Variables2d& p) const {
    if(q.empty() || q.contains(this->location())) {
        this->continuous_set().draw(c,projection(this->space(),p));
    }
}

Void HybridGridTreeSet::draw(CanvasInterface& canvas, const Set<DiscreteLocation>& locations, const Variables2d& axis_variables) const {
    for(LocationsConstIterator loc_iter=this->locations_begin(); loc_iter!=this->locations_end(); ++loc_iter) {
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

ValidatedSierpinskian HybridConstraintSet::overlaps(const HybridBoxType& bx) const {
    if(this->_sets.has_key(bx.location())) {
        return this->_sets[bx.location()].euclidean_set(bx.space()).overlaps(bx.continuous_set());
    } else {
        return false;
    }
}

ValidatedSierpinskian HybridConstraintSet::covers(const HybridBoxType& bx) const {
    if(this->_sets.has_key(bx.location())) {
        return this->_sets[bx.location()].euclidean_set(bx.space()).covers(bx.continuous_set());
    } else {
        return bx.continuous_set().is_empty();
    }
}

ValidatedSierpinskian HybridConstraintSet::separated(const HybridBoxType& bx) const {
    if(this->_sets.has_key(bx.location())) {
        return this->_sets[bx.location()].euclidean_set(bx.space()).separated(bx.continuous_set());
    } else {
        return true;
    }
}

OutputStream& HybridConstraintSet::write(OutputStream& os) const {
    return os << "HybridConstraintSet( "<< this->_sets << " )";
}



HybridBoundedConstraintSet::HybridBoundedConstraintSet()
    : _sets()
{
}

HybridBoundedConstraintSet::HybridBoundedConstraintSet(const DiscreteLocation& loc,
                                                               const InitializerList<RealVariableInterval>& bnd)
    : _sets()
{
    _sets.insert(loc,RealExpressionBoundedConstraintSet(bnd));
}

HybridBoundedConstraintSet::HybridBoundedConstraintSet(const DiscreteLocation& loc,
                                                               const InitializerList<RealVariableInterval>& bnd,
                                                               const InitializerList<ContinuousPredicate>& cnstr)
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

ValidatedSierpinskian HybridBoundedConstraintSet::overlaps(const HybridBoxType& bx) const {
    if(this->_sets.has_key(bx.location())) {
        return this->_sets[bx.location()].euclidean_set(bx.space()).overlaps(bx.continuous_set());
    } else {
        return false;
    }
}

ValidatedSierpinskian HybridBoundedConstraintSet::covers(const HybridBoxType& bx) const {
    if(this->_sets.has_key(bx.location())) {
        return this->_sets[bx.location()].euclidean_set(bx.space()).covers(bx.continuous_set());
    } else {
        return bx.continuous_set().is_empty();
    }
}

ValidatedSierpinskian HybridBoundedConstraintSet::separated(const HybridBoxType& bx) const {
    if(this->_sets.has_key(bx.location())) {
        return this->_sets[bx.location()].euclidean_set(bx.space()).separated(bx.continuous_set());
    } else {
        return true;
    }
}

ValidatedSierpinskian HybridBoundedConstraintSet::inside(const HybridBoxes& bxs) const {
    ValidatedSierpinskian result=true;
    for(Map<DiscreteLocation,RealExpressionBoundedConstraintSet>::ConstIterator iter=this->_sets.begin(); iter!=this->_sets.end(); ++iter) {
        DiscreteLocation const& loc=iter->first;
        RealExpressionBoundedConstraintSet const& set = iter->second;
        Set<RealVariable> vars=set.variables();
        RealSpace const& spc = bxs[loc].space();
        ExactBoxType const& bx = bxs[loc].continuous_set();
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
    for(Map<DiscreteLocation,RealExpressionBoundedConstraintSet>::ConstIterator iter=this->_sets.begin(); iter!=this->_sets.end(); ++iter) {
        result.insert(iter->first,over_approximation(RealVariablesBox(iter->second.bounds())));
    }
    return result;
}

OutputStream& HybridBoundedConstraintSet::write(OutputStream& os) const {
    return os << "HybridBoundedConstraintSet( "<< this->_sets << " )";
}

Void HybridBoundedConstraintSet::draw(CanvasInterface& c, const Set<DiscreteLocation>& q, const Variables2d& p) const {
    if(q.empty() || q.contains(this->location())) {
        Set<RealVariable> variables=this->variables(this->location());
        RealSpace space(List<RealVariable>(variables.begin(),variables.end()));
        this->euclidean_set(this->location(),space).draw(c,projection(space,p));
    }
}

} // namespace Ariadne
