/***************************************************************************
 *            hybrid_set.cpp
 *
 *  Copyright 2008--17  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "../function/functional.hpp"
#include "../config.hpp"

#include "../hybrid/hybrid_expression_set.hpp"
#include "../hybrid/hybrid_set.hpp"
#include "../hybrid/hybrid_paving.hpp"

#include "../numeric/real.hpp"

#include "../symbolic/expression_set.hpp"
#include "../geometry/function_set.hpp"

#include "../hybrid/hybrid_space.hpp"
#include "../hybrid/hybrid_time.hpp"
#include "../hybrid/hybrid_orbit.hpp"
#include "../hybrid/hybrid_automaton_interface.hpp"
#include "../output/graphics.hpp"
#include "../hybrid/hybrid_graphics.hpp"
#include "../numeric/rounding.hpp"
#include "../symbolic/assignment.hpp"
#include "../output/graphics_interface.hpp"
#include "../geometry/function_set.hpp"

namespace Ariadne {

ExactBoxType over_approximation(RealBox const&);
ExactBoxType under_approximation(RealBox const&);
ExactBoxType approximation(RealBox const&);
HybridExactBox under_approximation(HybridRealBox const& hbx);
HybridExactBox over_approximation(HybridRealBox const& hbx);

HybridExactBox under_approximation(HybridRealBox const& hbx) {
    return HybridExactBox(hbx.location(),hbx.space(),under_approximation(hbx.euclidean_set()));
}

HybridExactBox over_approximation(HybridRealBox const& hbx) {
    return HybridExactBox(hbx.location(),hbx.space(),over_approximation(hbx.euclidean_set()));
}


Orbit<HybridApproximatePoint>::Orbit(const HybridApproximatePoint& hpt)
    : _curves_ptr(new std::vector<HybridInterpolatedCurve>(1u,HybridInterpolatedCurve(hpt.location(),hpt.space(),InterpolatedCurve(0,hpt.point()))))
{ }

Nat
Orbit<HybridApproximatePoint>::size() const
{
    return this->_curves_ptr->size();
}

const InterpolatedCurve&
Orbit<HybridApproximatePoint>::curve(Nat m) const
{
    return (*this->_curves_ptr)[m].euclidean_set();
}

Void
Orbit<HybridApproximatePoint>::insert(HybridTime ht, const HybridApproximatePoint& hpt)
{
    ARIADNE_ASSERT(ht.discrete_time()<=this->size());
    // FIXME: Should allow non-exact times
    Real time=ht.continuous_time();
    FloatDPValue flt_time=cast_exact(time.get(dp));
    ARIADNE_ASSERT(decide(Real(flt_time)==time));
    if(this->size()==ht.discrete_time()) {
        this->_curves_ptr->push_back(HybridInterpolatedCurve(hpt.location(),hpt.space(),InterpolatedCurve(flt_time,hpt.point())));
    } else {
        (*this->_curves_ptr)[static_cast<unsigned int>(ht.discrete_time().get_si())].euclidean_set().insert(flt_time,hpt.point());
    }
}

Void Orbit<HybridApproximatePoint>::draw(CanvasInterface& canvas, const Set<DiscreteLocation>& locations, const Variables2d& axes) const {
    const Orbit<HybridApproximatePoint>& orbit=*this;
    for(Nat i=0; i!=orbit._curves_ptr->size(); ++i) {
        HybridInterpolatedCurve const& hcurve=this->_curves_ptr->at(i);
        if(locations.empty() || locations.contains(hcurve.location())) {
            RealSpace const& space=hcurve.space();
            if(valid_axis_variables(space,axes)) {
                hcurve.euclidean_set().draw(canvas,projection(space,axes));
            }
        }
    }
}

template<>
OutputStream&
operator<<(OutputStream& os, const Orbit< HybridApproximatePoint >& orb)
{
    return os << orb.curves();
}



struct Orbit<HybridGridCell>::Data {
    Data(const HybridGrid& grid)
        : initial(grid), reach(grid), intermediate(grid), final(grid) { }
    HybridGridTreePaving initial;
    HybridGridTreePaving reach;
    HybridGridTreePaving intermediate;
    HybridGridTreePaving final;
};

Orbit<HybridGridCell>::
Orbit(const HybridGridTreePaving& initial_set)
    : _data(new Data(initial_set.grid()))
{
    this->_data->initial=initial_set;
}

Orbit<HybridGridCell>::
Orbit(const HybridGridTreePaving& initial_set,
      const HybridGridTreePaving& reach_set,
      const HybridGridTreePaving& intermediate_set,
      const HybridGridTreePaving& final_set)
    : _data(new Data(initial_set.grid()))
{
    this->_data->initial=initial_set;
    this->_data->reach=reach_set;
    this->_data->intermediate=intermediate_set;
    this->_data->final=final_set;
}

HybridGridTreePaving const&
Orbit<HybridGridCell>::
initial() const
{
    return this->_data->initial;
}

HybridGridTreePaving const&
Orbit<HybridGridCell>::
reach() const
{
    return this->_data->reach;
}

HybridGridTreePaving const&
Orbit<HybridGridCell>::
intermediate() const
{
    return this->_data->intermediate;
}

HybridGridTreePaving const&
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


template<class X> HybridPoint<X>::HybridPoint(const DiscreteLocation& q, const Map<RealVariable,X>& x)
    : HybridBasicSet<Point<X>>(q,make_list(x.keys()),Point<X>(x.size()))
{
    SizeType i=0;
    for(auto iter=x.begin(); iter!=x.end(); ++iter, ++i) {
        this->point()[i]=iter->second;
    }
}

template<class X> HybridPoint<X>::HybridPoint(const DiscreteLocation& q, const List<Assignment<RealVariable,X>>& x)
    : HybridBasicSet<Point<X>>(q,left_hand_sides(x),Point<X>(x.size()))
{
    SizeType i=0;
    for(auto iter=x.begin(); iter!=x.end(); ++iter, ++i) {
        this->point()[i]=numeric_cast<X>(iter->right_hand_side());
    }
}

template<class X> HybridPoint<X>::HybridPoint(const DiscreteLocation& q, const InitializerList<Assignment<RealVariable,X>>& x)
    : HybridPoint<X>(q,List<Assignment<RealVariable,X>>(x))
{
}

template<class X> Map<RealVariable,X> HybridPoint<X>::values() const {
    Map<RealVariable,X> r;
    for(SizeType i=0; i!=this->space().dimension(); ++i) {
        r.insert(this->space()[i],this->point()[i]);
    }
    return r;
}



HybridBoxSet* HybridBoxSet::clone() const {
    return new HybridBoxSet(*this);
}
Set<DiscreteLocation> HybridBoxSet::locations() const {
    return { this->HybridVariablesBox<RealInterval>::location() };
}
Set<RealVariable> HybridBoxSet::variables(DiscreteLocation) const {
    return this->HybridVariablesBox<RealInterval>::variables();
}
RealSpace HybridBoxSet::space() const {
    return RealSpace( make_list(this->HybridVariablesBox<RealInterval>::variables()) );
}
RealSpace HybridBoxSet::space(DiscreteLocation loc) const {
    ARIADNE_ASSERT(this->location()==loc);
    return RealSpace( make_list(this->HybridVariablesBox<RealInterval>::variables()) );
}
SetInterface* HybridBoxSet::_euclidean_set(DiscreteLocation loc, RealSpace spc) const {
    ARIADNE_NOT_IMPLEMENTED; // FIXME: Box does not inherit from SetInterface...
}

using DyadicBox = Box<Interval<Dyadic>>;

LowerKleenean HybridBoxSet::is_empty() const {
    return this->euclidean_set(this->space()).is_empty();
}
LowerKleenean HybridBoxSet::inside(const HybridExactBoxes& hbxs) const {
    const DiscreteLocation& loc=this->location();
    if(hbxs.has_location(loc)) {
        const RealSpace& spc=hbxs.space(loc);
        return this->euclidean_set(spc).inside(static_cast<DyadicBox>(hbxs.euclidean_set(loc)));
    } else {
        return this->is_empty();
    }
}
LowerKleenean HybridBoxSet::inside(const HybridExactBox& hbx) const {
    return (hbx.location()==this->location() and this->euclidean_set(hbx.space()).inside(static_cast<DyadicBox>(hbx.euclidean_set()))) or this->is_empty();
}
LowerKleenean HybridBoxSet::overlaps(const HybridExactBox& hbx) const {
    return this->location()==hbx.location() and this->euclidean_set(hbx.space()).overlaps(static_cast<DyadicBox>(hbx.euclidean_set()));
}
LowerKleenean HybridBoxSet::separated(const HybridExactBox& hbx) const {
    return this->location()!=hbx.location() or this->euclidean_set(hbx.space()).separated(static_cast<DyadicBox>(hbx.euclidean_set()));
}
LowerKleenean HybridBoxSet::covers(const HybridExactBox& hbx) const {
    return (this->location()==hbx.location() and this->euclidean_set(hbx.space()).covers(static_cast<DyadicBox>(hbx.euclidean_set()))) or hbx.euclidean_set().is_empty();
}

HybridUpperBoxes HybridBoxSet::bounding_box() const {
    DiscreteLocation const& loc=this->location(); RealSpace spc(this->space()); RealBox bx=this->euclidean_set(spc);
    UpperBoxType bbx(bx,dp);
    HybridUpperBoxes res;
    ExactBoxType exbbx=reinterpret_cast<ExactBoxType const&>(bbx);  // FIXME: Should not need to convert to ExactBoxType here
    res.insert(loc,spc,exbbx);
    return res;
}

OutputStream& HybridBoxSet::write(OutputStream& os) const {
    return os << static_cast<HybridVariablesBox<RealInterval>const&>(*this);
}
Void HybridBoxSet::draw(CanvasInterface& c, const Set<DiscreteLocation>& qs, const Variables2d& vs) const {
    if(qs.empty() || qs.contains(this->location())) {
        RealSpace spc(List<RealVariable>(this->variables()));
        Projection2d prj(spc.dimension(),spc.index(vs.x_variable()),spc.index(vs.y_variable()));
        this->euclidean_set(spc).draw(c,prj);
    }
}




HybridConstraintSet::HybridConstraintSet()
    : _sets()
{
}

HybridConstraintSet::HybridConstraintSet(const DiscreteLocation& loc,
                                         const List<ContinuousPredicate>& cnstr)
{
    this->adjoin(loc,RealExpressionConstraintSet(cnstr));
}

HybridConstraintSet* HybridConstraintSet::clone() const {
    return new HybridConstraintSet(*this);
}

HybridConstraintSet& HybridConstraintSet::adjoin(const DiscreteLocation& loc, RealExpressionConstraintSet const& cs) {
    this->_sets.insert(loc,cs); return *this;
}

HybridConstraintSet& HybridConstraintSet::adjoin(const DiscreteLocation& loc, List<ContinuousPredicate> const& cnstr) {
    return this->adjoin(loc,RealExpressionConstraintSet(cnstr));
}

Set<DiscreteLocation> HybridConstraintSet::locations() const {
    return this->_sets.keys();
}

Bool HybridConstraintSet::has_location(DiscreteLocation loc) const {
    return this->_sets.has_key(loc);
}

Set<RealVariable> HybridConstraintSet::variables(DiscreteLocation loc) const {
    ARIADNE_ASSERT(this->_sets.has_key(loc));
    return _sets[loc].variables();
}

RealExpressionConstraintSet const& HybridConstraintSet::continuous_set(DiscreteLocation loc) const {
    ARIADNE_ASSERT(this->_sets.has_key(loc));
    return this->_sets[loc];
}

ConstraintSet const HybridConstraintSet::euclidean_set(DiscreteLocation loc, RealSpace spc) const {
    ARIADNE_ASSERT(this->_sets.has_key(loc));
    return ConstraintSet(this->_sets[loc].euclidean_set(spc));
}

RegularSetInterface* HybridConstraintSet::_euclidean_set(DiscreteLocation loc, RealSpace spc) const {
    return new ConstraintSet(this->euclidean_set(loc,spc));
}

LowerKleenean HybridConstraintSet::overlaps(const HybridExactBox& bx) const {
    if(this->_sets.has_key(bx.location())) {
        return this->_sets[bx.location()].euclidean_set(bx.space()).overlaps(bx.euclidean_set());
    } else {
        return false;
    }
}

LowerKleenean HybridConstraintSet::covers(const HybridExactBox& bx) const {
    if(this->_sets.has_key(bx.location())) {
        return this->_sets[bx.location()].euclidean_set(bx.space()).covers(bx.euclidean_set());
    } else {
        return bx.euclidean_set().is_empty();
    }
}

LowerKleenean HybridConstraintSet::separated(const HybridExactBox& bx) const {
    if(this->_sets.has_key(bx.location())) {
        return this->_sets[bx.location()].euclidean_set(bx.space()).separated(bx.euclidean_set());
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
                                                       const RealExpressionBoundedConstraintSet& set)
{
    this->adjoin(loc,set);
}

HybridBoundedConstraintSet::HybridBoundedConstraintSet(const DiscreteLocation& loc,
                                                       const InitializerList<RealVariableInterval>& bnd)
{
    this->adjoin(loc,RealExpressionBoundedConstraintSet(bnd));
}

HybridBoundedConstraintSet::HybridBoundedConstraintSet(const DiscreteLocation& loc,
                                                       const InitializerList<RealVariableInterval>& bnd,
                                                       const InitializerList<ContinuousPredicate>& cnstr)
{
    this->adjoin(loc,RealExpressionBoundedConstraintSet(bnd,cnstr));
}

HybridBoundedConstraintSet::HybridBoundedConstraintSet(const DiscreteLocation& loc,
                                                       const RealVariablesBox& bx)
{
    this->adjoin(loc,RealExpressionBoundedConstraintSet(bx));
}

HybridBoundedConstraintSet& HybridBoundedConstraintSet::adjoin(const DiscreteLocation& loc,
                                                               const RealExpressionBoundedConstraintSet& set)
{
    _sets.insert(loc,set); return *this;
}

HybridBoundedConstraintSet* HybridBoundedConstraintSet::clone() const {
    return new HybridBoundedConstraintSet(*this);
}

Set<RealVariable> HybridBoundedConstraintSet::variables(DiscreteLocation loc) const {
    ARIADNE_ASSERT(this->_sets.has_key(loc));
    return _sets[loc].variables();
}

BoundedConstraintSet const HybridBoundedConstraintSet::euclidean_set(DiscreteLocation loc, RealSpace spc) const {
    ARIADNE_ASSERT(this->_sets.has_key(loc));
    return BoundedConstraintSet(this->_sets[loc].euclidean_set(spc));
}

BoundedConstraintSet* HybridBoundedConstraintSet::_euclidean_set(DiscreteLocation loc, RealSpace spc) const {
    return new BoundedConstraintSet(this->euclidean_set(loc,spc));
}

LowerKleenean HybridBoundedConstraintSet::overlaps(const HybridExactBox& bx) const {
    if(this->_sets.has_key(bx.location())) {
        return this->_sets[bx.location()].euclidean_set(bx.space()).overlaps(bx.euclidean_set());
    } else {
        return false;
    }
}

LowerKleenean HybridBoundedConstraintSet::covers(const HybridExactBox& bx) const {
    if(this->_sets.has_key(bx.location())) {
        return this->_sets[bx.location()].euclidean_set(bx.space()).covers(bx.euclidean_set());
    } else {
        return bx.euclidean_set().is_empty();
    }
}

LowerKleenean HybridBoundedConstraintSet::separated(const HybridExactBox& bx) const {
    if(this->_sets.has_key(bx.location())) {
        return this->_sets[bx.location()].euclidean_set(bx.space()).separated(bx.euclidean_set());
    } else {
        return true;
    }
}

LowerKleenean HybridBoundedConstraintSet::inside(const HybridExactBoxes& bxs) const {
    LowerKleenean result=true;
    for(Map<DiscreteLocation,RealExpressionBoundedConstraintSet>::ConstIterator iter=this->_sets.begin(); iter!=this->_sets.end(); ++iter) {
        DiscreteLocation const& loc=iter->first;
        RealExpressionBoundedConstraintSet const& set = iter->second;
        Set<RealVariable> vars=set.variables();
        RealSpace const& spc = bxs[loc].space();
        ExactBoxType const& bx = bxs[loc].euclidean_set();
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

HybridExactBoxes HybridBoundedConstraintSet::bounding_box() const {
    HybridExactBoxes result;
    for(Map<DiscreteLocation,RealExpressionBoundedConstraintSet>::ConstIterator iter=this->_sets.begin(); iter!=this->_sets.end(); ++iter) {
        RealSpace spc=make_list(iter->second.variables());
        RealVariablesBox bnds=iter->second.bounds();
        RealBox rbx=bnds.euclidean_set(spc);
        ExactBoxType ebx=over_approximation(rbx);
        result.insert(iter->first,spc,ebx);
    }
    return result;
}

OutputStream& HybridBoundedConstraintSet::write(OutputStream& os) const {
    return os << "HybridBoundedConstraintSet( "<< this->_sets << " )";
}

Void HybridBoundedConstraintSet::draw(CanvasInterface& c, const Set<DiscreteLocation>& q, const Variables2d& p) const {
    for(auto loc : this->locations()) {
        if(q.empty() || q.contains(loc)) {
            Set<RealVariable> variables=this->variables(loc);
            RealSpace space(List<RealVariable>(variables.begin(),variables.end()));
            this->euclidean_set(loc,space).draw(c,projection(space,p));
        }
    }
}


HybridBoundedConstraintSet intersection(const HybridBoxesSet& hbxs, const HybridConstraintSet& hcs) {
    HybridBoundedConstraintSet res;
    for(auto loc : hbxs.locations()) {
        RealVariablesBox bx=hbxs.continuous_set(loc);
        RealExpressionConstraintSet cs=hcs.continuous_set(loc);
        res.adjoin(loc,intersection(bx,cs));
        //res.adjoin(loc,intersection(hbxs.continuous_set(loc),hcs.continuous_set(loc));
        if(hcs.has_location(loc)) { res.adjoin(loc,intersection(hbxs.continuous_set(loc),hcs.continuous_set(loc))); }
        else { res.adjoin(loc,hbxs.continuous_set(loc)); }
    }
    return res;
}

template<class EBS> Void HybridBasicSet<EBS>::adjoin_outer_approximation_to(HybridGridTreePaving& paving, Nat depth) const {
    if(this->space()==paving.space(this->location())) {
        paving[this->location()].adjoin_outer_approximation(this->euclidean_set(),depth);
    } else {
        ARIADNE_FAIL_MSG("HybridSet's state variables "<<this->space()<<
                         " do not match variables "<<paving.space()<<" of paving in location "<<this->location());
    }
}

template Void HybridBasicSet<Enclosure>::adjoin_outer_approximation_to(HybridGridTreePaving& paving, Nat depth) const;


template<class BS> Void draw_hybrid_basic_set(CanvasInterface& canvas, const DiscreteLocation& location, const Variables2d& axes, const HybridBasicSet<BS>& set)
{
    if(set.location()==location) {
        Projection2d projection(set.euclidean_set().dimension(),set.space().index(axes.x_variable()),set.space().index(axes.y_variable()));
        set.euclidean_set().draw(canvas,projection);
    }
}




template<class IVL> Void HybridBox<IVL>::draw(CanvasInterface& c, const Set<DiscreteLocation>& qs, const Variables2d& p) const {
    DiscreteLocation const& q=this->location();
    if(qs.empty() || qs.contains(q)) {
        Ariadne::draw_hybrid_basic_set(c,q,p,static_cast<HybridBasicSet<Box<IVL>>const&>(*this));
    }
}

template<class IVL> Void HybridBoxes<IVL>::draw(CanvasInterface& canvas, const Set<DiscreteLocation>& locations, const Variables2d& axis_variables) const {
    for(auto loc_iter=this->begin(); loc_iter!=this->end(); ++loc_iter) {
        if(locations.empty() || locations.contains(loc_iter->first)) {
            RealSpace const& space=this->space(loc_iter->first);
            Projection2d projection(space.dimension(),space.index(axis_variables.x_variable()),space.index(axis_variables.y_variable()));
            loc_iter->second.euclidean_set().draw(canvas,projection);
        }
    }
}



template<class DS, class HBS>
HybridSpaceSetConstIterator<DS,HBS>::
HybridSpaceSetConstIterator(const std::map<DiscreteLocation,DS>& map, const HybridSpace& spc, Bool end)
    : _loc_begin(map.begin()),
      _loc_end(map.end()),
      _loc_iter(end?_loc_end:_loc_begin),
      hspc(spc)
{
    if(_loc_iter!=_loc_end) {
        _bs_iter=_loc_iter->second.begin();
        this->increment_loc();
    }
}

template<class DS, class HBS> inline
Void
HybridSpaceSetConstIterator<DS,HBS>::increment_loc()
{
    while(_bs_iter==_loc_iter->second.end()) {
        ++_loc_iter;
        if(_loc_iter==_loc_end) { return; }
        _bs_iter=_loc_iter->second.begin();
    }
}


template<> String class_name<RealInterval>() { return "RealInterval"; }
template<> String class_name<InterpolatedCurve>() { return "InterpolatedCurve"; }
template<> String class_name<Box<RealInterval>>() { return "RealBox"; }
template<> String class_name<ExactBoxType>() { return "ExactFloatDPBox"; }

template<> String class_name<GridCell>() { return "GridCell"; }


// Instantiations
template class HybridPoint<ApproximateNumericType>;
template class HybridPoint<Real>;

template class HybridBasicSet<ExactBoxType>;
template class HybridBox<ExactIntervalType>;
template class HybridBoxes<ExactIntervalType>;

template class HybridBox<RealInterval>;


template class HybridSpaceSetConstIterator<GridTreePaving, HybridGridCell>;

} // namespace Ariadne
