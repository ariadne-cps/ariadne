/***************************************************************************
 *            hybrid/hybrid_set.cpp
 *
 *  Copyright  2008-20  Pieter Collins
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

#include "function/functional.hpp"
#include "config.hpp"
#include "hybrid/hybrid_expression_set.hpp"
#include "hybrid/hybrid_set.hpp"
#include "numeric/real.hpp"
#include "numeric/casts.hpp"
#include "numeric/rounding.hpp"
#include "expression-set/expression_set.hpp"
#include "set/function_set.hpp"
#include "hybrid/hybrid_space.hpp"
#include "symbolic/assignment.hpp"

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


/*
template<class X> HybridPoint<X>::HybridPoint(const DiscreteLocation& q, const Map<RealVariable,X>& x)
    : HybridPoint<X>(q,x.keys(),x.values())
{
}
*/

template<class X> HybridPoint<X>::HybridPoint(const DiscreteLocation& q, const List<Assignment<RealVariable,X>>& x)
    : HybridBasicSet<Point<X>>(q,left_hand_sides(x),Point<X>(x.size(),[&x](SizeType i){return x[i].right_hand_side();}))
{
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
    static_cast<void>(loc);
    ARIADNE_ASSERT(this->location()==loc);
    return RealSpace( make_list(this->HybridVariablesBox<RealInterval>::variables()) );
}
EuclideanSetInterfaceBase* HybridBoxSet::_euclidean_set(DiscreteLocation, RealSpace) const {
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
    ExactBoxType exbbx=cast_exact_box(bbx);
    res.insert(loc,spc,exbbx);
    return res;
}

OutputStream& HybridBoxSet::_write(OutputStream& os) const {
    return os << static_cast<HybridVariablesBox<RealInterval>const&>(*this);
}
HybridUpperBoxes HybridValidatedConstrainedImageSet::bounding_box() const {
    return this->Base::bounding_boxes();
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

HybridUpperBoxes HybridBoundedConstraintSet::bounding_box() const {
    HybridUpperBoxes result;
    for(Map<DiscreteLocation,RealExpressionBoundedConstraintSet>::ConstIterator iter=this->_sets.begin(); iter!=this->_sets.end(); ++iter) {
        RealSpace spc=make_list(iter->second.variables());
        RealVariablesBox bnds=iter->second.bounds();
        RealBox rbx=bnds.euclidean_set(spc);
        UpperBoxType ebx=over_approximation(rbx);
        result.insert(iter->first,spc,ebx);
    }
    return result;
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

template<class DS, class HBS>
HybridSpaceSetConstIterator<DS,HBS>::
HybridSpaceSetConstIterator(const Map<DiscreteLocation,DS>& map, const HybridSpace& spc, Bool end)
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

Vector<FloatDPApproximation> evaluate(ApproximateVectorMultivariateFunction const& f, Vector<FloatDPApproximation> const& v);

template<> String class_name<InterpolatedCurve>() { return "InterpolatedCurve"; }
template<> String class_name<Box<RealInterval>>() { return "RealBox"; }
template<> String class_name<ExactBoxType>() { return "ExactFloatDPBox"; }
template<> String class_name<UpperBoxType>() { return "UpperFloatDPBox"; }

template<> String class_name<GridCell>() { return "GridCell"; }


// Instantiations of mathematical hybrid points.
template class HybridPoint<FloatDPApproximation>;
template class HybridPoint<Real>;

template class HybridSpaceSetConstIterator<GridTreePaving, HybridGridCell>;




} // namespace Ariadne
